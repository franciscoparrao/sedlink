//! Index of Connectivity (Borselli 2008 / Cavalli 2013).
//!
//! IC = `log10(D_up / D_dn)`
//!
//! where:
//! - **`D_up`** = W̄ · S̄ · √A — upslope sediment delivery potential.
//!   `W̄` is the mean weighting factor over the upslope contributing area,
//!   `S̄` is the mean slope gradient (m/m) over the upslope contributing
//!   area, and `A` is the contributing area (flow accumulation × cell area).
//! - **`D_dn`** = Σ (`d_i` / (`W_i` · `S_i`)) — downslope impedance, summed
//!   along the flow path from the cell to the nearest stream (D8 single
//!   direction, or the D∞ primary-receiver path), using each path cell's
//!   local weight and slope gradient.
//!
//! The computation is generic over [`FlowNetwork`], so it accepts both
//! [`Network`](crate::Network) (D8) and [`DinfNetwork`](crate::DinfNetwork)
//! (D∞, the routing used by Cavalli et al. 2013).
//!
//! Following Cavalli et al. (2013) / `SedInConnect`, the slope gradient is
//! tan(θ) clamped to `[min_slope, 1.0]` to avoid division by zero and to
//! bound the impedance on very steep cells.
//!
//! Stream cells (flow accumulation ≥ threshold) have `D_dn` → 0, giving
//! IC → +clamp. (`SedInConnect` masks stream cells as `NoData` instead; keep
//! this in mind when comparing outputs.)
//!
//! When [`ConnectivityParams::targets`] is set, IC is computed relative to
//! an arbitrary target mask instead (Cavalli et al. 2013 / `SedInConnect`
//! "targets" version): `D_dn` is summed along the flow path to the nearest
//! target cell, the threshold-based stream network is not used, and cells
//! that never drain to a target get IC = NaN.
//!
//! D8 direction codes follow the `SurtGIS` convention: 1=E, 2=NE, 3=N, 4=NW,
//! 5=W, 6=SW, 7=S, 8=SE (counter-clockwise from East).

use ndarray::Array2;
use surtgis_core::raster::Raster;

use crate::{ConnectivityParams, FlowNetwork, SedlinkError};

/// Upper bound for the slope gradient (m/m), following Cavalli et al.
/// (2013) / `SedInConnect`: slopes steeper than 1.0 are capped so the
/// downslope impedance stays bounded.
const MAX_SLOPE: f64 = 1.0;

/// Result of a connectivity index computation.
#[derive(Debug, Clone)]
pub struct ConnectivityIndex {
    /// IC values per cell (row-major). NaN for `NoData` cells.
    pub ic: Array2<f64>,
    /// `D_up` component per cell.
    pub d_up: Array2<f64>,
    /// `D_dn` component per cell.
    pub d_dn: Array2<f64>,
    /// Destination mask: `true` for cells where `D_dn` = 0. Threshold-based
    /// stream cells by default, or the target mask when
    /// [`ConnectivityParams::targets`] is set.
    pub is_stream: Array2<bool>,
}

impl ConnectivityIndex {
    /// Compute the Index of Connectivity from a flow network and DEM.
    ///
    /// # Arguments
    ///
    /// * `net` - Pre-built flow network (D8 [`Network`](crate::Network) or
    ///   D∞ [`DinfNetwork`](crate::DinfNetwork)).
    /// * `dem` - DEM raster (must match the network's grid).
    /// * `params` - Connectivity parameters.
    ///
    /// # Errors
    ///
    /// Returns [`SedlinkError::InvalidParam`] if a parameter is outside its
    /// valid range, and [`SedlinkError::GridMismatch`] if the DEM or the
    /// weighting raster dimensions don't match the network.
    ///
    /// # Panics
    ///
    /// Panics only if internal grid-shape invariants are violated, which
    /// cannot happen for inputs that pass the validations above.
    pub fn compute<F: FlowNetwork>(
        net: &F,
        dem: &Raster<f32>,
        params: &ConnectivityParams,
    ) -> Result<Self, SedlinkError> {
        params.validate()?;

        let (rows, cols) = dem.shape();
        if rows != net.rows() || cols != net.cols() {
            return Err(SedlinkError::GridMismatch {
                expected_rows: net.rows(),
                expected_cols: net.cols(),
                got_rows: rows,
                got_cols: cols,
            });
        }

        if let Some(wt) = &params.weight.raster {
            let (wr, wc) = wt.dim();
            if wr != rows || wc != cols {
                return Err(SedlinkError::GridMismatch {
                    expected_rows: rows,
                    expected_cols: cols,
                    got_rows: wr,
                    got_cols: wc,
                });
            }
        }

        if let Some(tg) = &params.targets {
            let (tr, tc) = tg.dim();
            if tr != rows || tc != cols {
                return Err(SedlinkError::GridMismatch {
                    expected_rows: rows,
                    expected_cols: cols,
                    got_rows: tr,
                    got_cols: tc,
                });
            }
        }

        let n = rows * cols;
        let cellsize = net.cellsize();
        let threshold = params.stream_threshold;
        let min_slope = params.min_slope;
        let clamp = params.clamp;
        let weight = &params.weight;
        let min_weight = weight.min_value;
        let slope = net.slope_slice();
        let flow_acc = net.flow_acc_slice();

        // Per-cell weighting factor and slope gradient (tan θ, clamped to
        // [min_slope, MAX_SLOPE] following Cavalli et al. 2013).
        let mut w_local = vec![1.0f64; n];
        let mut s_local = vec![0.0f64; n];
        for idx in 0..n {
            s_local[idx] = f64::from(slope[idx]).tan().clamp(min_slope, MAX_SLOPE);
            if let Some(wt) = &weight.raster {
                w_local[idx] = wt[(idx / cols, idx % cols)].max(min_weight);
            }
        }

        // Destination mask: explicit targets if provided, else stream cells
        // by flow accumulation threshold.
        let is_stream: Vec<bool> = match &params.targets {
            Some(tg) => (0..n).map(|i| tg[(i / cols, i % cols)]).collect(),
            None => (0..n).map(|i| flow_acc[i] >= threshold).collect(),
        };

        // Compute D_dn by tracing downstream paths with memoization.
        let d_dn = Self::compute_d_dn(net, &is_stream, &w_local, &s_local);

        // Upslope sums of W and S for the D_up means (W̄, S̄).
        let sum_w = net.accumulate_upslope(&w_local);
        let sum_s = net.accumulate_upslope(&s_local);

        // Compute D_up and IC.
        let mut ic_data = Array2::<f64>::from_elem((rows, cols), f64::NAN);
        let mut d_up_data = Array2::<f64>::from_elem((rows, cols), f64::NAN);
        let mut d_dn_data = Array2::<f64>::from_elem((rows, cols), f64::NAN);

        for row in 0..rows {
            for col in 0..cols {
                let idx = row * cols + col;

                if net.is_solid_idx(idx) {
                    continue;
                }

                // D_up = W̄ · S̄ · √A, with means over the upslope
                // contributing area (flow_acc = cell count incl. self).
                let acc = flow_acc[idx];
                let w_mean = sum_w[idx] / acc;
                let s_mean = sum_s[idx] / acc;
                let a_up = (acc * cellsize * cellsize).sqrt();
                let d_up = w_mean * s_mean * a_up;

                let dn = d_dn[idx];

                d_up_data[(row, col)] = d_up;
                d_dn_data[(row, col)] = dn;

                if dn.is_nan() {
                    // Unreachable channel.
                    continue;
                }

                if dn < 1e-10 {
                    // At channel: IC → +clamp.
                    ic_data[(row, col)] = clamp;
                } else {
                    let ic = (d_up / dn).log10();
                    ic_data[(row, col)] = ic.clamp(-clamp, clamp);
                }
            }
        }

        let is_stream_arr = Array2::from_shape_vec((rows, cols), is_stream).unwrap();

        Ok(ConnectivityIndex {
            ic: ic_data,
            d_up: d_up_data,
            d_dn: d_dn_data,
            is_stream: is_stream_arr,
        })
    }

    /// Compute `D_dn` for each cell by tracing downstream paths.
    ///
    /// Uses memoization: once a cell's `D_dn` is computed, it is cached
    /// and reused by upstream cells that trace through it. Cells whose
    /// path ends in a pit (no channel reachable) are marked unreachable
    /// so they are never re-traced.
    fn compute_d_dn<F: FlowNetwork>(
        net: &F,
        is_stream: &[bool],
        w_local: &[f64],
        s_local: &[f64],
    ) -> Vec<f64> {
        let n = net.len();
        let mut d_dn = vec![f64::NAN; n];
        // Cells known to never reach a channel (path ends in a pit).
        let mut unreachable = vec![false; n];

        // Stream cells: D_dn = 0 (they are at the channel).
        for idx in 0..n {
            if is_stream[idx] {
                d_dn[idx] = 0.0;
            }
        }

        // Trace each non-stream cell downstream.
        for start in 0..n {
            if !d_dn[start].is_nan() || unreachable[start] {
                continue;
            }
            if net.is_solid_idx(start) {
                continue;
            }

            // Trace path, collecting contributions.
            let mut path: Vec<(usize, f64)> = Vec::new();
            let mut cur = start;
            let max_steps = n;
            let mut steps = 0;

            loop {
                // If we reached a cell with known D_dn (or one already known
                // to be unreachable), resolve the path.
                if !d_dn[cur].is_nan() || unreachable[cur] {
                    break;
                }

                // Pit: no downstream. Cell is unreachable.
                let Some((ds, step_dist)) = net.downstream_step(cur) else {
                    break;
                };

                let contribution = step_dist / (w_local[cur] * s_local[cur]);
                path.push((cur, contribution));

                cur = ds;
                steps += 1;
                if steps > max_steps {
                    break;
                }
            }

            if d_dn[cur].is_nan() {
                // Terminal cell never reaches a channel: mark the whole
                // path (and the terminal) unreachable to avoid re-tracing.
                unreachable[cur] = true;
                for &(idx, _) in &path {
                    unreachable[idx] = true;
                }
            } else {
                // Resolve the path backwards, accumulating D_dn.
                let mut cumulative = d_dn[cur];
                for &(idx, contribution) in path.iter().rev() {
                    cumulative += contribution;
                    d_dn[idx] = cumulative;
                }
            }
        }

        d_dn
    }

    /// Export IC as a `Raster<f64>` with the network's geotransform.
    #[must_use]
    pub fn ic_raster<F: FlowNetwork>(&self, net: &F) -> Raster<f64> {
        let mut r = Raster::from_array(self.ic.clone());
        r.set_transform(*net.transform());
        r.set_nodata(Some(f64::NAN));
        r
    }

    /// Export `D_up` as a `Raster<f64>`.
    #[must_use]
    pub fn d_up_raster<F: FlowNetwork>(&self, net: &F) -> Raster<f64> {
        let mut r = Raster::from_array(self.d_up.clone());
        r.set_transform(*net.transform());
        r
    }

    /// Export `D_dn` as a `Raster<f64>`.
    #[must_use]
    pub fn d_dn_raster<F: FlowNetwork>(&self, net: &F) -> Raster<f64> {
        let mut r = Raster::from_array(self.d_dn.clone());
        r.set_transform(*net.transform());
        r
    }
}
