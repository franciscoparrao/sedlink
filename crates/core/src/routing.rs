//! Sediment routing with a distance-decay sediment delivery ratio (SDR).
//!
//! Routes a per-cell sediment source (e.g. RUSLE soil loss) downslope to
//! the channel network and accumulates the delivered load along the
//! channels:
//!
//! 1. **Travel distance** `L_i`: length (m) of the flow path from cell
//!    `i` to the first stream cell, traced on the flow network.
//! 2. **Delivery ratio** — two models:
//!    - [`SedimentRouting::compute`]: `SDR_i = exp(−q · L_i / L_ref)`,
//!      the distance-decay family of hillslope delivery models (cf.
//!      SEDD, Ferro & Porto 2000). `q` is
//!      [`RoutingParams::sdr_exponent`] and `L_ref` is
//!      [`RoutingParams::reference_length`]. Paths longer than
//!      [`RoutingParams::max_distance`] deliver nothing (`SDR = 0`).
//!    - [`SedimentRouting::compute_from_ic`]:
//!      `SDR_i = sdr_max / (1 + exp((ic0 − IC_i) / k))`, the IC-based
//!      sigmoid of Vigiak et al. (2012) as used in `InVEST`, driven by a
//!      previously computed [`ConnectivityIndex`].
//!
//!    Stream cells have `SDR = 1` (already in the channel).
//! 3. **Hillslope delivery**: each cell's `source_i · SDR_i` is credited
//!    to the stream cell where its flow path enters the channel.
//! 4. **Channel flux**: deliveries are accumulated downstream along the
//!    channel network (in-channel transport is assumed conservative), so
//!    a channel cell's flux is the sediment yield of its whole upstream
//!    area.
//!
//! Units follow the source raster (e.g. t/cell/yr). With no source
//! raster, a unit source of 1.0 per cell is used, so the channel flux
//! counts "SDR-weighted contributing cells".
//!
//! Cells whose flow path ends in a pit before reaching a channel are
//! unreachable: their SDR and travel distance are NaN and they deliver
//! nothing.
//!
//! References:
//! - Ferro, V., & Porto, P. (2000). Sediment delivery distributed (SEDD)
//!   model. *Journal of Hydrologic Engineering*, 5(4).
//! - Vigiak, O., et al. (2012). Comparison of conceptual landscape
//!   metrics to define hillslope-scale sediment delivery ratio.
//!   *Geomorphology*, 138(1), 74–88.

use ndarray::Array2;

use crate::{ConnectivityIndex, FlowNetwork, IcSdrParams, RoutingParams, SedlinkError};

/// Result of a sediment routing computation.
#[derive(Debug, Clone)]
pub struct SedimentRouting {
    /// Sediment delivery ratio per cell in `[0, 1]`. 1.0 on stream cells,
    /// 0.0 beyond `max_distance`, NaN for `NoData` and unreachable cells.
    pub sdr: Array2<f64>,
    /// Flow-path travel distance to the channel (m). 0 on stream cells,
    /// NaN for `NoData` and unreachable cells.
    pub dist_to_stream: Array2<f64>,
    /// Sediment mass entering the channel at each stream cell (source
    /// units). 0 on non-stream cells, NaN for `NoData` cells.
    pub hillslope_delivery: Array2<f64>,
    /// Sediment flux accumulated along the channel network (source
    /// units). NaN off-channel and on `NoData` cells.
    pub channel_flux: Array2<f64>,
    /// Stream mask: `true` for cells identified as channels.
    pub is_stream: Array2<bool>,
}

impl SedimentRouting {
    /// Route sediment from a per-cell source to the channel network.
    ///
    /// # Arguments
    ///
    /// * `net` - Pre-built flow network (D8 [`Network`](crate::Network) or
    ///   D∞ [`DinfNetwork`](crate::DinfNetwork)).
    /// * `source` - Optional per-cell sediment source (e.g. RUSLE soil
    ///   loss, t/cell/yr). `None` uses a unit source of 1.0 per cell.
    /// * `stream_threshold` - Flow accumulation threshold (cell count)
    ///   delineating stream cells.
    /// * `params` - Routing parameters (SDR decay and truncation).
    ///
    /// # Errors
    ///
    /// Returns [`SedlinkError::InvalidParam`] if `stream_threshold` or a
    /// routing parameter is outside its valid range, and
    /// [`SedlinkError::GridMismatch`] if the source raster dimensions
    /// don't match the network.
    pub fn compute<F: FlowNetwork>(
        net: &F,
        source: Option<&Array2<f64>>,
        stream_threshold: f64,
        params: &RoutingParams,
    ) -> Result<Self, SedlinkError> {
        params.validate()?;
        if !stream_threshold.is_finite() || stream_threshold <= 0.0 {
            return Err(SedlinkError::InvalidParam {
                name: "stream_threshold",
                value: stream_threshold,
                constraint: "finite and > 0",
            });
        }

        let n = net.len();
        let flow_acc = net.flow_acc_slice();
        let is_stream: Vec<bool> = (0..n).map(|i| flow_acc[i] >= stream_threshold).collect();

        let (dist, entry) = Self::trace_to_channel(net, &is_stream);

        // SDR per cell: exp(−q · L / L_ref), truncated at max_distance.
        let q = params.sdr_exponent;
        let l_ref = params.reference_length;
        let max_dist = params.max_distance;

        let mut sdr = vec![f64::NAN; n];
        for idx in 0..n {
            if net.is_solid_idx(idx) {
                continue;
            }
            if is_stream[idx] {
                sdr[idx] = 1.0;
            } else if dist[idx].is_nan() {
                // Unreachable: stays NaN.
            } else if dist[idx] > max_dist {
                sdr[idx] = 0.0;
            } else {
                sdr[idx] = (-q * dist[idx] / l_ref).exp();
            }
        }

        Self::finish(net, is_stream, sdr, dist, &entry, source)
    }

    /// Route sediment using the IC-based SDR of Vigiak et al. (2012) /
    /// `InVEST`: `SDR = sdr_max / (1 + exp((ic0 − IC) / k))`.
    ///
    /// The stream mask is taken from the [`ConnectivityIndex`] (computed
    /// with its own `stream_threshold`), so routing and IC are always
    /// consistent. Cells with NaN IC (`NoData` or unreachable) get NaN
    /// SDR and deliver nothing; stream cells have SDR = 1.
    ///
    /// # Arguments
    ///
    /// * `net` - The same flow network the IC was computed on.
    /// * `ic` - Connectivity index result.
    /// * `source` - Optional per-cell sediment source; `None` uses 1.0.
    /// * `params` - IC-SDR calibration parameters.
    ///
    /// # Errors
    ///
    /// Returns [`SedlinkError::InvalidParam`] if a parameter is outside
    /// its valid range, and [`SedlinkError::GridMismatch`] if the IC or
    /// source raster dimensions don't match the network.
    pub fn compute_from_ic<F: FlowNetwork>(
        net: &F,
        ic: &ConnectivityIndex,
        source: Option<&Array2<f64>>,
        params: &IcSdrParams,
    ) -> Result<Self, SedlinkError> {
        params.validate()?;

        let rows = net.rows();
        let cols = net.cols();
        let (ir, ic_cols) = ic.ic.dim();
        if ir != rows || ic_cols != cols {
            return Err(SedlinkError::GridMismatch {
                expected_rows: rows,
                expected_cols: cols,
                got_rows: ir,
                got_cols: ic_cols,
            });
        }

        let n = rows * cols;
        let is_stream: Vec<bool> = ic.is_stream.iter().copied().collect();
        let (dist, entry) = Self::trace_to_channel(net, &is_stream);

        let mut sdr = vec![f64::NAN; n];
        for idx in 0..n {
            if net.is_solid_idx(idx) {
                continue;
            }
            if is_stream[idx] {
                sdr[idx] = 1.0;
            } else {
                let ic_val = ic.ic[(idx / cols, idx % cols)];
                if !ic_val.is_nan() {
                    sdr[idx] = params.sdr_max / (1.0 + ((params.ic0 - ic_val) / params.k).exp());
                }
            }
        }

        Self::finish(net, is_stream, sdr, dist, &entry, source)
    }

    /// Shared tail of both SDR models: hillslope delivery, channel flux,
    /// and output assembly. Validates the source raster dimensions.
    fn finish<F: FlowNetwork>(
        net: &F,
        is_stream: Vec<bool>,
        sdr: Vec<f64>,
        dist: Vec<f64>,
        entry: &[usize],
        source: Option<&Array2<f64>>,
    ) -> Result<Self, SedlinkError> {
        let rows = net.rows();
        let cols = net.cols();
        let n = rows * cols;

        if let Some(src) = source {
            let (sr, sc) = src.dim();
            if sr != rows || sc != cols {
                return Err(SedlinkError::GridMismatch {
                    expected_rows: rows,
                    expected_cols: cols,
                    got_rows: sr,
                    got_cols: sc,
                });
            }
        }

        // Hillslope delivery: credit source·SDR to the channel entry cell.
        let mut delivered = vec![0.0f64; n];
        for idx in 0..n {
            if net.is_solid_idx(idx) {
                continue;
            }
            let src = source.map_or(1.0, |s| s[(idx / cols, idx % cols)]);
            if is_stream[idx] {
                delivered[idx] += src;
            } else if entry[idx] != usize::MAX && sdr[idx] > 0.0 {
                delivered[entry[idx]] += src * sdr[idx];
            }
        }

        // Channel flux: accumulate deliveries downstream. Deliveries are
        // nonzero only at stream cells, so the upslope sum at a channel
        // cell is the total load entering the channel above it.
        let flux = net.accumulate_upslope(&delivered);

        // Assemble outputs.
        let shape = (rows, cols);
        let mut sdr_arr = Array2::from_shape_vec(shape, sdr).unwrap();
        let mut dist_arr = Array2::from_shape_vec(shape, dist).unwrap();
        let mut delivered_arr = Array2::from_shape_vec(shape, delivered).unwrap();
        let mut flux_arr = Array2::from_shape_vec(shape, flux).unwrap();
        let stream_arr = Array2::from_shape_vec(shape, is_stream).unwrap();

        for row in 0..rows {
            for col in 0..cols {
                let idx = row * cols + col;
                if net.is_solid_idx(idx) {
                    sdr_arr[(row, col)] = f64::NAN;
                    dist_arr[(row, col)] = f64::NAN;
                    delivered_arr[(row, col)] = f64::NAN;
                    flux_arr[(row, col)] = f64::NAN;
                } else if !stream_arr[(row, col)] {
                    // Flux is only meaningful along the channel.
                    flux_arr[(row, col)] = f64::NAN;
                }
            }
        }

        Ok(SedimentRouting {
            sdr: sdr_arr,
            dist_to_stream: dist_arr,
            hillslope_delivery: delivered_arr,
            channel_flux: flux_arr,
            is_stream: stream_arr,
        })
    }

    /// Trace each cell's flow path to the channel, returning the travel
    /// distance (m) and the stream cell where the path enters the channel.
    ///
    /// Memoized like `ConnectivityIndex::compute_d_dn`: resolved cells are
    /// reused by upstream cells, and paths ending in a pit are marked
    /// unreachable (distance NaN, entry `usize::MAX`) without re-tracing.
    fn trace_to_channel<F: FlowNetwork>(net: &F, is_stream: &[bool]) -> (Vec<f64>, Vec<usize>) {
        let n = net.len();
        let mut dist = vec![f64::NAN; n];
        let mut entry = vec![usize::MAX; n];
        let mut unreachable = vec![false; n];

        for idx in 0..n {
            if is_stream[idx] {
                dist[idx] = 0.0;
                entry[idx] = idx;
            }
        }

        for start in 0..n {
            if !dist[start].is_nan() || unreachable[start] || net.is_solid_idx(start) {
                continue;
            }

            let mut path: Vec<(usize, f64)> = Vec::new();
            let mut cur = start;
            let max_steps = n;
            let mut steps = 0;

            loop {
                if !dist[cur].is_nan() || unreachable[cur] {
                    break;
                }
                let Some((ds, step_dist)) = net.downstream_step(cur) else {
                    break;
                };
                path.push((cur, step_dist));
                cur = ds;
                steps += 1;
                if steps > max_steps {
                    break;
                }
            }

            if dist[cur].is_nan() {
                unreachable[cur] = true;
                for &(idx, _) in &path {
                    unreachable[idx] = true;
                }
            } else {
                let mut cumulative = dist[cur];
                let e = entry[cur];
                for &(idx, step) in path.iter().rev() {
                    cumulative += step;
                    dist[idx] = cumulative;
                    entry[idx] = e;
                }
            }
        }

        (dist, entry)
    }
}
