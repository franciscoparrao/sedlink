//! Multiple flow direction (MFD / FD8) routing.
//!
//! Flow from each cell is partitioned among **all** downslope neighbours
//! with fractions proportional to `(tan β)^p` (Freeman 1991), where `p`
//! is the Holmgren (1994) exponent: `p = 1` gives the most dispersive
//! partition, and as `p → ∞` the scheme converges to D8 single flow.
//! The default `p = 1.1` follows Freeman.
//!
//! Cells with no downslope neighbour drain through the flat-resolution
//! direction (single receiver, fraction 1); cells with neither are pits.
//! All edges descend (or follow the acyclic flat gradient), so the
//! network is acyclic by construction.
//!
//! Receivers are stored in CSR form (`offsets` + `entries`), keeping
//! memory linear in the number of edges for basin-scale DEMs.
//!
//! References:
//! - Freeman, T.G. (1991). Calculating catchment area with divergent
//!   flow based on a regular grid. *Computers & Geosciences*, 17(3).
//! - Quinn, P., et al. (1991). The prediction of hillslope flow paths…
//!   *Hydrological Processes*, 5(1).
//! - Holmgren, P. (1994). Multiple flow direction algorithms for runoff
//!   modelling in grid based elevation models. *Hydrological
//!   Processes*, 8(4).

use surtgis_core::raster::Raster;
use surtgis_core::{GeoTransform, Result as GstResult};

use crate::network::{D8_DISTANCE, D8_OFFSETS};
use crate::{Network, SedlinkError};

/// Default Holmgren exponent (Freeman 1991).
pub const MFD_DEFAULT_EXPONENT: f64 = 1.1;

/// A multiple-flow-direction (FD8) network derived from a DEM.
pub struct MfdNetwork {
    rows: usize,
    cols: usize,
    cellsize: f64,
    transform: GeoTransform,
    /// CSR: cell `i`'s receivers are `entries[offsets[i]..offsets[i+1]]`.
    offsets: Vec<usize>,
    /// `(receiver_cell, fraction)`; fractions of a cell sum to 1.
    entries: Vec<(u32, f32)>,
    /// Flow accumulation (fractional cell count, including self).
    /// NaN for `NoData` cells.
    flow_acc: Vec<f64>,
    /// Slope (radians) per cell, from the original DEM.
    slope: Vec<f32>,
    /// Solid mask: `true` cells are `NoData` walls.
    solid: Vec<bool>,
}

impl MfdNetwork {
    /// Build an MFD network with the default exponent (1.1).
    ///
    /// # Errors
    ///
    /// Returns [`SedlinkError::EmptyGrid`], [`SedlinkError::RotatedGrid`],
    /// or [`SedlinkError::NonSquareCells`] if the DEM grid is invalid.
    pub fn from_dem(dem: &Raster<f32>) -> Result<Self, SedlinkError> {
        Self::from_dem_with_exponent(dem, MFD_DEFAULT_EXPONENT)
    }

    /// Build an MFD network with an explicit Holmgren exponent.
    ///
    /// # Errors
    ///
    /// As [`from_dem`](Self::from_dem), plus
    /// [`SedlinkError::InvalidParam`] if `exponent` is not finite and
    /// positive.
    pub fn from_dem_with_exponent(dem: &Raster<f32>, exponent: f64) -> Result<Self, SedlinkError> {
        if !exponent.is_finite() || exponent <= 0.0 {
            return Err(SedlinkError::InvalidParam {
                name: "mfd_exponent",
                value: exponent,
                constraint: "finite and > 0",
            });
        }
        let (rows, cols, t, pw, z, solid) = crate::DinfNetwork::validate_dem(dem)?;
        let filled = Network::priority_fill(&z, &solid, rows, cols);
        let flats = crate::flats::FlatOffsets::compute(&filled, &solid, rows, cols);

        let (offsets, entries) =
            Self::compute_receivers(&filled, &solid, rows, cols, exponent, flats.as_ref());
        let flow_acc = Self::compute_flow_acc(&offsets, &entries, &solid, rows * cols);
        let slope = Network::compute_slope(&z, &solid, rows, cols, pw);

        Ok(Self {
            rows,
            cols,
            cellsize: pw,
            transform: t,
            offsets,
            entries,
            flow_acc,
            slope,
            solid,
        })
    }

    /// Downslope receivers with `(tan β)^p` fractions; flat cells get
    /// their resolved single receiver.
    fn compute_receivers(
        z: &[f32],
        solid: &[bool],
        rows: usize,
        cols: usize,
        exponent: f64,
        flats: Option<&crate::flats::FlatOffsets>,
    ) -> (Vec<usize>, Vec<(u32, f32)>) {
        let n = rows * cols;
        let mut offsets = Vec::with_capacity(n + 1);
        let mut entries: Vec<(u32, f32)> = Vec::new();
        offsets.push(0);

        let mut weights: Vec<(usize, f64)> = Vec::with_capacity(8);
        for idx in 0..n {
            if !solid[idx] {
                let (r, c) = ((idx / cols) as isize, (idx % cols) as isize);
                let here = f64::from(z[idx]);
                weights.clear();
                for (k, &(dr, dc)) in D8_OFFSETS.iter().enumerate() {
                    let (nr, nc) = (r + dr, c + dc);
                    if nr < 0 || nc < 0 || nr as usize >= rows || nc as usize >= cols {
                        continue;
                    }
                    let nidx = nr as usize * cols + nc as usize;
                    if solid[nidx] {
                        continue;
                    }
                    let drop = here - f64::from(z[nidx]);
                    if drop > 0.0 {
                        weights.push((nidx, (drop / D8_DISTANCE[k]).powf(exponent)));
                    }
                }

                if weights.is_empty() {
                    // Flat (or pit): single resolved receiver if available.
                    if let Some(f) = flats
                        && let Some(d) = f.direction(idx, z, solid)
                    {
                        let (dr, dc) = D8_OFFSETS[(d - 1) as usize];
                        let nidx = (r + dr) as usize * cols + (c + dc) as usize;
                        entries.push((nidx as u32, 1.0));
                    }
                } else {
                    let total: f64 = weights.iter().map(|&(_, w)| w).sum();
                    for &(nidx, w) in &weights {
                        entries.push((nidx as u32, (w / total) as f32));
                    }
                }
            }
            offsets.push(entries.len());
        }

        (offsets, entries)
    }

    /// Topological (in-degree) fractional accumulation, counting self.
    fn compute_flow_acc(
        offsets: &[usize],
        entries: &[(u32, f32)],
        solid: &[bool],
        n: usize,
    ) -> Vec<f64> {
        let mut acc = vec![1.0f64; n];
        let mut in_degree = vec![0u32; n];
        for &(recv, _) in entries {
            in_degree[recv as usize] += 1;
        }

        let mut queue: Vec<usize> = (0..n).filter(|&i| !solid[i] && in_degree[i] == 0).collect();
        while let Some(idx) = queue.pop() {
            let outflow = acc[idx];
            for &(recv, frac) in &entries[offsets[idx]..offsets[idx + 1]] {
                let recv = recv as usize;
                acc[recv] += outflow * f64::from(frac);
                in_degree[recv] -= 1;
                if in_degree[recv] == 0 {
                    queue.push(recv);
                }
            }
        }

        for idx in 0..n {
            if solid[idx] {
                acc[idx] = f64::NAN;
            }
        }
        acc
    }

    /// Receivers of a cell as `(receiver_index, fraction)` pairs.
    #[must_use]
    pub fn receivers(&self, idx: usize) -> &[(u32, f32)] {
        &self.entries[self.offsets[idx]..self.offsets[idx + 1]]
    }

    /// Flow accumulation (fractional cell count, incl. self) at (row, col).
    #[must_use]
    pub fn flow_acc(&self, row: usize, col: usize) -> f64 {
        self.flow_acc[row * self.cols + col]
    }

    /// Build a `Raster<f64>` of flow accumulation with the same geotransform.
    ///
    /// # Errors
    ///
    /// Propagates the underlying raster construction error (cannot occur
    /// for a well-formed network, whose grids always match its shape).
    pub fn flow_acc_raster(&self) -> GstResult<Raster<f64>> {
        let mut r = Raster::from_vec(self.flow_acc.clone(), self.rows, self.cols)?;
        r.set_transform(self.transform);
        r.set_nodata(Some(f64::NAN));
        Ok(r)
    }
}

impl crate::FlowNetwork for MfdNetwork {
    fn rows(&self) -> usize {
        self.rows
    }

    fn cols(&self) -> usize {
        self.cols
    }

    fn cellsize(&self) -> f64 {
        self.cellsize
    }

    fn transform(&self) -> &GeoTransform {
        &self.transform
    }

    fn is_solid_idx(&self, idx: usize) -> bool {
        self.solid[idx]
    }

    fn flow_acc_slice(&self) -> &[f64] {
        &self.flow_acc
    }

    fn slope_slice(&self) -> &[f32] {
        &self.slope
    }

    fn downstream_step(&self, idx: usize) -> Option<(usize, f64)> {
        let recv = self.receivers(idx);
        let (target, _) = recv
            .iter()
            .copied()
            .max_by(|a, b| a.1.total_cmp(&b.1))
            .map(|(t, f)| (t as usize, f))?;
        let (r, c) = (idx / self.cols, idx % self.cols);
        let (tr, tc) = (target / self.cols, target % self.cols);
        let dist = if tr != r && tc != c {
            std::f64::consts::SQRT_2 * self.cellsize
        } else {
            self.cellsize
        };
        Some((target, dist))
    }

    fn accumulate_upslope(&self, local: &[f64]) -> Vec<f64> {
        let n = self.rows * self.cols;
        let mut sum = local.to_vec();
        let mut in_degree = vec![0u32; n];
        for &(recv, _) in &self.entries {
            in_degree[recv as usize] += 1;
        }
        let mut queue: Vec<usize> = (0..n)
            .filter(|&i| !self.solid[i] && in_degree[i] == 0)
            .collect();
        while let Some(idx) = queue.pop() {
            let out = sum[idx];
            for &(recv, frac) in self.receivers(idx) {
                let recv = recv as usize;
                sum[recv] += out * f64::from(frac);
                in_degree[recv] -= 1;
                if in_degree[recv] == 0 {
                    queue.push(recv);
                }
            }
        }
        sum
    }
}
