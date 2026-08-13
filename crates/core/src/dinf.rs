//! D-infinity (D∞) flow routing.
//!
//! Computes continuous flow direction angles based on triangular facets
//! fitted to the 3×3 neighborhood (Tarboton 1997). Flow direction is the
//! steepest downslope angle, which can point in any direction (0–2π).
//!
//! For flow accumulation, the flow from each cell is partitioned between
//! the two D8 neighbors that bracket the D∞ angle, with proportions based
//! on angular proximity.
//!
//! ## D∞ convention
//!
//! Flow direction angles are in radians, measured counter-clockwise from
//! East (matching the `SurtGIS` D8 convention): 0 = E, π/2 = N, π = W,
//! 3π/2 = S. A value of -1.0 indicates a pit (no downslope neighbour).
//!
//! ## Integration with [`Network`](crate::Network)
//!
//! The D8 [`Network`] provides pit-filling and the D8 constants. `DinfNetwork`
//! reuses the priority-flood fill from `Network::priority_fill` and operates
//! on the same filled DEM, but computes continuous angles and fractional
//! flow accumulation instead of single-direction routing.
//!
//! Reference:
//! Tarboton, D.G. (1997). A new method for the determination of flow
//! directions and upslope areas in grid digital elevation models.
//! *Water Resources Research*, 33(2), 309–319.

use surtgis_core::raster::Raster;
use surtgis_core::{GeoTransform, Result as GstResult};

use crate::network::{D8_DISTANCE, D8_OFFSETS};
use crate::{Network, SedlinkError};

/// Tarboton (1997) facet decomposition.
///
/// Each facet is defined by a CARDINAL neighbour (e1 direction) and an
/// adjacent DIAGONAL neighbour (e2 direction). The `base_angle` is the
/// azimuth of the cardinal direction; `sign` controls whether θ is
/// added (+1) or subtracted (-1) to reach the diagonal.
///
/// Facet tuple: `(cardinal_idx, diagonal_idx, base_angle, sign)`.
const TARBOTON_FACETS: [(usize, usize, f64, f64); 8] = [
    (0, 1, 0.0, 1.0),
    (2, 1, std::f64::consts::FRAC_PI_2, -1.0),
    (2, 3, std::f64::consts::FRAC_PI_2, 1.0),
    (4, 3, std::f64::consts::PI, -1.0),
    (4, 5, std::f64::consts::PI, 1.0),
    (6, 5, 3.0 * std::f64::consts::FRAC_PI_2, -1.0),
    (6, 7, 3.0 * std::f64::consts::FRAC_PI_2, 1.0),
    (0, 7, 2.0 * std::f64::consts::PI, -1.0),
];

/// Sentinel value for pit cells (no downslope neighbour).
pub const DINF_PIT: f64 = -1.0;

/// Validated DEM parts: (rows, cols, transform, cellsize, elevations,
/// solid mask).
pub(crate) type DemParts = (usize, usize, GeoTransform, f64, Vec<f32>, Vec<bool>);

/// The two D8 receivers of a D∞ angle with their flow fractions.
///
/// Flow is split between the two D8 neighbours bracketing the angle,
/// proportional to the angular offset within the π/4 sector
/// (Tarboton 1997, eq. 6). On exact sector boundaries all flow goes
/// to a single neighbour (the second fraction is 0).
fn dinf_receivers(angle: f64) -> [(usize, f64); 2] {
    let pi4 = std::f64::consts::FRAC_PI_4;
    let sector = ((angle / pi4) % 8.0).floor() as usize;
    let lower_idx = sector.min(7);
    let upper_idx = (lower_idx + 1) % 8;
    let alpha = (angle - (lower_idx as f64) * pi4).clamp(0.0, pi4);
    let frac_upper = alpha / pi4;
    [(lower_idx, 1.0 - frac_upper), (upper_idx, frac_upper)]
}

/// A D-infinity flow network derived from a DEM.
///
/// Holds the D∞ flow direction angles, fractional flow accumulation,
/// and the two downstream receivers (with fractions) for each cell.
/// The network is acyclic by construction (pits are resolved at build
/// time via priority-flood, same as [`Network`]).
pub struct DinfNetwork {
    /// Number of rows.
    rows: usize,
    /// Number of columns.
    cols: usize,
    /// Uniform cell size in metres.
    cellsize: f64,
    /// Geotransform of the domain.
    transform: GeoTransform,
    /// D∞ flow direction angles in radians (0 = East, CCW).
    /// `DINF_PIT` (-1.0) for pit cells. NaN for `NoData`.
    angles: Vec<f64>,
    /// Flow accumulation (fractional cell count, including self) per cell.
    /// NaN for `NoData` cells.
    flow_acc: Vec<f64>,
    /// Downstream receiver indices and fractions.
    /// `downstream[i]` = `[(receiver0, frac0), (receiver1, frac1)]`.
    /// For pit cells, both receivers are `usize::MAX` with frac 0.
    downstream: Vec<[(usize, f64); 2]>,
    /// Slope (radians) per cell, computed from the DEM.
    slope: Vec<f32>,
    /// Solid mask: `true` cells are `NoData` walls.
    solid: Vec<bool>,
}

impl DinfNetwork {
    /// Build a D∞ flow network from a DEM raster.
    ///
    /// Performs pit-filling via priority-flood (reusing [`Network::priority_fill`]),
    /// then computes D∞ flow directions (Tarboton 1997) and fractional
    /// flow accumulation.
    ///
    /// # Errors
    ///
    /// Returns [`SedlinkError::EmptyGrid`], [`SedlinkError::RotatedGrid`],
    /// or [`SedlinkError::NonSquareCells`] if the DEM grid is invalid.
    pub fn from_dem(dem: &Raster<f32>) -> Result<Self, SedlinkError> {
        let (rows, cols, t, pw, z, solid) = Self::validate_dem(dem)?;
        let filled = Network::priority_fill(&z, &solid, rows, cols);
        Ok(Self::build(rows, cols, t, pw, &filled, &z, &solid))
    }

    /// Build a D∞ flow network from a pre-filled DEM (no priority-flood).
    ///
    /// The caller is responsible for hydrologically conditioning the DEM.
    /// This constructor is useful for testing against reference
    /// implementations that operate on already-conditioned DEMs.
    ///
    /// # Errors
    ///
    /// Returns [`SedlinkError::EmptyGrid`], [`SedlinkError::RotatedGrid`],
    /// or [`SedlinkError::NonSquareCells`] if the DEM grid is invalid.
    pub fn from_filled_dem(dem: &Raster<f32>) -> Result<Self, SedlinkError> {
        let (rows, cols, t, pw, z, solid) = Self::validate_dem(dem)?;
        Ok(Self::build(rows, cols, t, pw, &z, &z, &solid))
    }

    /// Validate DEM and extract metadata, elevations, and solid mask.
    pub(crate) fn validate_dem(dem: &Raster<f32>) -> Result<DemParts, SedlinkError> {
        let (rows, cols) = dem.shape();
        if rows == 0 || cols == 0 {
            return Err(SedlinkError::EmptyGrid);
        }

        let t = *dem.transform();
        if t.row_rotation != 0.0 || t.col_rotation != 0.0 {
            return Err(SedlinkError::RotatedGrid {
                row_rotation: t.row_rotation,
                col_rotation: t.col_rotation,
            });
        }

        let pw = t.pixel_width;
        let ph = t.pixel_height.abs();
        if !pw.is_finite() || !ph.is_finite() || pw <= 0.0 || ph <= 0.0 {
            return Err(SedlinkError::NonSquareCells {
                pixel_width: t.pixel_width,
                pixel_height: t.pixel_height,
            });
        }
        let rel_tol = 1e-9;
        if (pw - ph).abs() > rel_tol * pw.max(ph) {
            return Err(SedlinkError::NonSquareCells {
                pixel_width: t.pixel_width,
                pixel_height: t.pixel_height,
            });
        }

        let n = rows * cols;
        let mut z = vec![0.0f32; n];
        let mut solid = vec![false; n];
        for (i, &v) in dem.data().iter().enumerate() {
            if v.is_finite() && !dem.is_nodata(v) {
                z[i] = v;
            } else {
                solid[i] = true;
            }
        }

        Ok((rows, cols, t, pw, z, solid))
    }

    /// Build the network from a filled DEM and original elevations.
    fn build(
        rows: usize,
        cols: usize,
        t: GeoTransform,
        pw: f64,
        filled: &[f32],
        z_orig: &[f32],
        solid: &[bool],
    ) -> Self {
        let flats = crate::flats::FlatOffsets::compute(filled, solid, rows, cols);
        let angles = Self::compute_flow_dir(filled, solid, rows, cols, flats.as_ref());
        let flow_acc = Self::compute_flow_acc(&angles, solid, rows, cols);
        let downstream = Self::compute_downstream(&angles, solid, rows, cols);
        let slope = Network::compute_slope(z_orig, solid, rows, cols, pw);

        Self {
            rows,
            cols,
            cellsize: pw,
            transform: t,
            angles,
            flow_acc,
            downstream,
            slope,
            solid: solid.to_vec(),
        }
    }

    /// Compute D∞ flow direction angles from a (filled) DEM.
    ///
    /// Uses the Tarboton (1997) triangular facet decomposition: for each
    /// cell, fits 8 triangular facets to the 3×3 neighbourhood and selects
    /// the steepest downslope facet. The flow angle is continuous.
    fn compute_flow_dir(
        z: &[f32],
        solid: &[bool],
        rows: usize,
        cols: usize,
        flats: Option<&crate::flats::FlatOffsets>,
    ) -> Vec<f64> {
        let n = rows * cols;
        let cs = 1.0_f64; // We work in cell units; scale doesn't affect direction
        let diag = cs * std::f64::consts::SQRT_2;
        let mut angles = vec![f64::NAN; n];

        crate::par::for_each_row(&mut angles, cols, |r, out_row| {
            for (c, out) in out_row.iter_mut().enumerate() {
                let idx = r * cols + c;
                if solid[idx] {
                    continue;
                }

                let z0 = f64::from(z[idx]);

                // Gather 8 neighbours.
                let mut zn = [f64::NAN; 8];
                let mut all_valid = true;

                for (idx_n, &(dr, dc)) in D8_OFFSETS.iter().enumerate() {
                    let nr = r as isize + dr;
                    let nc = c as isize + dc;

                    if nr < 0 || nc < 0 || nr as usize >= rows || nc as usize >= cols {
                        all_valid = false;
                        continue;
                    }

                    let nidx = nr as usize * cols + nc as usize;
                    if solid[nidx] {
                        all_valid = false;
                        continue;
                    }
                    zn[idx_n] = f64::from(z[nidx]);
                }

                // Find steepest facet using Tarboton (1997) decomposition.
                let mut best_slope = 0.0_f64;
                let mut best_angle = DINF_PIT;

                for &(card_idx, diag_idx, base_angle, sign) in &TARBOTON_FACETS {
                    if zn[card_idx].is_nan() || zn[diag_idx].is_nan() {
                        continue;
                    }

                    // e1: slope from center toward cardinal neighbour (d1 = cs)
                    // e2: cross-slope from cardinal toward diagonal (d2 = cs)
                    let e1 = (z0 - zn[card_idx]) / cs;
                    let e2 = (zn[card_idx] - zn[diag_idx]) / cs;

                    if e1 == 0.0 && e2 == 0.0 {
                        continue;
                    }

                    let raw = e2.atan2(e1);
                    let ad = std::f64::consts::FRAC_PI_4;

                    let (theta, slope): (f64, f64);

                    if raw < 0.0 {
                        theta = 0.0;
                        slope = e1;
                    } else if raw > ad {
                        theta = ad;
                        slope = (z0 - zn[diag_idx]) / diag;
                    } else {
                        theta = raw;
                        slope = (e1 * e1 + e2 * e2).sqrt();
                    }

                    if slope > best_slope {
                        best_slope = slope;
                        let mut angle = base_angle + sign * theta;
                        let two_pi = 2.0 * std::f64::consts::PI;
                        if angle >= two_pi {
                            angle -= two_pi;
                        }
                        if angle < 0.0 {
                            angle += two_pi;
                        }
                        best_angle = angle;
                    }
                }

                // Fallback for edge/nodata cells: D8-style single-direction.
                if best_angle < 0.0 && !all_valid {
                    let cardinal_angles = [
                        0.0,
                        std::f64::consts::FRAC_PI_4,
                        std::f64::consts::FRAC_PI_2,
                        3.0 * std::f64::consts::FRAC_PI_4,
                        std::f64::consts::PI,
                        5.0 * std::f64::consts::FRAC_PI_4,
                        3.0 * std::f64::consts::FRAC_PI_2,
                        7.0 * std::f64::consts::FRAC_PI_4,
                    ];
                    for (idx_n, _) in D8_OFFSETS.iter().enumerate() {
                        if zn[idx_n].is_nan() {
                            continue;
                        }
                        let dist = D8_DISTANCE[idx_n] * cs;
                        let slope = (z0 - zn[idx_n]) / dist;
                        if slope > best_slope {
                            best_slope = slope;
                            best_angle = cardinal_angles[idx_n];
                        }
                    }
                }

                // Flat cells: use the resolved D8 flat direction's angle.
                if best_angle < 0.0
                    && let Some(f) = flats
                    && let Some(d) = f.direction(idx, z, solid)
                {
                    best_angle = f64::from(d - 1) * std::f64::consts::FRAC_PI_4;
                }

                *out = best_angle;
            }
        });

        angles
    }

    /// Compute D∞ flow accumulation from D∞ angles.
    ///
    /// Each cell counts itself plus the (fractional) count of all
    /// upstream cells that flow into it — headwater cells have
    /// accumulation = 1, matching the D8 [`Network`] convention.
    /// Processing order is a dependency-driven topological sort.
    fn compute_flow_acc(angles: &[f64], solid: &[bool], rows: usize, cols: usize) -> Vec<f64> {
        let n = rows * cols;
        let mut acc = vec![1.0_f64; n];

        // In-degree = number of incoming flow edges per cell.
        let mut in_degree = vec![0_u32; n];

        for r in 0..rows {
            for c in 0..cols {
                let idx = r * cols + c;
                if solid[idx] {
                    continue;
                }
                let a = angles[idx];
                if a.is_nan() || a < 0.0 {
                    continue;
                }

                let recv = dinf_receivers(a);
                for (nbr_idx, frac) in recv {
                    if frac <= 0.0 {
                        continue;
                    }
                    let (dr, dc) = D8_OFFSETS[nbr_idx];
                    let nr = r as isize + dr;
                    let nc = c as isize + dc;
                    if nr >= 0 && nc >= 0 && (nr as usize) < rows && (nc as usize) < cols {
                        let nidx = nr as usize * cols + nc as usize;
                        if !solid[nidx] {
                            in_degree[nidx] += 1;
                        }
                    }
                }
            }
        }

        // Topological propagation from headwaters.
        let mut queue: Vec<usize> = (0..n).filter(|&i| !solid[i] && in_degree[i] == 0).collect();

        let mut processed = 0_usize;
        while let Some(idx) = queue.pop() {
            processed += 1;
            let a = angles[idx];
            if a.is_nan() || a < 0.0 {
                continue;
            }

            let outflow = acc[idx];
            let r = idx / cols;
            let c = idx % cols;

            let recv = dinf_receivers(a);
            for (nbr_idx, frac) in recv {
                if frac <= 0.0 {
                    continue;
                }
                let (dr, dc) = D8_OFFSETS[nbr_idx];
                let nr = r as isize + dr;
                let nc = c as isize + dc;
                if nr < 0 || nc < 0 || nr as usize >= rows || nc as usize >= cols {
                    continue;
                }
                let nidx = nr as usize * cols + nc as usize;
                if solid[nidx] {
                    continue;
                }
                acc[nidx] += outflow * frac;
                in_degree[nidx] -= 1;
                if in_degree[nidx] == 0 {
                    queue.push(nidx);
                }
            }
        }

        // Forward remaining cells once (cycle guard).
        if processed < n {
            Self::forward_remaining(angles, solid, rows, cols, &in_degree, &mut acc);
        }

        // NaN out nodata cells.
        for idx in 0..n {
            if solid[idx] {
                acc[idx] = f64::NAN;
            }
        }

        acc
    }

    /// Cycle guard for [`Self::compute_flow_acc`]: forward the outflow of
    /// cells the topological pass never reached (should not happen after
    /// pit resolution) one single time.
    fn forward_remaining(
        angles: &[f64],
        solid: &[bool],
        rows: usize,
        cols: usize,
        in_degree: &[u32],
        acc: &mut [f64],
    ) {
        for r in 0..rows {
            for c in 0..cols {
                let idx = r * cols + c;
                if solid[idx] || in_degree[idx] == 0 {
                    continue;
                }
                let a = angles[idx];
                if a.is_nan() || a < 0.0 {
                    continue;
                }
                let outflow = acc[idx];
                let recv = dinf_receivers(a);
                for (nbr_idx, frac) in recv {
                    if frac <= 0.0 {
                        continue;
                    }
                    let (dr, dc) = D8_OFFSETS[nbr_idx];
                    let nr = r as isize + dr;
                    let nc = c as isize + dc;
                    if nr >= 0 && nc >= 0 && (nr as usize) < rows && (nc as usize) < cols {
                        let nidx = nr as usize * cols + nc as usize;
                        if !solid[nidx] {
                            acc[nidx] += outflow * frac;
                        }
                    }
                }
            }
        }
    }

    /// Compute downstream receivers (two D8 neighbours + fractions) for each cell.
    fn compute_downstream(
        angles: &[f64],
        solid: &[bool],
        rows: usize,
        cols: usize,
    ) -> Vec<[(usize, f64); 2]> {
        let n = rows * cols;
        let mut downstream = vec![[(usize::MAX, 0.0_f64), (usize::MAX, 0.0)]; n];

        crate::par::for_each_row(&mut downstream, cols, |r, out_row| {
            for (c, out) in out_row.iter_mut().enumerate() {
                let idx = r * cols + c;
                if solid[idx] {
                    continue;
                }
                let a = angles[idx];
                if a.is_nan() || a < 0.0 {
                    continue;
                }

                let recv = dinf_receivers(a);
                for (slot, (nbr_idx, frac)) in recv.iter().enumerate() {
                    if *frac <= 0.0 {
                        continue;
                    }
                    let (dr, dc) = D8_OFFSETS[*nbr_idx];
                    let nr = r as isize + dr;
                    let nc = c as isize + dc;
                    if nr >= 0 && nc >= 0 && (nr as usize) < rows && (nc as usize) < cols {
                        let nidx = nr as usize * cols + nc as usize;
                        if !solid[nidx] {
                            out[slot] = (nidx, *frac);
                        }
                    }
                }
            }
        });

        downstream
    }

    /// Primary receiver of a cell: the downstream neighbour carrying the
    /// largest flow fraction, or `None` for pits and `NoData` cells.
    fn primary_receiver(&self, idx: usize) -> Option<usize> {
        let [a, b] = self.downstream[idx];
        let (target, frac) = if a.1 >= b.1 { a } else { b };
        if target == usize::MAX || frac <= 0.0 {
            None
        } else {
            Some(target)
        }
    }

    /// Trace the downstream flow path from a cell, following the primary
    /// receiver (highest fraction). Returns cell indices along the path
    /// (excluding the starting cell, including the terminal cell).
    #[must_use]
    pub fn trace_downstream(&self, start: usize) -> Vec<usize> {
        let mut path = Vec::new();
        let mut cur = start;
        let max_steps = self.rows * self.cols;
        let mut steps = 0;

        while let Some(primary) = self.primary_receiver(cur) {
            path.push(primary);
            cur = primary;
            steps += 1;
            if steps > max_steps {
                break;
            }
        }

        path
    }

    /// Number of rows.
    #[must_use]
    pub fn rows(&self) -> usize {
        self.rows
    }

    /// Number of columns.
    #[must_use]
    pub fn cols(&self) -> usize {
        self.cols
    }

    /// Total number of cells.
    #[must_use]
    pub fn len(&self) -> usize {
        self.rows * self.cols
    }

    /// `true` if the grid holds no cells.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Uniform cell size in metres.
    #[must_use]
    pub fn cellsize(&self) -> f64 {
        self.cellsize
    }

    /// Geotransform of the domain.
    #[must_use]
    pub fn transform(&self) -> &GeoTransform {
        &self.transform
    }

    /// D∞ flow direction angle at (row, col) in radians.
    /// Returns `DINF_PIT` (-1.0) for pit cells, NaN for `NoData`.
    #[must_use]
    pub fn angle(&self, row: usize, col: usize) -> f64 {
        self.angles[row * self.cols + col]
    }

    /// Flow accumulation (cell count, fractional) at (row, col).
    #[must_use]
    pub fn flow_acc(&self, row: usize, col: usize) -> f64 {
        self.flow_acc[row * self.cols + col]
    }

    /// Slope (radians) at (row, col).
    #[must_use]
    pub fn slope(&self, row: usize, col: usize) -> f32 {
        self.slope[row * self.cols + col]
    }

    /// `true` if the cell is `NoData` (solid wall).
    #[must_use]
    pub fn is_solid(&self, row: usize, col: usize) -> bool {
        self.solid[row * self.cols + col]
    }

    /// Downstream receivers for a cell: `[(receiver_idx, fraction), ...]`.
    /// For pit cells, both entries have `usize::MAX` as the receiver.
    #[must_use]
    pub fn downstream(&self, idx: usize) -> &[(usize, f64); 2] {
        &self.downstream[idx]
    }

    /// Borrow the flow direction angle grid as a flat slice.
    #[must_use]
    pub fn angles_slice(&self) -> &[f64] {
        &self.angles
    }

    /// Borrow the flow accumulation grid as a flat slice.
    #[must_use]
    pub fn flow_acc_slice(&self) -> &[f64] {
        &self.flow_acc
    }

    /// Borrow the slope grid as a flat slice.
    #[must_use]
    pub fn slope_slice(&self) -> &[f32] {
        &self.slope
    }

    /// Build a `Raster<f64>` of D∞ flow direction angles with the same geotransform.
    ///
    /// # Errors
    ///
    /// Propagates the underlying raster construction error (cannot occur
    /// for a well-formed network, whose grids always match its shape).
    pub fn angles_raster(&self) -> GstResult<Raster<f64>> {
        let mut r = Raster::from_vec(self.angles.clone(), self.rows, self.cols)?;
        r.set_transform(self.transform);
        r.set_nodata(Some(f64::NAN));
        Ok(r)
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

    /// Build a `Raster<f32>` of slope (radians) with the same geotransform.
    ///
    /// # Errors
    ///
    /// Propagates the underlying raster construction error (cannot occur
    /// for a well-formed network, whose grids always match its shape).
    pub fn slope_raster(&self) -> GstResult<Raster<f32>> {
        let mut r = Raster::from_vec(self.slope.clone(), self.rows, self.cols)?;
        r.set_transform(self.transform);
        Ok(r)
    }
}

impl crate::FlowNetwork for DinfNetwork {
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
        let target = self.primary_receiver(idx)?;
        let (r, c) = (idx / self.cols, idx % self.cols);
        let (tr, tc) = (target / self.cols, target % self.cols);
        let dist = if tr != r && tc != c {
            std::f64::consts::SQRT_2 * self.cellsize
        } else {
            self.cellsize
        };
        Some((target, dist))
    }

    /// Fractional accumulation: each upstream cell's contribution is
    /// weighted by the flow fraction it sends through each receiver.
    fn accumulate_upslope(&self, local: &[f64]) -> Vec<f64> {
        let n = self.rows * self.cols;
        let mut sum = local.to_vec();

        let mut in_degree = vec![0_u32; n];
        for idx in 0..n {
            if self.solid[idx] {
                continue;
            }
            for &(ds, frac) in &self.downstream[idx] {
                if ds != usize::MAX && frac > 0.0 {
                    in_degree[ds] += 1;
                }
            }
        }

        let mut queue: Vec<usize> = (0..n)
            .filter(|&i| !self.solid[i] && in_degree[i] == 0)
            .collect();

        while let Some(idx) = queue.pop() {
            for &(ds, frac) in &self.downstream[idx] {
                if ds == usize::MAX || frac <= 0.0 {
                    continue;
                }
                sum[ds] += sum[idx] * frac;
                in_degree[ds] -= 1;
                if in_degree[ds] == 0 {
                    queue.push(ds);
                }
            }
        }

        sum
    }
}
