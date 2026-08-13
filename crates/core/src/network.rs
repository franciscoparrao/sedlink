//! D8 flow network construction and analysis.
//!
//! Builds a flow graph from a DEM using D8 (single-flow-direction) routing,
//! resolves pits via priority-flood, and provides downstream tracing with
//! memoization. Also computes Strahler stream order, stream magnitude, and
//! longitudinal profiles.
//!
//! ## D8 convention
//!
//! Flow directions follow the `SurtGIS` convention (counter-clockwise from
//! East):
//!
//! ```text
//!   4  3  2
//!   5  0  1
//!   6  7  8
//! ```
//!
//! Direction 0 = pit (no downstream). All other values index into
//! [`D8_OFFSETS`] and [`D8_DISTANCE`].

use surtgis_core::raster::Raster;
use surtgis_core::{GeoTransform, Result as GstResult};

use crate::SedlinkError;

/// D8 flow direction offsets (row, col) for directions 1–8.
/// Counter-clockwise from East (matching `SurtGIS` convention):
/// 1=E, 2=NE, 3=N, 4=NW, 5=W, 6=SW, 7=S, 8=SE.
pub const D8_OFFSETS: [(isize, isize); 8] = [
    (0, 1),   // 1: E
    (-1, 1),  // 2: NE
    (-1, 0),  // 3: N
    (-1, -1), // 4: NW
    (0, -1),  // 5: W
    (1, -1),  // 6: SW
    (1, 0),   // 7: S
    (1, 1),   // 8: SE
];

/// D8 flow path distances (in units of cell size) for directions 1–8.
/// Cardinal directions = 1.0, diagonal directions = √2.
pub const D8_DISTANCE: [f64; 8] = [
    1.0,
    std::f64::consts::SQRT_2,
    1.0,
    std::f64::consts::SQRT_2,
    1.0,
    std::f64::consts::SQRT_2,
    1.0,
    std::f64::consts::SQRT_2,
];

/// A D8 flow network derived from a DEM.
///
/// Holds the flow direction grid, flow accumulation, and the resolved
/// downstream cell index for each cell. The network is acyclic by
/// construction (pits are resolved at build time).
pub struct Network {
    /// Number of rows.
    rows: usize,
    /// Number of columns.
    cols: usize,
    /// Uniform cell size in metres.
    cellsize: f64,
    /// Geotransform of the domain.
    transform: GeoTransform,
    /// Flow direction per cell (0 = pit, 1–8 = D8 direction).
    pub(crate) flow_dir: Vec<u8>,
    /// Flow accumulation per cell (cell count, including self).
    pub(crate) flow_acc: Vec<f64>,
    /// Downstream cell index (row-major), or `usize::MAX` for pits.
    pub(crate) downstream: Vec<usize>,
    /// Precomputed list of upstream cell indices for each cell.
    pub(crate) upstream: Vec<Vec<usize>>,
    /// Slope (radians) per cell, computed from the DEM.
    pub(crate) slope: Vec<f32>,
    /// Solid mask: `true` cells are `NoData` walls.
    pub(crate) solid: Vec<bool>,
}

impl Network {
    /// Build a flow network from a DEM raster.
    ///
    /// Performs pit-filling via priority-flood, then computes D8 flow
    /// directions and flow accumulation.
    ///
    /// # Errors
    ///
    /// Returns [`SedlinkError::EmptyGrid`], [`SedlinkError::RotatedGrid`],
    /// or [`SedlinkError::NonSquareCells`] if the DEM grid is invalid.
    pub fn from_dem(dem: &Raster<f32>) -> Result<Self, SedlinkError> {
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

        // Pit-filling via priority-flood, then flat resolution.
        let filled = Self::priority_fill(&z, &solid, rows, cols);
        let flats = crate::flats::FlatOffsets::compute(&filled, &solid, rows, cols);

        // D8 flow direction from filled DEM (+ flat drainage).
        let flow_dir = Self::compute_flow_dir(&filled, &solid, rows, cols, flats.as_ref());

        // Flow accumulation.
        let flow_acc = Self::compute_flow_acc(&flow_dir, &solid, rows, cols);

        // Downstream index.
        let downstream = Self::compute_downstream(&flow_dir, rows, cols);

        // Upstream lists.
        let upstream = Self::compute_upstream(&downstream, n);

        // Slope from original (unfilled) DEM.
        let slope = Self::compute_slope(&z, &solid, rows, cols, pw);

        Ok(Self {
            rows,
            cols,
            cellsize: pw,
            transform: t,
            flow_dir,
            flow_acc,
            downstream,
            upstream,
            slope,
            solid,
        })
    }

    /// Priority-flood pit filling (Barnes et al. 2014).
    ///
    /// Uses a min-heap over true `f32` elevations seeded with all boundary
    /// cells. Cells below the spill elevation of their lowest processed
    /// neighbour are raised to exactly that spill elevation; the
    /// resulting flats are drained afterwards by
    /// [`FlatOffsets`](crate::flats::FlatOffsets) (two-gradient flat
    /// resolution, Garbrecht & Martz 1997).
    pub(crate) fn priority_fill(z: &[f32], solid: &[bool], rows: usize, cols: usize) -> Vec<f32> {
        use ordered_float::OrderedFloat;
        use std::cmp::Reverse;
        use std::collections::BinaryHeap;

        let n = rows * cols;
        let mut filled = z.to_vec();
        let mut visited = vec![false; n];

        // Min-heap over (elevation, index); OrderedFloat gives a total order.
        let mut heap: BinaryHeap<Reverse<(OrderedFloat<f32>, usize)>> = BinaryHeap::new();

        // Seed boundary cells.
        for c in 0..cols {
            for r_seed in [0usize, rows - 1] {
                let idx = r_seed * cols + c;
                if !visited[idx] && !solid[idx] {
                    visited[idx] = true;
                    heap.push(Reverse((OrderedFloat(filled[idx]), idx)));
                }
            }
        }
        for r in 1..rows.saturating_sub(1) {
            for c_seed in [0usize, cols - 1] {
                let idx = r * cols + c_seed;
                if !visited[idx] && !solid[idx] {
                    visited[idx] = true;
                    heap.push(Reverse((OrderedFloat(filled[idx]), idx)));
                }
            }
        }

        // 8-neighbour connectivity for flood.
        let nbrs: [(isize, isize); 8] = [
            (-1, 0),
            (1, 0),
            (0, -1),
            (0, 1),
            (-1, -1),
            (-1, 1),
            (1, -1),
            (1, 1),
        ];

        while let Some(Reverse((spill, idx))) = heap.pop() {
            let spill = spill.0;
            let r = idx / cols;
            let c = idx % cols;

            for &(dr, dc) in &nbrs {
                let nr = r as isize + dr;
                let nc = c as isize + dc;
                if nr < 0 || nc < 0 || nr as usize >= rows || nc as usize >= cols {
                    continue;
                }
                let nidx = nr as usize * cols + nc as usize;
                if visited[nidx] || solid[nidx] {
                    continue;
                }
                visited[nidx] = true;
                // Don't lower; only raise to the spill elevation.
                if filled[nidx] < spill {
                    filled[nidx] = spill;
                }
                heap.push(Reverse((OrderedFloat(filled[nidx]), nidx)));
            }
        }

        filled
    }

    /// Compute D8 flow direction from a (filled) DEM.
    ///
    /// Steepest descent: picks the neighbour with the largest drop per
    /// unit distance (diagonal steps are √2 longer), matching the
    /// standard D8 definition (O'Callaghan & Mark 1984; `TauDEM`).
    /// Cardinal directions are checked first so that gradient ties break
    /// towards cardinal flow. Cells with no lower neighbour consult the
    /// flat-resolution offsets.
    fn compute_flow_dir(
        z: &[f32],
        solid: &[bool],
        rows: usize,
        cols: usize,
        flats: Option<&crate::flats::FlatOffsets>,
    ) -> Vec<u8> {
        // Check order: cardinals first (1=E, 3=N, 5=W, 7=S), then diagonals
        // (2=NE, 4=NW, 6=SW, 8=SE), so gradient ties break towards cardinal
        // flow.
        const CHECK_ORDER: [(usize, (isize, isize)); 8] = [
            (0, (0, 1)),   // 1: E
            (2, (-1, 0)),  // 3: N
            (4, (0, -1)),  // 5: W
            (6, (1, 0)),   // 7: S
            (1, (-1, 1)),  // 2: NE
            (3, (-1, -1)), // 4: NW
            (5, (1, -1)),  // 6: SW
            (7, (1, 1)),   // 8: SE
        ];

        let n = rows * cols;
        let mut flow_dir = vec![0u8; n];

        crate::par::for_each_row(&mut flow_dir, cols, |r, out_row| {
            for (c, out) in out_row.iter_mut().enumerate() {
                let idx = r * cols + c;
                if solid[idx] {
                    continue;
                }

                let here = f64::from(z[idx]);
                let mut best_dir = 0u8;
                let mut best_grad = 0.0f64;

                for &(d8_index, (dr, dc)) in &CHECK_ORDER {
                    let nr = r as isize + dr;
                    let nc = c as isize + dc;
                    if nr < 0 || nc < 0 || nr as usize >= rows || nc as usize >= cols {
                        continue;
                    }
                    let nidx = nr as usize * cols + nc as usize;
                    if solid[nidx] {
                        continue;
                    }
                    let grad = (here - f64::from(z[nidx])) / D8_DISTANCE[d8_index];
                    if grad > best_grad {
                        best_grad = grad;
                        best_dir = (d8_index + 1) as u8;
                    }
                }

                if best_dir == 0
                    && let Some(f) = flats
                    && let Some(d) = f.direction(idx, z, solid)
                {
                    best_dir = d;
                }

                *out = best_dir;
            }
        });

        flow_dir
    }

    /// Compute flow accumulation (cell count) from flow directions.
    ///
    /// Processes cells in reverse topological order (downstream first)
    /// so each cell's accumulation is final before its upstream cells
    /// are processed.
    fn compute_flow_acc(flow_dir: &[u8], solid: &[bool], rows: usize, cols: usize) -> Vec<f64> {
        let n = rows * cols;
        let mut acc = vec![1.0f64; n];

        // Topological order: process from outlets upstream.
        // We use the iterative approach: process cells whose downstream
        // is already resolved (or is a pit).
        let downstream = Self::compute_downstream(flow_dir, rows, cols);

        // Count upstream dependencies for each cell.
        let mut pending = vec![0usize; n];
        for idx in 0..n {
            if solid[idx] {
                continue;
            }
            let ds = downstream[idx];
            if ds != usize::MAX {
                pending[ds] += 1;
            }
        }

        // Queue: headwater cells (no upstream contributors).
        let mut queue: Vec<usize> = (0..n).filter(|&i| !solid[i] && pending[i] == 0).collect();

        while let Some(idx) = queue.pop() {
            let ds = downstream[idx];
            if ds != usize::MAX {
                acc[ds] += acc[idx];
                pending[ds] -= 1;
                if pending[ds] == 0 {
                    queue.push(ds);
                }
            }
        }

        acc
    }

    /// Accumulate a per-cell quantity over each cell's upslope contributing
    /// area (including the cell itself), following the D8 flow network.
    ///
    /// `sum[i]` = `local[i]` + Σ `local[j]` over all cells `j` draining
    /// through `i`. Dividing by [`flow_acc`](Self::flow_acc) yields the
    /// upslope mean (e.g. W̄ and S̄ in the Borselli/Cavalli IC equations).
    pub(crate) fn accumulate_upslope(&self, local: &[f64]) -> Vec<f64> {
        let n = self.len();
        let mut sum = local.to_vec();

        let mut pending = vec![0usize; n];
        for idx in 0..n {
            if self.solid[idx] {
                continue;
            }
            let ds = self.downstream[idx];
            if ds != usize::MAX {
                pending[ds] += 1;
            }
        }

        let mut queue: Vec<usize> = (0..n)
            .filter(|&i| !self.solid[i] && pending[i] == 0)
            .collect();

        while let Some(idx) = queue.pop() {
            let ds = self.downstream[idx];
            if ds != usize::MAX {
                sum[ds] += sum[idx];
                pending[ds] -= 1;
                if pending[ds] == 0 {
                    queue.push(ds);
                }
            }
        }

        sum
    }

    /// Compute downstream cell index for each cell.
    fn compute_downstream(flow_dir: &[u8], rows: usize, cols: usize) -> Vec<usize> {
        let n = rows * cols;
        let mut downstream = vec![usize::MAX; n];

        for idx in 0..n {
            let dir = flow_dir[idx];
            if dir == 0 || dir > 8 {
                continue;
            }
            let (dr, dc) = D8_OFFSETS[(dir - 1) as usize];
            let r = idx / cols;
            let c = idx % cols;
            let nr = r as isize + dr;
            let nc = c as isize + dc;
            if nr >= 0 && nc >= 0 && (nr as usize) < rows && (nc as usize) < cols {
                downstream[idx] = nr as usize * cols + nc as usize;
            }
        }

        downstream
    }

    /// Compute upstream cell lists from downstream indices.
    fn compute_upstream(downstream: &[usize], n: usize) -> Vec<Vec<usize>> {
        let mut upstream: Vec<Vec<usize>> = vec![Vec::new(); n];
        for (idx, &ds) in downstream.iter().enumerate() {
            if ds != usize::MAX {
                upstream[ds].push(idx);
            }
        }
        upstream
    }

    /// Compute slope (radians) from the original DEM using central
    /// differences, one-sided at edges.
    pub(crate) fn compute_slope(
        z: &[f32],
        solid: &[bool],
        rows: usize,
        cols: usize,
        dx: f64,
    ) -> Vec<f32> {
        let n = rows * cols;
        let mut slope = vec![0.0f32; n];

        let sample = |r: isize, c: isize| -> Option<f64> {
            if r < 0 || c < 0 || r as usize >= rows || c as usize >= cols {
                return None;
            }
            let i = r as usize * cols + c as usize;
            if solid[i] {
                return None;
            }
            Some(f64::from(z[i]))
        };

        crate::par::for_each_row(&mut slope, cols, |r, out_row| {
            for (c, out) in out_row.iter_mut().enumerate() {
                let idx = r * cols + c;
                if solid[idx] {
                    *out = 0.0;
                    continue;
                }
                let here = sample(r as isize, c as isize).unwrap_or(0.0);
                let grad_x = match (
                    sample(r as isize, c as isize - 1),
                    sample(r as isize, c as isize + 1),
                ) {
                    (Some(west), Some(east)) => (east - west) / (2.0 * dx),
                    (Some(west), None) => (here - west) / dx,
                    (None, Some(east)) => (east - here) / dx,
                    (None, None) => 0.0,
                };
                let grad_y = match (
                    sample(r as isize - 1, c as isize),
                    sample(r as isize + 1, c as isize),
                ) {
                    (Some(north), Some(south)) => (south - north) / (2.0 * dx),
                    (Some(north), None) => (here - north) / dx,
                    (None, Some(south)) => (south - here) / dx,
                    (None, None) => 0.0,
                };
                let mag = (grad_x * grad_x + grad_y * grad_y).sqrt();
                *out = (mag.atan() as f32).clamp(0.0, std::f32::consts::FRAC_PI_2);
            }
        });

        slope
    }

    /// Delineate watersheds for a set of pour points.
    ///
    /// Returns a per-cell label: `k + 1` when the cell drains through
    /// `pour_points[k]` (the **nearest downstream** pour point wins for
    /// nested basins), `0` when the cell drains elsewhere or is `NoData`.
    /// A pour point cell belongs to its own basin.
    ///
    /// Pour points are flat (row-major) cell indices; use
    /// [`snap_to_stream`](Self::snap_to_stream) to place them on the
    /// channel first.
    #[must_use]
    pub fn watersheds(&self, pour_points: &[usize]) -> Vec<u32> {
        const UNRESOLVED: u32 = u32::MAX;

        let n = self.len();
        let mut label = vec![UNRESOLVED; n];

        for (k, &p) in pour_points.iter().enumerate() {
            if p < n {
                label[p] = (k + 1) as u32;
            }
        }

        for start in 0..n {
            if label[start] != UNRESOLVED || self.solid[start] {
                if self.solid[start] {
                    label[start] = 0;
                }
                continue;
            }
            // Walk downstream until a resolved cell or a pit.
            let mut path = vec![start];
            let mut cur = start;
            let resolved = loop {
                let ds = self.downstream[cur];
                if ds == usize::MAX {
                    break 0; // pit/outlet without pour point
                }
                if label[ds] != UNRESOLVED {
                    break label[ds];
                }
                path.push(ds);
                cur = ds;
            };
            for idx in path {
                label[idx] = resolved;
            }
        }

        label
    }

    /// Snap a pour point to the highest flow accumulation cell within a
    /// Chebyshev `radius`, returning the flat index of the snapped cell.
    /// Useful to place outlet coordinates (from GPS or maps) onto the
    /// modelled channel.
    #[must_use]
    pub fn snap_to_stream(&self, row: usize, col: usize, radius: usize) -> usize {
        let mut best = row * self.cols + col;
        let mut best_acc = f64::MIN;
        let (r0, c0) = (row as isize, col as isize);
        let rad = radius as isize;
        for dr in -rad..=rad {
            for dc in -rad..=rad {
                let (r, c) = (r0 + dr, c0 + dc);
                if r < 0 || c < 0 || r as usize >= self.rows || c as usize >= self.cols {
                    continue;
                }
                let idx = r as usize * self.cols + c as usize;
                if !self.solid[idx] && self.flow_acc[idx] > best_acc {
                    best_acc = self.flow_acc[idx];
                    best = idx;
                }
            }
        }
        best
    }

    /// Trace the downstream flow path from a cell to the nearest stream
    /// or outlet, returning the list of cell indices along the path
    /// (excluding the starting cell, including the terminal cell).
    #[must_use]
    pub fn trace_downstream(&self, start: usize) -> Vec<usize> {
        let mut path = Vec::new();
        let mut cur = start;
        let max_steps = self.rows * self.cols;
        let mut steps = 0;

        loop {
            let ds = self.downstream[cur];
            if ds == usize::MAX {
                break;
            }
            path.push(ds);
            cur = ds;
            steps += 1;
            if steps > max_steps {
                break;
            }
        }

        path
    }

    /// Trace the downstream path from a cell to the nearest stream cell,
    /// returning the list of cell indices (excluding start, including
    /// the stream cell).
    #[must_use]
    pub fn trace_to_stream(&self, start: usize, threshold: f64) -> Vec<usize> {
        let mut path = Vec::new();
        let mut cur = start;
        let max_steps = self.rows * self.cols;
        let mut steps = 0;

        loop {
            if self.flow_acc[cur] >= threshold && cur != start {
                break;
            }
            let ds = self.downstream[cur];
            if ds == usize::MAX {
                break;
            }
            path.push(ds);
            cur = ds;
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

    /// Flow direction at (row, col). 0 = pit, 1–8 = D8 direction.
    #[must_use]
    pub fn flow_dir(&self, row: usize, col: usize) -> u8 {
        self.flow_dir[row * self.cols + col]
    }

    /// Flow accumulation (cell count) at (row, col).
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

    /// Downstream cell index (row-major), or `usize::MAX` for pits.
    #[must_use]
    pub fn downstream(&self, idx: usize) -> usize {
        self.downstream[idx]
    }

    /// Upstream cell indices for a given cell.
    #[must_use]
    pub fn upstream(&self, idx: usize) -> &[usize] {
        &self.upstream[idx]
    }

    /// Borrow the flow direction grid as a flat slice.
    #[must_use]
    pub fn flow_dir_slice(&self) -> &[u8] {
        &self.flow_dir
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

    /// Build a `Raster<f64>` of flow accumulation with the same geotransform.
    ///
    /// # Errors
    ///
    /// Propagates the underlying raster construction error (cannot occur
    /// for a well-formed network, whose grids always match its shape).
    pub fn flow_acc_raster(&self) -> GstResult<Raster<f64>> {
        let mut r = Raster::from_vec(self.flow_acc.clone(), self.rows, self.cols)?;
        r.set_transform(self.transform);
        Ok(r)
    }

    /// Build a `Raster<u8>` of flow directions with the same geotransform.
    ///
    /// # Errors
    ///
    /// Propagates the underlying raster construction error (cannot occur
    /// for a well-formed network, whose grids always match its shape).
    pub fn flow_dir_raster(&self) -> GstResult<Raster<u8>> {
        let mut r = Raster::from_vec(self.flow_dir.clone(), self.rows, self.cols)?;
        r.set_transform(self.transform);
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

impl crate::FlowNetwork for Network {
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
        let dir = self.flow_dir[idx];
        if dir == 0 || dir > 8 {
            return None;
        }
        let ds = self.downstream[idx];
        if ds == usize::MAX {
            return None;
        }
        Some((ds, D8_DISTANCE[(dir - 1) as usize] * self.cellsize))
    }

    fn accumulate_upslope(&self, local: &[f64]) -> Vec<f64> {
        Network::accumulate_upslope(self, local)
    }
}

/// Strahler stream order for each cell.
#[derive(Debug, Clone)]
pub struct StrahlerOrder {
    /// Order value per cell (0 for non-stream cells).
    pub values: Vec<u8>,
    /// Maximum order in the network.
    pub max_order: u8,
}

/// Stream magnitude (number of upstream sources) for each cell.
#[derive(Debug, Clone)]
pub struct StreamMagnitude {
    /// Magnitude per cell (0 for non-stream cells).
    pub values: Vec<f64>,
    /// Maximum magnitude in the network.
    pub max_magnitude: f64,
}

/// Longitudinal profile: distance upstream (m) and elevation (m) along
/// the flow path from each stream cell to the outlet.
#[derive(Debug, Clone)]
pub struct StreamProfile {
    /// Distance upstream from the outlet, in metres.
    pub distance: Vec<f64>,
    /// Elevation at each profile point, in metres.
    pub elevation: Vec<f64>,
}

/// Network analysis: Strahler ordering, magnitude, and profiles.
pub struct NetworkAnalysis<'a> {
    net: &'a Network,
}

impl<'a> NetworkAnalysis<'a> {
    /// Create a network analysis handle.
    #[must_use]
    pub fn new(net: &'a Network) -> Self {
        Self { net }
    }

    /// Compute Strahler stream order.
    ///
    /// Stream cells are defined by `flow_acc >= threshold`. The Strahler
    /// order is computed by processing cells in topological order
    /// (downstream first):
    /// - Order 1 for stream cells with no upstream stream cells.
    /// - If exactly one upstream stream cell has order *n*, the cell gets
    ///   order *n*.
    /// - If two or more upstream stream cells have the same maximum order
    ///   *n*, the cell gets order *n* + 1.
    #[must_use]
    pub fn strahler_order(&self, threshold: f64) -> StrahlerOrder {
        let n = self.net.len();
        let mut order = vec![0u8; n];

        // Identify stream cells.
        // Solid (NoData) cells keep acc = 1.0, so exclude them explicitly
        // before thresholding.
        let is_stream: Vec<bool> = (0..n)
            .map(|i| !self.net.solid[i] && self.net.flow_acc[i] >= threshold)
            .collect();

        // topological_order() yields downstream cells first; Strahler
        // needs each cell processed AFTER its upstream cells, so walk it
        // in reverse.
        let topo = self.topological_order();

        for &idx in topo.iter().rev() {
            if !is_stream[idx] {
                continue;
            }

            let mut max_up = 0u8;
            let mut count_max = 0;

            for &up in self.net.upstream(idx) {
                if is_stream[up] {
                    let up_order = order[up];
                    if up_order > max_up {
                        max_up = up_order;
                        count_max = 1;
                    } else if up_order == max_up {
                        count_max += 1;
                    }
                }
            }

            order[idx] = if max_up == 0 {
                1
            } else if count_max >= 2 {
                max_up + 1
            } else {
                max_up
            };
        }

        let max_order = order.iter().copied().max().unwrap_or(0);

        StrahlerOrder {
            values: order,
            max_order,
        }
    }

    /// Compute stream magnitude (number of upstream sources).
    ///
    /// For each stream cell, magnitude = 1 + sum of magnitudes of upstream
    /// stream cells.
    pub fn stream_magnitude(&self, threshold: f64) -> StreamMagnitude {
        let n = self.net.len();
        let mut mag = vec![0.0f64; n];

        // Solid (NoData) cells keep acc = 1.0, so exclude them explicitly
        // before thresholding.
        let is_stream: Vec<bool> = (0..n)
            .map(|i| !self.net.solid[i] && self.net.flow_acc[i] >= threshold)
            .collect();

        // Reverse for the same reason as strahler_order: magnitudes of
        // upstream cells must be final before their downstream cell.
        let topo = self.topological_order();

        for &idx in topo.iter().rev() {
            if !is_stream[idx] {
                continue;
            }

            // Shreve magnitude: sources are 1; a junction is the SUM of
            // its tributaries (no per-cell increment).
            let mut sum = 0.0;
            for &up in self.net.upstream(idx) {
                if is_stream[up] {
                    sum += mag[up];
                }
            }
            mag[idx] = if sum == 0.0 { 1.0 } else { sum };
        }

        let max_magnitude = mag.iter().copied().fold(0.0f64, f64::max);

        StreamMagnitude {
            values: mag,
            max_magnitude,
        }
    }

    /// Extract the longitudinal profile for a given stream cell: trace
    /// downstream to the outlet, recording distance and elevation.
    #[must_use]
    pub fn longitudinal_profile(&self, start: usize, dem: &[f32]) -> StreamProfile {
        let mut distance = Vec::new();
        let mut elevation = Vec::new();

        let mut cur = start;
        let mut dist = 0.0f64;
        let cellsize = self.net.cellsize();

        // Include the starting cell.
        distance.push(dist);
        elevation.push(f64::from(dem[cur]));

        loop {
            let ds = self.net.downstream[cur];
            if ds == usize::MAX {
                break;
            }
            let dir = self.net.flow_dir[cur];
            if dir == 0 || dir > 8 {
                break;
            }
            let step = D8_DISTANCE[(dir - 1) as usize] * cellsize;
            dist += step;
            distance.push(dist);
            elevation.push(f64::from(dem[ds]));
            cur = ds;
        }

        StreamProfile {
            distance,
            elevation,
        }
    }

    /// Compute topological order (downstream cells first) using Kahn's
    /// algorithm on the downstream→upstream graph.
    ///
    /// Processes cells from outlets/pits upstream, so that when a cell
    /// is processed, all of its downstream cells have already been
    /// processed.
    fn topological_order(&self) -> Vec<usize> {
        let n = self.net.len();
        let downstream = &self.net.downstream;

        // pending[i] = 1 if cell i has a downstream cell (not a pit),
        // 0 if it's a pit/outlet. When the downstream cell is processed,
        // we decrement pending[i] to 0 and enqueue i.
        let mut pending: Vec<u8> = (0..n)
            .map(|i| u8::from(downstream[i] != usize::MAX))
            .collect();

        // Start with cells that have no downstream (pits/outlets).
        let mut queue: Vec<usize> = (0..n).filter(|&i| pending[i] == 0).collect();

        let mut order = Vec::with_capacity(n);
        while let Some(idx) = queue.pop() {
            order.push(idx);
            for &up in self.net.upstream(idx) {
                pending[up] -= 1;
                if pending[up] == 0 {
                    queue.push(up);
                }
            }
        }

        order
    }
}
