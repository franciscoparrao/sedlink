//! Flow-network abstraction shared by D8 and D∞ routing.
//!
//! [`FlowNetwork`] exposes the grid geometry, flow accumulation, slope,
//! and downstream traversal that [`ConnectivityIndex`](crate::ConnectivityIndex)
//! needs, so the IC can be computed over either [`Network`](crate::Network)
//! (D8) or [`DinfNetwork`](crate::DinfNetwork) (D∞) interchangeably.
//!
//! ## Conventions
//!
//! - Flow accumulation is a cell count **including the cell itself**
//!   (headwater cells have accumulation 1) for both implementations.
//! - The downstream traversal follows a single path; for D∞ this is the
//!   receiver with the largest flow fraction.

use surtgis_core::GeoTransform;

/// A flow network derived from a DEM, usable by the connectivity index.
///
/// Implemented by [`Network`](crate::Network) (D8 single-flow-direction)
/// and [`DinfNetwork`](crate::DinfNetwork) (D∞, Tarboton 1997).
pub trait FlowNetwork {
    /// Number of rows.
    fn rows(&self) -> usize;

    /// Number of columns.
    fn cols(&self) -> usize;

    /// Total number of cells.
    fn len(&self) -> usize {
        self.rows() * self.cols()
    }

    /// `true` if the grid holds no cells.
    fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Uniform cell size in metres.
    fn cellsize(&self) -> f64;

    /// Geotransform of the domain.
    fn transform(&self) -> &GeoTransform;

    /// `true` if the cell at flat index `idx` is `NoData` (solid wall).
    fn is_solid_idx(&self, idx: usize) -> bool;

    /// Flow accumulation per cell (cell count including self).
    fn flow_acc_slice(&self) -> &[f64];

    /// Slope (radians) per cell, computed from the DEM.
    fn slope_slice(&self) -> &[f32];

    /// Primary downstream neighbour of `idx` and the step distance in
    /// metres, or `None` for pits and `NoData` cells.
    ///
    /// For D∞ networks the primary receiver is the one carrying the
    /// largest flow fraction.
    fn downstream_step(&self, idx: usize) -> Option<(usize, f64)>;

    /// Accumulate a per-cell quantity over each cell's upslope
    /// contributing area (including the cell itself).
    ///
    /// `sum[i]` = `local[i]` + Σ `local[j]` over all cells `j` draining
    /// through `i` (fractionally weighted for D∞). Dividing by
    /// [`flow_acc_slice`](Self::flow_acc_slice) yields the upslope mean
    /// (e.g. W̄ and S̄ in the Borselli/Cavalli IC equations).
    fn accumulate_upslope(&self, local: &[f64]) -> Vec<f64>;
}
