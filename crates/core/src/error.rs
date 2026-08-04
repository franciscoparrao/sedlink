//! Error types for the sedlink crate.

use thiserror::Error;

/// Errors returned by the sedlink public API.
#[derive(Debug, Error)]
#[non_exhaustive]
pub enum SedlinkError {
    /// The DEM raster is empty (zero rows or columns).
    #[error("DEM raster is empty")]
    EmptyGrid,

    /// The DEM has rotation terms in its geotransform; sedlink requires an
    /// axis-aligned, north-up grid.
    #[error(
        "rotated grids are not supported (row_rotation={row_rotation}, col_rotation={col_rotation})"
    )]
    RotatedGrid {
        /// Rotation about the X axis from the geotransform.
        row_rotation: f64,
        /// Rotation about the Y axis from the geotransform.
        col_rotation: f64,
    },

    /// The DEM cells are not square.
    #[error("non-square cells: pixel_width={pixel_width}, pixel_height={pixel_height}")]
    NonSquareCells {
        /// Pixel width from the geotransform.
        pixel_width: f64,
        /// Pixel height from the geotransform (sign included).
        pixel_height: f64,
    },

    /// A secondary raster does not share the DEM's dimensions.
    #[error(
        "grid mismatch: DEM is {expected_rows}x{expected_cols}, other raster is {got_rows}x{got_cols}"
    )]
    GridMismatch {
        /// DEM rows.
        expected_rows: usize,
        /// DEM columns.
        expected_cols: usize,
        /// Rows of the offending raster.
        got_rows: usize,
        /// Columns of the offending raster.
        got_cols: usize,
    },

    /// A secondary raster shares dimensions but not the DEM's geotransform.
    #[error("geotransform mismatch between DEM and secondary raster")]
    TransformMismatch,

    /// A parameter is outside its valid range.
    #[error("invalid parameter {name}={value}: must satisfy {constraint}")]
    InvalidParam {
        /// Parameter name.
        name: &'static str,
        /// Offending value.
        value: f64,
        /// Human-readable constraint that was violated.
        constraint: &'static str,
    },

    /// The DEM contains no valid (non-NoData) cells.
    #[error("DEM contains no valid cells")]
    NoValidCells,

    /// A required raster was not provided.
    #[error("required raster '{name}' was not provided")]
    MissingRaster {
        /// Name of the missing raster.
        name: &'static str,
    },

    /// The flow network contains a cycle (should be impossible after pit
    /// resolution, but detected as a safety guard).
    #[error("flow network contains a cycle at cell (row={row}, col={col})")]
    CycleDetected {
        /// Row of the offending cell.
        row: usize,
        /// Column of the offending cell.
        col: usize,
    },
}
