//! Terrain pre-processing for hydraulic solvers.
//!
//! A 2D shallow-water run needs, besides the DEM, a handful of terrain
//! facts that are normally extracted by hand or with external scripts:
//! where on the channel the inflow hydrograph enters, the mean channel
//! slope of the modelled reach, and which cells belong to the basin.
//! [`ChannelSetup`] derives all of them from the flow network so the
//! solver script carries no hard-coded row/column constants.
//!
//! The values map directly onto the constants a
//! [Hydroflux](https://github.com/franciscoparrao/hydroflux) `solver-2d`
//! run declares (`INFLOW_ROW`/`INFLOW_COL`, `SLOPE_MEAN`,
//! `ACC_THRESHOLD`); see `docs/hydroflux-coupling.md`.

use crate::{Network, NetworkAnalysis, SedlinkError};

/// Terrain facts describing the modelled reach and its basin.
#[derive(Debug, Clone)]
pub struct ChannelSetup {
    /// Inflow cell `(row, col)` — the pour point after snapping to the
    /// highest-accumulation cell in its neighbourhood.
    pub inflow: (usize, usize),
    /// Flow accumulation (cell count) at the inflow cell.
    pub inflow_accumulation: f64,
    /// Mean slope (m/m) of the flow path from the inflow to its outlet,
    /// i.e. total drop over total path length.
    pub mean_channel_slope: f64,
    /// Length (m) of that flow path.
    pub channel_length: f64,
    /// Elevation drop (m) along the path (start minus end).
    pub channel_drop: f64,
    /// Number of cells draining through the inflow cell.
    pub basin_cells: usize,
    /// Area (m²) of those cells.
    pub basin_area: f64,
    /// Number of channel cells in the basin at `stream_threshold`.
    pub stream_cells: usize,
}

impl ChannelSetup {
    /// Derive the setup for a pour point.
    ///
    /// `pour_point` is `(row, col)`; it is snapped to the highest
    /// accumulation cell within `snap_radius` (0 disables snapping).
    /// `stream_threshold` is the accumulation threshold delineating
    /// channels, used only for the `stream_cells` count.
    ///
    /// # Errors
    ///
    /// Returns [`SedlinkError::InvalidParam`] if the pour point lies
    /// outside the grid or on a `NoData` cell, or if `stream_threshold`
    /// is not finite and positive.
    pub fn derive(
        net: &Network,
        dem: &[f32],
        pour_point: (usize, usize),
        snap_radius: usize,
        stream_threshold: f64,
    ) -> Result<Self, SedlinkError> {
        let (rows, cols) = (net.rows(), net.cols());
        if pour_point.0 >= rows || pour_point.1 >= cols {
            return Err(SedlinkError::InvalidParam {
                name: "pour_point",
                value: pour_point.0 as f64,
                constraint: "inside the DEM grid",
            });
        }
        if !stream_threshold.is_finite() || stream_threshold <= 0.0 {
            return Err(SedlinkError::InvalidParam {
                name: "stream_threshold",
                value: stream_threshold,
                constraint: "finite and > 0",
            });
        }
        if dem.len() != rows * cols {
            return Err(SedlinkError::GridMismatch {
                expected_rows: rows,
                expected_cols: cols,
                got_rows: dem.len(),
                got_cols: 1,
            });
        }

        let idx = net.snap_to_stream(pour_point.0, pour_point.1, snap_radius);
        if net.is_solid(idx / cols, idx % cols) {
            return Err(SedlinkError::InvalidParam {
                name: "pour_point",
                value: pour_point.0 as f64,
                constraint: "on a valid (non-NoData) cell",
            });
        }
        let inflow = (idx / cols, idx % cols);

        // Longitudinal profile from the inflow to its outlet.
        let analysis = NetworkAnalysis::new(net);
        let profile = analysis.longitudinal_profile(idx, dem);
        let channel_length = profile.distance.last().copied().unwrap_or(0.0);
        let channel_drop = match (profile.elevation.first(), profile.elevation.last()) {
            (Some(&top), Some(&bottom)) => top - bottom,
            _ => 0.0,
        };
        let mean_channel_slope = if channel_length > 0.0 {
            channel_drop / channel_length
        } else {
            0.0
        };

        // Basin draining through the inflow cell.
        let labels = net.watersheds(&[idx]);
        let basin_cells = labels.iter().filter(|&&l| l == 1).count();
        let cellsize = net.cellsize();
        let basin_area = basin_cells as f64 * cellsize * cellsize;

        let flow_acc = net.flow_acc_slice();
        let stream_cells = (0..rows * cols)
            .filter(|&i| labels[i] == 1 && flow_acc[i] >= stream_threshold)
            .count();

        Ok(Self {
            inflow,
            inflow_accumulation: net.flow_acc(inflow.0, inflow.1),
            mean_channel_slope,
            channel_length,
            channel_drop,
            basin_cells,
            basin_area,
            stream_cells,
        })
    }

    /// Serialise as a flat JSON object (no external dependency).
    #[must_use]
    pub fn to_json(&self) -> String {
        format!(
            concat!(
                "{{\n",
                "  \"inflow_row\": {},\n",
                "  \"inflow_col\": {},\n",
                "  \"inflow_accumulation\": {},\n",
                "  \"mean_channel_slope\": {},\n",
                "  \"channel_length_m\": {},\n",
                "  \"channel_drop_m\": {},\n",
                "  \"basin_cells\": {},\n",
                "  \"basin_area_m2\": {},\n",
                "  \"stream_cells\": {}\n",
                "}}\n"
            ),
            self.inflow.0,
            self.inflow.1,
            self.inflow_accumulation,
            self.mean_channel_slope,
            self.channel_length,
            self.channel_drop,
            self.basin_cells,
            self.basin_area,
            self.stream_cells,
        )
    }
}
