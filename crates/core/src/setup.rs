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
    /// Mean slope (m/m) of the flow path from the inflow to its outlet
    /// (or to `max_reach_length`), i.e. total drop over total length.
    pub mean_channel_slope: f64,
    /// Least-squares slope (m/m) fitted over the same profile. More
    /// robust than the endpoint mean when the reach mixes flats and
    /// steps; the two can differ by ~10% on real reaches.
    pub fitted_channel_slope: f64,
    /// Length (m) of the profiled flow path.
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
    /// `max_reach_length` truncates the profiled reach (m); `None`
    /// profiles down to the outlet.
    ///
    /// # Errors
    ///
    /// Returns [`SedlinkError::InvalidParam`] if the pour point lies
    /// outside the grid or on a `NoData` cell, if `stream_threshold`
    /// is not finite and positive, or if `max_reach_length` is given
    /// but not finite and positive.
    pub fn derive(
        net: &Network,
        dem: &[f32],
        pour_point: (usize, usize),
        snap_radius: usize,
        stream_threshold: f64,
        max_reach_length: Option<f64>,
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
        if let Some(len) = max_reach_length
            && (!len.is_finite() || len <= 0.0)
        {
            return Err(SedlinkError::InvalidParam {
                name: "max_reach_length",
                value: len,
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

        // Longitudinal profile from the inflow to its outlet, optionally
        // truncated to the requested reach length.
        let analysis = NetworkAnalysis::new(net);
        let profile = analysis.longitudinal_profile(idx, dem);
        let n_pts = match max_reach_length {
            Some(len) => profile
                .distance
                .iter()
                .take_while(|&&d| d <= len)
                .count()
                .max(1),
            None => profile.distance.len(),
        };
        let (dist, elev) = (&profile.distance[..n_pts], &profile.elevation[..n_pts]);

        let channel_length = dist.last().copied().unwrap_or(0.0);
        let channel_drop = match (elev.first(), elev.last()) {
            (Some(&top), Some(&bottom)) => top - bottom,
            _ => 0.0,
        };
        let mean_channel_slope = if channel_length > 0.0 {
            channel_drop / channel_length
        } else {
            0.0
        };
        let fitted_channel_slope = fit_slope(dist, elev);

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
            fitted_channel_slope,
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
                "  \"fitted_channel_slope\": {},\n",
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
            self.fitted_channel_slope,
            self.channel_length,
            self.channel_drop,
            self.basin_cells,
            self.basin_area,
            self.stream_cells,
        )
    }
}

/// Least-squares slope of elevation vs distance, sign-flipped so
/// descending reaches give a positive value. Returns 0 for fewer than
/// two points or zero distance variance.
fn fit_slope(dist: &[f64], elev: &[f64]) -> f64 {
    let n = dist.len();
    if n < 2 {
        return 0.0;
    }
    let nf = n as f64;
    let mean_d = dist.iter().sum::<f64>() / nf;
    let mean_z = elev.iter().sum::<f64>() / nf;
    let (mut cov, mut var) = (0.0, 0.0);
    for (&d, &z) in dist.iter().zip(elev) {
        cov += (d - mean_d) * (z - mean_z);
        var += (d - mean_d) * (d - mean_d);
    }
    if var > 0.0 { -cov / var } else { 0.0 }
}
