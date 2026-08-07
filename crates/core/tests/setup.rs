//! Tests for solver setup extraction (`ChannelSetup`).

use approx::assert_relative_eq;
use ndarray::Array2;
use sedlink_core::{ChannelSetup, Network, SedlinkError};
use surtgis_core::GeoTransform;
use surtgis_core::raster::Raster;

const N: usize = 5;

/// Uniform south slope: 2 m drop per 5 m cell → gradient 0.4.
fn simple_slope_dem() -> Raster<f32> {
    let mut data = Array2::<f32>::zeros((N, N));
    for r in 0..N {
        for c in 0..N {
            data[(r, c)] = 10.0 - r as f32 * 2.0;
        }
    }
    let mut dem = Raster::from_array(data);
    dem.set_transform(GeoTransform::new(0.0, 0.0, 5.0, -5.0));
    dem
}

fn flat_dem_data(dem: &Raster<f32>) -> Vec<f32> {
    dem.data().as_slice().unwrap().to_vec()
}

#[test]
fn test_setup_reach_slope_and_basin() {
    let dem = simple_slope_dem();
    let net = Network::from_dem(&dem).unwrap();
    let z = flat_dem_data(&dem);

    // Pour point at the top of column 2, no snapping: the reach runs the
    // full column (4 steps × 5 m = 20 m, dropping 10 - 2 = 8 m).
    let setup = ChannelSetup::derive(&net, &z, (0, 2), 0, 5.0).unwrap();

    assert_eq!(setup.inflow, (0, 2));
    assert_relative_eq!(setup.channel_length, 20.0, epsilon = 1e-10);
    assert_relative_eq!(setup.channel_drop, 8.0, epsilon = 1e-6);
    assert_relative_eq!(setup.mean_channel_slope, 0.4, epsilon = 1e-6);

    // Only the pour point itself drains through it (headwater cell).
    assert_eq!(setup.basin_cells, 1);
    assert_relative_eq!(setup.basin_area, 25.0, epsilon = 1e-10);
    assert_relative_eq!(setup.inflow_accumulation, 1.0, epsilon = 1e-10);
}

#[test]
fn test_setup_basin_at_outlet() {
    let dem = simple_slope_dem();
    let net = Network::from_dem(&dem).unwrap();
    let z = flat_dem_data(&dem);

    // Outlet of column 2: the whole column drains through it.
    let setup = ChannelSetup::derive(&net, &z, (4, 2), 0, 5.0).unwrap();

    assert_eq!(setup.inflow, (4, 2));
    assert_eq!(setup.basin_cells, N);
    assert_relative_eq!(setup.inflow_accumulation, 5.0, epsilon = 1e-10);
    // At the outlet there is no downstream path left.
    assert_relative_eq!(setup.channel_length, 0.0, epsilon = 1e-10);
    assert_relative_eq!(setup.mean_channel_slope, 0.0, epsilon = 1e-10);
    // With threshold 5 only the outlet row qualifies as channel, and the
    // basin holds exactly one such cell.
    assert_eq!(setup.stream_cells, 1);
}

#[test]
fn test_setup_snaps_to_channel() {
    let dem = simple_slope_dem();
    let net = Network::from_dem(&dem).unwrap();
    let z = flat_dem_data(&dem);

    // A hillslope pour point with radius 3 snaps onto the bottom row,
    // where accumulation peaks.
    let setup = ChannelSetup::derive(&net, &z, (1, 1), 3, 5.0).unwrap();
    assert_eq!(setup.inflow.0, 4, "should snap to the outlet row");
    assert_relative_eq!(setup.inflow_accumulation, 5.0, epsilon = 1e-10);
}

#[test]
fn test_setup_json_is_parseable() {
    let dem = simple_slope_dem();
    let net = Network::from_dem(&dem).unwrap();
    let z = flat_dem_data(&dem);
    let setup = ChannelSetup::derive(&net, &z, (0, 2), 0, 5.0).unwrap();

    let json = setup.to_json();
    for key in [
        "inflow_row",
        "inflow_col",
        "mean_channel_slope",
        "basin_area_m2",
        "stream_cells",
    ] {
        assert!(json.contains(key), "JSON missing key {key}:\n{json}");
    }
    assert!(json.trim_start().starts_with('{') && json.trim_end().ends_with('}'));

    // Every value must be a finite number: NaN and inf are not valid
    // JSON literals and would break any downstream parser.
    const SEP: &str = "\": ";
    for line in json.lines().filter(|l| l.contains(SEP)) {
        let value = line.split_once(SEP).unwrap().1.trim().trim_end_matches(',');
        let parsed: f64 = value
            .parse()
            .unwrap_or_else(|_| panic!("value {value:?} is not a number in:\n{json}"));
        assert!(parsed.is_finite(), "non-finite value in:\n{json}");
    }
}

#[test]
fn test_setup_rejects_bad_inputs() {
    let dem = simple_slope_dem();
    let net = Network::from_dem(&dem).unwrap();
    let z = flat_dem_data(&dem);

    // Pour point outside the grid.
    assert!(matches!(
        ChannelSetup::derive(&net, &z, (99, 0), 0, 5.0),
        Err(SedlinkError::InvalidParam { .. })
    ));
    // Invalid threshold.
    assert!(matches!(
        ChannelSetup::derive(&net, &z, (0, 2), 0, 0.0),
        Err(SedlinkError::InvalidParam { .. })
    ));
    // DEM slice of the wrong length.
    assert!(matches!(
        ChannelSetup::derive(&net, &z[..10], (0, 2), 0, 5.0),
        Err(SedlinkError::GridMismatch { .. })
    ));
}
