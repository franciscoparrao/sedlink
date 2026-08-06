//! Tests for watershed delineation on the D8 network.

use ndarray::Array2;
use sedlink_core::Network;
use surtgis_core::GeoTransform;
use surtgis_core::raster::Raster;

const N: usize = 5;

const fn cell(r: usize, c: usize) -> usize {
    r * N + c
}

/// Uniform south slope: every column drains straight south.
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

#[test]
fn test_watershed_single_column() {
    let dem = simple_slope_dem();
    let net = Network::from_dem(&dem).unwrap();

    // Pour point at the bottom of column 2: its basin is exactly column 2.
    let labels = net.watersheds(&[cell(4, 2)]);
    for r in 0..N {
        for c in 0..N {
            let expected = u32::from(c == 2);
            assert_eq!(
                labels[cell(r, c)],
                expected,
                "cell ({r}, {c}) label mismatch"
            );
        }
    }
}

#[test]
fn test_watershed_nested_pour_points() {
    let dem = simple_slope_dem();
    let net = Network::from_dem(&dem).unwrap();

    // Pour point 1 mid-column, pour point 2 at the bottom of the same
    // column: the nearest downstream pour point wins, so rows 0-2 belong
    // to basin 1 and rows 3-4 to basin 2.
    let labels = net.watersheds(&[cell(2, 2), cell(4, 2)]);
    for r in 0..N {
        let expected = if r <= 2 { 1 } else { 2 };
        assert_eq!(labels[cell(r, 2)], expected, "row {r} of column 2");
    }
    // Other columns drain past no pour point.
    assert_eq!(labels[cell(2, 0)], 0);
}

#[test]
fn test_snap_to_stream_finds_max_accumulation() {
    let dem = simple_slope_dem();
    let net = Network::from_dem(&dem).unwrap();

    // From a hillslope cell (1, 1), radius 3 reaches the bottom row where
    // accumulation peaks (acc = 5).
    let snapped = net.snap_to_stream(1, 1, 3);
    assert_eq!(snapped / N, 4, "should snap to the bottom row");

    // Radius 0 stays in place.
    assert_eq!(net.snap_to_stream(1, 1, 0), cell(1, 1));
}
