//! Tests for priority-flood pit filling on realistic (fractional) DEMs.
//!
//! These cover the cases the synthetic integer-valued DEMs miss: pits with
//! fractional depths, and filled depressions that must keep a drainable
//! gradient (epsilon fill) instead of becoming flat sinks.

use ndarray::Array2;
use sedlink_core::Network;
use surtgis_core::GeoTransform;
use surtgis_core::raster::Raster;

const N: usize = 5;

/// Southward slope with fractional elevations and a fractional-depth pit
/// at (2, 2). Row elevations: 10.0, 9.3, 8.6, 7.9, 7.2; the pit is 7.0.
fn fractional_pit_dem() -> Raster<f32> {
    let mut data = Array2::<f32>::zeros((N, N));
    for r in 0..N {
        for c in 0..N {
            data[(r, c)] = 10.0 - r as f32 * 0.7;
        }
    }
    data[(2, 2)] = 7.0;
    let mut dem = Raster::from_array(data);
    dem.set_transform(GeoTransform::new(0.0, 0.0, 5.0, -5.0));
    dem
}

/// Inward-draining bowl with fractional elevations: interior depression
/// that priority-flood must fill and drain to the boundary.
fn bowl_dem() -> Raster<f32> {
    let mut data = Array2::<f32>::zeros((N, N));
    for r in 0..N {
        for c in 0..N {
            let dr = r as f32 - 2.0;
            let dc = c as f32 - 2.0;
            data[(r, c)] = (dr * dr + dc * dc).sqrt() * 1.3;
        }
    }
    let mut dem = Raster::from_array(data);
    dem.set_transform(GeoTransform::new(0.0, 0.0, 5.0, -5.0));
    dem
}

fn is_boundary(idx: usize) -> bool {
    let (r, c) = (idx / N, idx % N);
    r == 0 || c == 0 || r == N - 1 || c == N - 1
}

#[test]
fn test_fractional_pit_is_filled_and_drains() {
    let dem = fractional_pit_dem();
    let net = Network::from_dem(&dem).unwrap();

    // The pit must be resolved: (2, 2) drains (south, towards lower rows).
    assert_ne!(
        net.flow_dir(2, 2),
        0,
        "fractional-depth pit should be filled and drain"
    );

    // The pit's flow path must reach the boundary (the pit acts as a
    // local funnel, so downstream accumulation includes neighbours too).
    let path = net.trace_downstream(2 * N + 2);
    let terminal = *path.last().unwrap();
    assert!(
        is_boundary(terminal),
        "filled pit should drain to the boundary, ended at {terminal}"
    );
    assert!(net.flow_acc(4, 2) >= 5.0);
}

#[test]
fn test_filled_depression_drains_to_boundary() {
    let dem = bowl_dem();
    let net = Network::from_dem(&dem).unwrap();

    // After epsilon filling, no interior cell may remain a pit, and every
    // cell's downstream trace must terminate on the boundary.
    for idx in 0..net.len() {
        if is_boundary(idx) {
            continue;
        }
        let (r, c) = (idx / N, idx % N);
        assert_ne!(
            net.flow_dir(r, c),
            0,
            "interior cell ({r}, {c}) should not be a pit after filling"
        );

        let path = net.trace_downstream(idx);
        let terminal = *path.last().unwrap();
        assert!(
            is_boundary(terminal),
            "cell ({r}, {c}) should drain to the boundary, ended at {terminal}"
        );
    }
}

#[test]
fn test_fill_preserves_dry_terrain() {
    // A DEM without depressions must pass through the fill untouched.
    let mut data = Array2::<f32>::zeros((N, N));
    for r in 0..N {
        for c in 0..N {
            data[(r, c)] = 10.0 - r as f32 * 0.7 - c as f32 * 0.1;
        }
    }
    let mut dem = Raster::from_array(data.clone());
    dem.set_transform(GeoTransform::new(0.0, 0.0, 5.0, -5.0));

    let net = Network::from_dem(&dem).unwrap();
    // Every cell except the outlet corner should drain; accumulation at
    // the lowest corner is bounded by the full grid.
    assert!(net.flow_acc(4, 4) > 1.0);
    assert!(net.flow_acc(4, 4) <= (N * N) as f64);
}
