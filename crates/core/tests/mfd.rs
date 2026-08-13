//! Tests for multiple flow direction (FD8) routing.

use approx::assert_relative_eq;
use ndarray::Array2;
use sedlink_core::{MfdNetwork, Network, SedlinkError};
use surtgis_core::GeoTransform;
use surtgis_core::raster::Raster;

fn dem_from(data: Array2<f32>) -> Raster<f32> {
    let mut dem = Raster::from_array(data);
    dem.set_transform(GeoTransform::new(0.0, 0.0, 5.0, -5.0));
    dem
}

/// Uniform south slope (SW/S/SE all descend: MFD disperses here).
fn simple_slope_dem() -> Raster<f32> {
    let mut data = Array2::<f32>::zeros((5, 5));
    for r in 0..5 {
        for c in 0..5 {
            data[(r, c)] = 10.0 - r as f32 * 2.0;
        }
    }
    dem_from(data)
}

#[test]
fn test_mfd_fractions_sum_to_one() {
    // Cone: divergent flow everywhere.
    let mut data = Array2::<f32>::zeros((7, 7));
    for r in 0..7 {
        for c in 0..7 {
            let (dr, dc) = (r as f32 - 3.0, c as f32 - 3.0);
            data[(r, c)] = 20.0 - (dr * dr + dc * dc).sqrt();
        }
    }
    let net = MfdNetwork::from_dem(&dem_from(data)).unwrap();

    let mut any_divergent = false;
    for idx in 0..49 {
        let recv = net.receivers(idx);
        if !recv.is_empty() {
            let total: f64 = recv.iter().map(|&(_, f)| f64::from(f)).sum();
            assert_relative_eq!(total, 1.0, epsilon = 1e-6);
        }
        if recv.len() > 1 {
            any_divergent = true;
        }
    }
    // The cone summit must spread flow over multiple receivers.
    assert!(any_divergent, "a cone should produce divergent flow");
    assert!(net.receivers(3 * 7 + 3).len() > 1, "summit should diverge");
}

/// One-cell-wide descending channel between NoData walls: each channel
/// cell has a single valid downslope neighbour, so MFD must equal D8.
fn walled_channel_dem() -> Raster<f32> {
    let mut data = Array2::<f32>::from_elem((5, 3), f32::NAN);
    for r in 0..5 {
        data[(r, 1)] = 10.0 - r as f32 * 2.0;
    }
    dem_from(data)
}

#[test]
fn test_mfd_matches_d8_on_walled_channel() {
    let dem = walled_channel_dem();
    let d8 = Network::from_dem(&dem).unwrap();
    let mfd = MfdNetwork::from_dem(&dem).unwrap();
    for r in 0..5 {
        assert_relative_eq!(mfd.flow_acc(r, 1), d8.flow_acc(r, 1), epsilon = 1e-9);
        assert_relative_eq!(mfd.flow_acc(r, 1), (r + 1) as f64, epsilon = 1e-9);
    }
}

#[test]
fn test_mfd_disperses_on_open_slope() {
    // On an open south slope SW/S/SE all descend: interior cells spread
    // over three receivers and headwaters still hold exactly 1.
    let dem = simple_slope_dem();
    let mfd = MfdNetwork::from_dem(&dem).unwrap();
    assert_eq!(mfd.receivers(2 * 5 + 2).len(), 3);
    for c in 0..5 {
        assert_relative_eq!(mfd.flow_acc(0, c), 1.0, epsilon = 1e-12);
    }
}

#[test]
fn test_mfd_large_exponent_converges_to_d8() {
    // Tilted plane with a diagonal bias: multiple downslope neighbours,
    // but with p = 50 nearly all flow follows the steepest one, so the
    // total accumulation along the D8 path approaches the D8 value.
    let mut data = Array2::<f32>::zeros((6, 6));
    for r in 0..6 {
        for c in 0..6 {
            data[(r, c)] = 30.0 - r as f32 * 2.0 - c as f32 * 0.6;
        }
    }
    let dem = dem_from(data);
    let d8 = Network::from_dem(&dem).unwrap();
    let sharp = MfdNetwork::from_dem_with_exponent(&dem, 50.0).unwrap();
    let soft = MfdNetwork::from_dem_with_exponent(&dem, 1.1).unwrap();

    // Interior cell deep in the slope: p=50 within 1% of D8; p=1.1 differs.
    let (r, c) = (4, 2);
    let d8_acc = d8.flow_acc(r, c);
    assert!(
        (sharp.flow_acc(r, c) - d8_acc).abs() / d8_acc < 0.01,
        "p=50 should converge to D8: mfd={}, d8={}",
        sharp.flow_acc(r, c),
        d8_acc
    );
    assert!(
        (soft.flow_acc(r, c) - d8_acc).abs() / d8_acc > 0.01,
        "p=1.1 should disperse: mfd={}, d8={}",
        soft.flow_acc(r, c),
        d8_acc
    );
}

#[test]
fn test_mfd_mass_conservation() {
    // Funnel draining south: total accumulation on the bottom row must
    // equal the full grid (every cell counted once, fractions included).
    let mut data = Array2::<f32>::zeros((5, 5));
    for r in 0..5 {
        for c in 0..5 {
            data[(r, c)] = 10.0 - r as f32 - (c as f32 - 2.0).abs() * 0.3;
        }
    }
    let net = MfdNetwork::from_dem(&dem_from(data)).unwrap();
    // Conservation is measured at the pits (cells with no receivers):
    // every cell's unit input must end up in exactly one of them. The
    // tolerance absorbs the f32 quantisation of the stored fractions
    // (relative error ~1e-7 per split).
    let total: f64 = (0..25)
        .filter(|&i| net.receivers(i).is_empty())
        .map(|i| net.flow_acc(i / 5, i % 5))
        .sum();
    assert_relative_eq!(total, 25.0, epsilon = 1e-5);
}

#[test]
fn test_mfd_invalid_exponent_rejected() {
    let dem = simple_slope_dem();
    for bad in [0.0, -1.0, f64::NAN] {
        assert!(matches!(
            MfdNetwork::from_dem_with_exponent(&dem, bad),
            Err(SedlinkError::InvalidParam { .. })
        ));
    }
}

#[test]
fn test_mfd_ic_runs_and_matches_d8_on_walled_channel() {
    use sedlink_core::{ConnectivityIndex, ConnectivityParams};
    let dem = walled_channel_dem();
    let d8 = Network::from_dem(&dem).unwrap();
    let mfd = MfdNetwork::from_dem(&dem).unwrap();
    let params = ConnectivityParams {
        stream_threshold: 5.0,
        ..Default::default()
    };
    let ic_d8 = ConnectivityIndex::compute(&d8, &dem, &params).unwrap();
    let ic_mfd = ConnectivityIndex::compute(&mfd, &dem, &params).unwrap();
    for r in 0..5 {
        let (a, b) = (ic_d8.ic[(r, 1)], ic_mfd.ic[(r, 1)]);
        assert!((a - b).abs() < 1e-9, "IC mismatch at row {r}: {a} vs {b}");
    }
}
