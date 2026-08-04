//! Tests for sediment routing with distance-decay SDR.
//!
//! Uses the 5×5 uniform south slope (cellsize 5 m, threshold 5 → bottom
//! row is the channel, straight south flow), where travel distances and
//! deliveries have closed forms: row r is L = (4 − r) · 5 m from the
//! channel.

use approx::assert_relative_eq;
use ndarray::Array2;
use sedlink_core::{
    ConnectivityIndex, ConnectivityParams, DinfNetwork, IcSdrParams, Network, RoutingParams,
    SedimentRouting, SedlinkError,
};
use surtgis_core::GeoTransform;
use surtgis_core::raster::Raster;

fn simple_slope_dem() -> Raster<f32> {
    let mut data = Array2::<f32>::zeros((5, 5));
    for r in 0..5 {
        for c in 0..5 {
            data[(r, c)] = 10.0 - r as f32 * 2.0;
        }
    }
    let mut dem = Raster::from_array(data);
    dem.set_transform(GeoTransform::new(0.0, 0.0, 5.0, -5.0));
    dem
}

/// q=1, L_ref=5 → SDR at row r is exp(−(4−r)).
fn test_params() -> RoutingParams {
    RoutingParams {
        sdr_exponent: 1.0,
        reference_length: 5.0,
        max_distance: 100_000.0,
    }
}

#[test]
fn test_sdr_distance_decay() {
    let dem = simple_slope_dem();
    let net = Network::from_dem(&dem).unwrap();

    let routing = SedimentRouting::compute(&net, None, 5.0, &test_params()).unwrap();

    for r in 0..4 {
        let expected_dist = (4 - r) as f64 * 5.0;
        let expected_sdr = (-((4 - r) as f64)).exp();
        for c in 0..5 {
            assert_relative_eq!(
                routing.dist_to_stream[(r, c)],
                expected_dist,
                epsilon = 1e-10
            );
            assert_relative_eq!(routing.sdr[(r, c)], expected_sdr, epsilon = 1e-12);
        }
    }

    // Stream cells: in the channel already.
    for c in 0..5 {
        assert_relative_eq!(routing.sdr[(4, c)], 1.0, epsilon = 1e-12);
        assert_relative_eq!(routing.dist_to_stream[(4, c)], 0.0, epsilon = 1e-12);
        assert!(routing.is_stream[(4, c)]);
    }
}

#[test]
fn test_hillslope_delivery_and_channel_flux() {
    let dem = simple_slope_dem();
    let net = Network::from_dem(&dem).unwrap();

    let routing = SedimentRouting::compute(&net, None, 5.0, &test_params()).unwrap();

    // Unit source: each channel cell receives its own 1.0 plus the
    // column's hillslope cells damped by exp(−k), k = rows above.
    let expected: f64 = 1.0 + (1..=4).map(|k| (-(k as f64)).exp()).sum::<f64>();
    for c in 0..5 {
        assert_relative_eq!(
            routing.hillslope_delivery[(4, c)],
            expected,
            epsilon = 1e-10
        );
        // Bottom-row channel cells are isolated pits, so the flux equals
        // the local delivery.
        assert_relative_eq!(routing.channel_flux[(4, c)], expected, epsilon = 1e-10);
    }

    // Hillslope cells deliver nothing locally and have no channel flux.
    assert_relative_eq!(routing.hillslope_delivery[(1, 2)], 0.0, epsilon = 1e-12);
    assert!(routing.channel_flux[(1, 2)].is_nan());
}

#[test]
fn test_max_distance_truncates_delivery() {
    let dem = simple_slope_dem();
    let net = Network::from_dem(&dem).unwrap();

    let params = RoutingParams {
        sdr_exponent: 1.0,
        reference_length: 5.0,
        max_distance: 7.0, // only row 3 (L = 5 m) can deliver
    };
    let routing = SedimentRouting::compute(&net, None, 5.0, &params).unwrap();

    for c in 0..5 {
        assert_relative_eq!(routing.sdr[(3, c)], (-1.0f64).exp(), epsilon = 1e-12);
        for r in 0..3 {
            assert_relative_eq!(routing.sdr[(r, c)], 0.0, epsilon = 1e-12);
        }
        let expected = 1.0 + (-1.0f64).exp();
        assert_relative_eq!(
            routing.hillslope_delivery[(4, c)],
            expected,
            epsilon = 1e-10
        );
    }
}

#[test]
fn test_source_raster_scales_delivery() {
    let dem = simple_slope_dem();
    let net = Network::from_dem(&dem).unwrap();

    let unit = SedimentRouting::compute(&net, None, 5.0, &test_params()).unwrap();

    let source = Array2::<f64>::from_elem((5, 5), 2.0);
    let doubled = SedimentRouting::compute(&net, Some(&source), 5.0, &test_params()).unwrap();

    for c in 0..5 {
        assert_relative_eq!(
            doubled.hillslope_delivery[(4, c)],
            2.0 * unit.hillslope_delivery[(4, c)],
            epsilon = 1e-10
        );
    }
}

#[test]
fn test_no_channel_means_unreachable() {
    let dem = simple_slope_dem();
    let net = Network::from_dem(&dem).unwrap();

    // Threshold higher than any accumulation → no stream cells at all.
    let routing = SedimentRouting::compute(&net, None, 100.0, &test_params()).unwrap();

    for r in 0..5 {
        for c in 0..5 {
            assert!(
                routing.sdr[(r, c)].is_nan(),
                "cell ({r}, {c}) should be unreachable without a channel"
            );
            assert!(routing.dist_to_stream[(r, c)].is_nan());
        }
    }
}

#[test]
fn test_routing_dinf_matches_d8_on_south_slope() {
    // Straight south slope: D∞ routing is identical to D8, so SDR and
    // deliveries must match exactly.
    let dem = simple_slope_dem();
    let d8 = Network::from_dem(&dem).unwrap();
    let dinf = DinfNetwork::from_dem(&dem).unwrap();

    let r8 = SedimentRouting::compute(&d8, None, 5.0, &test_params()).unwrap();
    let rf = SedimentRouting::compute(&dinf, None, 5.0, &test_params()).unwrap();

    for r in 0..5 {
        for c in 0..5 {
            let (a, b) = (r8.sdr[(r, c)], rf.sdr[(r, c)]);
            assert!(
                (a - b).abs() < 1e-9,
                "SDR mismatch at ({r}, {c}): D8={a}, D∞={b}"
            );
        }
    }
    for c in 0..5 {
        assert_relative_eq!(
            r8.hillslope_delivery[(4, c)],
            rf.hillslope_delivery[(4, c)],
            epsilon = 1e-9
        );
    }
}

#[test]
fn test_ic_sdr_follows_sigmoid() {
    let dem = simple_slope_dem();
    let net = Network::from_dem(&dem).unwrap();

    let cparams = ConnectivityParams {
        stream_threshold: 5.0,
        ..Default::default()
    };
    let ic = ConnectivityIndex::compute(&net, &dem, &cparams).unwrap();

    let params = IcSdrParams::default();
    let routing = SedimentRouting::compute_from_ic(&net, &ic, None, &params).unwrap();

    // Hillslope cells: SDR must equal the sigmoid of their IC value.
    for r in 0..4 {
        for c in 0..5 {
            let ic_val = ic.ic[(r, c)];
            let expected = params.sdr_max / (1.0 + ((params.ic0 - ic_val) / params.k).exp());
            assert_relative_eq!(routing.sdr[(r, c)], expected, epsilon = 1e-12);
            assert!(routing.sdr[(r, c)] > 0.0 && routing.sdr[(r, c)] < 1.0);
        }
    }

    // Stream cells (bottom row): SDR = 1, consistent with ic.is_stream.
    for c in 0..5 {
        assert!(ic.is_stream[(4, c)]);
        assert_relative_eq!(routing.sdr[(4, c)], 1.0, epsilon = 1e-12);
    }

    // IC increases towards the channel, so the sigmoid SDR must too.
    for c in 0..5 {
        for r in 0..3 {
            assert!(
                routing.sdr[(r + 1, c)] > routing.sdr[(r, c)],
                "SDR should increase downslope at col {c}, row {r}"
            );
        }
    }
}

#[test]
fn test_ic_sdr_invalid_inputs_rejected() {
    let dem = simple_slope_dem();
    let net = Network::from_dem(&dem).unwrap();
    let ic = ConnectivityIndex::compute(
        &net,
        &dem,
        &ConnectivityParams {
            stream_threshold: 5.0,
            ..Default::default()
        },
    )
    .unwrap();

    // Bad params.
    for bad in [
        IcSdrParams {
            sdr_max: 0.0,
            ..Default::default()
        },
        IcSdrParams {
            sdr_max: 1.5,
            ..Default::default()
        },
        IcSdrParams {
            k: 0.0,
            ..Default::default()
        },
        IcSdrParams {
            ic0: f64::NAN,
            ..Default::default()
        },
    ] {
        assert!(matches!(
            SedimentRouting::compute_from_ic(&net, &ic, None, &bad),
            Err(SedlinkError::InvalidParam { .. })
        ));
    }

    // IC computed on a different grid.
    let mut small = Array2::<f32>::zeros((3, 3));
    for r in 0..3 {
        for c in 0..3 {
            small[(r, c)] = 10.0 - r as f32;
        }
    }
    let mut small_dem = Raster::from_array(small);
    small_dem.set_transform(GeoTransform::new(0.0, 0.0, 5.0, -5.0));
    let small_net = Network::from_dem(&small_dem).unwrap();
    let small_ic = ConnectivityIndex::compute(
        &small_net,
        &small_dem,
        &ConnectivityParams {
            stream_threshold: 2.0,
            ..Default::default()
        },
    )
    .unwrap();
    assert!(matches!(
        SedimentRouting::compute_from_ic(&net, &small_ic, None, &IcSdrParams::default()),
        Err(SedlinkError::GridMismatch { .. })
    ));
}

#[test]
fn test_routing_invalid_inputs_rejected() {
    let dem = simple_slope_dem();
    let net = Network::from_dem(&dem).unwrap();

    // Bad params.
    let bad = RoutingParams {
        sdr_exponent: -1.0,
        ..RoutingParams::default()
    };
    assert!(matches!(
        SedimentRouting::compute(&net, None, 5.0, &bad),
        Err(SedlinkError::InvalidParam { .. })
    ));

    let bad = RoutingParams {
        reference_length: 0.0,
        ..RoutingParams::default()
    };
    assert!(matches!(
        SedimentRouting::compute(&net, None, 5.0, &bad),
        Err(SedlinkError::InvalidParam { .. })
    ));

    // Bad threshold.
    assert!(matches!(
        SedimentRouting::compute(&net, None, 0.0, &RoutingParams::default()),
        Err(SedlinkError::InvalidParam { .. })
    ));

    // Mismatched source raster.
    let source = Array2::<f64>::from_elem((3, 3), 1.0);
    assert!(matches!(
        SedimentRouting::compute(&net, Some(&source), 5.0, &RoutingParams::default()),
        Err(SedlinkError::GridMismatch { .. })
    ));
}
