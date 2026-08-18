//! Tests for the sediment connectivity index and network analysis.
//!
//! These tests verify numerical correctness of the IC computation against
//! hand-computed values and known properties of the Borselli index.

use approx::assert_relative_eq;
use ndarray::Array2;
use sedlink_core::{ConnectivityIndex, ConnectivityParams, Network, NetworkAnalysis, SedlinkError};
use surtgis_core::GeoTransform;
use surtgis_core::raster::Raster;

/// Row-major index in the 5×5 test grids.
const fn cell(r: usize, c: usize) -> usize {
    r * 5 + c
}

/// Build a simple synthetic DEM: a uniform southward slope with a channel
/// at the bottom row.
///
/// ```text
///  10  10  10  10  10
///   8   8   8   8   8
///   6   6   6   6   6
///   4   4   4   4   4
///   2   2   2   2   2
/// ```
///
/// All cells flow south (D8 dir 7). Flow accumulation increases southward.
fn simple_slope_dem() -> Raster<f32> {
    let rows = 5;
    let cols = 5;
    let mut data = Array2::<f32>::zeros((rows, cols));
    for r in 0..rows {
        for c in 0..cols {
            data[(r, c)] = 10.0 - r as f32 * 2.0;
        }
    }
    let mut dem = Raster::from_array(data);
    dem.set_transform(GeoTransform::new(0.0, 0.0, 5.0, -5.0));
    dem
}

/// Build a DEM with a convergent flow pattern (funnel) to test D8 routing.
fn funnel_dem() -> Raster<f32> {
    let rows = 5;
    let cols = 5;
    let mut data = Array2::<f32>::zeros((rows, cols));
    // Ridge at top, valley at bottom center
    for r in 0..rows {
        for c in 0..cols {
            let dist_from_center = ((c as f32) - 2.0).abs();
            data[(r, c)] = 10.0 - r as f32 - dist_from_center * 0.5;
        }
    }
    let mut dem = Raster::from_array(data);
    dem.set_transform(GeoTransform::new(0.0, 0.0, 5.0, -5.0));
    dem
}

#[test]
fn test_network_build_simple_slope() {
    let dem = simple_slope_dem();
    let net = Network::from_dem(&dem).unwrap();

    assert_eq!(net.rows(), 5);
    assert_eq!(net.cols(), 5);
    assert_eq!(net.cellsize(), 5.0);

    // All interior cells should flow south (D8 dir 7).
    for r in 0..4 {
        for c in 0..5 {
            assert_eq!(
                net.flow_dir(r, c),
                7,
                "cell ({}, {}) should flow south",
                r,
                c
            );
        }
    }

    // Bottom row cells are pits (dir 0).
    for c in 0..5 {
        assert_eq!(net.flow_dir(4, c), 0, "bottom row should be pit");
    }

    // Flow accumulation: in D8 southward flow, each cell receives flow from
    // only the cell directly above it. So row r has acc = r + 1.
    for r in 0..5 {
        for c in 0..5 {
            let expected = (r + 1) as f64;
            assert_relative_eq!(net.flow_acc(r, c), expected, epsilon = 1e-10);
        }
    }
}

#[test]
fn test_network_slope_uniform() {
    let dem = simple_slope_dem();
    let net = Network::from_dem(&dem).unwrap();

    // Slope should be uniform: arctan(2.0 / 5.0) ≈ 0.3805 rad
    let expected_slope = (2.0_f64 / 5.0).atan() as f32;
    for r in 0..5 {
        for c in 0..5 {
            let slope = net.slope(r, c);
            assert_relative_eq!(slope, expected_slope, epsilon = 1e-4);
        }
    }
}

#[test]
fn test_ic_simple_slope() {
    let dem = simple_slope_dem();
    let net = Network::from_dem(&dem).unwrap();

    let params = ConnectivityParams {
        stream_threshold: 5.0, // Bottom row becomes stream (acc=5)
        ..Default::default()
    };

    let ic = ConnectivityIndex::compute(&net, &dem, &params).unwrap();

    // Bottom row cells are stream cells → IC should be clamped to +10.
    for c in 0..5 {
        let val = ic.ic[(4, c)];
        assert!(!val.is_nan(), "stream cell IC should not be NaN");
        assert_relative_eq!(val, 10.0, epsilon = 1e-10);
    }

    // Top row cells should have lower IC (farther from channel).
    let ic_top = ic.ic[(0, 2)];
    assert!(!ic_top.is_nan(), "top cell IC should not be NaN");
    assert!(
        ic_top < 10.0,
        "top cell IC should be < 10 (not a stream cell)"
    );
    assert!(ic_top >= -10.0, "top cell IC should be >= -10 (clamped)");
}

#[test]
fn test_ic_d_up_increases_upslope() {
    let dem = simple_slope_dem();
    let net = Network::from_dem(&dem).unwrap();

    let params = ConnectivityParams {
        stream_threshold: 5.0,
        ..Default::default()
    };

    let ic = ConnectivityIndex::compute(&net, &dem, &params).unwrap();

    // D_up = W * sin(θ) * sqrt(A * Δ²)
    // Since slope is uniform and W=1, D_up ∝ sqrt(A) ∝ sqrt(flow_acc)
    // Flow acc increases southward, so D_up should increase southward.
    let d_up_top = ic.d_up[(0, 2)];
    let d_up_mid = ic.d_up[(2, 2)];
    let d_up_bot = ic.d_up[(4, 2)];

    assert!(d_up_top > 0.0, "D_up should be positive");
    assert!(
        d_up_mid > d_up_top,
        "D_up should increase downstream (more accumulation)"
    );
    assert!(d_up_bot > d_up_mid, "D_up should increase downstream");
}

#[test]
fn test_ic_d_dn_increases_upslope() {
    let dem = simple_slope_dem();
    let net = Network::from_dem(&dem).unwrap();

    let params = ConnectivityParams {
        stream_threshold: 5.0,
        ..Default::default()
    };

    let ic = ConnectivityIndex::compute(&net, &dem, &params).unwrap();

    // D_dn = sum(d_i / (W_i * sin(θ_i))) along path to channel
    // Since slope is uniform and W=1, D_dn ∝ path_length
    // Path length increases northward, so D_dn should increase northward.
    let d_dn_top = ic.d_dn[(0, 2)];
    let d_dn_mid = ic.d_dn[(2, 2)];
    let d_dn_bot = ic.d_dn[(4, 2)];

    assert!(d_dn_bot < 1e-10, "stream cell D_dn should be ~0");
    assert!(d_dn_mid > d_dn_bot, "D_dn should increase upstream");
    assert!(d_dn_top > d_dn_mid, "D_dn should increase upstream");
}

#[test]
fn test_ic_weighting_factor() {
    let dem = simple_slope_dem();
    let net = Network::from_dem(&dem).unwrap();

    // Without weight.
    let params_no_weight = ConnectivityParams {
        stream_threshold: 5.0,
        ..Default::default()
    };
    let ic_no_weight = ConnectivityIndex::compute(&net, &dem, &params_no_weight).unwrap();

    // With uniform weight = 2.0.
    let weight = ndarray::Array2::<f64>::from_elem((5, 5), 2.0);
    let params_weight = ConnectivityParams {
        stream_threshold: 5.0,
        weight: sedlink_core::WeightingFactor {
            raster: Some(weight),
            min_value: 0.001,
        },
        ..Default::default()
    };
    let ic_weight = ConnectivityIndex::compute(&net, &dem, &params_weight).unwrap();

    // With W=2 everywhere:
    // D_up ∝ W → doubles
    // D_dn ∝ 1/W → halves
    // IC = log10(D_up / D_dn) → log10(2 * D_up / (D_dn / 2)) = log10(4 * D_up / D_dn) = IC_no_weight + log10(4)
    let ic_top_no = ic_no_weight.ic[(0, 2)];
    let ic_top_w = ic_weight.ic[(0, 2)];

    assert_relative_eq!(ic_top_w, ic_top_no + 4.0_f64.log10(), epsilon = 1e-4);
}

#[test]
fn test_ic_dimension_mismatch() {
    let dem = simple_slope_dem();
    let net = Network::from_dem(&dem).unwrap();

    // Create a mismatched DEM.
    let mut bad_data = Array2::<f32>::zeros((3, 3));
    bad_data.fill(10.0);
    let mut bad_dem = Raster::from_array(bad_data);
    bad_dem.set_transform(GeoTransform::new(0.0, 5.0, 10.0, -10.0));

    let result = ConnectivityIndex::compute(&net, &bad_dem, &ConnectivityParams::default());
    assert!(result.is_err());
}

#[test]
fn test_strahler_order_simple() {
    let dem = simple_slope_dem();
    let net = Network::from_dem(&dem).unwrap();
    let analysis = NetworkAnalysis::new(&net);

    // With threshold=5, bottom row is stream (acc=5 ≥ 5).
    let order = analysis.strahler_order(5.0);

    // Bottom row cells: order 1 (no upstream stream cells).
    for c in 0..5 {
        assert_eq!(order.values[4 * 5 + c], 1, "bottom row should be order 1");
    }

    // Row 3 cells: not stream cells (acc=4 < 5), so order = 0.
    for c in 0..5 {
        assert_eq!(
            order.values[3 * 5 + c],
            0,
            "row 3 should be order 0 (not stream)"
        );
    }

    // Max order should be 1 (single-channel network).
    assert_eq!(order.max_order, 1);
}

#[test]
fn test_strahler_order_convergent() {
    let dem = funnel_dem();
    let net = Network::from_dem(&dem).unwrap();
    let analysis = NetworkAnalysis::new(&net);

    // With threshold=3, cells with acc >= 3 are streams.
    let order = analysis.strahler_order(3.0);

    // All stream cells should have order >= 1.
    for i in 0..net.len() {
        if net.flow_acc(i / net.cols(), i % net.cols()) >= 3.0 {
            assert!(order.values[i] >= 1, "stream cell should have order >= 1");
        }
    }

    // Max order should be at least 1 (there are stream cells).
    assert!(order.max_order >= 1);
}

#[test]
fn test_stream_magnitude() {
    let dem = simple_slope_dem();
    let net = Network::from_dem(&dem).unwrap();
    let analysis = NetworkAnalysis::new(&net);

    let mag = analysis.stream_magnitude(5.0);

    // Bottom row cells: magnitude 1 (no upstream stream cells).
    for c in 0..5 {
        assert_relative_eq!(mag.values[cell(4, c)], 1.0, epsilon = 1e-10);
    }

    // Row 0 cells: not stream cells (acc=1 < 5), so magnitude = 0.
    for c in 0..5 {
        assert_relative_eq!(mag.values[cell(0, c)], 0.0, epsilon = 1e-10);
    }
}

#[test]
fn test_longitudinal_profile() {
    let dem = simple_slope_dem();
    let net = Network::from_dem(&dem).unwrap();
    let analysis = NetworkAnalysis::new(&net);

    // Profile from top-center cell.
    let profile = analysis.longitudinal_profile(cell(0, 2), dem.data().as_slice().unwrap());

    // Should have 5 points (5 rows).
    assert_eq!(profile.distance.len(), 5);
    assert_eq!(profile.elevation.len(), 5);

    // Distance should increase downstream.
    for i in 1..profile.distance.len() {
        assert!(profile.distance[i] > profile.distance[i - 1]);
    }

    // Elevation should decrease downstream (southward slope).
    for i in 1..profile.elevation.len() {
        assert!(profile.elevation[i] < profile.elevation[i - 1]);
    }
}

#[test]
fn test_trace_downstream() {
    let dem = simple_slope_dem();
    let net = Network::from_dem(&dem).unwrap();

    // From top-center (0, 2), trace should go straight down to bottom.
    let path = net.trace_downstream(cell(0, 2));
    assert_eq!(path.len(), 4); // 4 downstream cells (rows 1-4)
    assert_eq!(path[0], cell(1, 2));
    assert_eq!(path[3], cell(4, 2));
}

#[test]
fn test_trace_to_stream() {
    let dem = simple_slope_dem();
    let net = Network::from_dem(&dem).unwrap();

    // From top-center (0, 2), trace to stream (threshold=5, bottom row is stream).
    let path = net.trace_to_stream(cell(0, 2), 5.0);
    // Should stop at the stream cell (row 4).
    assert!(!path.is_empty());
    // Last cell should be a stream cell.
    let last = *path.last().unwrap();
    let last_row = last / 5;
    let last_col = last % 5;
    assert!(net.flow_acc(last_row, last_col) >= 5.0);
}

#[test]
fn test_rotated_grid_rejected() {
    let mut dem = simple_slope_dem();
    // Set a rotation term.
    let t = dem.transform();
    dem.set_transform(GeoTransform {
        origin_x: t.origin_x,
        origin_y: t.origin_y,
        pixel_width: t.pixel_width,
        pixel_height: t.pixel_height,
        row_rotation: 0.1,
        col_rotation: 0.0,
    });

    let result = Network::from_dem(&dem);
    assert!(result.is_err());
}

#[test]
fn test_non_square_cells_rejected() {
    let mut dem = simple_slope_dem();
    let t = dem.transform();
    dem.set_transform(GeoTransform::new(t.origin_x, t.origin_y, 5.0, -10.0));

    let result = Network::from_dem(&dem);
    assert!(result.is_err());
}

#[test]
fn test_empty_grid_rejected() {
    let dem = Raster::<f32>::new(0, 0);
    let result = Network::from_dem(&dem);
    assert!(result.is_err());
}

#[test]
fn test_d_up_uses_upslope_means() {
    // South slope, cellsize 5, drop 2 per cell → slope gradient tan θ = 0.4.
    // Weight varies by row: w(r) = r + 1. For cell (2, 2) the upslope
    // contributing area is {(0,2), (1,2), (2,2)} with weights {1, 2, 3}:
    //   W̄ = 2.0, S̄ = 0.4, A = 3 · 25 m²
    //   D_up = W̄ · S̄ · √A = 2.0 · 0.4 · √75
    let dem = simple_slope_dem();
    let net = Network::from_dem(&dem).unwrap();

    let mut weight = Array2::<f64>::zeros((5, 5));
    for r in 0..5 {
        for c in 0..5 {
            weight[(r, c)] = (r + 1) as f64;
        }
    }

    let params = ConnectivityParams {
        stream_threshold: 5.0,
        weight: sedlink_core::WeightingFactor {
            raster: Some(weight),
            min_value: 0.001,
        },
        ..Default::default()
    };
    let ic = ConnectivityIndex::compute(&net, &dem, &params).unwrap();

    let expected_d_up = 2.0 * 0.4 * 75.0_f64.sqrt();
    assert_relative_eq!(ic.d_up[(2, 2)], expected_d_up, epsilon = 1e-5);

    // D_dn for (2, 2) uses LOCAL weights and slopes along the path to the
    // stream (row 4): cells (2,2) w=3 and (3,2) w=4, step 5 m, S = 0.4:
    //   D_dn = 5/(3·0.4) + 5/(4·0.4)
    let expected_d_dn = 5.0 / (3.0 * 0.4) + 5.0 / (4.0 * 0.4);
    assert_relative_eq!(ic.d_dn[(2, 2)], expected_d_dn, epsilon = 1e-5);
}

#[test]
fn test_ic_dinf_matches_d8_on_south_slope() {
    // On a uniform south slope the D∞ angle is exactly 3π/2, so all flow
    // goes to the south neighbour — routing, accumulation, and downslope
    // paths are identical to D8, and so must be the IC.
    let dem = simple_slope_dem();
    let d8 = Network::from_dem(&dem).unwrap();
    let dinf = sedlink_core::DinfNetwork::from_dem(&dem).unwrap();

    let params = ConnectivityParams {
        stream_threshold: 5.0,
        ..Default::default()
    };

    let ic_d8 = ConnectivityIndex::compute(&d8, &dem, &params).unwrap();
    let ic_dinf = ConnectivityIndex::compute(&dinf, &dem, &params).unwrap();

    for r in 0..5 {
        for c in 0..5 {
            let a = ic_d8.ic[(r, c)];
            let b = ic_dinf.ic[(r, c)];
            assert!(
                (a - b).abs() < 1e-9,
                "IC mismatch at ({r}, {c}): D8={a}, D∞={b}"
            );
        }
    }
}

#[test]
fn test_weight_raster_mismatch_rejected() {
    let dem = simple_slope_dem();
    let net = Network::from_dem(&dem).unwrap();

    let params = ConnectivityParams {
        stream_threshold: 5.0,
        weight: sedlink_core::WeightingFactor {
            raster: Some(Array2::<f64>::from_elem((3, 3), 1.0)),
            min_value: 0.001,
        },
        ..Default::default()
    };

    let result = ConnectivityIndex::compute(&net, &dem, &params);
    assert!(matches!(result, Err(SedlinkError::GridMismatch { .. })));
}

#[test]
fn test_ic_targets_basic() {
    // Single target at (2, 2) on the uniform south slope. Only column 2
    // above the target drains to it; everything else (other columns, and
    // the cells downstream of the target) never reaches a target.
    let dem = simple_slope_dem();
    let net = Network::from_dem(&dem).unwrap();

    let mut targets = Array2::<bool>::from_elem((5, 5), false);
    targets[(2, 2)] = true;

    let params = ConnectivityParams {
        targets: Some(targets),
        ..Default::default()
    };
    let ic = ConnectivityIndex::compute(&net, &dem, &params).unwrap();

    // Target cell: D_dn = 0 → IC = +clamp.
    assert_relative_eq!(ic.d_dn[(2, 2)], 0.0, epsilon = 1e-12);
    assert_relative_eq!(ic.ic[(2, 2)], params.clamp, epsilon = 1e-12);
    assert!(ic.is_stream[(2, 2)]);

    // Upstream of the target (column 2): D_dn accumulates 5 m steps with
    // W = 1 and S = 0.4 → 12.5 per cell.
    assert_relative_eq!(ic.d_dn[(1, 2)], 12.5, epsilon = 1e-5);
    assert_relative_eq!(ic.d_dn[(0, 2)], 25.0, epsilon = 1e-5);
    assert!(ic.ic[(1, 2)].is_finite());
    assert!(ic.ic[(0, 2)].is_finite());

    // Downstream of the target and all other columns: unreachable → NaN.
    assert!(ic.ic[(3, 2)].is_nan());
    assert!(ic.ic[(4, 2)].is_nan());
    for r in 0..5 {
        for c in [0usize, 1, 3, 4] {
            assert!(
                ic.ic[(r, c)].is_nan(),
                "cell ({r}, {c}) does not drain to the target"
            );
        }
    }
}

#[test]
fn test_ic_targets_equal_stream_mask_reproduces_classic_ic() {
    // Passing the threshold-based stream mask as the target mask must
    // reproduce the classic IC cell by cell, on both a parallel and a
    // convergent flow pattern.
    for dem in [simple_slope_dem(), funnel_dem()] {
        let net = Network::from_dem(&dem).unwrap();

        let classic = ConnectivityIndex::compute(
            &net,
            &dem,
            &ConnectivityParams {
                stream_threshold: 5.0,
                ..Default::default()
            },
        )
        .unwrap();

        let targeted = ConnectivityIndex::compute(
            &net,
            &dem,
            &ConnectivityParams {
                // Deliberately different threshold: it must be ignored.
                stream_threshold: 9999.0,
                targets: Some(classic.is_stream.clone()),
                ..Default::default()
            },
        )
        .unwrap();

        for r in 0..5 {
            for c in 0..5 {
                let a = classic.ic[(r, c)];
                let b = targeted.ic[(r, c)];
                assert!(
                    (a.is_nan() && b.is_nan()) || (a - b).abs() < 1e-12,
                    "IC mismatch at ({r}, {c}): classic={a}, targeted={b}"
                );
            }
        }
    }
}

#[test]
fn test_ic_targets_shape_mismatch_rejected() {
    let dem = simple_slope_dem();
    let net = Network::from_dem(&dem).unwrap();

    let params = ConnectivityParams {
        targets: Some(Array2::<bool>::from_elem((3, 3), true)),
        ..Default::default()
    };

    let result = ConnectivityIndex::compute(&net, &dem, &params);
    assert!(matches!(result, Err(SedlinkError::GridMismatch { .. })));
}

#[test]
fn test_invalid_params_rejected() {
    let dem = simple_slope_dem();
    let net = Network::from_dem(&dem).unwrap();

    for params in [
        ConnectivityParams {
            stream_threshold: -1.0,
            ..Default::default()
        },
        ConnectivityParams {
            clamp: 0.0,
            ..Default::default()
        },
        ConnectivityParams {
            min_slope: 0.0,
            ..Default::default()
        },
        ConnectivityParams {
            min_slope: 1.5,
            ..Default::default()
        },
    ] {
        let result = ConnectivityIndex::compute(&net, &dem, &params);
        assert!(matches!(result, Err(SedlinkError::InvalidParam { .. })));
    }
}
