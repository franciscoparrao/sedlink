//! Tests for the combined flood–sediment hazard classification.

use ndarray::Array2;
use sedlink_core::{HazardParams, SedlinkError, combined_hazard};

/// 3×3 case with hand-picked values covering the full matrix:
/// depth rows = {0.2 (low), 1.0 (medium), 3.0 (high)},
/// IC cols = {-4 (low), 0 (medium), 4 (high)} with explicit breaks.
#[test]
fn test_hazard_matrix_classes() {
    let depth =
        Array2::from_shape_vec((3, 3), vec![0.2, 0.2, 0.2, 1.0, 1.0, 1.0, 3.0, 3.0, 3.0]).unwrap();
    let ic = Array2::from_shape_vec((3, 3), vec![-4.0, 0.0, 4.0, -4.0, 0.0, 4.0, -4.0, 0.0, 4.0])
        .unwrap();

    let params = HazardParams {
        ic_breaks: Some([-2.0, 2.0]),
        ..Default::default()
    };
    let hazard = combined_hazard(&depth, &ic, &params).unwrap();

    // Full matrix, row-major: 1..9.
    for r in 0..3 {
        for c in 0..3 {
            let expected = (r * 3 + c + 1) as u8;
            assert_eq!(hazard.class[(r, c)], expected, "cell ({r}, {c})");
        }
    }
    assert_eq!(hazard.class_counts[0], 0);
    assert_eq!(hazard.ic_breaks, [-2.0, 2.0]);
}

#[test]
fn test_hazard_dry_and_nodata_are_zero() {
    let depth = Array2::from_shape_vec((2, 2), vec![0.0, 0.01, f64::NAN, 1.0]).unwrap();
    let ic = Array2::from_shape_vec((2, 2), vec![3.0, 3.0, 3.0, f64::NAN]).unwrap();

    let params = HazardParams {
        ic_breaks: Some([-2.0, 2.0]),
        ..Default::default()
    };
    let hazard = combined_hazard(&depth, &ic, &params).unwrap();

    // Dry (0.0 and 0.01 < wet_depth 0.05), NaN depth, NaN IC → all 0.
    assert!(hazard.class.iter().all(|&c| c == 0));
    assert_eq!(hazard.class_counts[0], 4);
}

#[test]
fn test_hazard_percentile_breaks_from_wet_area() {
    // 1×5: one dry cell (IC = 100 must NOT influence the percentiles)
    // and four wet cells with IC 1..4 → p33 ≈ 2, p66 ≈ 3.
    let depth = Array2::from_shape_vec((1, 5), vec![0.0, 1.0, 1.0, 1.0, 1.0]).unwrap();
    let ic = Array2::from_shape_vec((1, 5), vec![100.0, 1.0, 2.0, 3.0, 4.0]).unwrap();

    let hazard = combined_hazard(&depth, &ic, &HazardParams::default()).unwrap();

    assert!(hazard.ic_breaks[0] >= 1.0 && hazard.ic_breaks[0] <= 3.0);
    assert!(hazard.ic_breaks[1] > hazard.ic_breaks[0] || hazard.ic_breaks[1] <= 4.0);
    // Lowest-IC wet cell is class low (4 with medium depth 1.0 m),
    // highest-IC wet cell is class high (6).
    assert_eq!(hazard.class[(0, 1)], 4);
    assert_eq!(hazard.class[(0, 4)], 6);
    assert_eq!(hazard.class[(0, 0)], 0);
}

#[test]
fn test_hazard_all_dry_gives_nan_breaks() {
    let depth = Array2::from_shape_vec((1, 2), vec![0.0, 0.0]).unwrap();
    let ic = Array2::from_shape_vec((1, 2), vec![1.0, 2.0]).unwrap();
    let hazard = combined_hazard(&depth, &ic, &HazardParams::default()).unwrap();
    assert!(hazard.ic_breaks[0].is_nan());
    assert!(hazard.class.iter().all(|&c| c == 0));
}

#[test]
fn test_hazard_invalid_inputs_rejected() {
    let depth = Array2::<f64>::zeros((2, 2));
    let ic = Array2::<f64>::zeros((3, 3));
    assert!(matches!(
        combined_hazard(&depth, &ic, &HazardParams::default()),
        Err(SedlinkError::GridMismatch { .. })
    ));

    let ic = Array2::<f64>::zeros((2, 2));
    for params in [
        HazardParams {
            wet_depth: -1.0,
            ..Default::default()
        },
        HazardParams {
            depth_breaks: [2.0, 0.5],
            ..Default::default()
        },
        HazardParams {
            ic_breaks: Some([1.0, 1.0]),
            ..Default::default()
        },
    ] {
        assert!(matches!(
            combined_hazard(&depth, &ic, &params),
            Err(SedlinkError::InvalidParam { .. })
        ));
    }
}
