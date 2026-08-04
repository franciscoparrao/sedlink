//! Tests for D-infinity flow routing.
//!
//! These tests verify numerical correctness of the D∞ computation against
//! hand-computed values, known properties of the Tarboton (1997) algorithm,
//! and parity with the D8 implementation in [`Network`].

use approx::assert_relative_eq;
use ndarray::Array2;
use sedlink_core::{DINF_PIT, DinfNetwork, Network};
use surtgis_core::GeoTransform;
use surtgis_core::raster::Raster;

/// Row-major index in the 5×5 test grids.
const fn cell(r: usize, c: usize) -> usize {
    r * 5 + c
}

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

fn funnel_dem() -> Raster<f32> {
    let rows = 5;
    let cols = 5;
    let mut data = Array2::<f32>::zeros((rows, cols));
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
fn test_dinf_south_slope_angle() {
    let dem = simple_slope_dem();
    let net = DinfNetwork::from_dem(&dem).unwrap();

    let target = 3.0 * std::f64::consts::FRAC_PI_2;
    for r in 1..4 {
        for c in 1..4 {
            let angle = net.angle(r, c);
            assert!(
                (angle - target).abs() < 0.5,
                "cell ({}, {}) should flow south (~{:.3} rad), got {:.3}",
                r,
                c,
                target,
                angle
            );
        }
    }
}

#[test]
fn test_dinf_bottom_row_is_pit() {
    let dem = simple_slope_dem();
    let net = DinfNetwork::from_dem(&dem).unwrap();

    for c in 0..5 {
        let angle = net.angle(4, c);
        assert!(
            angle < 0.0,
            "bottom row cell ({}, {}) should be a pit (angle < 0), got {:.3}",
            4,
            c,
            angle
        );
    }
}

#[test]
fn test_dinf_flow_acc_south_slope() {
    let dem = simple_slope_dem();
    let net = DinfNetwork::from_dem(&dem).unwrap();

    // Accumulation counts the cell itself (same convention as D8): row r
    // receives r upstream cells plus itself, so acc = r + 1.
    for r in 0..5 {
        for c in 0..5 {
            let acc = net.flow_acc(r, c);
            let expected = (r + 1) as f64;
            assert_relative_eq!(acc, expected, epsilon = 1e-6);
        }
    }
}

#[test]
fn test_dinf_flow_acc_headwater_one() {
    let dem = simple_slope_dem();
    let net = DinfNetwork::from_dem(&dem).unwrap();

    for c in 0..5 {
        let acc = net.flow_acc(0, c);
        assert_relative_eq!(acc, 1.0, epsilon = 1e-10);
    }
}

#[test]
fn test_dinf_angle_range() {
    let dem = funnel_dem();
    let net = DinfNetwork::from_dem(&dem).unwrap();

    let two_pi = 2.0 * std::f64::consts::PI;
    for r in 0..net.rows() {
        for c in 0..net.cols() {
            let a = net.angle(r, c);
            if a.is_nan() {
                continue;
            }
            assert!(
                a < 0.0 || (a >= 0.0 && a <= two_pi + 0.01),
                "angle at ({},{}) should be in [-1, 2π], got {}",
                r,
                c,
                a
            );
        }
    }
}

#[test]
fn test_dinf_trace_downstream_south_slope() {
    let dem = simple_slope_dem();
    let net = DinfNetwork::from_dem(&dem).unwrap();

    let path = net.trace_downstream(cell(0, 2));
    assert_eq!(path.len(), 4);
    assert_eq!(path[0], cell(1, 2));
    assert_eq!(path[3], cell(4, 2));
}

#[test]
fn test_dinf_pit_returns_dinf_pit() {
    let dem = simple_slope_dem();
    let net = DinfNetwork::from_dem(&dem).unwrap();

    for c in 0..5 {
        let angle = net.angle(4, c);
        assert_eq!(angle, DINF_PIT);
    }
}

#[test]
fn test_dinf_downstream_receivers() {
    let dem = simple_slope_dem();
    let net = DinfNetwork::from_dem(&dem).unwrap();

    let ds = net.downstream(cell(1, 2));
    let total_frac: f64 = ds.iter().map(|(_, f)| f).sum();
    assert_relative_eq!(total_frac, 1.0, epsilon = 1e-10);
    let primary = ds[0].0;
    assert_eq!(primary, cell(2, 2));
}

#[test]
fn test_dinf_mass_conservation_bowl() {
    // Inward-draining bowl: every valid cell's flow eventually
    // reaches the central pit → pit accumulation = N = 25 (incl. self).
    // Uses `from_filled_dem` to avoid priority-flood filling the bowl.
    let rows = 5;
    let cols = 5;
    let mut data = Array2::<f32>::zeros((rows, cols));
    for r in 0..rows {
        for c in 0..cols {
            let dr = r as f32 - 2.0;
            let dc = c as f32 - 2.0;
            data[(r, c)] = (dr * dr + dc * dc).sqrt() * 10.0;
        }
    }
    let mut dem = Raster::from_array(data);
    dem.set_transform(GeoTransform::new(0.0, 0.0, 5.0, -5.0));

    let net = DinfNetwork::from_filled_dem(&dem).unwrap();
    let center_acc = net.flow_acc(2, 2);

    assert!(
        (center_acc - 25.0).abs() < 1e-6,
        "bowl pit should accumulate all 25 cells (incl. self), got {}",
        center_acc
    );
}

#[test]
fn test_dinf_parity_with_d8_south_slope() {
    // Both conventions count the cell itself, so on a straight south
    // slope D∞ accumulation equals D8 exactly.
    let dem = simple_slope_dem();
    let d8_net = Network::from_dem(&dem).unwrap();
    let dinf_net = DinfNetwork::from_dem(&dem).unwrap();

    for r in 0..5 {
        for c in 0..5 {
            let d8_acc = d8_net.flow_acc(r, c);
            let dinf_acc = dinf_net.flow_acc(r, c);
            assert_relative_eq!(d8_acc, dinf_acc, epsilon = 1e-6);
        }
    }
}

#[test]
fn test_dinf_convergent_flow() {
    let dem = funnel_dem();
    let net = DinfNetwork::from_dem(&dem).unwrap();

    for r in 1..4 {
        for c in 1..4 {
            let angle = net.angle(r, c);
            assert!(
                !angle.is_nan() && angle >= 0.0,
                "interior cell ({}, {}) should have valid angle, got {}",
                r,
                c,
                angle
            );
        }
    }

    let acc_top = net.flow_acc(0, 2);
    let acc_bot = net.flow_acc(4, 2);
    assert!(
        acc_bot > acc_top,
        "accumulation should increase downstream: top={}, bot={}",
        acc_top,
        acc_bot
    );
}

#[test]
fn test_dinf_raster_output() {
    let dem = simple_slope_dem();
    let net = DinfNetwork::from_dem(&dem).unwrap();

    let acc_raster = net.flow_acc_raster().unwrap();
    let (rows, cols) = acc_raster.shape();
    assert_eq!(rows, 5);
    assert_eq!(cols, 5);

    // Accumulation includes the cell itself: row r has acc = r + 1.
    let acc = acc_raster.get(2, 2).unwrap();
    assert_relative_eq!(acc, 3.0, epsilon = 1e-6);
}

#[test]
fn test_dinf_rotated_grid_rejected() {
    let mut dem = simple_slope_dem();
    let t = dem.transform();
    dem.set_transform(GeoTransform {
        origin_x: t.origin_x,
        origin_y: t.origin_y,
        pixel_width: t.pixel_width,
        pixel_height: t.pixel_height,
        row_rotation: 0.1,
        col_rotation: 0.0,
    });

    let result = DinfNetwork::from_dem(&dem);
    assert!(result.is_err());
}

#[test]
fn test_dinf_non_square_cells_rejected() {
    let mut dem = simple_slope_dem();
    let t = dem.transform();
    dem.set_transform(GeoTransform::new(t.origin_x, t.origin_y, 5.0, -10.0));

    let result = DinfNetwork::from_dem(&dem);
    assert!(result.is_err());
}

#[test]
fn test_dinf_empty_grid_rejected() {
    let dem = Raster::<f32>::new(0, 0);
    let result = DinfNetwork::from_dem(&dem);
    assert!(result.is_err());
}
