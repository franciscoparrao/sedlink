//! Cross-implementation parity: sedlink vs the independent NumPy
//! reference in `tools/reference_ic.py`.
//!
//! The fixture DEM (40×40, 10 m cells) has fractional elevations, an
//! interior depression, and a NoData patch; the weighting raster varies
//! per cell. The reference implements the same documented conventions
//! (priority-flood-epsilon in f32, D8 by steepest gradient, tan-slope
//! clamped to [0.005, 1], upslope means for `D_up`) from the Borselli /
//! Cavalli equations.
//!
//! This is NOT SedInConnect parity (separate pending task): it guards
//! against transcription and indexing bugs by comparing two independent
//! implementations cell by cell.
//!
//! Regenerate the fixture with: `python3 tools/reference_ic.py`.

use ndarray::Array2;
use sedlink_core::{ConnectivityIndex, ConnectivityParams, Network, WeightingFactor};
use surtgis_core::GeoTransform;
use surtgis_core::raster::Raster;

const DEM_TXT: &str = include_str!("data/parity_dem.txt");
const WEIGHT_TXT: &str = include_str!("data/parity_weight.txt");
const ACC_TXT: &str = include_str!("data/parity_acc.txt");
const D_UP_TXT: &str = include_str!("data/parity_d_up.txt");
const D_DN_TXT: &str = include_str!("data/parity_d_dn.txt");
const IC_TXT: &str = include_str!("data/parity_ic.txt");

/// Parse a fixture grid: header `rows cols cellsize`, then row-major values.
fn parse_grid(text: &str) -> (usize, usize, f64, Vec<f64>) {
    let mut lines = text.lines();
    let header: Vec<&str> = lines.next().unwrap().split_whitespace().collect();
    let rows: usize = header[0].parse().unwrap();
    let cols: usize = header[1].parse().unwrap();
    let cellsize: f64 = header[2].parse().unwrap();

    let values: Vec<f64> = lines
        .flat_map(str::split_whitespace)
        .map(|v| v.parse().unwrap())
        .collect();
    assert_eq!(values.len(), rows * cols, "fixture size mismatch");
    (rows, cols, cellsize, values)
}

fn fixture_dem() -> (Raster<f32>, Array2<f64>) {
    let (rows, cols, cellsize, dem_vals) = parse_grid(DEM_TXT);
    let (wr, wc, _, weight_vals) = parse_grid(WEIGHT_TXT);
    assert_eq!((rows, cols), (wr, wc));

    let dem_f32: Vec<f32> = dem_vals.iter().map(|&v| v as f32).collect();
    let mut dem = Raster::from_array(Array2::from_shape_vec((rows, cols), dem_f32).unwrap());
    dem.set_transform(GeoTransform::new(0.0, 0.0, cellsize, -cellsize));

    let weight = Array2::from_shape_vec((rows, cols), weight_vals).unwrap();
    (dem, weight)
}

/// Both NaN, or numerically close (relative tolerance with an absolute
/// floor to absorb f32 slope rounding and summation-order differences).
fn check(name: &str, idx: (usize, usize), got: f64, want: f64, tol: f64) {
    if want.is_nan() {
        assert!(got.is_nan(), "{name} at {idx:?}: expected NaN, got {got}");
        return;
    }
    assert!(!got.is_nan(), "{name} at {idx:?}: expected {want}, got NaN");
    let diff = (got - want).abs();
    let scale = want.abs().max(1.0);
    assert!(
        diff <= tol * scale,
        "{name} at {idx:?}: got {got}, want {want} (diff {diff:.3e})"
    );
}

#[test]
fn test_parity_with_reference_implementation() {
    let (dem, weight) = fixture_dem();
    let net = Network::from_dem(&dem).unwrap();

    let (rows, cols, _, acc_want) = parse_grid(ACC_TXT);
    let (_, _, _, d_up_want) = parse_grid(D_UP_TXT);
    let (_, _, _, d_dn_want) = parse_grid(D_DN_TXT);
    let (_, _, _, ic_want) = parse_grid(IC_TXT);

    let params = ConnectivityParams {
        stream_threshold: 40.0,
        min_slope: 0.005,
        clamp: 10.0,
        weight: WeightingFactor {
            raster: Some(weight),
            min_value: 0.001,
        },
    };
    let ic = ConnectivityIndex::compute(&net, &dem, &params).unwrap();

    let mut checked = 0usize;
    for r in 0..rows {
        for c in 0..cols {
            let i = r * cols + c;

            // Flow accumulation: NoData cells keep acc = 1.0 in sedlink's
            // D8 network but are NaN in the reference outputs' D_up/D_dn;
            // compare accumulation only on valid cells.
            if !net.is_solid(r, c) {
                check("flow_acc", (r, c), net.flow_acc(r, c), acc_want[i], 1e-9);
            }

            check("d_up", (r, c), ic.d_up[(r, c)], d_up_want[i], 1e-6);
            check("d_dn", (r, c), ic.d_dn[(r, c)], d_dn_want[i], 1e-6);
            check("ic", (r, c), ic.ic[(r, c)], ic_want[i], 1e-6);
            checked += 1;
        }
    }
    assert_eq!(checked, rows * cols);

    // Sanity on the fixture itself: streams exist and NoData is confined
    // to the 4×4 patch.
    let stream_count = ic.is_stream.iter().filter(|&&s| s).count();
    assert_eq!(stream_count, 147, "fixture stream count changed");
    let nan_ic = ic.ic.iter().filter(|v| v.is_nan()).count();
    assert_eq!(nan_ic, 16, "only the NoData patch should be NaN");
}
