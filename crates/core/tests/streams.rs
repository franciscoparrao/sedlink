//! Tests for stream network vector export.

use ndarray::Array2;
use sedlink_core::{Network, stream_links, streams_geojson};
use surtgis_core::GeoTransform;
use surtgis_core::raster::Raster;

/// Two walled channels (cols 1 and 3, NoData elsewhere) converging
/// diagonally into a single outlet cell (4, 2): a Y-shaped network.
fn y_network_dem() -> Raster<f32> {
    let mut data = Array2::<f32>::from_elem((5, 5), f32::NAN);
    for r in 0..4 {
        data[(r, 1)] = 10.0 - r as f32 * 2.0;
        data[(r, 3)] = 10.0 - r as f32 * 2.0;
    }
    data[(4, 2)] = 1.0;
    let mut dem = Raster::from_array(data);
    dem.set_transform(GeoTransform::new(0.0, 0.0, 5.0, -5.0));
    dem
}

#[test]
fn test_stream_links_y_network() {
    let dem = y_network_dem();
    let net = Network::from_dem(&dem).unwrap();

    // Threshold 1: every valid cell is a stream cell.
    let links = stream_links(&net, 1.0);

    // Two head links (4 cells each + the junction cell) and the junction
    // link (single cell, also the outlet).
    assert_eq!(links.len(), 3, "expected 2 head links + 1 junction link");

    let mut head_links = 0;
    for link in &links {
        match link.cells.len() {
            5 => {
                head_links += 1;
                assert_eq!(link.order, 1);
                // 3 cardinal steps + 1 diagonal step of 5 m cells.
                let expected = 3.0 * 5.0 + 5.0 * std::f64::consts::SQRT_2;
                assert!((link.length - expected).abs() < 1e-9);
            }
            1 => {
                assert_eq!(link.order, 2, "junction cell should be order 2");
                assert!((link.outlet_acc - 9.0).abs() < 1e-9);
            }
            other => panic!("unexpected link of {other} cells"),
        }
    }
    assert_eq!(head_links, 2);

    // Same-family regression: stream magnitude at the junction must sum
    // its two order-1 tributaries (1 + 1 + 1 = 3) — this and the order-2
    // check both fail if the topological walk runs downstream-first.
    let mag = sedlink_core::NetworkAnalysis::new(&net).stream_magnitude(1.0);
    assert!((mag.values[4 * 5 + 2] - 3.0).abs() < 1e-9);
}

#[test]
fn test_streams_geojson_shape() {
    let dem = y_network_dem();
    let net = Network::from_dem(&dem).unwrap();
    let links = stream_links(&net, 1.0);
    let geojson = streams_geojson(&net, &links);

    assert!(geojson.starts_with("{\"type\":\"FeatureCollection\""));
    assert_eq!(geojson.matches("{\"type\":\"Feature\"").count(), 3);
    assert_eq!(geojson.matches("LineString").count(), 2);
    assert_eq!(geojson.matches("Point").count(), 1);
    // Cell centre of (4, 2) with 5 m cells from origin (0, 0):
    // x = (2 + 0.5) * 5 = 12.5, y = -(4 + 0.5) * 5 = -22.5.
    assert!(geojson.contains("[12.5,-22.5]"));
    // Valid JSON (parseable by any consumer).
    assert!(geojson.trim_end().ends_with("]}"));
}
