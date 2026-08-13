//! Vector export of the stream network.
//!
//! Splits the D8 stream network (cells with accumulation ≥ threshold)
//! into links at junctions: each link runs from a head or junction cell
//! downstream to the next junction or outlet (inclusive), so the
//! Strahler order is constant along a link. Links are exported as a
//! `GeoJSON` `FeatureCollection` of `LineString` features (single-cell
//! links become `Point` features) with order, length, and accumulation
//! attributes.
//!
//! Coordinates are cell centres in the DEM's CRS — note that strict
//! `GeoJSON` (RFC 7946) mandates `WGS84`; `QGIS` and most GIS tools read
//! projected-coordinate `GeoJSON` regardless.

use crate::network::D8_DISTANCE;
use crate::{Network, NetworkAnalysis};

/// One stream link (junction-to-junction reach).
#[derive(Debug, Clone)]
pub struct StreamLink {
    /// Flat cell indices along the link (upstream → downstream).
    pub cells: Vec<usize>,
    /// Strahler order (constant along the link).
    pub order: u8,
    /// Along-path length (m); 0 for single-cell links.
    pub length: f64,
    /// Flow accumulation at the downstream end.
    pub outlet_acc: f64,
}

/// Split the stream network into links at junctions.
///
/// # Panics
///
/// Panics only if internal invariants are violated (links always hold
/// at least their starting cell).
#[must_use]
pub fn stream_links(net: &Network, threshold: f64) -> Vec<StreamLink> {
    let n = net.len();
    // Solid (NoData) cells keep acc = 1.0 in the D8 network, so they
    // must be excluded explicitly before thresholding.
    let is_stream: Vec<bool> = (0..n)
        .map(|i| {
            !net.is_solid(i / net.cols(), i % net.cols()) && net.flow_acc_slice()[i] >= threshold
        })
        .collect();
    let up_streams: Vec<usize> = (0..n)
        .map(|i| net.upstream(i).iter().filter(|&&u| is_stream[u]).count())
        .collect();
    let order = NetworkAnalysis::new(net).strahler_order(threshold);
    let cellsize = net.cellsize();

    let mut links = Vec::new();
    for start in 0..n {
        if !is_stream[start] || (up_streams[start] != 0 && up_streams[start] < 2) {
            continue; // only heads (0 upstream) and junctions (≥2) start links
        }
        let mut cells = vec![start];
        let mut length = 0.0;
        let mut cur = start;
        loop {
            let dir = net.flow_dir(cur / net.cols(), cur % net.cols());
            let ds = net.downstream(cur);
            if dir == 0 || ds == usize::MAX || !is_stream[ds] {
                break;
            }
            length += D8_DISTANCE[(dir - 1) as usize] * cellsize;
            cells.push(ds);
            if up_streams[ds] >= 2 {
                break; // junction terminates the link (inclusive)
            }
            cur = ds;
        }
        let outlet_acc = net.flow_acc_slice()[*cells.last().unwrap()];
        links.push(StreamLink {
            order: order.values[start],
            cells,
            length,
            outlet_acc,
        });
    }
    links
}

/// Serialise links as a `GeoJSON` `FeatureCollection` (coordinates in the
/// DEM's CRS).
#[must_use]
pub fn streams_geojson(net: &Network, links: &[StreamLink]) -> String {
    let t = net.transform();
    let coord = |idx: usize| {
        let (r, c) = (idx / net.cols(), idx % net.cols());
        (
            t.origin_x + (c as f64 + 0.5) * t.pixel_width,
            t.origin_y + (r as f64 + 0.5) * t.pixel_height,
        )
    };

    let mut features = Vec::with_capacity(links.len());
    for (id, link) in links.iter().enumerate() {
        let coords: Vec<String> = link
            .cells
            .iter()
            .map(|&i| {
                let (x, y) = coord(i);
                format!("[{x},{y}]")
            })
            .collect();
        let geometry = if coords.len() == 1 {
            format!("{{\"type\":\"Point\",\"coordinates\":{}}}", coords[0])
        } else {
            format!(
                "{{\"type\":\"LineString\",\"coordinates\":[{}]}}",
                coords.join(",")
            )
        };
        features.push(format!(
            "{{\"type\":\"Feature\",\"geometry\":{geometry},\"properties\":{{\
             \"id\":{id},\"strahler\":{},\"length_m\":{},\"n_cells\":{},\
             \"outlet_accumulation\":{}}}}}",
            link.order,
            link.length,
            link.cells.len(),
            link.outlet_acc,
        ));
    }
    format!(
        "{{\"type\":\"FeatureCollection\",\"features\":[{}]}}\n",
        features.join(",")
    )
}
