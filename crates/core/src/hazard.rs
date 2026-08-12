//! Combined flood–sediment hazard classification.
//!
//! Bivariate hazard matrix over two co-registered rasters: simulated
//! flood depth (e.g. from a shallow-water solver such as
//! [Hydroflux](https://github.com/franciscoparrao/hydroflux)) and the
//! Index of Connectivity computed by sedlink. Each wet cell receives a
//! class from a 3×3 matrix:
//!
//! ```text
//!                 connectivity (IC)
//!                 low   med   high
//! depth  low       1     2     3
//!        medium    4     5     6
//!        high      7     8     9
//! ```
//!
//! so class = `(depth_class − 1) · 3 + ic_class`. Dry cells (depth below
//! [`HazardParams::wet_depth`]) and `NoData` cells are class 0.
//!
//! This is a **classification**, not new physics: the same
//! intensity-matrix practice used in Alpine flood hazard mapping
//! (e.g. Swiss federal guidelines, Loat & Petrascheck 1997, whose depth
//! thresholds 0.5 m / 2.0 m are the defaults here), with sediment
//! connectivity as the second axis. High classes mark cells that are
//! both deeply flooded and well connected to sediment sources.
//!
//! IC thresholds default to the 33rd/66th percentiles of the IC over
//! the wet area, i.e. a relative classification within the flooded
//! zone; pass explicit breaks for absolute comparisons across events.

use ndarray::Array2;

use crate::SedlinkError;

/// Parameters for the combined hazard classification.
#[derive(Debug, Clone)]
pub struct HazardParams {
    /// Minimum depth (m) for a cell to count as wet. Default 0.05.
    pub wet_depth: f64,
    /// Depth class boundaries (m): `< breaks[0]` low, `< breaks[1]`
    /// medium, else high. Default `[0.5, 2.0]` (Swiss intensity
    /// thresholds).
    pub depth_breaks: [f64; 2],
    /// IC class boundaries. `None` uses the 33rd/66th percentiles of
    /// the IC over the wet cells.
    pub ic_breaks: Option<[f64; 2]>,
}

impl Default for HazardParams {
    fn default() -> Self {
        Self {
            wet_depth: 0.05,
            depth_breaks: [0.5, 2.0],
            ic_breaks: None,
        }
    }
}

impl HazardParams {
    /// Validate parameter ranges.
    ///
    /// # Errors
    ///
    /// Returns [`SedlinkError::InvalidParam`] if `wet_depth` is not
    /// finite and ≥ 0, or if a break pair is not finite and strictly
    /// increasing.
    pub fn validate(&self) -> Result<(), SedlinkError> {
        if !self.wet_depth.is_finite() || self.wet_depth < 0.0 {
            return Err(SedlinkError::InvalidParam {
                name: "wet_depth",
                value: self.wet_depth,
                constraint: "finite and >= 0",
            });
        }
        let increasing = |b: &[f64; 2]| b[0].is_finite() && b[1].is_finite() && b[0] < b[1];
        if !increasing(&self.depth_breaks) {
            return Err(SedlinkError::InvalidParam {
                name: "depth_breaks",
                value: self.depth_breaks[0],
                constraint: "finite and strictly increasing",
            });
        }
        if let Some(b) = &self.ic_breaks
            && !increasing(b)
        {
            return Err(SedlinkError::InvalidParam {
                name: "ic_breaks",
                value: b[0],
                constraint: "finite and strictly increasing",
            });
        }
        Ok(())
    }
}

/// Result of a combined hazard classification.
#[derive(Debug, Clone)]
pub struct CombinedHazard {
    /// Class per cell: 0 = dry or `NoData`, 1–9 = hazard matrix
    /// (`(depth_class − 1) · 3 + ic_class`).
    pub class: Array2<u8>,
    /// IC breaks actually used (given, or wet-area percentiles). NaN
    /// when there are no wet cells with valid IC.
    pub ic_breaks: [f64; 2],
    /// Number of cells per class (index 0–9).
    pub class_counts: [usize; 10],
}

/// Classify combined flood–sediment hazard from co-registered depth and
/// IC rasters.
///
/// # Errors
///
/// Returns [`SedlinkError::InvalidParam`] on invalid parameters and
/// [`SedlinkError::GridMismatch`] when the rasters' dimensions differ.
pub fn combined_hazard(
    depth: &Array2<f64>,
    ic: &Array2<f64>,
    params: &HazardParams,
) -> Result<CombinedHazard, SedlinkError> {
    params.validate()?;

    let dim = depth.dim();
    if ic.dim() != dim {
        return Err(SedlinkError::GridMismatch {
            expected_rows: dim.0,
            expected_cols: dim.1,
            got_rows: ic.dim().0,
            got_cols: ic.dim().1,
        });
    }

    let wet = |h: f64| h.is_finite() && h >= params.wet_depth;

    // Resolve IC breaks: given, or 33rd/66th percentiles over wet cells.
    let ic_breaks = if let Some(b) = params.ic_breaks {
        b
    } else {
        let mut vals: Vec<f64> = depth
            .iter()
            .zip(ic.iter())
            .filter(|&(&h, &i)| wet(h) && i.is_finite())
            .map(|(_, &i)| i)
            .collect();
        if vals.is_empty() {
            [f64::NAN, f64::NAN]
        } else {
            vals.sort_by(f64::total_cmp);
            let pick = |p: f64| vals[((vals.len() - 1) as f64 * p).round() as usize];
            [pick(0.33), pick(0.66)]
        }
    };

    let class_of = |v: f64, breaks: &[f64; 2]| -> u8 {
        if v < breaks[0] {
            1
        } else if v < breaks[1] {
            2
        } else {
            3
        }
    };

    let mut class_counts = [0usize; 10];
    let class = Array2::from_shape_fn(dim, |rc| {
        let (h, i) = (depth[rc], ic[rc]);
        let c = if !wet(h) || !i.is_finite() || !ic_breaks[0].is_finite() {
            0
        } else {
            (class_of(h, &params.depth_breaks) - 1) * 3 + class_of(i, &ic_breaks)
        };
        class_counts[usize::from(c)] += 1;
        c
    });

    Ok(CombinedHazard {
        class,
        ic_breaks,
        class_counts,
    })
}
