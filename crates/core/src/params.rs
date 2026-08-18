//! Configuration parameters for connectivity and routing.

use crate::SedlinkError;

/// Weighting factor for sediment transport (W in the IC equations).
///
/// The weighting factor modulates both `D_up` and `D_dn`. Common choices:
/// - **C-factor** (soil erodibility, RUSLE): higher W → more sediment delivered.
/// - **NDVI** or **vegetation cover**: higher W → more resistance to transport.
/// - **Lithology / geological unit**: constant per unit.
///
/// If no weighting raster is provided, W defaults to 1.0 everywhere.
#[derive(Debug, Clone)]
pub struct WeightingFactor {
    /// Optional raster of per-cell weighting values. If `None`, W=1.0.
    pub raster: Option<ndarray::Array2<f64>>,
    /// Minimum weight to avoid division by zero in `D_dn`.
    pub min_value: f64,
}

impl Default for WeightingFactor {
    fn default() -> Self {
        Self {
            raster: None,
            min_value: 0.001,
        }
    }
}

/// Parameters for the Index of Connectivity (Borselli 2008 / Cavalli 2013).
#[derive(Debug, Clone)]
pub struct ConnectivityParams {
    /// Flow accumulation threshold (cell count) to delineate stream cells.
    /// Cells with flow accumulation ≥ this threshold are treated as channel
    /// heads; `D_dn` = 0 at stream cells. Ignored when [`targets`] is set.
    ///
    /// [`targets`]: ConnectivityParams::targets
    pub stream_threshold: f64,
    /// Optional target mask (`true` = target cell). When set, IC is computed
    /// relative to these targets (Cavalli et al. 2013 / `SedInConnect`
    /// "targets" version): `D_dn` is the impedance along the flow path to
    /// the nearest target, the stream network is not used as a destination,
    /// and cells that never drain to a target get IC = NaN. Typical targets:
    /// reservoirs, road networks, check dams, catchment outlets.
    pub targets: Option<ndarray::Array2<bool>>,
    /// Weighting factor (W). See [`WeightingFactor`].
    pub weight: WeightingFactor,
    /// Minimum slope gradient (tan θ, m/m) to avoid division by zero in
    /// `D_dn`. Cavalli et al. (2013) use 0.005. Must be in (0, 1].
    pub min_slope: f64,
    /// Clamp IC to [-clamp, +clamp]. Stream cells get +clamp.
    pub clamp: f64,
}

impl ConnectivityParams {
    /// Validate parameter ranges.
    ///
    /// # Errors
    ///
    /// Returns [`SedlinkError::InvalidParam`] if `stream_threshold` or
    /// `clamp` is not finite and positive, if `min_slope` is outside
    /// (0, 1], or if the weighting factor's `min_value` is not finite
    /// and positive.
    pub fn validate(&self) -> Result<(), SedlinkError> {
        if !self.stream_threshold.is_finite() || self.stream_threshold <= 0.0 {
            return Err(SedlinkError::InvalidParam {
                name: "stream_threshold",
                value: self.stream_threshold,
                constraint: "finite and > 0",
            });
        }
        if !self.min_slope.is_finite() || self.min_slope <= 0.0 || self.min_slope > 1.0 {
            return Err(SedlinkError::InvalidParam {
                name: "min_slope",
                value: self.min_slope,
                constraint: "in (0, 1]",
            });
        }
        if !self.clamp.is_finite() || self.clamp <= 0.0 {
            return Err(SedlinkError::InvalidParam {
                name: "clamp",
                value: self.clamp,
                constraint: "finite and > 0",
            });
        }
        if !self.weight.min_value.is_finite() || self.weight.min_value <= 0.0 {
            return Err(SedlinkError::InvalidParam {
                name: "weight.min_value",
                value: self.weight.min_value,
                constraint: "finite and > 0",
            });
        }
        Ok(())
    }
}

impl Default for ConnectivityParams {
    fn default() -> Self {
        Self {
            stream_threshold: 1000.0,
            targets: None,
            weight: WeightingFactor::default(),
            min_slope: 0.005,
            clamp: 10.0,
        }
    }
}

/// Parameters for sediment routing with a distance-decay SDR:
/// `SDR = exp(−sdr_exponent · L / reference_length)`.
#[derive(Debug, Clone)]
pub struct RoutingParams {
    /// Sediment delivery ratio (SDR) decay exponent `q`. Controls how
    /// fast delivery decreases with flow path length. Must be ≥ 0.
    pub sdr_exponent: f64,
    /// Reference length (m) for SDR scaling. Must be > 0.
    pub reference_length: f64,
    /// Maximum routing distance (m). Cells farther than this from the
    /// channel deliver nothing (SDR = 0). Must be > 0.
    pub max_distance: f64,
}

impl Default for RoutingParams {
    fn default() -> Self {
        Self {
            sdr_exponent: 0.5,
            reference_length: 1000.0,
            max_distance: 100_000.0,
        }
    }
}

/// Parameters for the IC-based SDR (Vigiak et al. 2012 / InVEST):
/// `SDR = sdr_max / (1 + exp((ic0 − IC) / k))`.
///
/// Defaults follow `InVEST`: `sdr_max = 0.8`, `ic0 = 0.5`, `k = 2.0`.
#[derive(Debug, Clone)]
pub struct IcSdrParams {
    /// Maximum SDR a hillslope cell can reach. Must be in (0, 1].
    pub sdr_max: f64,
    /// Sigmoid midpoint: IC value at which SDR = `sdr_max` / 2.
    pub ic0: f64,
    /// Sigmoid width: smaller `k` → sharper transition. Must be > 0.
    pub k: f64,
}

impl Default for IcSdrParams {
    fn default() -> Self {
        Self {
            sdr_max: 0.8,
            ic0: 0.5,
            k: 2.0,
        }
    }
}

impl IcSdrParams {
    /// Validate parameter ranges.
    ///
    /// # Errors
    ///
    /// Returns [`SedlinkError::InvalidParam`] if `sdr_max` is outside
    /// (0, 1], `ic0` is not finite, or `k` is not finite and > 0.
    pub fn validate(&self) -> Result<(), SedlinkError> {
        if !self.sdr_max.is_finite() || self.sdr_max <= 0.0 || self.sdr_max > 1.0 {
            return Err(SedlinkError::InvalidParam {
                name: "sdr_max",
                value: self.sdr_max,
                constraint: "in (0, 1]",
            });
        }
        if !self.ic0.is_finite() {
            return Err(SedlinkError::InvalidParam {
                name: "ic0",
                value: self.ic0,
                constraint: "finite",
            });
        }
        if !self.k.is_finite() || self.k <= 0.0 {
            return Err(SedlinkError::InvalidParam {
                name: "k",
                value: self.k,
                constraint: "finite and > 0",
            });
        }
        Ok(())
    }
}

impl RoutingParams {
    /// Validate parameter ranges.
    ///
    /// # Errors
    ///
    /// Returns [`SedlinkError::InvalidParam`] if `sdr_exponent` is not
    /// finite and ≥ 0, or if `reference_length` / `max_distance` is not
    /// finite and > 0.
    pub fn validate(&self) -> Result<(), SedlinkError> {
        if !self.sdr_exponent.is_finite() || self.sdr_exponent < 0.0 {
            return Err(SedlinkError::InvalidParam {
                name: "sdr_exponent",
                value: self.sdr_exponent,
                constraint: "finite and >= 0",
            });
        }
        if !self.reference_length.is_finite() || self.reference_length <= 0.0 {
            return Err(SedlinkError::InvalidParam {
                name: "reference_length",
                value: self.reference_length,
                constraint: "finite and > 0",
            });
        }
        if !self.max_distance.is_finite() || self.max_distance <= 0.0 {
            return Err(SedlinkError::InvalidParam {
                name: "max_distance",
                value: self.max_distance,
                constraint: "finite and > 0",
            });
        }
        Ok(())
    }
}
