//! Sedlink — Sediment connectivity index and fluvial network analysis.
//!
//! Implements the Index of Connectivity (IC) of Borselli et al. (2008) and
//! Cavalli et al. (2013), river network analysis (Strahler ordering,
//! magnitude, longitudinal profiles), and sediment routing with a
//! distance-decay sediment delivery ratio ([`SedimentRouting`]) over D8
//! or D∞ flow networks.
//!
//! ## Architecture
//!
//! The crate is organised into three layers:
//!
//! 1. **Flow network** ([`Network`]) — builds a D8 flow graph from a DEM,
//!    resolves pits, and provides downstream tracing with caching.
//!    [`DinfNetwork`] provides the same using D∞ (Tarboton 1997)
//!    continuous flow directions with fractional accumulation. Both
//!    implement the [`FlowNetwork`] trait.
//! 2. **Connectivity** ([`connectivity::ConnectivityIndex`]) — computes
//!    `D_up` (upslope sediment delivery) and `D_dn` (downslope impedance)
//!    over any [`FlowNetwork`], yielding IC = `log10(D_up / D_dn)`.
//! 3. **Network analysis** ([`network::NetworkAnalysis`]) — Strahler
//!    stream order, stream magnitude, and longitudinal profile extraction
//!    (D8 only).
//!
//! ## Entry point
//!
//! ```no_run
//! use sedlink_core::{ConnectivityIndex, ConnectivityParams, DinfNetwork, Network, NetworkAnalysis};
//! use surtgis_core::Raster;
//! # fn run(dem: Raster<f32>) -> Result<(), sedlink_core::SedlinkError> {
//! let net = Network::from_dem(&dem)?;
//! let ic = ConnectivityIndex::compute(&net, &dem, &ConnectivityParams::default())?;
//! let analysis = NetworkAnalysis::new(&net);
//! let order = analysis.strahler_order(1000.0);
//! # Ok(()) }
//! ```
//!
//! ## References
//!
//! - Borselli, L., Cassi, P., & Torri, D. (2008). Prolegomena to sediment
//!   and flow connectivity in the landscape. *CATENA*, 75(3), 268–277.
//! - Cavalli, M., et al. (2013). The Index of Connectivity: a valuable
//!   tool to understand the role of sediment connectivity in hydrology.
//!   *Journal of Hydrology*, 489, 145–156.

#![deny(missing_docs)]
#![warn(clippy::pedantic)]
#![allow(
    clippy::cast_precision_loss,
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::cast_possible_wrap,
    clippy::similar_names,
    clippy::module_name_repetitions
)]

mod connectivity;
mod dinf;
mod error;
mod flats;
mod flow;
mod network;
mod par;
mod params;
mod routing;
mod setup;

pub use connectivity::ConnectivityIndex;
pub use dinf::{DINF_PIT, DinfNetwork};
pub use error::SedlinkError;
pub use flow::FlowNetwork;
pub use network::{Network, NetworkAnalysis, StrahlerOrder, StreamProfile};
pub use params::{ConnectivityParams, IcSdrParams, RoutingParams, WeightingFactor};
pub use routing::SedimentRouting;
pub use setup::ChannelSetup;
