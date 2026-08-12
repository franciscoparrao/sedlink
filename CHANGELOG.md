# Changelog

All notable changes to this project will be documented in this file.
The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.3.0] - 2026-08-12

### Added
- Two-gradient flat resolution (Garbrecht & Martz 1997; Barnes et al.
  2014): BFS distance to flat outlets (primary) plus distance from
  higher terrain (tie-break), applied to both D8 and D∞. Natural plains
  and filled depressions now drain deterministically; flats touching the
  map edge drain off-grid.
- Watershed delineation (`Network::watersheds`) from pour points with
  nested-basin support, plus `Network::snap_to_stream` and a
  `sedlink watershed` CLI command.
- `parallel` feature (default): per-cell kernels (D8/D∞ flow directions,
  slope, D∞ receivers) run on the rayon thread pool with results
  identical to the sequential build.
- Solver setup extraction (`ChannelSetup` + `sedlink prep`): derives the
  inflow cell, mean reach slope, basin extent and channel-cell count for
  a pour point, emitting flat JSON plus optional accumulation and basin
  rasters. Replaces the terrain constants a 2D hydraulic run otherwise
  hard-codes; verified against Hydroflux's Huasco DEM
  (`docs/hydroflux-coupling.md`).

### Changed
- Dependencies: `surtgis-core` 1.0 → 1.2 (TIFF decoder hardening, NoData
  fixes) and `tiff` 0.10 → 0.11, keeping both crates on the same TIFF
  stack.
- **Numeric**: priority-flood now fills to the exact spill elevation
  (previously `next_up` epsilon); flat drainage comes from the flat
  resolver instead of ulp-sized gradients, giving more natural paths
  away from high terrain. Flow directions inside flats can differ from
  v0.2.0. The NumPy reference and parity fixture were updated to match.

## [0.2.0] - 2026-08-04

### Added
- D∞ flow routing (Tarboton 1997) with fractional flow accumulation
  (`DinfNetwork`).
- `FlowNetwork` trait: the connectivity index and sediment routing are
  generic over D8 and D∞ networks; CLI commands accept `--flow d8|dinf`.
- Sediment routing (`SedimentRouting`) with two SDR models: distance
  decay `exp(−q·L/L_ref)` (SEDD family) and the IC-based sigmoid of
  Vigiak et al. (2012) / InVEST (`compute_from_ic`), including
  hillslope-to-channel delivery and accumulated channel flux. New
  `sedlink route` CLI command.
- Parameter validation (`ConnectivityParams::validate`,
  `RoutingParams::validate`, `IcSdrParams::validate`) and dimension
  checks for weighting/source rasters — invalid inputs now return typed
  errors instead of panicking.
- Cross-implementation parity harness: independent NumPy reference
  (`tools/reference_ic.py`) and a cell-by-cell fixture test over a
  fractional DEM with an interior depression and NoData.

### Changed
- **Breaking**: `ConnectivityIndex::compute` now takes
  `&ConnectivityParams` (was by value) and is generic over
  `FlowNetwork`. Migration: pass `&params`.
- **Breaking / numeric**: the IC now follows Borselli 2008 /
  Cavalli 2013 exactly — `D_up` uses upslope means W̄·S̄ (previously
  local cell values) and the slope gradient is tan θ clamped to
  [0.005, 1.0] (previously sin θ, no upper cap). IC values change on
  any non-uniform terrain; the new values are the reference-correct ones.
- **Breaking / numeric**: D8 flow directions now use steepest gradient
  (drop/distance, the standard definition) instead of lowest neighbour
  elevation; diagonal-vs-cardinal choices can differ.
- D∞ flow accumulation now counts the cell itself (headwater = 1),
  matching the D8 convention.

### Fixed
- Priority-flood pit filling truncated `f32` elevations through an `i64`
  heap key, leaving fractional-depth pits unfilled and breaking drainage
  on real DEMs; the heap now orders true elevations (`ordered-float`).
- Filled depressions became flat sinks with no outlet; the fill is now
  the epsilon variant (Barnes et al. 2014), so filled areas drain.
- `D_dn` paths that end in a pit are cached as unreachable instead of
  being re-traced quadratically.

### Removed
- Phantom dependencies (`anyhow`, `serde`, `num-traits`, `rayon` in core;
  `indicatif`, `tracing` in the CLI) and the `parallel` feature flag that
  enabled no parallelism.

## [0.1.0] - 2026-06-10

### Added
- D8 flow network from a DEM: priority-flood pit filling, flow
  accumulation, downstream tracing (`Network`).
- Index of Connectivity (Borselli 2008 / Cavalli 2013) with configurable
  weighting factor (`ConnectivityIndex`).
- Network analysis: Strahler order, stream magnitude, longitudinal
  profiles (`NetworkAnalysis`).
- CLI (`sedlink`) with `ic`, `order`, `acc`, and `slope` commands over
  GeoTIFF via `surtgis-core`.
