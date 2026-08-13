# sedlink

**Sediment connectivity index and fluvial network analysis in Rust.**

sedlink computes the Index of Connectivity (IC) of Borselli et al. (2008) /
Cavalli et al. (2013), fluvial network metrics (Strahler order, stream
magnitude, longitudinal profiles), and sediment routing with a sediment
delivery ratio (SDR) — all from a DEM, in a single native binary with no
Python or GIS-suite dependency.

Part of the [SurtGIS](https://github.com/franciscoparrao/surtgis) ecosystem:
raster I/O and grid types come from `surtgis-core`.

## Features

- **Flow networks**: D8 (steepest gradient), D∞ (Tarboton 1997) and
  MFD/FD8 (Freeman 1991 / Holmgren 1994), all behind a common
  `FlowNetwork` trait. Priority-flood pit filling with
  two-gradient flat resolution (Garbrecht & Martz 1997; Barnes et al.
  2014), so filled depressions and natural plains drain deterministically.
- **Watersheds**: basin delineation from pour points (nested basins
  supported) with accumulation-snapping of outlet coordinates.
- **Parallel**: per-cell kernels run on all cores (rayon, default
  `parallel` feature) with results identical to the sequential build.
- **Index of Connectivity**: `IC = log10(D_up / D_dn)` with upslope means
  W̄·S̄ over the contributing area and slope gradient tan θ clamped to
  [0.005, 1.0], following Cavalli et al. (2013).
- **Network analysis**: Strahler stream order, stream magnitude,
  longitudinal profiles per tributary.
- **Sediment routing**: two SDR models — distance decay
  `SDR = exp(−q·L/L_ref)` (SEDD family, Ferro & Porto 2000) and the
  IC-based sigmoid `SDR = SDR_max / (1 + exp((IC0 − IC)/k))`
  (Vigiak et al. 2012 / InVEST) — with hillslope-to-channel delivery and
  accumulated channel flux.
- **Validated**: 50-test suite including hand-computed fixtures and a
  cell-by-cell parity harness against an independent NumPy reference
  implementation (`tools/reference_ic.py`).

## Install

```bash
cargo install sedlink-cli
```

Or build from source:

```bash
git clone https://github.com/franciscoparrao/sedlink
cd sedlink && cargo build --release
```

## Usage

```bash
# Index of Connectivity (D8; use --flow dinf for D∞)
sedlink ic --dem dem.tif --output ic.tif --threshold 1000 [--weight cfactor.tif]

# Strahler stream order
sedlink order --dem dem.tif --output order.tif --threshold 1000

# Flow accumulation / slope
sedlink acc --dem dem.tif --output acc.tif [--flow dinf]
sedlink slope --dem dem.tif --output slope.tif

# Solver setup for a pour point: inflow cell, reach slope, basin
sedlink prep --dem dem.tif --pour-point "135,66" --snap 5 \
    --output setup.json --basin basin.tif

# Watersheds from pour points (snapped to the channel within 5 cells)
sedlink watershed --dem dem.tif --output basins.tif \
    --pour-points "120,45;300,200" --snap 5

# Sediment routing: distance-decay SDR
sedlink route --dem dem.tif --output sdr.tif --source rusle.tif \
    --flux flux.tif --exponent 0.5 --ref-length 1000

# Sediment routing: IC-based SDR (InVEST-style)
sedlink route --dem dem.tif --output sdr.tif --sdr-model ic \
    --sdr-max 0.8 --ic0 0.5 --k 2
```

As a library:

```rust
use sedlink_core::{ConnectivityIndex, ConnectivityParams, Network};

let net = Network::from_dem(&dem)?;
let ic = ConnectivityIndex::compute(&net, &dem, &ConnectivityParams::default())?;
```

## Conventions

- D8 direction codes follow the SurtGIS convention: 1=E, 2=NE, 3=N, 4=NW,
  5=W, 6=SW, 7=S, 8=SE (counter-clockwise from East); 0 = pit.
- Flow accumulation counts the cell itself (headwater = 1) for both D8
  and D∞.
- DEMs must be axis-aligned with square cells, in a projected CRS with
  metric units.
- Stream cells receive `IC = +clamp`; SedInConnect masks them as NoData —
  keep this in mind when comparing outputs.

## References

- Borselli, L., Cassi, P., & Torri, D. (2008). Prolegomena to sediment and
  flow connectivity in the landscape. *CATENA*, 75(3), 268–277.
- Cavalli, M., Trevisani, S., Comiti, F., & Marchi, L. (2013). Geomorphometric
  assessment of spatial sediment connectivity in small Alpine catchments.
  *Geomorphology*, 188, 31–41.
- Tarboton, D. G. (1997). A new method for the determination of flow
  directions and upslope areas in grid DEMs. *WRR*, 33(2), 309–319.
- Barnes, R., Lehman, C., & Mulla, D. (2014). Priority-flood: An optimal
  depression-filling and watershed-labeling algorithm. *C&G*, 62, 117–127.
- Ferro, V., & Porto, P. (2000). Sediment delivery distributed (SEDD) model.
  *Journal of Hydrologic Engineering*, 5(4), 411–422.
- Vigiak, O., Borselli, L., Newham, L.T.H., McInnes, J., & Roberts, A.M.
  (2012). Comparison of conceptual landscape metrics to define hillslope-scale
  sediment delivery ratio. *Geomorphology*, 138(1), 74–88.

## License

Licensed under either of [Apache License 2.0](LICENSE-APACHE) or
[MIT License](LICENSE-MIT) at your option.
