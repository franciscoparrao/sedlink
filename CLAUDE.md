# sedlink — Conectividad de sedimentos y redes fluviales en Rust

> **Estado:** MVP v0.1 (implementado). Creado 2026-06-10.
> Familia de motores Rust del autor: SurtGIS, Hydroflux, Smelt, Anvil, Cantus, Criterium.
> Doc madre: `~/proyectos/ideas-motores-rust.md` (idea B2).

## Qué es

Motor para índice de conectividad de sedimentos (Borselli IC), análisis de
redes de drenaje y routing de sedimentos a partir de un DEM.

## El gap que llena

SurtGIS hace flow/streams pero no conectividad de sedimentos. El IC vive hoy en
plugins SAGA/scripts dispersos. Une **terrain (SurtGIS)** con **hazards
(Hydroflux)** en el eje sedimentológico.

## Arquitectura

- **sedlink-core**: crate de núcleo con tipos, errores, algoritmos de conectividad
  (IC Borselli) y análisis de redes (Strahler order, magnitude, perfiles longitudinales).
- **sedlink-cli**: CLI con clap para `ic`, `order`, `acc`, `slope`.
- Dependencia de `surtgis-core` para tipos `Raster`, `GeoTransform`, I/O GeoTIFF.

## Alcance MVP (v0.1) — IMPLEMENTADO

- [x] Índice de conectividad IC (Borselli 2008 / Cavalli 2013).
- [x] Componentes upslope/downslope; weighting factor configurable (C-factor, roughness).
- [x] Análisis de red: orden Strahler, magnitud, perfiles longitudinales por tributario.
- [x] D∞ flow routing (Tarboton 1997) como módulo (`dinf.rs`).
- [x] Trait `FlowNetwork` (D8 + D∞): `ConnectivityIndex::compute` es genérico
      y la CLI acepta `--flow d8|dinf` en `ic` y `acc`. Convención unificada:
      flow_acc incluye la propia celda (headwater = 1) en ambos routings.
- [x] (v0.7, unreleased) IC hacia targets arbitrarios (Cavalli 2013 /
      SedInConnect "targets"): `ConnectivityParams::targets` (máscara bool)
      reemplaza la red de streams como destino de D_dn; celdas que no drenan
      a un target → IC NaN. CLI: `sedlink ic --targets targets.tif`.
- [x] (v0.2) Sediment routing con dos modelos de SDR (`SedimentRouting`):
      distancia-decay SDR = exp(−q·L/L_ref) (familia SEDD, `compute`) y
      sigmoide basado en IC SDR = sdr_max/(1+exp((IC0−IC)/k)) (Vigiak 2012 /
      InVEST, `compute_from_ic` — usa `ic.is_stream` para consistencia).
      Entrega de ladera al canal y flujo acumulado en la red. CLI:
      `sedlink route` (--sdr-model distance|ic, --source, --flux,
      --delivery, --flow d8|dinf).

## Validación / paridad numérica

Objetivo: cross-check contra **SedInConnect** y casos publicados de cuencas
alpinas. Estado actual (2026-07-31):

- Tests de **propiedades** (monotonía, conservación de masa, valores
  calculados a mano) en `crates/core/tests/{connectivity,dinf,fill,routing}.rs`.
- La fórmula fue alineada con Borselli 2008 / Cavalli 2013: D_up usa
  promedios upslope W̄·S̄, slope gradient tan θ clampeado a [0.005, 1.0],
  D8 por gradiente (drop/distancia), priority-flood-epsilon (Barnes 2014).
- **Paridad cross-implementation**: `tools/reference_ic.py` (NumPy,
  independiente, desde las ecuaciones) genera el fixture de
  `crates/core/tests/data/` (DEM 40×40 fraccionario + depresión + NoData
  + peso variable); `tests/parity.rs` compara celda a celda (acc 1e-9;
  D_up/D_dn/IC 1e-6). Regenerar con `python3 tools/reference_ic.py`.
- **Pendiente**: paridad contra SedInConnect real (app Windows) con el
  mismo DEM fixture.
- Diferencia conocida: sedlink asigna IC=+clamp en celdas de stream;
  SedInConnect las enmascara como NoData.

## Venue objetivo

**Earth Surface Processes and Landforms (ESPL)** o **Geomorphology**.

## Conexiones con tu ecosistema

- **SurtGIS**: depende de `surtgis-core` para I/O y tipos raster. Usa D8 flow routing
  (convención SurtGIS: 1=E, 2=NE, 3=N, 4=NW, 5=W, 6=SW, 7=S, 8=SE).
- **Hydroflux**: acoplamiento con propagación granular/sedimentos.

## Próximos pasos

1. Cross-validation numérica contra SedInConnect: requiere fixture externo
   (SedInConnect es app Windows; alternativa: script Python de referencia
   independiente). No hay herramienta IC en el MCP Gateway (verificado).
3. Evaluar publicación en ESPL o Geomorphology.
4. Coordinar release de surtgis-core (hoy es path dependency) para poder
   publicar en crates.io / Zenodo.

## Uso

```bash
# Build
cargo build --release

# Compute IC
./target/release/sedlink-cli ic --dem input.tif --output ic.tif --threshold 1000

# Compute Strahler order
./target/release/sedlink-cli order --dem input.tif --output order.tif --threshold 1000

# Compute flow accumulation
./target/release/sedlink-cli acc --dem input.tif --output acc.tif

# Compute slope
./target/release/sedlink-cli slope --dem input.tif --output slope.tif
```

## Tests

```bash
cargo test
```
