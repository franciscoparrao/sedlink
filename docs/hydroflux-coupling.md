# Acoplamiento sedlink ↔ Hydroflux

> Estado: **v1 implementado** (`sedlink prep`). Verificado contra el DEM real
> de Huasco de Hydroflux el 2026-08-07.
> Hydroflux: `~/proyectos/postdoc/hydroflux` (workspace Rust: `autograd`,
> `solver-1d`, `solver-2d`).

## Qué acopla — y qué no

**Alcance real, no el que suponíamos.** El roadmap de Hydroflux
(`outline.md`) coloca **transporte de sedimento en "Años 7+ / 2032+"**. El
solver de hoy es shallow-water puro: `Conserved2D { h, hu, hv }`,
`PointSource { row, col, q_mass }` con `q_mass` en **m³/s de agua**. Por lo
tanto **no** existe (ni corresponde forzar) un acople "routing de sedimento →
propagación granular": eso es trabajo de 2032.

Lo que sí acopla hoy, y es inmediatamente útil, es la **preparación de
terreno** que un run de Hydroflux hace a mano. En
`solver-2d/examples/huasco_2d_event.rs` conviven hoy:

```rust
const ACC_THRESHOLD: f64 = 1_000_000.0;
const SLOPE_MEAN: f64 = 0.0074; // from 1D longitudinal-profile extraction
const INFLOW_ROW: usize = 135;
const INFLOW_COL: usize = 66;
```

Constantes derivadas fuera del código, sin trazabilidad ni forma de
recalcularlas al cambiar de cuenca. `sedlink prep` las produce de forma
reproducible desde el DEM.

## Sustrato común (por qué el acople es barato)

Ambos motores ya comparten base: `surtgis-core::Raster` para I/O GeoTIFF
(Hydroflux lo usa en `mesh_from_geotiff` y en sus ejemplos), `ndarray` 0.16,
edición 2024, convención de grilla row-major con fila 0 al norte. **No hace
falta capa de conversión**: los GeoTIFF que emite sedlink los lee
`mesh_from_geotiff` directamente.

## Interfaz v1: `sedlink prep`

```bash
sedlink prep --dem huasco_subset_dem.tif --pour-point "135,66" \
    --snap 5 --threshold 1000 \
    --output setup.json --acc acc.tif --basin basin.tif
```

Emite JSON plano (sin dependencias) con lo que el run necesita:

| Campo | Reemplaza en el ejemplo | Origen en sedlink |
|---|---|---|
| `inflow_row`, `inflow_col` | `INFLOW_ROW`, `INFLOW_COL` | `Network::snap_to_stream` |
| `mean_channel_slope` | `SLOPE_MEAN` | `NetworkAnalysis::longitudinal_profile` (drop/largo) |
| `channel_length_m`, `channel_drop_m` | — (contexto del tramo) | idem |
| `basin_cells`, `basin_area_m2` | — (dominio) | `Network::watersheds` |
| `stream_cells` | calibra `ACC_THRESHOLD` | máscara `flow_acc ≥ threshold` |

Y opcionalmente los rasters `--acc` (acumulación) y `--basin` (máscara de
cuenca, 1 = dentro) co-registrados con el DEM.

## Verificación sobre Huasco (2026-08-07)

DEM real de Hydroflux (`huasco_subset_dem.tif`, 67×200, celdas de 30 m),
pour point en el inflow hardcodeado:

```
inflow_row 135, inflow_col 66      → coincide con las constantes del ejemplo
mean_channel_slope 0.01007          → vs SLOPE_MEAN = 0.0074 hardcodeado
channel_length_m 4518.3, drop 45.5 m
basin_cells 63, stream_cells 0      → ver caveat
```

**Dos discrepancias que hay que entender antes de usar los valores:**

1. **Pendiente 0.0101 vs 0.0074 (+36 %).** Mismo orden, no idénticas. La del
   ejemplo viene de "1D longitudinal-profile extraction", probablemente sobre
   un tramo distinto (más largo, cuenca completa) o con ajuste por regresión
   en vez de drop/largo entre extremos. **Sin reconciliar**: antes de
   sustituir la constante hay que fijar qué tramo y qué estimador se quiere.

2. **`basin_cells` = 63 y `stream_cells` = 0 no significan que el punto esté
   fuera del canal.** El raster de acumulación de cuenca completa da 8 579 576
   en (135, 66) — 99,7 % del máximo del dominio: el punto **sí** está sobre el
   cauce principal. sedlink reporta 63 porque calcula la acumulación **solo
   dentro del subset**, y el área aportante real de Huasco queda fuera de esa
   ventana de 67×200.

   **Regla general**: cuando el DEM es una ventana recortada de una cuenca
   mayor, la acumulación local subestima el área aportante en los puntos de
   borde. Por eso mismo el ejemplo inyecta un hidrograma medido (DGA Santa
   Juana 03820003) en vez de derivar Q del terreno. Para obtener acumulación
   correcta en un punto de borde, correr `sedlink prep` sobre el **DEM
   completo** de la cuenca y luego recortar.

## Roadmap del acople

- **v1 (hecho)**: preparación de terreno — inflow, pendiente de tramo,
  dominio, umbral de canal.
- **v2 (cuando haya caso de uso)**: post-proceso de peligro combinado —
  cruzar la profundidad simulada (`write_depth_geotiff`) con el IC/SDR de
  sedlink para mapas de "dónde hay agua **y** conectividad de sedimento
  alta". No requiere física nueva en el solver; es álgebra de rasters
  co-registrados. Encaja con la narrativa de peligro acoplado del postdoc.
- **v3 (2032+, si Hydroflux entra en sedimento)**: suministro de sedimento
  (`hillslope_delivery`, `channel_flux`) como condición inicial de volumen
  movilizable, y `FlowNetwork` como grafo de propagación. Recién ahí aplica
  el diseño de "trait compartido" que se había supuesto.

## Nota de higiene

Este acople **no modifica el repo de Hydroflux**: todo vive en sedlink y el
intercambio es por archivos (GeoTIFF + JSON). Hydroflux tiene trabajo sin
commitear (papers, ejemplos de autograd) y su propia convención de crates;
cualquier crate `coupling/` allá es decisión de ese proyecto.
