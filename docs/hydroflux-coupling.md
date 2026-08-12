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

1. **Pendiente 0.0101 vs 0.0074 — RECONCILIADA (2026-08-12).** El 0.0074 es
   el `fitted_slope_linear` de
   `examples/huasco_channel/extract_longitudinal_profile.py`: regresión
   lineal sobre 50 celdas (1805.5 m) aguas abajo del gauge Santa Juana,
   snapeado al cauce principal con media-ventana de **150 celdas** (el gauge
   crudo cae en ladera: acc = 34 celdas, z = 655 m), sobre el DEM de cuenca
   completa de paper1. La diferencia con el 0.0101 se descompone en:
   **(a) tramo distinto** — el subset camina otro sector del río, más
   empinado (factor dominante); **(b) estimador** — regresión vs drop/largo
   difieren ~10 % en este tramo dominado por un flat con dos escalones
   (mediana de pendiente por segmento ≈ 3e-6); **(c) snap** — radio 3 vs 150.

   Con los mismos parámetros, sedlink **reproduce la extracción 1D al 0.5 %**:

   ```bash
   sedlink prep --dem factors/06_rio_huasco/hydrology/filled.tif \
       --pour-point "995,1941" --snap 150 --max-reach 1805.5
   # → mean 0.006902 (vs 0.006740), fitted 0.007481 (vs 0.007443),
   #   drop 12.1690 m (exacto); DEM 5068×4975 en ~18 s
   ```

   `ChannelSetup` expone ambos estimadores (`mean_channel_slope`,
   `fitted_channel_slope`) y `--max-reach` acota el tramo, así que la
   constante deja de ser un número sin procedencia: se regenera con un
   comando. Cuál usar es decisión de modelación de Hydroflux — el fitted
   del tramo del gauge (0.0074) o el local del punto de inyección en el
   subset (~0.010).

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
- **v2 (hecho, 2026-08-12)**: peligro combinado — `sedlink hazard` cruza la
  profundidad simulada (`write_depth_geotiff`) con el IC en una matriz
  bivariada 3×3 (clases de profundidad con umbrales suizos 0.5/2.0 m ×
  clases de IC por percentiles del área inundada o umbrales explícitos).
  Clase 0 = seco; 1–9 = `(clase_prof − 1)·3 + clase_IC`. Verificado
  end-to-end con la salida real del evento Huasco 2017 (día 19): 273
  celdas inundadas, 15 en clase 9 (agua profunda + conectividad alta).

  *Nota de interpretación*: sobre una mancha de inundación la mayoría de
  las celdas mojadas son cauce (IC = +clamp), así que los percentiles se
  comprimen hacia arriba; para discriminar fino dentro del canal usar
  `--ic-breaks` explícitos o clasificar contra el SDR en vez del IC.
- **v3 (2032+, si Hydroflux entra en sedimento)**: suministro de sedimento
  (`hillslope_delivery`, `channel_flux`) como condición inicial de volumen
  movilizable, y `FlowNetwork` como grafo de propagación. Recién ahí aplica
  el diseño de "trait compartido" que se había supuesto.

## Nota de higiene

Este acople **no modifica el repo de Hydroflux**: todo vive en sedlink y el
intercambio es por archivos (GeoTIFF + JSON). Hydroflux tiene trabajo sin
commitear (papers, ejemplos de autograd) y su propia convención de crates;
cualquier crate `coupling/` allá es decisión de ese proyecto.
