# Acoplamiento sedlink → Hydroflux (diseño, v0.3)

> Estado: **diseño**. Hydroflux (solver de hazards / propagación granular de la
> familia de motores) aún no existe como crate; este documento fija el contrato
> que sedlink expone para que el acoplamiento sea directo cuando se cree.
> Actualizado 2026-08-06.

## Rol de cada motor

- **sedlink** responde *cuánto sedimento llega al canal y dónde*: IC, SDR,
  entrega de ladera (`hillslope_delivery`) y flujo acumulado en la red
  (`channel_flux`), sobre redes D8/D∞ (`FlowNetwork`).
- **Hydroflux** responderá *qué hace ese sedimento durante un evento*:
  propagación granular/hiperconcentrada (debris flows, lahares) con física de
  onda/reología, paso de tiempo explícito.

La frontera natural: sedlink es **pre-evento y estacionario** (suministro,
susceptibilidad); Hydroflux es **evento y dinámico** (propagación). El
acoplamiento es one-way en v1 (sedlink → Hydroflux), sin retroalimentación de
erosión del evento hacia el IC.

## Contrato de datos (lo que sedlink ya provee hoy)

Todas las salidas son rasters `f64` co-registrados con el DEM de entrada
(misma grilla, `GeoTransform` compartido, NaN = NoData), más la topología:

| Insumo Hydroflux | API sedlink | Semántica |
|---|---|---|
| Suministro de sedimento por celda de canal (condición inicial de volumen movilizable) | `SedimentRouting::hillslope_delivery` | masa (unidades de la fuente) que ENTRA al canal en cada celda |
| Carga acumulada aguas arriba (para hidrogramas de sedimento en el punto de inicio) | `SedimentRouting::channel_flux` | yield total upstream por celda de canal |
| Susceptibilidad de ladera (dónde inicializar aportes laterales) | `ConnectivityIndex::ic` | IC de Borselli/Cavalli |
| Red de propagación (grafo por el que corre el flujo) | `Network::{flow_dir_slice, downstream, upstream}` + `FlowNetwork::downstream_step` | D8: 1 receptor; D∞: receptor primario + fracciones |
| Jerarquía para segmentar tramos | `NetworkAnalysis::strahler_order` | orden por celda |
| Dominio del evento | `Network::watersheds(&[outlet])` | máscara de la cuenca que drena al punto de inicio |
| Geometría | `FlowNetwork::{cellsize, transform}` | resolución y georreferencia |

## Interfaz propuesta (lado Hydroflux, cuando exista)

```rust
// hydroflux-core
pub struct SedimentSupply<'a> {
    /// Masa movilizable por celda (canal), co-registrada con el DEM.
    pub supply: &'a Array2<f64>,
    /// Máscara de dominio (cuenca del evento).
    pub domain: &'a [u32],
    /// Red de flujo para la fase de transporte en canal.
    pub network: &'a dyn FlowNetwork,   // re-export de sedlink_core
}
```

Decisiones:
1. **`FlowNetwork` es el punto de acople**, no `Network` concreto — Hydroflux
   consume el trait (re-exportado o duplicado mínimo) y funciona con D8 y D∞.
2. **Sin serialización intermedia** en proceso: ambos crates comparten
   `surtgis-core::Raster` y `ndarray`, el paso es por referencia. Para acople
   entre procesos (CLI a CLI), los GeoTIFF que ya emite `sedlink route`
   (`--flux`, `--delivery`) + `sedlink watershed` son el formato de intercambio.
3. **Unidades**: sedlink no impone unidades (siguen a la fuente). El contrato
   exige declararlas en el metadato del GeoTIFF (`units=` tag) — pendiente
   menor en `save_raster`.

## Qué falta (en orden, cuando Hydroflux parta)

1. Crear `hydroflux` con `sedlink-core = { version = "0.3" }` como dep.
2. Decidir si `FlowNetwork` se re-exporta o si Hydroflux define su propio trait
   con blanket impl (evita dependencia dura si Hydroflux quiere otras redes).
3. Caso de validación conjunto: cuenca alpina con IC + evento documentado
   (candidato natural para el paper ESPL: "del índice estático al evento").
