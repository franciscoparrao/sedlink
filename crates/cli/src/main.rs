//! sedlink-cli — CLI for sediment connectivity and fluvial network analysis.
//!
//! Entry point for computing the Index of Connectivity (Borselli 2008),
//! Strahler stream order, and related analyses from a DEM.

use std::path::PathBuf;

use clap::{Parser, Subcommand, ValueEnum};
use sedlink_core::{
    ConnectivityIndex, ConnectivityParams, DinfNetwork, FlowNetwork, IcSdrParams, Network,
    NetworkAnalysis, RoutingParams, SedimentRouting,
};
use surtgis_core::raster::Raster;

/// Flow routing model for accumulation and downstream tracing.
#[derive(Clone, Copy, Debug, ValueEnum)]
enum FlowModel {
    /// D8 single flow direction (O'Callaghan & Mark 1984).
    D8,
    /// D-infinity continuous flow direction (Tarboton 1997).
    Dinf,
}

/// Sediment delivery ratio model.
#[derive(Clone, Copy, Debug, ValueEnum)]
enum SdrModel {
    /// Distance decay: SDR = exp(-q L / L_ref) (SEDD family).
    Distance,
    /// IC-based sigmoid: SDR = sdr_max / (1 + exp((ic0 - IC)/k))
    /// (Vigiak 2012 / InVEST).
    Ic,
}

#[derive(Parser)]
#[command(
    name = "sedlink",
    version,
    about = "Sediment connectivity index and fluvial network analysis"
)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Compute the Index of Connectivity (Borselli 2008)
    Ic {
        /// Input DEM (GeoTIFF)
        #[arg(short, long)]
        dem: PathBuf,
        /// Output IC raster (GeoTIFF)
        #[arg(short, long)]
        output: PathBuf,
        /// Stream threshold (cell count)
        #[arg(short, long, default_value_t = 1000.0)]
        threshold: f64,
        /// Optional weighting raster (GeoTIFF)
        #[arg(short, long)]
        weight: Option<PathBuf>,
        /// Clamp IC to [-value, +value]
        #[arg(short, long, default_value_t = 10.0)]
        clamp: f64,
        /// Flow routing model
        #[arg(long, value_enum, default_value_t = FlowModel::D8)]
        flow: FlowModel,
    },
    /// Compute Strahler stream order
    Order {
        /// Input DEM (GeoTIFF)
        #[arg(short, long)]
        dem: PathBuf,
        /// Output order raster (GeoTIFF)
        #[arg(short, long)]
        output: PathBuf,
        /// Stream threshold (cell count)
        #[arg(short, long, default_value_t = 1000.0)]
        threshold: f64,
    },
    /// Compute flow accumulation
    Acc {
        /// Input DEM (GeoTIFF)
        #[arg(short, long)]
        dem: PathBuf,
        /// Output accumulation raster (GeoTIFF)
        #[arg(short, long)]
        output: PathBuf,
        /// Flow routing model
        #[arg(long, value_enum, default_value_t = FlowModel::D8)]
        flow: FlowModel,
    },
    /// Compute slope (radians)
    Slope {
        /// Input DEM (GeoTIFF)
        #[arg(short, long)]
        dem: PathBuf,
        /// Output slope raster (GeoTIFF)
        #[arg(short, long)]
        output: PathBuf,
    },
    /// Delineate watersheds from pour points (D8)
    Watershed {
        /// Input DEM (GeoTIFF)
        #[arg(short, long)]
        dem: PathBuf,
        /// Output basin-label raster (GeoTIFF; 1-based labels, 0 = none)
        #[arg(short, long)]
        output: PathBuf,
        /// Pour points as "row,col" pairs separated by ';' (e.g. "120,45;300,200")
        #[arg(short, long)]
        pour_points: String,
        /// Snap each pour point to the max-accumulation cell within this
        /// Chebyshev radius (cells); 0 = no snapping
        #[arg(long, default_value_t = 0)]
        snap: usize,
    },
    /// Route sediment to the channel network (distance-decay SDR)
    Route {
        /// Input DEM (GeoTIFF)
        #[arg(short, long)]
        dem: PathBuf,
        /// Output SDR raster (GeoTIFF)
        #[arg(short, long)]
        output: PathBuf,
        /// Optional sediment source raster (e.g. RUSLE soil loss, GeoTIFF).
        /// Without it a unit source of 1.0 per cell is used.
        #[arg(short, long)]
        source: Option<PathBuf>,
        /// Stream threshold (cell count)
        #[arg(short, long, default_value_t = 1000.0)]
        threshold: f64,
        /// SDR model
        #[arg(long, value_enum, default_value_t = SdrModel::Distance)]
        sdr_model: SdrModel,
        /// [distance] SDR decay exponent q in SDR = exp(-q L / L_ref)
        #[arg(long, default_value_t = 0.5)]
        exponent: f64,
        /// [distance] Reference length L_ref (m)
        #[arg(long, default_value_t = 1000.0)]
        ref_length: f64,
        /// [distance] Maximum routing distance (m); farther cells deliver nothing
        #[arg(long, default_value_t = 100_000.0)]
        max_distance: f64,
        /// [ic] Maximum SDR of a hillslope cell
        #[arg(long, default_value_t = 0.8)]
        sdr_max: f64,
        /// [ic] Sigmoid midpoint IC0
        #[arg(long, default_value_t = 0.5)]
        ic0: f64,
        /// [ic] Sigmoid width k
        #[arg(long, default_value_t = 2.0)]
        k: f64,
        /// [ic] Optional weighting raster for the IC (GeoTIFF)
        #[arg(short, long)]
        weight: Option<PathBuf>,
        /// Optional output for channel sediment flux (GeoTIFF)
        #[arg(long)]
        flux: Option<PathBuf>,
        /// Optional output for hillslope delivery at channel cells (GeoTIFF)
        #[arg(long)]
        delivery: Option<PathBuf>,
        /// Flow routing model
        #[arg(long, value_enum, default_value_t = FlowModel::D8)]
        flow: FlowModel,
    },
}

fn main() -> anyhow::Result<()> {
    tracing_subscriber::fmt::init();
    let cli = Cli::parse();

    match cli.command {
        Commands::Ic {
            dem,
            output,
            threshold,
            weight,
            clamp,
            flow,
        } => {
            let dem_raster = load_dem(&dem)?;

            let weight_arr = if let Some(wpath) = weight {
                let w = load_raster_f64(&wpath)?;
                Some(w.data().clone())
            } else {
                None
            };

            let params = ConnectivityParams {
                stream_threshold: threshold,
                clamp,
                weight: sedlink_core::WeightingFactor {
                    raster: weight_arr,
                    min_value: 0.001,
                },
                ..Default::default()
            };

            let ic_raster = match flow {
                FlowModel::D8 => {
                    let net = Network::from_dem(&dem_raster)?;
                    let ic = ConnectivityIndex::compute(&net, &dem_raster, &params)?;
                    ic.ic_raster(&net)
                }
                FlowModel::Dinf => {
                    let net = DinfNetwork::from_dem(&dem_raster)?;
                    let ic = ConnectivityIndex::compute(&net, &dem_raster, &params)?;
                    ic.ic_raster(&net)
                }
            };
            save_raster(&ic_raster, &output)?;
            eprintln!("IC raster ({flow:?}) saved to {}", output.display());
        }
        Commands::Order {
            dem,
            output,
            threshold,
        } => {
            let dem_raster = load_dem(&dem)?;
            let net = Network::from_dem(&dem_raster)?;
            let analysis = NetworkAnalysis::new(&net);
            let order = analysis.strahler_order(threshold);

            let mut out = Raster::from_array(
                ndarray::Array2::from_shape_vec((net.rows(), net.cols()), order.values).unwrap(),
            );
            out.set_transform(*net.transform());
            save_raster(&out, &output)?;
            eprintln!(
                "Strahler order (max={}) saved to {}",
                order.max_order,
                output.display()
            );
        }
        Commands::Acc { dem, output, flow } => {
            let dem_raster = load_dem(&dem)?;
            let acc = match flow {
                FlowModel::D8 => Network::from_dem(&dem_raster)?.flow_acc_raster()?,
                FlowModel::Dinf => DinfNetwork::from_dem(&dem_raster)?.flow_acc_raster()?,
            };
            save_raster(&acc, &output)?;
            eprintln!("Flow accumulation ({flow:?}) saved to {}", output.display());
        }
        Commands::Slope { dem, output } => {
            let dem_raster = load_dem(&dem)?;
            let net = Network::from_dem(&dem_raster)?;
            let slope = net.slope_raster()?;
            save_raster(&slope, &output)?;
            eprintln!("Slope saved to {}", output.display());
        }
        Commands::Watershed {
            dem,
            output,
            pour_points,
            snap,
        } => {
            let dem_raster = load_dem(&dem)?;
            let net = Network::from_dem(&dem_raster)?;

            let mut points = Vec::new();
            for part in pour_points.split(';') {
                let (r, c) = part
                    .split_once(',')
                    .ok_or_else(|| anyhow::anyhow!("pour point '{part}' is not 'row,col'"))?;
                let (r, c): (usize, usize) = (r.trim().parse()?, c.trim().parse()?);
                if r >= net.rows() || c >= net.cols() {
                    anyhow::bail!(
                        "pour point ({r}, {c}) outside the {}x{} grid",
                        net.rows(),
                        net.cols()
                    );
                }
                points.push(if snap > 0 {
                    net.snap_to_stream(r, c, snap)
                } else {
                    r * net.cols() + c
                });
            }

            let labels = net.watersheds(&points);
            let data: Vec<f64> = labels.iter().map(|&l| f64::from(l)).collect();
            let mut out = Raster::from_vec(data, net.rows(), net.cols())?;
            out.set_transform(*net.transform());
            save_raster(&out, &output)?;
            eprintln!(
                "Watersheds for {} pour point(s) saved to {}",
                points.len(),
                output.display()
            );
        }
        Commands::Route {
            dem,
            output,
            source,
            threshold,
            sdr_model,
            exponent,
            ref_length,
            max_distance,
            sdr_max,
            ic0,
            k,
            weight,
            flux,
            delivery,
            flow,
        } => {
            let dem_raster = load_dem(&dem)?;
            let source_arr = if let Some(spath) = source {
                Some(load_raster_f64(&spath)?.data().clone())
            } else {
                None
            };
            let weight_arr = if let Some(wpath) = weight {
                Some(load_raster_f64(&wpath)?.data().clone())
            } else {
                None
            };
            let dist_params = RoutingParams {
                sdr_exponent: exponent,
                reference_length: ref_length,
                max_distance,
            };
            let ic_params = IcSdrParams { sdr_max, ic0, k };

            let cfg = RouteConfig {
                source: source_arr,
                weight: weight_arr,
                threshold,
                model: sdr_model,
                dist_params,
                ic_params,
            };
            match flow {
                FlowModel::D8 => {
                    let net = Network::from_dem(&dem_raster)?;
                    run_route(
                        &net,
                        &dem_raster,
                        &cfg,
                        &output,
                        flux.as_deref(),
                        delivery.as_deref(),
                    )?;
                }
                FlowModel::Dinf => {
                    let net = DinfNetwork::from_dem(&dem_raster)?;
                    run_route(
                        &net,
                        &dem_raster,
                        &cfg,
                        &output,
                        flux.as_deref(),
                        delivery.as_deref(),
                    )?;
                }
            }
            eprintln!(
                "SDR raster ({flow:?}, {sdr_model:?}) saved to {}",
                output.display()
            );
        }
    }

    Ok(())
}

/// Inputs shared by both SDR models of the `route` command.
struct RouteConfig {
    source: Option<ndarray::Array2<f64>>,
    weight: Option<ndarray::Array2<f64>>,
    threshold: f64,
    model: SdrModel,
    dist_params: RoutingParams,
    ic_params: IcSdrParams,
}

/// Compute sediment routing on a flow network and write the requested rasters.
fn run_route<F: FlowNetwork>(
    net: &F,
    dem: &Raster<f32>,
    cfg: &RouteConfig,
    sdr_path: &std::path::Path,
    flux_path: Option<&std::path::Path>,
    delivery_path: Option<&std::path::Path>,
) -> anyhow::Result<()> {
    let routing = match cfg.model {
        SdrModel::Distance => {
            SedimentRouting::compute(net, cfg.source.as_ref(), cfg.threshold, &cfg.dist_params)?
        }
        SdrModel::Ic => {
            let cparams = ConnectivityParams {
                stream_threshold: cfg.threshold,
                weight: sedlink_core::WeightingFactor {
                    raster: cfg.weight.clone(),
                    min_value: 0.001,
                },
                ..Default::default()
            };
            let ic = ConnectivityIndex::compute(net, dem, &cparams)?;
            SedimentRouting::compute_from_ic(net, &ic, cfg.source.as_ref(), &cfg.ic_params)?
        }
    };

    let write = |data: &ndarray::Array2<f64>, path: &std::path::Path| -> anyhow::Result<()> {
        let mut r = Raster::from_array(data.clone());
        r.set_transform(*net.transform());
        r.set_nodata(Some(f64::NAN));
        surtgis_core::io::write_geotiff(&r, path, None)?;
        Ok(())
    };

    write(&routing.sdr, sdr_path)?;
    if let Some(p) = flux_path {
        write(&routing.channel_flux, p)?;
        eprintln!("Channel flux saved to {}", p.display());
    }
    if let Some(p) = delivery_path {
        write(&routing.hillslope_delivery, p)?;
        eprintln!("Hillslope delivery saved to {}", p.display());
    }
    Ok(())
}

fn load_dem(path: &PathBuf) -> anyhow::Result<Raster<f32>> {
    let raster = surtgis_core::io::read_geotiff(path, None)?;
    Ok(raster)
}

fn load_raster_f64(path: &PathBuf) -> anyhow::Result<Raster<f64>> {
    let raster = surtgis_core::io::read_geotiff(path, None)?;
    Ok(raster)
}

fn save_raster<T>(raster: &Raster<T>, path: &PathBuf) -> anyhow::Result<()>
where
    T: surtgis_core::raster::RasterElement + surtgis_core::io::NativeGraySample,
    [T]: tiff::encoder::TiffValue,
{
    surtgis_core::io::write_geotiff(raster, path, None)?;
    Ok(())
}
