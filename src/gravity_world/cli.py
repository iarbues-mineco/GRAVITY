from __future__ import annotations

import argparse
from pathlib import Path

from .cartogram import build_latent_country_maps
from .charts import build_all_country_inflow_comparison
from .dyadic import normalize_cepii_controls
from .flows import FLOW_METHOD_COLUMNS, normalize_abel_cohen_flows
from .harmonize import build_country_reference
from .latent_geodesic import estimate_latent_geodesic_ppml
from .latent_geodesic_fe import estimate_latent_geodesic_fe_ppml
from .model import (
    estimate_cepii_gravity_model,
    estimate_cepii_stock_gravity_model,
    estimate_cepii_stock_ppml_model,
    estimate_cepii_stock_unemployment_gravity_model,
    estimate_minimal_gravity_model,
)
from .panel import (
    assemble_cepii_bilateral_panel,
    assemble_cepii_stock_bilateral_panel,
    assemble_minimal_bilateral_panel,
)
from .pipeline import inventory, run_downloads
from .settings import Settings
from .stock import normalize_un_desa_stock
from .transform import assemble_country_year_covariates


MODEL_PREFIX_CHOICES = [
    "minimal_gravity_ols",
    "cepii_gravity_ols",
    "cepii_stock_gravity_ols",
    "cepii_stock_unemployment_gravity_ols",
    "cepii_stock_ppml",
    "latent_geodesic_ppml",
    "latent_geodesic_penalized_ppml",
    "latent_geodesic_fe_ppml",
    "latent_geodesic_fe_penalized_ppml",
]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="World migration gravity data acquisition CLI.")
    subparsers = parser.add_subparsers(dest="command", required=True)

    inventory_parser = subparsers.add_parser("inventory", help="List configured data sources.")
    inventory_parser.set_defaults(command_name="inventory")

    download_parser = subparsers.add_parser("download", help="Download one or more data sources.")
    download_group = download_parser.add_mutually_exclusive_group(required=True)
    download_group.add_argument("--all", action="store_true", help="Download all configured sources.")
    download_group.add_argument("--source", action="append", dest="sources", help="Source id to download.")
    download_parser.set_defaults(command_name="download")

    harmonize_parser = subparsers.add_parser(
        "build-country-reference",
        help="Build the canonical country reference and source override tables.",
    )
    harmonize_parser.set_defaults(command_name="build-country-reference")

    flows_parser = subparsers.add_parser(
        "normalize-abel-cohen",
        help="Normalize the Abel-Cohen bilateral flow file into a standard 5-year dyadic table.",
    )
    flows_parser.add_argument(
        "--preferred-flow",
        choices=FLOW_METHOD_COLUMNS,
        default="da_pb_closed",
        help="Flow estimate column to expose as the default `flow` variable.",
    )
    flows_parser.set_defaults(command_name="normalize-abel-cohen")

    stock_parser = subparsers.add_parser(
        "normalize-un-desa-stock",
        help="Normalize the UN DESA destination-origin migrant stock workbook into a long bilateral stock table.",
    )
    stock_parser.set_defaults(command_name="normalize-un-desa-stock")

    cepii_parser = subparsers.add_parser(
        "normalize-cepii",
        help="Normalize CEPII GeoDist bilateral controls into a harmonized dyadic table.",
    )
    cepii_parser.set_defaults(command_name="normalize-cepii")

    panel_parser = subparsers.add_parser(
        "assemble-minimal-panel",
        help="Assemble a minimal bilateral estimation panel with start-year origin and destination covariates.",
    )
    panel_parser.set_defaults(command_name="assemble-minimal-panel")

    cepii_panel_parser = subparsers.add_parser(
        "assemble-cepii-panel",
        help="Assemble an extended bilateral panel that adds CEPII distance and dyadic controls.",
    )
    cepii_panel_parser.set_defaults(command_name="assemble-cepii-panel")

    cepii_stock_panel_parser = subparsers.add_parser(
        "assemble-cepii-stock-panel",
        help="Assemble the CEPII panel with lagged UN DESA bilateral migrant stock added at the period start year.",
    )
    cepii_stock_panel_parser.set_defaults(command_name="assemble-cepii-stock-panel")

    model_parser = subparsers.add_parser(
        "estimate-minimal-model",
        help="Estimate a baseline log-linear gravity model and write summary outputs.",
    )
    model_parser.set_defaults(command_name="estimate-minimal-model")

    cepii_model_parser = subparsers.add_parser(
        "estimate-cepii-model",
        help="Estimate an extended log-linear gravity model with CEPII dyadic controls.",
    )
    cepii_model_parser.set_defaults(command_name="estimate-cepii-model")

    cepii_stock_model_parser = subparsers.add_parser(
        "estimate-cepii-stock-model",
        help="Estimate the CEPII model augmented with lagged bilateral migrant stock.",
    )
    cepii_stock_model_parser.set_defaults(command_name="estimate-cepii-stock-model")

    cepii_stock_unemployment_model_parser = subparsers.add_parser(
        "estimate-cepii-stock-unemployment-model",
        help="Estimate the CEPII stock model augmented with origin and destination unemployment rates.",
    )
    cepii_stock_unemployment_model_parser.set_defaults(command_name="estimate-cepii-stock-unemployment-model")

    cepii_stock_ppml_model_parser = subparsers.add_parser(
        "estimate-cepii-stock-ppml-model",
        help="Estimate the CEPII stock model with PPML on the full nonnegative flow sample.",
    )
    cepii_stock_ppml_model_parser.set_defaults(command_name="estimate-cepii-stock-ppml-model")

    latent_geodesic_parser = subparsers.add_parser(
        "estimate-latent-geodesic-model",
        help="Estimate an experimental PPML with latent geodesic country coordinates anchored on a user-specified set of fixed countries.",
    )
    latent_geodesic_parser.add_argument(
        "--penalty-weight",
        type=float,
        default=0.0,
        help="Quadratic displacement penalty multiplier toward the real-world reference coordinates. Use 0 for the unpenalized model.",
    )
    latent_geodesic_parser.add_argument(
        "--anchors",
        nargs="+",
        default=["ESP", "FRA"],
        help="ISO3 country codes to keep fixed in the latent map. Provide at least two, for example `--anchors USA CHN`.",
    )
    latent_geodesic_parser.set_defaults(command_name="estimate-latent-geodesic-model")

    latent_geodesic_fe_parser = subparsers.add_parser(
        "estimate-latent-geodesic-fe-model",
        help="Estimate an experimental latent geodesic PPML with static origin and destination fixed effects.",
    )
    latent_geodesic_fe_parser.add_argument(
        "--penalty-weight",
        type=float,
        default=0.0,
        help="Quadratic displacement penalty multiplier toward the real-world reference coordinates. Use 0 for the unpenalized model.",
    )
    latent_geodesic_fe_parser.add_argument(
        "--anchors",
        nargs="+",
        default=["ESP", "FRA"],
        help="ISO3 country codes to keep fixed in the latent map. Provide at least two, for example `--anchors USA CHN`.",
    )
    latent_geodesic_fe_parser.add_argument(
        "--progress-every",
        type=int,
        default=1,
        help="Print one progress line every N outer coordinate iterations.",
    )
    latent_geodesic_fe_parser.set_defaults(command_name="estimate-latent-geodesic-fe-model")

    latent_map_parser = subparsers.add_parser(
        "plot-latent-country-maps",
        help="Render country-surface SVG maps for a latent geodesic model: one with countries fixed in place and one with shifted country shapes.",
    )
    latent_map_parser.add_argument(
        "--model-prefix",
        default="latent_geodesic_ppml",
        choices=["latent_geodesic_ppml", "latent_geodesic_penalized_ppml", "latent_geodesic_fe_ppml", "latent_geodesic_fe_penalized_ppml"],
        help="Latent model output prefix to visualize.",
    )
    latent_map_parser.set_defaults(command_name="plot-latent-country-maps")

    chart_parser = subparsers.add_parser(
        "plot-inflow-comparison",
        help="Build an all-country observed vs fitted inflow comparison chart for a model output.",
    )
    chart_parser.add_argument(
        "--model-prefix",
        default="cepii_gravity_ols",
        choices=MODEL_PREFIX_CHOICES,
        help="Model output prefix to plot.",
    )
    chart_parser.add_argument(
        "--scale",
        default="linear",
        choices=["linear", "log"],
        help="Horizontal axis scaling for the SVG chart.",
    )
    chart_parser.set_defaults(command_name="plot-inflow-comparison")

    assemble_parser = subparsers.add_parser(
        "assemble-country-year",
        help="Merge downloaded World Bank indicators into a country-year covariate table.",
    )
    assemble_parser.set_defaults(command_name="assemble-country-year")

    return parser


def _print_paths(paths: list[Path]) -> None:
    for path in paths:
        print(path)


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    settings = Settings.discover()

    if args.command == "inventory":
        for record in inventory(settings):
            print(f"{record['id']}\t{record['kind']}\t{record['automation']}\t{record['description']}")
        return

    if args.command == "download":
        selected = None if args.all else args.sources
        outputs = run_downloads(settings, selected)
        for source_id, paths in outputs.items():
            print(f"[{source_id}]")
            _print_paths(paths)
        return

    if args.command == "build-country-reference":
        output_paths = build_country_reference(settings)
        _print_paths(output_paths)
        return

    if args.command == "normalize-abel-cohen":
        output_paths = normalize_abel_cohen_flows(settings, preferred_flow=args.preferred_flow)
        _print_paths(output_paths)
        return

    if args.command == "normalize-un-desa-stock":
        output_paths = normalize_un_desa_stock(settings)
        _print_paths(output_paths)
        return

    if args.command == "normalize-cepii":
        output_paths = normalize_cepii_controls(settings)
        _print_paths(output_paths)
        return

    if args.command == "assemble-minimal-panel":
        output_paths = assemble_minimal_bilateral_panel(settings)
        _print_paths(output_paths)
        return

    if args.command == "assemble-cepii-panel":
        output_paths = assemble_cepii_bilateral_panel(settings)
        _print_paths(output_paths)
        return

    if args.command == "assemble-cepii-stock-panel":
        output_paths = assemble_cepii_stock_bilateral_panel(settings)
        _print_paths(output_paths)
        return

    if args.command == "estimate-minimal-model":
        output_paths = estimate_minimal_gravity_model(settings)
        _print_paths(output_paths)
        return

    if args.command == "estimate-cepii-model":
        output_paths = estimate_cepii_gravity_model(settings)
        _print_paths(output_paths)
        return

    if args.command == "estimate-cepii-stock-model":
        output_paths = estimate_cepii_stock_gravity_model(settings)
        _print_paths(output_paths)
        return

    if args.command == "estimate-cepii-stock-unemployment-model":
        output_paths = estimate_cepii_stock_unemployment_gravity_model(settings)
        _print_paths(output_paths)
        return

    if args.command == "estimate-cepii-stock-ppml-model":
        output_paths = estimate_cepii_stock_ppml_model(settings)
        _print_paths(output_paths)
        return

    if args.command == "estimate-latent-geodesic-model":
        output_paths = estimate_latent_geodesic_ppml(settings, penalty_weight=args.penalty_weight, anchors=args.anchors)
        _print_paths(output_paths)
        return

    if args.command == "estimate-latent-geodesic-fe-model":
        output_paths = estimate_latent_geodesic_fe_ppml(settings, penalty_weight=args.penalty_weight, anchors=args.anchors, progress_every=args.progress_every)
        _print_paths(output_paths)
        return

    if args.command == "plot-inflow-comparison":
        output_paths = build_all_country_inflow_comparison(settings, model_prefix=args.model_prefix, scale=args.scale)
        _print_paths(output_paths)
        return

    if args.command == "plot-latent-country-maps":
        output_paths = build_latent_country_maps(settings, model_prefix=args.model_prefix)
        _print_paths(output_paths)
        return

    if args.command == "assemble-country-year":
        output_path = assemble_country_year_covariates(settings.raw_dir, settings.processed_dir)
        print(output_path)
        return

    parser.error(f"Unsupported command: {args.command}")


if __name__ == "__main__":
    main()
