from __future__ import annotations

import argparse
from pathlib import Path

from .flows import FLOW_METHOD_COLUMNS, normalize_abel_cohen_flows
from .harmonize import build_country_reference
from .panel import assemble_minimal_bilateral_panel
from .pipeline import inventory, run_downloads
from .settings import Settings
from .transform import assemble_country_year_covariates


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

    panel_parser = subparsers.add_parser(
        "assemble-minimal-panel",
        help="Assemble a minimal bilateral estimation panel with start-year origin and destination covariates.",
    )
    panel_parser.set_defaults(command_name="assemble-minimal-panel")

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

    if args.command == "assemble-minimal-panel":
        output_paths = assemble_minimal_bilateral_panel(settings)
        _print_paths(output_paths)
        return

    if args.command == "assemble-country-year":
        output_path = assemble_country_year_covariates(settings.raw_dir, settings.processed_dir)
        print(output_path)
        return

    parser.error(f"Unsupported command: {args.command}")


if __name__ == "__main__":
    main()