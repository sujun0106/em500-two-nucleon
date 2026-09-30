"""Command-line entry point for the parser vertical slice."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

from .index import build_index
from .parser import parse_file
from .query import execute_query
from .render import RenderError, plot_ddx_by_incident_energy_and_angle, plot_differential_spectrum, plot_excitation_function, plot_product_isotope_maps, write_table


def main() -> int:
    parser = argparse.ArgumentParser(description="Parse one EXFOR .x4 entry as JSON.")
    commands = parser.add_subparsers(dest="command")
    parse_command = commands.add_parser("parse", help="parse one .x4 entry as JSON")
    parse_command.add_argument("path", type=Path, help="path to an EXFOR .x4 entry")
    parse_command.add_argument("--subentry", help="emit one subentry by its eight-character ID")
    index_command = commands.add_parser("build-index", help="build a local SQLite metadata index")
    index_command.add_argument("source", type=Path, help="path to the exforall directory")
    index_command.add_argument("destination", type=Path, help="SQLite index destination")
    index_command.add_argument("--batch-size", type=int, default=500)
    query_command = commands.add_parser("query", help="execute a QuerySpec JSON file")
    query_command.add_argument("spec", type=Path, help="path to QuerySpec JSON")
    query_command.add_argument("--index", required=True, type=Path, help="SQLite metadata index")
    query_command.add_argument("--source", required=True, type=Path, help="path to the exforall directory")
    query_command.add_argument("--output", type=Path, help="optional file to receive the QueryOutcome JSON")
    table_command = commands.add_parser("render-table", help="write CSV and JSON from a QueryOutcome JSON file")
    table_command.add_argument("outcome", type=Path, help="path to QueryOutcome JSON")
    table_command.add_argument("--output", required=True, type=Path, help="output filename stem (without extension)")
    plot_command = commands.add_parser("plot-excitation", help="plot experimental EXFOR excitation-function points")
    plot_command.add_argument("outcome", type=Path, help="path to QueryOutcome JSON")
    plot_command.add_argument("--output", required=True, type=Path, help="PNG output path")
    map_command = commands.add_parser("plot-nuclide-map", help="plot separate product-isotope maps per EXFOR accession")
    map_command.add_argument("outcome", type=Path, help="path to QueryOutcome JSON")
    map_command.add_argument("--output", required=True, type=Path, help="output directory")
    differential_command = commands.add_parser("plot-differential", help="plot selected double-differential EXFOR spectra")
    differential_command.add_argument("outcome", type=Path, help="path to QueryOutcome JSON")
    differential_command.add_argument("--output", required=True, type=Path, help="PNG output path")
    ddx_grid_command = commands.add_parser("plot-ddx-grid", help="plot DDX grouped by incident energy and emission angle")
    ddx_grid_command.add_argument("outcome", type=Path, help="path to QueryOutcome JSON")
    ddx_grid_command.add_argument("--output", required=True, type=Path, help="output directory")

    arguments = sys.argv[1:]
    if arguments and arguments[0] not in {"parse", "build-index", "query", "render-table", "plot-excitation", "plot-nuclide-map", "plot-differential", "plot-ddx-grid", "-h", "--help"}:
        arguments.insert(0, "parse")
    args = parser.parse_args(arguments)

    if args.command == "build-index":
        def report(statistics: dict[str, int]) -> None:
            print(json.dumps({"progress": statistics}, sort_keys=True), flush=True)

        statistics = build_index(
            args.source,
            args.destination,
            batch_size=args.batch_size,
            progress=report,
        )
        print(json.dumps({"complete": statistics}, sort_keys=True))
        return 0

    if args.command == "query":
        try:
            spec = json.loads(args.spec.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            parser.error(f"cannot read QuerySpec JSON: {exc}")
        outcome = execute_query(spec, index_path=args.index, source_root=args.source)
        serialized = json.dumps(outcome, indent=2, sort_keys=True) + "\n"
        if args.output:
            args.output.parent.mkdir(parents=True, exist_ok=True)
            args.output.write_text(serialized, encoding="utf-8")
        print(serialized, end="")
        return 0

    if args.command in {"render-table", "plot-excitation", "plot-nuclide-map", "plot-differential", "plot-ddx-grid"}:
        try:
            outcome = json.loads(args.outcome.read_text(encoding="utf-8"))
            rendered = write_table(outcome, args.output) if args.command == "render-table" else plot_excitation_function(outcome, args.output) if args.command == "plot-excitation" else plot_product_isotope_maps(outcome, args.output) if args.command == "plot-nuclide-map" else plot_differential_spectrum(outcome, args.output) if args.command == "plot-differential" else plot_ddx_by_incident_energy_and_angle(outcome, args.output)
        except (OSError, json.JSONDecodeError, RenderError) as exc:
            parser.error(str(exc))
        print(json.dumps(rendered, indent=2, sort_keys=True))
        return 0

    entry = parse_file(args.path)
    if args.subentry:
        for subentry in entry["subentries"]:
            if subentry["id"] == args.subentry:
                print(json.dumps(subentry, indent=2, sort_keys=True))
                return 0
        parser.error(f"subentry {args.subentry!r} was not found")

    print(json.dumps(entry, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
