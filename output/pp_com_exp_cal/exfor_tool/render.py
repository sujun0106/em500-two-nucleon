"""Traceable table and experimental-data-only figure rendering for QueryOutcome."""

from __future__ import annotations

import csv
import json
import os
from pathlib import Path
import tempfile
from typing import Any
import re


_ENERGY_TO_MEV = {"EV": 1.0e-6, "KEV": 1.0e-3, "MEV": 1.0, "GEV": 1.0e3}
_AREA_TO_BARN = {"B": 1.0, "MB": 1.0e-3, "UB": 1.0e-6, "NB": 1.0e-9, "PB": 1.0e-12}
_DDX_TO_MB = {
    "B/SR/MEV": 1.0e3,
    "MB/SR/MEV": 1.0,
    "UB/SR/MEV": 1.0e-3,
    "MUB/SR/MEV": 1.0e-3,
    "NB/SR/MEV": 1.0e-6,
    # Changing both barn prefixes and energy denominators: 1 b/GeV =
    # 1 mb/MeV.  These are explicit dimensional conversions, not rescaling.
    "B/SR/GEV": 1.0,
    "MB/SR/GEV": 1.0e-3,
    "UB/SR/GEV": 1.0e-6,
    "MUB/SR/GEV": 1.0e-6,
    "NB/SR/GEV": 1.0e-9,
}
_DDX_MOMENTUM_TO_MB = {
    "B/SR/MEVC": 1.0e3,
    "MB/SR/MEVC": 1.0,
    "MUB/SR/MEVC": 1.0e-3,
    "B/SR/GEVC": 1.0,
    "MB/SR/GEVC": 1.0e-3,
    "MUB/SR/GEVC": 1.0e-6,
}
_MOMENTUM_TO_MEV_C = {"MEV/C": 1.0, "GEV/C": 1.0e3}
_OUTGOING_LABELS = {
    "0-NN-1": "neutron",
    "1-H-1": "proton",
    "1-H-2": "deuteron",
    "1-H-3": "triton",
    "2-HE-3": "helium-3",
    "2-HE-4": "alpha",
    "0-G-0": "gamma",
}


class RenderError(ValueError):
    """Raised when a terminal query result cannot be faithfully rendered."""


def flatten_table(outcome: dict[str, Any]) -> dict[str, Any]:
    """Return a lossless long-form representation of every selected EXFOR cell."""

    _require_result(outcome)
    rows: list[dict[str, Any]] = []
    for record in outcome["records"]:
        for table_index, table in enumerate(record["data"]):
            for row_index, cells in enumerate(table["rows"]):
                selected_energies = table.get("incident_energy_mev") or [None] * len(table["rows"])
                energy = selected_energies[row_index]
                for column_index, cell in enumerate(cells):
                    rows.append(
                        {
                            "accession": record["accession"],
                            "subentry": record["subentry"],
                            "reaction": record["reaction"]["raw_code"],
                            "relative_path": record["source"]["relative_path"],
                            "table_index": table_index,
                            "row_index": row_index,
                            "incident_energy_mev": energy,
                            "field": table["headers"][column_index],
                            "unit": table["units"][column_index],
                            "raw": cell["raw"],
                            "value": cell["value"],
                        }
                    )
    return {
        "schema": "exfor-tool/long-table/v1",
        "normalized_query": outcome["normalized_query"],
        "selection": outcome.get("selection", {}),
        "provenance": outcome["provenance"],
        "columns": ["accession", "subentry", "reaction", "relative_path", "table_index", "row_index", "incident_energy_mev", "field", "unit", "raw", "value"],
        "rows": rows,
    }


def write_table(outcome: dict[str, Any], output_stem: str | Path) -> dict[str, str | int]:
    """Write the long-form selected data as UTF-8 JSON and CSV, without rounding values."""

    table = flatten_table(outcome)
    stem = Path(output_stem)
    stem.parent.mkdir(parents=True, exist_ok=True)
    json_path = stem.with_suffix(".json")
    csv_path = stem.with_suffix(".csv")
    json_path.write_text(json.dumps(table, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=table["columns"])
        writer.writeheader()
        writer.writerows(table["rows"])
    return {"json": str(json_path), "csv": str(csv_path), "rows": len(table["rows"])}


def plot_excitation_function(outcome: dict[str, Any], output_path: str | Path) -> dict[str, Any]:
    """Plot one marker/errorbar series per EXFOR subentry; never fit or merge datasets."""

    _require_result(outcome)
    series, excluded = _excitation_series(outcome)
    if not series:
        raise RenderError("No selected table has numeric EN and DATA columns with supported energy and area units.")
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "exfor-tool-matplotlib"))
    import matplotlib.pyplot as plt

    figure, axis = plt.subplots(figsize=(7.2, 4.8), constrained_layout=True)
    for item in series:
        error = item["error"]
        axis.errorbar(item["energy_mev"], item["value_barn"], yerr=error, fmt=item["marker"], linestyle="none", capsize=2.5, label=item["label"])
    query = outcome["normalized_query"]
    reaction = query["reaction"]
    channel = f" ({reaction['channel']})" if reaction["channel"] else ""
    axis.set_title(f"EXFOR experimental data: {reaction['projectile']} + {reaction['target']}{channel}")
    axis.set_xlabel("Incident energy (MeV)")
    axis.set_ylabel("Cross section (barn)")
    axis.grid(True, alpha=0.25)
    if len(series) <= 12:
        axis.legend(fontsize="small", title="EXFOR subentry")
    figure.savefig(output, dpi=180)
    plt.close(figure)
    return {"path": str(output), "series": len(series), "excluded": excluded, "data_only": True, "y_unit": "B"}


def plot_product_isotope_maps(outcome: dict[str, Any], output_directory: str | Path) -> dict[str, Any]:
    """Create one independent nuclide map per accession; never merge campaigns."""

    _require_result(outcome)
    directory = Path(output_directory)
    directory.mkdir(parents=True, exist_ok=True)
    grouped: dict[str, list[tuple[int, int, float]]] = {}
    excluded: list[dict[str, str]] = []
    for record in outcome["records"]:
        product = _product_coordinates(record["reaction"].get("product_code"))
        value = _single_data_value(record)
        if product is None or value is None:
            excluded.append({"subentry": record["subentry"], "reason": "requires one nuclide product and one numeric DATA value"})
            continue
        grouped.setdefault(record["accession"], []).append((*product, value))
    if not grouped:
        raise RenderError("No independently plottable product-isotope points were found.")
    os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "exfor-tool-matplotlib"))
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    paths: list[str] = []
    for accession, points in grouped.items():
        positive = [value for _, _, value in points if value > 0]
        if not positive:
            excluded.append({"subentry": accession, "reason": "all selected DATA values are non-positive"})
            continue
        figure, axis = plt.subplots(figsize=(7.2, 5.2), constrained_layout=True)
        scatter = axis.scatter([neutrons for neutrons, _, _ in points], [protons for _, protons, _ in points], c=[value for _, _, value in points], marker="s", s=190, cmap="viridis", norm=LogNorm(vmin=min(positive), vmax=max(positive) if max(positive) > min(positive) else min(positive) * 10))
        figure.colorbar(scatter, ax=axis, label="Cross section (EXFOR source unit)")
        axis.set_xlabel("Neutron number N")
        axis.set_ylabel("Proton number Z")
        axis.set_title(f"EXFOR product-isotope distribution — accession {accession}")
        axis.grid(True, alpha=0.2)
        path = directory / f"nuclide-map-{accession}.png"
        figure.savefig(path, dpi=180)
        plt.close(figure)
        paths.append(str(path))
    if not paths:
        raise RenderError("No product-isotope map could be rendered from positive selected DATA values.")
    return {"paths": paths, "accessions": len(paths), "excluded": excluded, "data_only": True, "separate_by": "accession"}


def plot_differential_spectrum(outcome: dict[str, Any], output_path: str | Path) -> dict[str, Any]:
    """Plot selected double-differential spectra with bin-centre energy and reported errors."""

    _require_result(outcome)
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    series: list[dict[str, Any]] = []
    excluded: list[dict[str, str]] = []
    for record in outcome["records"]:
        for table_index, table in enumerate(record["data"]):
            points = _differential_points(table)
            if points is None:
                excluded.append({"subentry": record["subentry"], "table": str(table_index), "reason": "requires E-CM-MIN/E-CM-MAX and DATA-CM numeric columns"})
                continue
            series.append({"label": f"{record['accession']}/{record['subentry']}", **points})
    if not series:
        raise RenderError("No selected table has a supported double-differential spectrum layout.")
    os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "exfor-tool-matplotlib"))
    import matplotlib.pyplot as plt

    figure, axis = plt.subplots(figsize=(7.2, 4.8), constrained_layout=True)
    for item in series:
        axis.errorbar(item["energy"], item["data"], xerr=item["xerr"], yerr=item["yerr"], fmt="o", linestyle="none", capsize=2.0, markersize=3.0, label=item["label"])
    query = outcome["normalized_query"]
    axis.set_title(f"EXFOR experimental DDX: {query['reaction']['projectile']} + {query['reaction']['target']}")
    axis.set_xlabel("Outgoing-particle energy (MeV, bin centre)")
    axis.set_ylabel("Double-differential cross section (EXFOR source unit)")
    axis.set_yscale("log")
    axis.grid(True, which="both", alpha=0.25)
    if len(series) <= 12:
        axis.legend(fontsize="small", title="EXFOR subentry")
    figure.savefig(output, dpi=180)
    plt.close(figure)
    return {"path": str(output), "series": len(series), "excluded": excluded, "data_only": True, "x": "outgoing_energy_bin_centre"}


def plot_ddx_by_incident_energy_and_angle(outcome: dict[str, Any], output_directory: str | Path) -> dict[str, Any]:
    """Render one experimental DDX figure per incident energy and reference frame."""

    _require_result(outcome)
    directory = Path(output_directory)
    directory.mkdir(parents=True, exist_ok=True)
    groups: dict[tuple[str, float, str, str, str], list[dict[str, Any]]] = {}
    excluded: list[dict[str, str]] = []
    for record in outcome["records"]:
        for table_index, table in enumerate(record["data"]):
            extracted = _ddx_angle_series(record, table)
            if not extracted:
                excluded.append({"subentry": record["subentry"], "table": str(table_index), "reason": "requires numeric incident energy, ANG, outgoing-energy and DDX columns in the same table"})
                continue
            for item in extracted:
                key = (item["incident_kind"], item["incident_value"], item["incident_unit"], item["reference_frame"], item["outgoing_quantity"])
                groups.setdefault(key, []).append(item)
    if not groups:
        raise RenderError("No selected table has a supported incident-energy/angle DDX layout.")
    os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "exfor-tool-matplotlib"))
    import matplotlib.pyplot as plt

    paths: list[str] = []
    for (incident_kind, incident, incident_unit, reference_frame, outgoing_quantity), series in sorted(groups.items()):
        figure, axis = plt.subplots(figsize=(10.0, 5.4), constrained_layout=True)
        markers = ("o", "s", "^", "D", "v", "P", "X", "<", ">", "h", "*", "8", "p", "H", "d", "|", "_", "+", "x", "1")
        colors = tuple(f"C{index}" for index in range(10)) + tuple(f"tab:{name}" for name in ("blue", "orange", "green", "red", "purple", "brown", "pink", "gray", "olive", "cyan"))
        for ordinal, item in enumerate(sorted(series, key=lambda value: (value["angle_deg"], value["label"]))):
            axis.errorbar(item["energy"], item["data"], xerr=item["xerr"], yerr=item["yerr"], fmt=f"{markers[ordinal % len(markers)]}-", color=colors[ordinal % len(colors)], linewidth=1.0, markersize=2.5, capsize=1.5, label=f"{item['angle_deg']:g}° — {item['label']}")
        outgoing_code = outcome["normalized_query"].get("reaction", {}).get("outgoing_particle")
        outgoing = _OUTGOING_LABELS.get(outgoing_code, "outgoing-particle")
        frame_label = "centre-of-mass frame" if reference_frame == "cm" else "laboratory frame"
        incident_label = f"incident energy {incident:g} {incident_unit}" if incident_kind == "energy" else f"incident momentum {incident:g} {incident_unit}"
        axis.set_title(f"EXFOR experimental {outgoing} DDX — {frame_label}, {incident_label}", fontsize=14)
        if outgoing_quantity == "momentum":
            axis.set_xlabel(f"Outgoing-{outgoing} momentum ({frame_label}, MeV/c)")
            axis.set_ylabel("Double-differential cross section (mb/SR/(MeV/c))")
        else:
            axis.set_xlabel(f"Outgoing-{outgoing} energy ({frame_label}, MeV)")
            axis.set_ylabel("Double-differential cross section (mb/SR/MeV)")
        axis.set_yscale("log")
        axis.grid(True, which="both", alpha=0.25)
        if len(series) <= 16:
            axis.legend(fontsize="x-small", title="Angle — EXFOR subentry", ncol=2)
        prefix = "Ep" if incident_kind == "energy" else "Pp"
        unit_suffix = "MeV" if incident_unit == "MEV" else "GeV_c" if incident_unit == "GEV/C" else incident_unit.replace("/", "_")
        outgoing_suffix = "" if outgoing_quantity == "energy" else "_outgoing_momentum"
        path = directory / f"ddx_{prefix}_{incident:g}_{unit_suffix}_{reference_frame}{outgoing_suffix}.png"
        figure.savefig(path, dpi=180)
        plt.close(figure)
        paths.append(str(path))
    return {
        "paths": paths,
        "incident_energies_mev": sorted({incident for kind, incident, _, _, _ in groups if kind == "energy"}),
        "incident_momenta": [{"value": incident, "unit": unit} for kind, incident, unit, _, _ in sorted(groups) if kind == "momentum"],
        "reference_frames": sorted({frame for _, _, _, frame, _ in groups}),
        "outgoing_quantities": sorted({quantity for _, _, _, _, quantity in groups}),
        "series": sum(len(value) for value in groups.values()),
        "excluded": excluded,
        "data_only": True,
        "ddx_unit": "MB/SR/MEV",
        "separate_by": "incident_condition_and_reference_frame",
    }


def _ddx_angle_series(record: dict[str, Any], table: dict[str, Any]) -> list[dict[str, Any]]:
    headers = table["headers"]
    angle_header = "ANG-CM" if "ANG-CM" in headers else "ANG" if "ANG" in headers else None
    angle_index = headers.index(angle_header) if angle_header else None
    if angle_index is not None and table["units"][angle_index].upper() not in {"DEG", "ADEG"}:
        return []
    common_angle = _common_scalar(table, "ANG-CM", {"DEG", "ADEG"}) if angle_index is None else None
    if common_angle is None and angle_index is None:
        common_angle = _common_scalar(table, "ANG", {"DEG", "ADEG"})
    if angle_index is None and common_angle is None:
        return []
    data_header = "DATA-CM" if "DATA-CM" in headers else "DATA" if "DATA" in headers else None
    if "E-CM-MIN" in headers and "E-CM-MAX" in headers:
        outgoing_headers, outgoing_quantity, x_factor = ("E-CM-MIN", "E-CM-MAX"), "energy", 1.0
    elif "E-CM" in headers:
        outgoing_headers, outgoing_quantity, x_factor = ("E-CM",), "energy", 1.0
    elif "E-DGD" in headers:
        outgoing_headers, outgoing_quantity, x_factor = ("E-DGD",), "energy", 1.0
    elif "E" in headers:
        outgoing_headers, outgoing_quantity, x_factor = ("E",), "energy", 1.0
    elif "MOM-SEC" in headers:
        momentum_unit = table["units"][headers.index("MOM-SEC")].upper()
        x_factor = _MOMENTUM_TO_MEV_C.get(momentum_unit)
        outgoing_headers, outgoing_quantity = ("MOM-SEC",), "momentum"
    else:
        outgoing_headers, outgoing_quantity, x_factor = (), "", None
    if data_header is None or x_factor is None or any(header not in headers for header in outgoing_headers):
        return []
    data_index = headers.index(data_header)
    data_unit = table["units"][data_index]
    factor = (_DDX_TO_MB if outgoing_quantity == "energy" else _DDX_MOMENTUM_TO_MB).get(data_unit.upper())
    if factor is None:
        return []
    outgoing_indices = [headers.index(header) for header in outgoing_headers]
    incident_energies = table.get("incident_energy_mev")
    if incident_energies is not None and len(incident_energies) != len(table["rows"]):
        return []
    momentum = _common_quantity(table, "MOM", {"EV/C", "KEV/C", "MEV/C", "GEV/C"})
    error_index = next((headers.index(name) for name in ("DATA-ERR", "ERR-T", "ERR-S") if name in headers), None)
    collected: dict[tuple[str, float, str, float], dict[str, list[float]]] = {}
    for row_number, row in enumerate(table["rows"]):
        energy = incident_energies[row_number] if incident_energies else None
        if energy is not None:
            incident_kind, incident_value, incident_unit = "energy", energy, "MEV"
        elif momentum is not None:
            incident_kind, (incident_value, incident_unit) = "momentum", momentum
        else:
            continue
        angle = row[angle_index]["value"] if angle_index is not None else common_angle
        value = row[data_index]["value"]
        if angle is None or value is None or value <= 0:
            continue
        coordinates = [row[index]["value"] for index in outgoing_indices]
        if any(value is None for value in coordinates):
            continue
        centre = sum(coordinates) / len(coordinates) * x_factor
        half_width = (coordinates[-1] - coordinates[0]) / 2 * x_factor if len(coordinates) == 2 else 0.0
        error = row[error_index]["value"] if error_index is not None else None
        if error is not None and table["units"][error_index].upper() == "PER-CENT":
            error = value * error / 100
        elif error is not None and table["units"][error_index].upper() != data_unit.upper():
            error = None
        value *= factor
        if error is not None:
            # EXFOR occasionally encodes a lower uncertainty with a negative
            # sign.  The source table remains lossless; a symmetric plotted
            # error-bar length must be non-negative.
            error = abs(error * factor)
        bucket = collected.setdefault((incident_kind, incident_value, incident_unit, angle), {"energy": [], "xerr": [], "data": [], "yerr": []})
        bucket["energy"].append(centre)
        bucket["xerr"].append(half_width)
        bucket["data"].append(value)
        bucket["yerr"].append(error or 0.0)
    reference_frame = _ddx_reference_frame(headers, table.get("common", []))
    return [{"incident_kind": kind, "incident_value": value, "incident_unit": unit, "reference_frame": reference_frame, "outgoing_quantity": outgoing_quantity, "angle_deg": angle, "label": f"{record['accession']}/{record['subentry']}", **points} for (kind, value, unit, angle), points in collected.items()]


def _ddx_reference_frame(headers: list[str], common_tables: list[dict[str, Any]] | None = None) -> str:
    """Infer the frame from EXFOR quantity/kinematic field codes.

    An explicit ``-CM`` field means that the reported energy or cross section is
    in the centre-of-mass frame.  Unqualified DDX fields are laboratory-frame
    values under the EXFOR dictionary convention.  They must never share a
    rendered axis without an explicit transformation.
    """

    common_headers = [header for table in common_tables or [] for header in table["headers"]]
    return "cm" if any("-CM" in header for header in [*headers, *common_headers]) else "lab"


def _common_quantity(table: dict[str, Any], header: str, permitted_units: set[str]) -> tuple[float, str] | None:
    for common in table.get("common", []):
        if header not in common["headers"] or len(common["rows"]) != 1:
            continue
        column = common["headers"].index(header)
        unit = common["units"][column].upper()
        value = common["rows"][0][column]["value"]
        if unit in permitted_units and value is not None:
            return value, unit
    return None


def _common_scalar(table: dict[str, Any], header: str, permitted_units: set[str]) -> float | None:
    for common in table.get("common", []):
        if header not in common["headers"] or len(common["rows"]) != 1:
            continue
        column = common["headers"].index(header)
        if common["units"][column].upper() not in permitted_units:
            continue
        return common["rows"][0][column]["value"]
    return None


def _differential_points(table: dict[str, Any]) -> dict[str, Any] | None:
    headers = table["headers"]
    required = ("E-CM-MIN", "E-CM-MAX", "DATA-CM")
    if any(header not in headers for header in required):
        return None
    lower, upper, data = (headers.index(header) for header in required)
    if any(table["units"][index].upper() != "MEV" for index in (lower, upper)):
        return None
    error = headers.index("ERR-T") if "ERR-T" in headers else None
    percent_error = error is not None and table["units"][error].upper() == "PER-CENT"
    x: list[float] = []
    xerr: list[float] = []
    y: list[float] = []
    yerr: list[float] = []
    for row in table["rows"]:
        lo, hi, value = (row[index]["value"] for index in (lower, upper, data))
        if lo is None or hi is None or value is None or value <= 0:
            continue
        x.append((lo + hi) / 2)
        xerr.append((hi - lo) / 2)
        y.append(value)
        yerr.append(value * row[error]["value"] / 100 if percent_error and row[error]["value"] is not None else 0.0)
    return {"energy": x, "xerr": xerr, "data": y, "yerr": yerr} if x else None


def _product_coordinates(code: Any) -> tuple[int, int] | None:
    match = re.fullmatch(r"(\d+)-[A-Z]{1,2}-(\d+)(?:-[A-Z/+]+)?", str(code or ""))
    if not match:
        return None
    protons, mass = (int(value) for value in match.groups())
    return mass - protons, protons


def _single_data_value(record: dict[str, Any]) -> float | None:
    values: list[float] = []
    for table in record["data"]:
        if "DATA" not in table["headers"]:
            continue
        column = table["headers"].index("DATA")
        values.extend(cell[column]["value"] for cell in table["rows"] if cell[column]["value"] is not None)
    return values[0] if len(values) == 1 else None


def _excitation_series(outcome: dict[str, Any]) -> tuple[list[dict[str, Any]], list[dict[str, str]]]:
    series: list[dict[str, Any]] = []
    excluded: list[dict[str, str]] = []
    markers = ("o", "s", "^", "D", "v", "P", "X", "<", ">", "h")
    for record in outcome["records"]:
        for table_index, table in enumerate(record["data"]):
            try:
                energy_column = table["headers"].index("EN")
                data_column = table["headers"].index("DATA")
                energy_factor = _ENERGY_TO_MEV[table["units"][energy_column].upper()]
                data_factor = _AREA_TO_BARN[table["units"][data_column].upper()]
            except (ValueError, KeyError):
                excluded.append({"subentry": record["subentry"], "table": str(table_index), "reason": "requires EN and DATA columns with supported units"})
                continue
            errors = _error_columns(table, table["units"][data_column])
            points = _series_points(table["rows"], energy_column, data_column, energy_factor, data_factor, errors)
            if not points[0]:
                excluded.append({"subentry": record["subentry"], "table": str(table_index), "reason": "no numeric EN/DATA points"})
                continue
            ordinal = len(series)
            series.append(
                {
                    "label": f"{record['accession']}/{record['subentry']}",
                    "marker": markers[ordinal % len(markers)],
                    "energy_mev": points[0],
                    "value_barn": points[1],
                    "error": points[2],
                }
            )
    return series, excluded


def _error_columns(table: dict[str, Any], data_unit: str) -> tuple[int | None, int | None]:
    headers = table["headers"]
    if "DATA-ERR" in headers:
        index = headers.index("DATA-ERR")
        return (index, index) if table["units"][index].upper() == data_unit.upper() else (None, None)
    lower = headers.index("DATA-ERR-") if "DATA-ERR-" in headers else None
    upper = headers.index("DATA-ERR+") if "DATA-ERR+" in headers else None
    if any(index is not None and table["units"][index].upper() != data_unit.upper() for index in (lower, upper)):
        return None, None
    return lower, upper


def _series_points(rows: list[list[dict[str, Any]]], energy_column: int, data_column: int, energy_factor: float, data_factor: float, errors: tuple[int | None, int | None]) -> tuple[list[float], list[float], list[float] | list[list[float]] | None]:
    energies: list[float] = []
    values: list[float] = []
    lower_errors: list[float] = []
    upper_errors: list[float] = []
    for row in rows:
        energy, value = row[energy_column]["value"], row[data_column]["value"]
        if energy is None or value is None:
            continue
        energies.append(energy * energy_factor)
        values.append(value * data_factor)
        lower, upper = errors
        lower_errors.append((row[lower]["value"] or 0.0) * data_factor if lower is not None else 0.0)
        upper_errors.append((row[upper]["value"] or 0.0) * data_factor if upper is not None else 0.0)
    if errors == (None, None):
        return energies, values, None
    return energies, values, lower_errors if errors[0] == errors[1] else [lower_errors, upper_errors]


def _require_result(outcome: Any) -> None:
    if not isinstance(outcome, dict) or outcome.get("status") != "result":
        status = outcome.get("status") if isinstance(outcome, dict) else None
        raise RenderError(f"Rendering requires a QueryOutcome with status 'result', got {status!r}.")
