"""Deterministic QuerySpec execution against a local EXFOR index and source tree."""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

from .index import index_metadata, search_reactions
from .parser import ExforParseError, parse_file


_ELEMENTS = (
    "N H HE LI BE B C N O F NE NA MG AL SI P S CL AR K CA SC TI V CR MN FE CO NI CU ZN GA GE AS SE BR KR "
    "RB SR Y ZR NB MO TC RU RH PD AG CD IN SN SB TE I XE CS BA LA CE PR ND PM SM EU GD TB DY HO ER TM YB "
    "LU HF TA W RE OS IR PT AU HG TL PB BI PO AT RN FR RA AC TH PA U NP PU AM CM BK CF ES FM MD NO LR RF DB SG BH HS MT DS RG CN NH FL MC LV TS OG"
).split()
_SYMBOL_TO_Z = {symbol: number for number, symbol in enumerate(_ELEMENTS) if number}
_PROJECTILES = {"N": "N", "NEUTRON": "N", "P": "P", "PROTON": "P", "D": "D", "DEUTERON": "D", "T": "T", "TRITON": "T", "HE3": "HE3", "A": "A", "ALPHA": "A", "G": "G", "GAMMA": "G", "PHOTON": "G"}
_OUTGOING_PARTICLES = {"N": "0-NN-1", "NEUTRON": "0-NN-1", "P": "1-H-1", "PROTON": "1-H-1", "D": "1-H-2", "DEUTERON": "1-H-2", "T": "1-H-3", "TRITON": "1-H-3", "HE3": "2-HE-3", "A": "2-HE-4", "ALPHA": "2-HE-4", "G": "0-G-0", "GAMMA": "0-G-0", "PHOTON": "0-G-0"}
_OBSERVABLES = {
    "cross_section": "SIG",
    "production_cross_section": "SIG",
    "differential_cross_section": "DA",
    "double_differential_cross_section": "DA/DE",
    "product_isotope_distribution": "SIG",
}
_ENERGY_FACTORS = {"EV": 1.0e-6, "KEV": 1.0e-3, "MEV": 1.0, "GEV": 1.0e3}
_KINEMATIC_HEADERS = {"angle": ("ANG", "ANG-CM"), "outgoing_energy": ("E", "E-OUT", "E-SEC")}
_ANGLE_FACTORS = {"DEG": 1.0, "ADEG": 1.0, "RAD": 180.0 / math.pi}


def execute_query(spec: dict[str, Any], *, index_path: str | Path, source_root: str | Path) -> dict[str, Any]:
    """Validate and execute a QuerySpec without guessing missing physics intent."""

    normalized_or_outcome = _normalize_spec(spec)
    if "status" in normalized_or_outcome:
        return normalized_or_outcome
    normalized = normalized_or_outcome
    provenance = _provenance(index_path, source_root)
    candidates = search_reactions(index_path, **normalized["criteria"])
    if normalized["yield_code"]:
        candidates = [candidate for candidate in candidates if _yield_code(candidate["raw_code"]) == normalized["yield_code"]]
    if not candidates:
        return {"status": "no_match", "normalized_query": normalized["public"], "search_scope": normalized["criteria"], "provenance": provenance}

    records: list[dict[str, Any]] = []
    skipped: list[dict[str, str]] = []
    for candidate in candidates:
        try:
            subentry = _read_subentry(source_root, candidate)
        except (OSError, UnicodeDecodeError, ExforParseError, ValueError) as exc:
            skipped.append({"subentry": candidate["subentry_id"], "reason": str(exc)})
            continue
        data, energy_status = _select_data(subentry, normalized["energy"], normalized["kinematics"])
        if energy_status == "no_energy_column":
            skipped.append({"subentry": candidate["subentry_id"], "reason": "incident energy was requested but EXFOR has no EN column"})
            continue
        if data:
            records.append(_record(candidate, subentry, data))

    if not records:
        return {
            "status": "no_match",
            "normalized_query": normalized["public"],
            "search_scope": normalized["criteria"],
            "provenance": provenance,
            "details": {"matched_subentries": len(candidates), "excluded": skipped},
        }
    return {
        "status": "result",
        "normalized_query": normalized["public"],
        "selection": normalized["selection"],
        "records": records,
        "table_data": records,
        "plot_data": None,
        "provenance": provenance,
        "details": {"matched_subentries": len(candidates), "returned_subentries": len(records), "excluded": skipped},
    }


def _normalize_spec(spec: Any) -> dict[str, Any]:
    if not isinstance(spec, dict):
        return _invalid("QuerySpec must be a JSON object.")
    reaction = spec.get("reaction")
    if not isinstance(reaction, dict):
        return _clarify(["reaction.projectile", "reaction.target"], "A reaction object is required.")
    missing = [f"reaction.{field}" for field in ("projectile", "target") if not reaction.get(field)]
    observable = spec.get("observable")
    if not observable:
        missing.append("observable")
    if missing:
        return _clarify(missing, "Specify the reaction participants and requested observable.")
    projectile = _PROJECTILES.get(str(reaction["projectile"]).strip().upper())
    target = _nuclide_code(reaction["target"])
    if projectile is None:
        return _unsupported("projectile", str(reaction["projectile"]), "This projectile alias is not implemented.")
    if target is None:
        return _invalid("reaction.target must identify a nuclide, for example '56Fe', 'Fe-56', or '26-FE-56'.")
    quantity = _OBSERVABLES.get(str(observable).strip().lower())
    if quantity is None:
        return _unsupported("observable", str(observable), "Use a supported observable name or extend the tool mapping.")
    process = _process_code(reaction.get("channel"), projectile)
    if reaction.get("channel") is not None and process is None:
        return _invalid("reaction.channel must be an EXFOR process code such as '2N' or a channel such as '(n,2n)'.")
    product_value = reaction.get("product")
    product = _nuclide_code(product_value) if product_value else None
    if product_value and product is None:
        return _invalid("reaction.product must identify a nuclide, for example '99Mo'.")
    outgoing_value = reaction.get("outgoing_particle")
    outgoing = _OUTGOING_PARTICLES.get(str(outgoing_value).strip().upper()) if outgoing_value else None
    if outgoing_value and outgoing is None:
        return _unsupported("reaction.outgoing_particle", str(outgoing_value), "Supported outgoing-particle aliases include neutron, proton, deuteron, triton, helium-3, alpha, and gamma.")
    if product and outgoing and product != outgoing:
        return _invalid("reaction.product and reaction.outgoing_particle select different EXFOR products.")
    if outgoing:
        product = outgoing
    is_distribution = str(observable).strip().lower() == "product_isotope_distribution"
    is_double_differential = str(observable).strip().lower() == "double_differential_cross_section"
    if is_double_differential and not product:
        return _clarify(["reaction.outgoing_particle"], "A double-differential query requires the outgoing particle, for example neutron or proton; the tool will not mix particle species.")
    if is_distribution and not reaction.get("channel"):
        return _clarify(["reaction.channel"], "A product-isotope distribution requires an explicit reaction channel, for example '(p,x)' or '(p,f)'.")
    if is_distribution and product:
        return _invalid("reaction.product must be omitted for a product-isotope distribution.")
    energy_or_outcome = _normalize_energy(spec.get("incident_energy"))
    if isinstance(energy_or_outcome, dict) and "status" in energy_or_outcome:
        return energy_or_outcome
    if is_distribution and energy_or_outcome is None:
        return _clarify(["incident_energy"], "A product-isotope distribution requires a specified incident energy.")
    kinematics_or_outcome = _normalize_kinematics(spec.get("kinematics", {}))
    if isinstance(kinematics_or_outcome, dict) and "status" in kinematics_or_outcome:
        return kinematics_or_outcome
    criteria = {"projectile_code": projectile, "target_code": target, "quantity_code": quantity}
    if process:
        criteria["process_code"] = process
    if product:
        criteria["product_code"] = product
    selection = spec.get("selection", {})
    if not isinstance(selection, dict):
        return _invalid("selection must be an object when supplied.")
    dataset_policy = selection.get("dataset_policy", "all")
    if dataset_policy != "all":
        return _unsupported("selection.dataset_policy", str(dataset_policy), "Only 'all' is implemented; the tool never silently deduplicates datasets.")
    yield_code = None
    if is_distribution:
        yield_definition = selection.get("yield_definition")
        aliases = {"independent": "IND", "ind": "IND", "cumulative": "CUM", "cum": "CUM"}
        yield_code = aliases.get(str(yield_definition).strip().lower()) if yield_definition is not None else None
        if yield_code is None:
            return _clarify(["selection.yield_definition"], "Choose independent or cumulative yield; the tool will not mix them.")
    return {
        "criteria": criteria,
        "energy": energy_or_outcome,
        "kinematics": kinematics_or_outcome,
        "selection": {"dataset_policy": "all", "yield_definition": selection.get("yield_definition")},
        "yield_code": yield_code,
        "public": {
            "reaction": {"projectile": projectile, "target": target, "channel": process, "product": product, "outgoing_particle": outgoing},
            "observable": str(observable).strip().lower(),
            "incident_energy": energy_or_outcome,
            "kinematics": kinematics_or_outcome,
        },
    }


def _normalize_energy(energy: Any) -> dict[str, Any] | None:
    if energy is None:
        return None
    if not isinstance(energy, dict):
        return _invalid("incident_energy must be an object.")
    mode = energy.get("mode", "exact")
    unit = str(energy.get("unit", "")).upper()
    if mode not in {"exact", "nearest", "range"}:
        return _invalid("incident_energy.mode must be exact, nearest, or range.")
    if unit not in _ENERGY_FACTORS:
        return _unsupported("incident_energy.unit", unit, "Supported incident-energy units are eV, keV, MeV, and GeV.")
    keys = ("minimum", "maximum") if mode == "range" else ("value",)
    try:
        values = {key: float(energy[key]) for key in keys}
    except (KeyError, TypeError, ValueError):
        return _invalid(f"incident_energy for mode {mode!r} requires numeric {', '.join(keys)}.")
    if any(not math.isfinite(value) or value < 0 for value in values.values()):
        return _invalid("incident energy must be a finite non-negative value.")
    if mode == "range" and values["minimum"] > values["maximum"]:
        return _invalid("incident_energy.minimum cannot exceed incident_energy.maximum.")
    result = {"mode": mode, "unit": "MeV", **{key: value * _ENERGY_FACTORS[unit] for key, value in values.items()}}
    return result


def _normalize_kinematics(kinematics: Any) -> dict[str, dict[str, Any]]:
    if not isinstance(kinematics, dict):
        return _invalid("kinematics must be an object when supplied.")
    normalized: dict[str, dict[str, Any]] = {}
    for name, selector in kinematics.items():
        if name not in _KINEMATIC_HEADERS:
            return _unsupported("kinematics", name, "Supported kinematic selections are angle and outgoing_energy.")
        factors = _ANGLE_FACTORS if name == "angle" else _ENERGY_FACTORS
        normalized_selector = _normalize_selector(selector, factors, "degree" if name == "angle" else "MeV")
        if "status" in normalized_selector:
            return normalized_selector
        normalized[name] = normalized_selector
    return normalized


def _normalize_selector(selector: Any, factors: dict[str, float], canonical_unit: str) -> dict[str, Any]:
    if not isinstance(selector, dict):
        return _invalid("A kinematic selection must be an object with value/unit/mode.")
    mode = selector.get("mode", "exact")
    unit = str(selector.get("unit", "")).upper()
    if mode not in {"exact", "nearest", "range"}:
        return _invalid("A kinematic selection mode must be exact, nearest, or range.")
    if unit not in factors:
        return _unsupported("kinematics.unit", unit, f"Unsupported unit; expected one of {', '.join(factors)}.")
    keys = ("minimum", "maximum") if mode == "range" else ("value",)
    try:
        values = {key: float(selector[key]) for key in keys}
    except (KeyError, TypeError, ValueError):
        return _invalid(f"Kinematic selection for mode {mode!r} requires numeric {', '.join(keys)}.")
    if any(not math.isfinite(value) for value in values.values()) or (mode == "range" and values["minimum"] > values["maximum"]):
        return _invalid("Kinematic values must be finite and have an ordered range.")
    return {"mode": mode, "unit": canonical_unit, **{key: value * factors[unit] for key, value in values.items()}}


def _nuclide_code(value: Any) -> str | None:
    compact = str(value).strip().upper().replace(" ", "").replace("_", "-")
    if compact.count("-") == 2:
        z, symbol, mass = compact.split("-")
        if z.isdigit() and mass.isdigit() and _SYMBOL_TO_Z.get(symbol) == int(z):
            return f"{int(z)}-{symbol}-{int(mass)}"
    import re
    mass_symbol = re.fullmatch(r"(\d+)([A-Z]{1,2})", compact)
    if mass_symbol:
        mass_text, symbol = mass_symbol.groups()
        z = _SYMBOL_TO_Z.get(symbol)
        return f"{z}-{symbol}-{int(mass_text)}" if z is not None else None
    match = re.fullmatch(r"(?:(\d+)-?)?([A-Z]{1,2})-?(\d+)", compact)
    if not match:
        return None
    z_text, symbol, mass_text = match.groups()
    z = _SYMBOL_TO_Z.get(symbol)
    if z is None or (z_text and int(z_text) != z):
        return None
    return f"{z}-{symbol}-{int(mass_text)}"


def _process_code(channel: Any, projectile: str) -> str | None:
    if channel is None or str(channel).strip() == "":
        return None
    text = str(channel).upper().replace(" ", "")
    if text.startswith("(") and text.endswith(")"):
        text = text[1:-1]
    if "," in text:
        incoming, outgoing = text.split(",", 1)
        if _PROJECTILES.get(incoming) != projectile or not outgoing:
            return None
        text = outgoing
    aliases = {"GAMMA": "G", "A": "A", "ALPHA": "A", "FISSION": "F", "ELASTIC": "EL", "INELASTIC": "INL"}
    return aliases.get(text, text) if text.replace("/", "").replace("+", "").replace("-", "").isalnum() else None


def _read_subentry(source_root: str | Path, candidate: dict[str, Any]) -> dict[str, Any]:
    path = Path(source_root) / candidate["relative_path"]
    entry = parse_file(path)
    for index, subentry in enumerate(entry["subentries"]):
        if subentry["id"] == candidate["subentry_id"]:
            return _with_inherited_common(entry["subentries"], index)
    raise ValueError("indexed subentry is absent from its source entry")


def _with_inherited_common(subentries: list[dict[str, Any]], selected_index: int) -> dict[str, Any]:
    """Apply the entry control subentry's COMMON tables to later data subentries.

    EXFOR entries conventionally place shared experimental conditions in a
    leading ``...001`` control subentry.  Keep those raw tables separate but
    make their values effective for every later selected subentry.
    """

    selected = dict(subentries[selected_index])
    inherited = [table for subentry in subentries[:selected_index] if subentry["id"].endswith("001") for table in subentry["common"]]
    selected["common"] = [*selected["common"], *inherited]
    selected["inherited_common"] = inherited
    selected["entry_bibliography"] = {key: values for control in subentries[:selected_index] if control["id"].endswith("001") for key, values in control.get("bib", {}).items()}
    return selected


def _select_data(subentry: dict[str, Any], energy: dict[str, Any] | None, kinematics: dict[str, dict[str, Any]]) -> tuple[list[dict[str, Any]], str | None]:
    selected: list[dict[str, Any]] = []
    found_energy = energy is None
    for table in subentry["data"]:
        row_indices = list(range(len(table["rows"])))
        energy_points = _energy_points(table, subentry["common"])
        if energy is not None:
            if energy_points is None:
                continue
            found_energy = True
            row_indices = _select_indices(energy_points, energy)
        for name, selector in kinematics.items():
            points = _kinematic_points(table, subentry["common"], name)
            if points is None:
                row_indices = []
                break
            allowed = set(_select_indices(points, selector))
            row_indices = [index for index in row_indices if index in allowed]
        rows = [table["rows"][index] for index in row_indices]
        if energy_points is not None:
            energy_points = [energy_points[index] for index in row_indices]
        if rows:
            selected.append({"headers": table["headers"], "units": table["units"], "rows": rows, "incident_energy_mev": energy_points, "common": subentry["common"]})
    return selected, None if found_energy else "no_energy_column"


def _kinematic_points(table: dict[str, Any], common_tables: list[dict[str, Any]], name: str) -> list[float | None] | None:
    headers = _KINEMATIC_HEADERS[name]
    factors = _ANGLE_FACTORS if name == "angle" else _ENERGY_FACTORS
    for header in headers:
        if header in table["headers"]:
            column = table["headers"].index(header)
            factor = factors.get(table["units"][column].upper())
            return [row[column]["value"] * factor if factor and row[column]["value"] is not None else None for row in table["rows"]]
    for common in common_tables:
        for header in headers:
            if header in common["headers"] and len(common["rows"]) == 1:
                column = common["headers"].index(header)
                factor = factors.get(common["units"][column].upper())
                value = common["rows"][0][column]["value"]
                if factor and value is not None:
                    return [value * factor for _ in table["rows"]]
    return None


def _energy_points(table: dict[str, Any], common_tables: list[dict[str, Any]]) -> list[float | None] | None:
    try:
        column = table["headers"].index("EN")
    except ValueError:
        column = None
    if column is not None:
        factor = _ENERGY_FACTORS.get(table["units"][column].upper())
        return [row[column]["value"] * factor if factor and row[column]["value"] is not None else None for row in table["rows"]]
    for common in common_tables:
        if "EN" not in common["headers"] or len(common["rows"]) != 1:
            continue
        column = common["headers"].index("EN")
        factor = _ENERGY_FACTORS.get(common["units"][column].upper())
        value = common["rows"][0][column]["value"]
        if factor and value is not None:
            return [value * factor for _ in table["rows"]]
    return None


def _select_indices(points: list[float | None], energy: dict[str, Any]) -> list[int]:
    valid = [(index, point) for index, point in enumerate(points) if point is not None]
    if energy["mode"] == "range":
        return [index for index, point in valid if energy["minimum"] <= point <= energy["maximum"]]
    if energy["mode"] == "exact":
        tolerance = max(abs(energy["value"]) * 1e-9, 1e-12)
        return [index for index, point in valid if abs(point - energy["value"]) <= tolerance]
    if not valid:
        return []
    nearest = min(abs(point - energy["value"]) for _, point in valid)
    tolerance = max(abs(energy["value"]) * 1e-12, 1e-15)
    return [index for index, point in valid if abs(abs(point - energy["value"]) - nearest) <= tolerance]


def _record(candidate: dict[str, Any], subentry: dict[str, Any], data: list[dict[str, Any]]) -> dict[str, Any]:
    parsed = _yield_code(candidate["raw_code"])
    return {
        "accession": candidate["entry_id"],
        "subentry": candidate["subentry_id"],
        "reaction": {**{key: candidate[key] for key in ("raw_code", "target_code", "projectile_code", "process_code", "product_code", "quantity_code")}, "yield_code": parsed},
        "bibliography": subentry["bib"],
        "entry_bibliography": subentry.get("entry_bibliography", {}),
        "source": {"relative_path": candidate["relative_path"], "subentry_status": candidate["subentry_status"]},
        "data": data,
    }


def _provenance(index_path: str | Path, source_root: str | Path) -> dict[str, str]:
    metadata = index_metadata(index_path)
    return {
        "index_source_root": metadata.get("source_root", ""),
        "queried_source_root": str(Path(source_root).resolve()),
        "git_remote": metadata.get("git_remote", ""),
        "git_branch": metadata.get("git_branch", ""),
        "git_commit": metadata.get("git_commit", ""),
    }


def _clarify(fields: list[str], message: str) -> dict[str, Any]:
    return {"status": "clarification_needed", "missing_fields": fields, "message": message}


def _invalid(message: str) -> dict[str, Any]:
    return {"status": "invalid_query", "message": message}


def _unsupported(field: str, value: str, message: str) -> dict[str, Any]:
    return {"status": "unsupported", "field": field, "value": value, "message": message}


def _yield_code(raw_reaction: str) -> str | None:
    """Read the EXFOR yield/state field without changing the portable M2 index schema."""

    from .model import parse_reaction_code

    return parse_reaction_code(raw_reaction).get("yield_code")
