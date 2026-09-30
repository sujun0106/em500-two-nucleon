"""A small, dependency-free parser for EXFOR master-file entries.

This module deliberately exposes the raw record structure before later stages
interpret reaction codes or normalize physical quantities.
"""

from __future__ import annotations

from pathlib import Path
import re
from typing import Any


FIELD_WIDTH = 11
FIELDS_PER_LINE = 6
_NUMBER = re.compile(
    r"^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+|[+-]\d+)?$"
)


class ExforParseError(ValueError):
    """Raised when a record cannot be interpreted as the supported EXFOR subset."""


def parse_file(path: str | Path) -> dict[str, Any]:
    """Parse an EXFOR entry file and return a JSON-serializable dictionary."""

    text = Path(path).read_text(encoding="ascii")
    return parse_entry(text)


def parse_metadata_file(path: str | Path) -> dict[str, Any]:
    """Read only ENTRY/SUBENT identity and BIB metadata from an EXFOR file."""

    with Path(path).open(encoding="ascii") as source:
        return _parse_metadata_lines(source)


def parse_entry(text: str) -> dict[str, Any]:
    """Parse one EXFOR ENTRY containing zero or more SUBENT records.

    The vertical slice supports BIB, COMMON, DATA, NOCOMMON, NOSUBENT, and
    fixed-width (11-character) table fields. Unsupported malformed records
    raise :class:`ExforParseError` rather than silently discarding content.
    """

    lines = text.splitlines()
    entry: dict[str, Any] | None = None
    current: dict[str, Any] | None = None
    index = 0

    while index < len(lines):
        line = lines[index]
        keyword = _keyword(line)

        if keyword == "ENTRY":
            if entry is not None:
                raise ExforParseError("an input string must contain one ENTRY")
            entry = {"id": _record_id(line, "ENTRY"), "subentries": []}
        elif keyword in {"SUBENT", "NOSUBENT"}:
            if entry is None:
                raise ExforParseError("SUBENT occurred before ENTRY")
            current = {
                "id": _record_id(line, keyword),
                "status": "no_subentry" if keyword == "NOSUBENT" else "active",
                "bib": {},
                "common": [],
                "data": [],
            }
            entry["subentries"].append(current)
        elif keyword == "BIB":
            _require_subentry(current, "BIB")
            bib, index = _parse_bib(lines, index + 1)
            current["bib"] = bib
        elif keyword in {"COMMON", "DATA"}:
            _require_subentry(current, keyword)
            table, index = _parse_table(lines, index, keyword)
            current[keyword.lower()].append(table)
        elif keyword in {"ENDENTRY", "ENDSUBENT", "ENDBIB", "ENDCOMMON", "ENDDATA", "NOCOMMON"}:
            pass
        index += 1

    if entry is None:
        raise ExforParseError("ENTRY record not found")
    return entry


def parse_metadata(text: str) -> dict[str, Any]:
    """Parse only the metadata needed by the local index.

    DATA and COMMON values can be very large and are intentionally skipped.
    Full table parsing remains available through :func:`parse_entry`.
    """

    return _parse_metadata_lines(text.splitlines())


def _parse_metadata_lines(lines: Any) -> dict[str, Any]:
    """Parse index metadata from an iterable of EXFOR lines."""

    entry: dict[str, Any] | None = None
    current: dict[str, Any] | None = None
    line_iterator = iter(lines)
    for raw_line in line_iterator:
        line = raw_line.rstrip("\r\n")
        keyword = _keyword(line)
        if keyword == "ENTRY":
            if entry is not None:
                raise ExforParseError("an input string must contain one ENTRY")
            entry = {"id": _record_id(line, "ENTRY"), "subentries": []}
        elif keyword in {"SUBENT", "NOSUBENT"}:
            if entry is None:
                raise ExforParseError("SUBENT occurred before ENTRY")
            current = {
                "id": _record_id(line, keyword),
                "status": "no_subentry" if keyword == "NOSUBENT" else "active",
                "bib": {},
                "common_table_count": 0,
                "data_table_count": 0,
            }
            entry["subentries"].append(current)
        elif keyword == "BIB":
            _require_subentry(current, "BIB")
            current["bib"] = _parse_bib_stream(line_iterator)
        elif keyword == "COMMON":
            _require_subentry(current, "COMMON")
            current["common_table_count"] += 1
        elif keyword == "DATA":
            _require_subentry(current, "DATA")
            current["data_table_count"] += 1
    if entry is None:
        raise ExforParseError("ENTRY record not found")
    return entry


def _parse_bib(lines: list[str], index: int) -> tuple[dict[str, list[str]], int]:
    bib: dict[str, list[str]] = {}
    last_key: str | None = None
    while index < len(lines) and _keyword(lines[index]) != "ENDBIB":
        line = lines[index]
        key = _keyword(line)
        value = line[10:].strip()
        if key:
            bib.setdefault(key, []).append(value)
            last_key = key
        elif value:
            if last_key is None:
                raise ExforParseError("BIB continuation occurred before a BIB keyword")
            bib[last_key][-1] = f"{bib[last_key][-1]}\n{value}"
        index += 1
    if index == len(lines):
        raise ExforParseError("BIB block has no ENDBIB record")
    return bib, index


def _parse_bib_stream(lines: Any) -> dict[str, list[str]]:
    bib: dict[str, list[str]] = {}
    last_key: str | None = None
    for raw_line in lines:
        line = raw_line.rstrip("\r\n")
        if _keyword(line) == "ENDBIB":
            return bib
        key = _keyword(line)
        value = line[10:].strip()
        if key:
            bib.setdefault(key, []).append(value)
            last_key = key
        elif value:
            if last_key is None:
                raise ExforParseError("BIB continuation occurred before a BIB keyword")
            bib[last_key][-1] = f"{bib[last_key][-1]}\n{value}"
    raise ExforParseError("BIB block has no ENDBIB record")


def _parse_table(
    lines: list[str], index: int, table_type: str
) -> tuple[dict[str, Any], int]:
    count_values = lines[index][10:].split()
    if len(count_values) < 2:
        raise ExforParseError(f"{table_type} record has no column/row counts")
    try:
        column_count, declared_count = (int(count_values[0]), int(count_values[1]))
    except ValueError as exc:
        raise ExforParseError(f"invalid {table_type} dimensions: {count_values!r}") from exc
    if column_count <= 0 or declared_count < 0:
        raise ExforParseError(f"invalid {table_type} dimensions")

    line_count = (column_count + FIELDS_PER_LINE - 1) // FIELDS_PER_LINE
    cursor = index + 1
    headers, cursor = _read_fields(lines, cursor, column_count, line_count)
    units, cursor = _read_fields(lines, cursor, column_count, line_count)
    if table_type == "COMMON":
        remaining_line_count = declared_count - 2 * line_count
        if remaining_line_count < 0 or remaining_line_count % line_count:
            raise ExforParseError(f"invalid COMMON line count: {declared_count}")
        row_count = remaining_line_count // line_count
    else:
        row_count = declared_count
    rows: list[list[dict[str, Any]]] = []
    for _ in range(row_count):
        raw_values, cursor = _read_fields(lines, cursor, column_count, line_count)
        rows.append([_cell(value) for value in raw_values])

    end_keyword = f"END{table_type}"
    if cursor >= len(lines) or _keyword(lines[cursor]) != end_keyword:
        raise ExforParseError(f"{table_type} block has no {end_keyword} record")
    return {
        "headers": headers,
        "units": units,
        "rows": rows,
    }, cursor


def _read_fields(
    lines: list[str], index: int, field_count: int, line_count: int
) -> tuple[list[str], int]:
    values: list[str] = []
    for _ in range(line_count):
        if index >= len(lines):
            raise ExforParseError("unexpected end of table")
        values.extend(_fields(lines[index]))
        index += 1
    return values[:field_count], index


def _cell(raw: str) -> dict[str, Any]:
    value = _parse_number(raw)
    return {"raw": raw, "value": value}


def _parse_number(raw: str) -> float | None:
    # EXFOR fixed-width numbers may pad the mantissa/exponent with blanks.
    # Only recognize that spacing before a signed exponent; preserve raw text.
    value = re.sub(r"(?<=[0-9.])\s+(?=[+-]\d+$)", "", raw.strip())
    if not value or not _NUMBER.fullmatch(value):
        return None
    normalized = value.replace("D", "E").replace("d", "e")
    if "E" not in normalized.upper():
        tail = re.search(r"(?<=[0-9.])([+-]\d+)$", normalized)
        if tail:
            normalized = f"{normalized[:tail.start()]}E{tail.group(1)}"
    return float(normalized)


def _fields(line: str) -> list[str]:
    return [line[offset : offset + FIELD_WIDTH].strip() for offset in range(0, 66, FIELD_WIDTH)]


def _keyword(line: str) -> str:
    return line[:10].strip()


def _record_id(line: str, keyword: str) -> str:
    record_id = line[10:22].strip()
    if not record_id:
        raise ExforParseError(f"{keyword} record has no identifier")
    return record_id


def _require_subentry(current: dict[str, Any] | None, keyword: str) -> None:
    if current is None:
        raise ExforParseError(f"{keyword} occurred outside SUBENT")
