"""Build and query a portable SQLite metadata index over EXFOR master files."""

from __future__ import annotations

import json
from pathlib import Path
import sqlite3
import subprocess
from collections.abc import Callable
from typing import Any, Iterator

from .model import parse_reaction_code
from .parser import ExforParseError, parse_metadata_file


SCHEMA = """
PRAGMA foreign_keys = ON;
CREATE TABLE metadata (
    key TEXT PRIMARY KEY,
    value TEXT NOT NULL
);
CREATE TABLE entries (
    entry_id TEXT PRIMARY KEY,
    relative_path TEXT NOT NULL,
    byte_size INTEGER NOT NULL,
    mtime_ns INTEGER NOT NULL,
    parse_status TEXT NOT NULL,
    parse_error TEXT
);
CREATE TABLE subentries (
    subentry_id TEXT PRIMARY KEY,
    entry_id TEXT NOT NULL REFERENCES entries(entry_id),
    status TEXT NOT NULL,
    bib_json TEXT NOT NULL,
    common_table_count INTEGER NOT NULL,
    data_table_count INTEGER NOT NULL
);
CREATE TABLE reactions (
    subentry_id TEXT NOT NULL REFERENCES subentries(subentry_id),
    ordinal INTEGER NOT NULL,
    raw_code TEXT NOT NULL,
    parse_status TEXT NOT NULL,
    target_code TEXT,
    projectile_code TEXT,
    process_code TEXT,
    product_code TEXT,
    quantity_code TEXT,
    qualifiers_json TEXT NOT NULL,
    PRIMARY KEY (subentry_id, ordinal)
);
CREATE INDEX reactions_lookup ON reactions (
    projectile_code, target_code, process_code, product_code, quantity_code
);
"""


def build_index(
    source_root: str | Path,
    destination: str | Path,
    *,
    batch_size: int = 500,
    progress: Callable[[dict[str, int]], None] | None = None,
) -> dict[str, int]:
    """Build an atomically replaced SQLite index from an ``exforall`` tree."""

    source_root = Path(source_root)
    destination = Path(destination)
    if not source_root.is_dir():
        raise ValueError(f"EXFOR source directory does not exist: {source_root}")
    if batch_size <= 0:
        raise ValueError("batch_size must be positive")
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_suffix(f"{destination.suffix}.tmp")
    temporary.unlink(missing_ok=True)

    statistics = {"entries": 0, "parsed_entries": 0, "parse_errors": 0, "subentries": 0, "reactions": 0}
    connection = sqlite3.connect(temporary)
    try:
        connection.executescript(SCHEMA)
        _set_metadata(connection, "source_root", str(source_root.resolve()))
        _set_metadata(connection, "git_commit", _git_value(source_root, "rev-parse", "HEAD"))
        _set_metadata(connection, "git_branch", _git_value(source_root, "branch", "--show-current"))
        _set_metadata(connection, "git_remote", _git_value(source_root, "remote", "get-url", "origin"))

        for path in sorted(source_root.rglob("*.x4")):
            statistics["entries"] += 1
            relative_path = str(path.relative_to(source_root))
            file_stat = path.stat()
            try:
                entry = parse_metadata_file(path)
            except (ExforParseError, UnicodeDecodeError) as exc:
                connection.execute(
                    "INSERT INTO entries VALUES (?, ?, ?, ?, ?, ?)",
                    (path.stem, relative_path, file_stat.st_size, file_stat.st_mtime_ns, "parse_error", str(exc)),
                )
                statistics["parse_errors"] += 1
            else:
                connection.execute(
                    "INSERT INTO entries VALUES (?, ?, ?, ?, ?, ?)",
                    (entry["id"], relative_path, file_stat.st_size, file_stat.st_mtime_ns, "parsed", None),
                )
                statistics["parsed_entries"] += 1
                for subentry in entry["subentries"]:
                    _insert_subentry(connection, entry["id"], subentry, statistics)

            if statistics["entries"] % batch_size == 0:
                connection.commit()
                if progress is not None:
                    progress(dict(statistics))

        for key, value in statistics.items():
            _set_metadata(connection, key, str(value))
        connection.commit()
    finally:
        connection.close()
    temporary.replace(destination)
    return statistics


def search_reactions(index_path: str | Path, **criteria: str) -> list[dict[str, Any]]:
    """Return indexed reaction records matching exact normalized metadata fields."""

    allowed = {"target_code", "projectile_code", "process_code", "product_code", "quantity_code"}
    unknown = set(criteria) - allowed
    if unknown:
        raise ValueError(f"unsupported search criteria: {sorted(unknown)!r}")
    clauses = [f"r.{field} = ?" for field in criteria]
    values = [criteria[field].upper() for field in criteria]
    query = """
        SELECT r.subentry_id, s.entry_id, e.relative_path, r.raw_code,
               r.parse_status, r.target_code, r.projectile_code,
               r.process_code, r.product_code, r.quantity_code, r.qualifiers_json,
               s.status AS subentry_status, s.bib_json,
               s.common_table_count, s.data_table_count
        FROM reactions AS r
        JOIN subentries AS s ON s.subentry_id = r.subentry_id
        JOIN entries AS e ON e.entry_id = s.entry_id
    """
    if clauses:
        query += " WHERE " + " AND ".join(clauses)
    query += " ORDER BY r.subentry_id, r.ordinal"
    connection = sqlite3.connect(index_path)
    try:
        connection.row_factory = sqlite3.Row
        rows = [dict(row) for row in connection.execute(query, values)]
    finally:
        connection.close()
    return rows


def index_metadata(index_path: str | Path) -> dict[str, str]:
    """Return the provenance and build metadata stored in an index."""

    connection = sqlite3.connect(index_path)
    try:
        return {key: value for key, value in connection.execute("SELECT key, value FROM metadata")}
    finally:
        connection.close()


def iter_parse_errors(index_path: str | Path) -> Iterator[dict[str, str]]:
    """Yield entries not covered by the current raw parser."""

    connection = sqlite3.connect(index_path)
    try:
        connection.row_factory = sqlite3.Row
        rows = [
            dict(row)
            for row in connection.execute(
                "SELECT entry_id, relative_path, parse_error FROM entries WHERE parse_status = 'parse_error'"
            )
        ]
    finally:
        connection.close()
    yield from rows


def _insert_subentry(
    connection: sqlite3.Connection, entry_id: str, subentry: dict[str, Any], statistics: dict[str, int]
) -> None:
    connection.execute(
        "INSERT INTO subentries VALUES (?, ?, ?, ?, ?, ?)",
        (
            subentry["id"],
            entry_id,
            subentry["status"],
            json.dumps(subentry["bib"], sort_keys=True),
            subentry["common_table_count"],
            subentry["data_table_count"],
        ),
    )
    statistics["subentries"] += 1
    for ordinal, raw_code in enumerate(subentry["bib"].get("REACTION", [])):
        reaction = parse_reaction_code(raw_code)
        connection.execute(
            "INSERT INTO reactions VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
            (
                subentry["id"],
                ordinal,
                reaction["raw"],
                reaction["parse_status"],
                reaction["target_code"],
                reaction["projectile_code"],
                reaction["process_code"],
                reaction["product_code"],
                reaction["quantity_code"],
                json.dumps(reaction["qualifiers"]),
            ),
        )
        statistics["reactions"] += 1


def _git_value(source_root: Path, *arguments: str) -> str | None:
    repository_root = source_root.parent
    completed = subprocess.run(
        ["git", "-C", str(repository_root), *arguments],
        check=False,
        capture_output=True,
        text=True,
    )
    return completed.stdout.strip() if completed.returncode == 0 else None


def _set_metadata(connection: sqlite3.Connection, key: str, value: str | None) -> None:
    connection.execute("INSERT INTO metadata VALUES (?, ?)", (key, value or ""))
