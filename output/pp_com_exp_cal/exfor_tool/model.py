"""Normalized, lossless-enough metadata used by the local EXFOR index."""

from __future__ import annotations

from typing import Any


def parse_reaction_code(raw: str) -> dict[str, Any]:
    """Extract the simple EXFOR REACTION form while retaining its raw code.

    Complex algebraic, vector, and multi-reaction expressions intentionally
    return ``parse_status='unparsed'``. They remain indexed by raw code and
    can be added without changing this public metadata shape.
    """

    result: dict[str, Any] = {
        "raw": raw,
        "parse_status": "unparsed",
        "target_code": None,
        "projectile_code": None,
        "process_code": None,
        "product_code": None,
        "yield_code": None,
        "quantity_code": None,
        "qualifiers": [],
    }
    start = raw.find("(")
    if start < 0:
        return result
    outer_end = _matching_parenthesis(raw, start)
    if outer_end is None or raw[outer_end + 1 :].strip():
        return result

    body = raw[start + 1 : outer_end]
    inner_start = body.find("(")
    if inner_start <= 0:
        return result
    inner_end = _matching_parenthesis(body, inner_start)
    if inner_end is None:
        return result

    target = body[:inner_start].strip()
    projectile_process = body[inner_start + 1 : inner_end].split(",")
    suffix = body[inner_end + 1 :]
    if len(projectile_process) != 2:
        return result
    fields = suffix.split(",")
    if len(fields) < 3:
        return result

    projectile, process = (item.strip() for item in projectile_process)
    product, yield_code, quantity = (item.strip() for item in fields[:3])
    if not target or not projectile or not process or not quantity:
        return result
    result.update(
        {
            "parse_status": "simple",
            "target_code": target.upper(),
            "projectile_code": projectile.upper(),
            "process_code": process.upper(),
            "product_code": product.upper() or None,
            "yield_code": yield_code.upper() or None,
            "quantity_code": quantity.upper(),
            "qualifiers": [item.strip().upper() for item in fields[3:] if item.strip()],
        }
    )
    return result


def _matching_parenthesis(text: str, start: int) -> int | None:
    depth = 0
    for index, character in enumerate(text[start:], start=start):
        if character == "(":
            depth += 1
        elif character == ")":
            depth -= 1
            if depth == 0:
                return index
    return None
