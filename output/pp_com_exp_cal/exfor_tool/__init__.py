"""Deterministic parsing utilities for the EXFOR master-file format."""

from .index import build_index, index_metadata, iter_parse_errors, search_reactions
from .parser import ExforParseError, parse_entry, parse_file, parse_metadata, parse_metadata_file
from .query import execute_query
from .render import RenderError, flatten_table, plot_ddx_by_incident_energy_and_angle, plot_differential_spectrum, plot_excitation_function, plot_product_isotope_maps, write_table

__all__ = [
    "ExforParseError",
    "build_index",
    "execute_query",
    "index_metadata",
    "iter_parse_errors",
    "parse_entry",
    "parse_file",
    "parse_metadata",
    "parse_metadata_file",
    "RenderError",
    "search_reactions",
    "flatten_table",
    "plot_excitation_function",
    "plot_differential_spectrum",
    "plot_ddx_by_incident_energy_and_angle",
    "plot_product_isotope_maps",
    "write_table",
]
