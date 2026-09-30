#!/usr/bin/env python3
"""Replot frozen EXFOR pp experiments; no database or network needed."""
from pathlib import Path
import json
from exfor_tool.pp_angular_grid import plot_pp_angular_grid
root=Path(__file__).resolve().parent
config=json.loads((root/'panel_queries.json').read_text())
outcomes=[json.loads((root/path).read_text()) for path in config['outcomes']]
result=plot_pp_angular_grid(outcomes,root)
print(json.dumps({k:result[k] for k in ('png','pdf','energies_MeV','points')},indent=2))
