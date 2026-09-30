#!/usr/bin/env python3
"""Replot the frozen tool query results, without EXFOR database or network."""
from pathlib import Path
import json
from exfor_tool.angular_grid import plot_np_angular_grid

def main():
    root=Path(__file__).resolve().parent
    config=json.loads((root/'panel_queries.json').read_text())
    outcomes=[json.loads((root/f).read_text()) for f in config['outcomes']]
    prepared={float(e):json.loads((root/path).read_text()) for e,path in config.get('prepared_panels',{}).items()}
    overrides=json.loads((root/'unit_overrides.json').read_text())
    result=plot_np_angular_grid(outcomes,root,unit_overrides=overrides,prepared_panels=prepared)
    print(json.dumps({k:result[k] for k in ['status','png','pdf','energies_MeV','points']},indent=2))
if __name__=='__main__':main()
