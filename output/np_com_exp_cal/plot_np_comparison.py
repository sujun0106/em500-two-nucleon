#!/usr/bin/env python3
"""Overlay actual Fortran output with the frozen EXFOR data, no fitting."""
from pathlib import Path
import argparse,csv,json
from exfor_tool.angular_grid import plot_np_angular_grid

def main():
    root=Path(__file__).resolve().parent
    parser=argparse.ArgumentParser();parser.add_argument('--batch',type=Path);args=parser.parse_args()
    batch=args.batch or Path((root/'latest_fortran_batch.txt').read_text().strip())
    if not batch.is_absolute():batch=root/batch
    if not (batch/'completed_utc.txt').exists():raise RuntimeError('Fortran batch is incomplete')
    config=json.loads((root/'panel_queries.json').read_text())
    curves={}
    for energy in config['energies_MeV']:
        token=f'{energy:g}'.replace('.','p');run=batch/(token+'MeV')
        with (run/'np_angular_distribution.csv').open() as f:rows=list(csv.DictReader(f))
        curves[energy]={'theta_cm_deg':[float(r['theta_cm_deg']) for r in rows],
                       'dsigma_domega_mb_sr':[float(r['dsigma_domega_mb_sr']) for r in rows],
                       'source':str((run/'np_angular_distribution.csv').relative_to(root)),
                       'model':'fixed EM500 NN-only; relativistic kinematics; N=96; pmax=1000 MeV/c; Jmax=32'}
    outcomes=[json.loads((root/f).read_text()) for f in config['outcomes']]
    prepared={float(e):json.loads((root/path).read_text()) for e,path in config.get('prepared_panels',{}).items()}
    overrides=json.loads((root/'unit_overrides.json').read_text())
    result=plot_np_angular_grid(outcomes,root,unit_overrides=overrides,theory_curves=curves,prepared_panels=prepared)
    print(json.dumps({k:result[k] for k in ['png','pdf','energies_MeV','points']},indent=2))
if __name__=='__main__':main()
