#!/usr/bin/env python3
"""Overlay frozen pp experiments with actual Fortran output; no fitting."""
from pathlib import Path
import argparse,csv,json
from exfor_tool.pp_angular_grid import plot_pp_angular_grid

def main():
    root=Path(__file__).resolve().parent
    parser=argparse.ArgumentParser()
    parser.add_argument('--batch',type=Path,help='Folder containing <energy>MeV run subfolders')
    args=parser.parse_args()
    batch=args.batch or Path((root/'latest_fortran_batch.txt').read_text().strip())
    if not batch.is_absolute():batch=root/batch
    config=json.loads((root/'panel_queries.json').read_text())
    curves={}
    for energy in config['energies_MeV']:
        token=f'{energy:g}'.replace('.','p');run=batch/(token+'MeV')
        source=run/'pp_angular_distribution.csv'
        with source.open() as f:rows=list(csv.DictReader(f))
        summary_file=run/'pp_summary.txt'
        summary={}
        for line in summary_file.read_text().splitlines():
            if '=' in line:
                key,value=line.split('=',1);summary[key.strip()]=value.strip()
        curves[energy]={'theta_cm_deg':[float(r['theta_cm_deg']) for r in rows],
                       'dsigma_domega_mb_sr':[float(r['dsigma_domega_mb_sr']) for r in rows],
                       'source':str(source.relative_to(root)) if source.is_relative_to(root) else str(source),
                       'model':'EM500 + Coulomb (Fortran)', 'run_summary':summary}
    outcomes=[json.loads((root/f).read_text()) for f in config['outcomes']]
    result=plot_pp_angular_grid(outcomes,root,theory_curves=curves)
    print(json.dumps({k:result[k] for k in ('png','pdf','energies_MeV','points')},indent=2))
if __name__=='__main__':main()
