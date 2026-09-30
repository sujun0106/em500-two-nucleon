"""Traceable CM np angular-distribution extraction and multi-panel rendering.

Consumes QueryOutcome only. No .x4 parsing, evaluation curves, interpolation,
lab/CM kinematic transformation, normalization fit, or dataset averaging.
"""
from __future__ import annotations
import csv
import copy
import json
import math
from pathlib import Path
from .render import RenderError, _require_result, _ddx_reference_frame

_AREA = {'B/SR':1000., 'MB/SR':1., 'UB/SR':.001, 'MUB/SR':.001, 'NB/SR':1.e-6}
_ERROR_ORDER = ('ERR-T','DATA-ERR','DATA-ERR1','ERR-S')

def _value(table, row, header):
    if header in table['headers']:
        i=table['headers'].index(header)
        return row[i]['value'],table['units'][i],row[i]['raw']
    for common in table.get('common',[]):
        if header in common['headers'] and len(common['rows'])==1:
            i=common['headers'].index(header)
            return common['rows'][0][i]['value'],common['units'][i],common['rows'][0][i]['raw']
    return None,None,None

def extract_np_cm(outcome, unit_overrides=None):
    _require_result(outcome)
    points=[]; excluded=[]
    for record in outcome['records']:
        code=record['reaction']['raw_code']
        # Retain every compatible series; exclusions are recorded, never silent.
        if code != '(1-H-1(N,EL)1-H-1,,DA)':
            excluded.append({'subentry':record['subentry'],'reason':'Requires plain np DA; other outgoing-particle, relative, coefficient, residual or spectrum-averaged definitions are not mixed.','reaction':code});continue
        status=' '.join(record.get('bibliography',{}).get('STATUS',[]))
        if 'SPSDD' in status:
            excluded.append({'subentry':record['subentry'],'reason':'EXFOR explicitly marks data superseded (SPSDD).','status':status});continue
        for ti,table in enumerate(record['data']):
            if _ddx_reference_frame(table['headers'],table.get('common',[]))!='cm':
                excluded.append({'subentry':record['subentry'],'table':ti,'reason':'No explicit CM field; no lab-to-CM transformation attempted.'});continue
            angle_header=next((h for h in ('ANG-CM','COS-CM') if h in table['headers'] or any(h in c['headers'] for c in table.get('common',[]))),None)
            data_header=next((h for h in ('DATA-CM','DATA') if h in table['headers']),None)
            if not angle_header or not data_header:
                excluded.append({'subentry':record['subentry'],'table':ti,'reason':'No pointwise CM angle/data columns.'});continue
            for ri,row in enumerate(table['rows']):
                energy=(table.get('incident_energy_mev') or [None]*len(table['rows']))[ri]
                angle,au,ar=_value(table,row,angle_header);value,vu,vr=_value(table,row,data_header)
                if energy is None or angle is None or value is None or vu.upper() not in _AREA:
                    excluded.append({'subentry':record['subentry'],'table':ti,'row':ri,'reason':'Missing numeric energy/angle/cross section or unsupported cross-section unit.'});continue
                if angle_header=='COS-CM':
                    if au.upper()!='NO-DIM' or abs(angle)>1+1.e-12:raise RenderError('Invalid COS-CM value/unit')
                    theta=math.degrees(math.acos(max(-1.,min(1.,angle))))
                elif au.upper() in ('ADEG','DEG'):theta=angle
                elif au.upper()=='RAD':theta=math.degrees(angle)
                else:raise RenderError('Unsupported angle unit')
                applicable=[o for o in (unit_overrides or []) if o['subentry']==record['subentry'] and math.isclose(o['energy_MeV'],energy,abs_tol=1.e-10)]
                if len(applicable)>1:raise RenderError('Multiple unit overrides match one point')
                override=applicable[0] if applicable else None
                if override and (vu.upper()!=override['from_unit'] or override['to_unit'] not in _AREA):raise RenderError('Unit override does not match source unit')
                effective_unit=override['to_unit'] if override else vu.upper()
                y=value*_AREA[effective_unit]
                if not (0<=theta<=180) or y<0:raise RenderError('Invalid physical angle/cross section')
                err=None;error_field='';errors={}
                for h in ('ERR-T','DATA-ERR','DATA-ERR1','ERR-S','ERR-SYS','ERR-1','ERR-2'):
                    ev,eu,er=_value(table,row,h)
                    if ev is not None:errors[h]={'value':ev,'unit':eu,'raw':er}
                for h in _ERROR_ORDER:
                    if h not in errors:continue
                    ev=errors[h]['value'];eu=errors[h]['unit'].upper()
                    if eu=='PER-CENT':err=abs(y*ev/100.)
                    elif eu in _AREA:err=abs(ev*_AREA[override['to_unit'] if override and eu==override['from_unit'] else eu])
                    else:continue
                    error_field=h;break
                flag,_,_=_value(table,row,'FLAG')
                points.append({'energy_MeV':energy,'accession':record['accession'],'subentry':record['subentry'],'table_index':ti,'row_index':ri,'theta_cm_deg':theta,'dsigma_cm_mb_sr':y,'error_mb_sr':err,'error_field':error_field,'flag':flag,'angle_header':angle_header,'angle_original_value':angle,'angle_original_unit':au,'data_header':data_header,'data_original_value':value,'data_original_unit':vu,'errors_original_json':json.dumps(errors,sort_keys=True),'source_relative_path':record['source']['relative_path'],'data_effective_unit':effective_unit,'unit_override_json':json.dumps(override,ensure_ascii=False) if override else ''})
    return {'schema':'exfor-tool/np-cm-angular/v1','points':points,'excluded':excluded,'provenance':outcome['provenance'],'frame_policy':'Existing EXFOR tool policy: explicit -CM field denotes CM; angle must explicitly be ANG-CM or COS-CM. No kinematic transformation.','unit_overrides':unit_overrides or [],'error_policy':'ERR-T, otherwise DATA-ERR/DATA-ERR1, otherwise ERR-S. Reported components preserved separately; no invented or quadrature-combined errors.'}

def plot_np_angular_grid(outcomes, output_directory, unit_overrides=None, theory_curves=None, prepared_panels=None):
    """Render 9 exact-energy outcomes and export both raw and plotted data."""
    if len(outcomes)!=9:raise RenderError('A 3x3 figure requires nine outcomes')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    output=Path(output_directory);output.mkdir(parents=True,exist_ok=True)
    allpoints=[];panels=[]
    for outcome in outcomes:
        selector=outcome['normalized_query']['incident_energy']
        if selector['mode']!='exact':raise RenderError('Panels require exact energy selectors')
        prepared=(prepared_panels or {}).get(selector['value'])
        result=copy.deepcopy(prepared) if prepared is not None else extract_np_cm(outcome,unit_overrides=unit_overrides)
        pts=result['points']
        if not pts:raise RenderError(f'No compatible CM angular data at {selector["value"]} MeV')
        if any(not math.isclose(x['energy_MeV'],selector['value'],abs_tol=1.e-10) for x in pts):raise RenderError('Mixed energies in panel')
        energy=selector['value'];stem=f'np_exp_{energy:g}MeV'.replace('.','p')
        with (output/(stem+'.csv')).open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(pts[0]));w.writeheader();w.writerows(pts)
        result['energy_MeV']=energy;result['csv']=stem+'.csv';panels.append(result);allpoints.extend(pts)
    energies=[x['energy_MeV'] for x in panels]
    if len(set(energies))!=9 or energies!=sorted(energies):raise RenderError('Nine distinct increasing energies required')
    if theory_curves is not None and set(theory_curves)!=set(energies):raise RenderError('Theory energies must exactly match experimental panels')
    with (output/'experimental_points.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(allpoints[0]));w.writeheader();w.writerows(allpoints)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.labelsize':11,'axes.titlesize':12,'axes.linewidth':.8,'savefig.facecolor':'white'})
    fig,axs=plt.subplots(3,3,figsize=(13.2,10.8),sharex=True)
    colors=['#1765a6','#cc4938','#399667','#8859a5'];markers=['o','s','^','D']
    for number,(ax,panel) in enumerate(zip(axs.flat,panels)):
        series={}
        for x in panel['points']:series.setdefault(x['subentry'],[]).append(x)
        ymax=0.
        for ordinal,(subentry,ps) in enumerate(series.items()):
            ps=sorted(ps,key=lambda x:x['theta_cm_deg']);xs=[x['theta_cm_deg'] for x in ps];ys=[x['dsigma_cm_mb_sr'] for x in ps]
            es=[x['error_mb_sr'] if x['error_mb_sr'] is not None else float('nan') for x in ps]
            ax.errorbar(xs,ys,yerr=es,fmt=markers[ordinal%4],linestyle='none',color=colors[ordinal%4],markersize=4.7,markerfacecolor='white',markeredgewidth=1.1,elinewidth=.9,capsize=2.2,label=subentry,zorder=3)
            ymax=max(ymax,max(x['dsigma_cm_mb_sr']+(x['error_mb_sr'] or 0.) for x in ps))
        if theory_curves is not None:
            curve=theory_curves[panel['energy_MeV']]
            tx=curve['theta_cm_deg'];ty=curve['dsigma_domega_mb_sr']
            if len(tx)!=len(ty) or len(tx)<2 or any(not math.isfinite(v) for v in tx+ty):raise RenderError('Invalid theory curve')
            if any(b<=a for a,b in zip(tx,tx[1:])) or min(tx)<0 or max(tx)>180 or min(ty)<0:raise RenderError('Invalid theory angles/cross sections')
            ax.plot(tx,ty,color='#242424',linewidth=1.5,label='EM500 (Fortran)',zorder=2)
            ymax=max(ymax,max(ty));panel['theory']=curve
        ax.set_title(f'({chr(97+number)})  '+r'$E_n^{\mathrm{lab}}$'+f' = {panel["energy_MeV"]:g} MeV',loc='left',pad=9)
        if any(x['unit_override_json'] for x in panel['points']):
            ax.text(.04,.78,'Provisional unit correction*',transform=ax.transAxes,fontsize=9,color='#a43a2f')
        ax.set_xlim(0,180);ax.set_xticks([0,45,90,135,180]);ax.set_ylim(0,ymax*1.26)
        ax.grid(alpha=.18,linewidth=.6);ax.set_axisbelow(True)
        ax.legend(loc='upper center',frameon=False,fontsize=8.0,ncol=min(2,len(series)+(theory_curves is not None)),handletextpad=.4,columnspacing=1.0)
        ax.tick_params(direction='in',top=True,right=True)
        if number//3==2:ax.set_xlabel(r'$\theta_{n,\mathrm{CM}}$ (deg)')
        if number%3==0:ax.set_ylabel(r'$d\sigma/d\Omega$ (mb/sr)')
        panel['n_points']=len(panel['points']);panel['subentries']=list(series)
    title='Neutron–proton elastic scattering: experiment and EM500 (Fortran)' if theory_curves is not None else 'Neutron–proton elastic scattering: experimental angular distributions'
    fig.suptitle(title,fontsize=15,y=.975)
    fig.text(.5,.026,'CM outgoing-neutron angles; independent vertical scales. Symbols: experimental datasets. Bars: reported point uncertainties.',ha='center',fontsize=9.5,color='#444444')
    note='Nominal incident energies; beam widths are retained in the source data.'
    for rec in outcomes[-1].get('records',[]):
        for table in rec['data']:
            if table['rows']:
                width, unit, _ = _value(table,table['rows'][0],'EN-RSL')
                if width is not None:
                    note=f'Nominal incident energies. Last panel: reported EN-RSL = {width:g} {unit}.'
    if prepared_panels:
        note='Franz 2000: published absolute normalization; shared 3.1% scale uncertainty not included in point error bars.'
    if unit_overrides:
        note='* '+ '; '.join(f"{o['energy_MeV']:g} MeV ({o['subentry']}): {o['from_unit']} interpreted as {o['to_unit']} for values and absolute errors; unverified." for o in unit_overrides)
    fig.text(.5,.009,note,ha='center',fontsize=9.5,color='#444444')
    fig.subplots_adjust(left=.075,right=.982,bottom=.095,top=.92,wspace=.24,hspace=.32)
    stem='np_exp_cal_angular_3x3' if theory_curves is not None else 'np_experimental_angular_3x3'
    png=output/(stem+'.png');pdf=output/(stem+'.pdf')
    fig.savefig(png,dpi=240);fig.savefig(pdf);plt.close(fig)
    manifest={'status':'result','data_only':theory_curves is None,'png':png.name,'pdf':pdf.name,'energies_MeV':energies,'points':len(allpoints),'panels':panels,'unit_overrides':unit_overrides or []}
    (output/('comparison_manifest.json' if theory_curves is not None else 'plot_manifest.json')).write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n')
    return manifest
