"""Plot source-level transfer gains from the frozen v302 aggregates."""
import argparse
import json
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np


def draw(data,output):
    roads=json.loads((data/'v302_height_only.json').read_text(encoding='utf8'))['primary']
    rail=json.loads((data/'v302_rail_to_rail.json').read_text(encoding='utf8'))['height_only/primary']
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(1,2,figsize=(11.4,4.8),gridspec_kw={'width_ratios':[1.15,1]},layout='constrained')
    labels={'AUT':'Austria','CHE':'Switzerland','CLS':'Chile (CLS)','IND':'India','ITA':'Italy','NOR':'Norway'}
    for i,(region,row) in enumerate(roads['regions']['roadside'].items()):
        for anchor,offset,color,marker in [('lowest',-.12,'#216b78','o'),('median',.12,'#b96530','s')]:
            axes[0].scatter(row['scores'][anchor]['gain'],i+offset,c=color,marker=marker,s=48,
                            label=('Lowest-elevation anchor' if anchor=='lowest' else 'Median-elevation anchor') if i==0 else None)
    axes[0].set_yticks(range(6),[labels[r] for r in roads['regions']['roadside']])
    axes[0].invert_yaxis();axes[0].axvline(0,color='#888888',lw=1)
    axes[0].set_xlabel('Prediction improvement (baseline RMSE − model RMSE)')
    axes[0].set_title('Whole-region holdout on mountain roads',loc='left',fontweight='bold')
    axes[0].legend(loc='lower right',frameon=False,fontsize=8)
    labels2=['Parma railway','Croatia / Bosnia railway','Kashmir railway 2014','Kashmir railway 2017']
    series=[rail['observed_slopes']['IT28'],rail['observed_slopes']['HR_BA15']]
    series += [rail['Kashmir'][f'station/{year}/Track+Embankment']['result'] for year in (2014,2017)]
    for i,row in enumerate(series):
        value=row.get('slope',row.get('observed_slope'))
        low,high=row.get('delete_one_range',row.get('delete_one_site_slope_range'))
        axes[1].plot([low,high],[i,i],color='#777777',lw=2)
        axes[1].scatter(value,i,s=60,c='#216b78' if i<2 else '#b96530',zorder=3)
    axes[1].axvline(0,color='#888888',lw=1)
    axes[1].set_yticks(range(4),labels2);axes[1].invert_yaxis()
    axes[1].set_xlabel('Community stature slope per 1 km elevation')
    axes[1].set_title('Railway gradients differ among sources',loc='left',fontweight='bold')
    fig.suptitle('Height-only LEDA matches · Species-equal mean z(log canopy height)',fontsize=12)
    fig.supxlabel('Road predictions use one measured anchor per held-out road. Railway lines show leave-one-site-out slope ranges, not confidence intervals.',fontsize=8)
    output.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(output,dpi=180)
    plt.close(fig)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--data',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();draw(a.data,a.output)
