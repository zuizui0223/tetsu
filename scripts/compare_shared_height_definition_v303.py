"""Same AusTraits stature definition across Europe/Australia; area remains unmatched."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from analyze_europe_rail_to_Kashmir_elevation_v302 import european_tables
from audit_SI_CZ_original_continuous_trait_coverage_v220 import binomial,read_leda
from analyze_region_rail_elevation_transfer_v302 import fit_foreign,predict_score
from analyze_australian_stature_v303 import run as australian_run


def slope(x,y):
    return float(np.dot(x-x.mean(),y-y.mean())/np.dot(x-x.mean(),x-x.mean()))


def run(folder):
    folder=Path(folder);data=folder/'nonEurope_sources'
    traits=pd.read_csv(data/'austraits_height.csv',dtype=str)
    traits=traits[(traits.trait_name=='plant_height')&(traits.unit=='m')&(traits.life_stage=='adult')&
        (traits.entity_type=='species')&(traits.value_type=='maximum')].copy()
    traits['value']=pd.to_numeric(traits.value,errors='coerce');traits=traits[traits.value>0]
    traits['taxon_name']=traits.taxon_name.map(lambda s:' '.join(s.split()))
    values=traits.groupby(['taxon_name','dataset_id']).value.median().groupby(level=0).median().to_dict()
    assert hashlib.sha256((folder/'canopy_height.txt').read_bytes()).hexdigest()=='cbf5e761390d78a843512dfc1139acfc1154d46d6c1995a827354840a9f78a92'
    leda,_=read_leda(folder/'canopy_height.txt','height_m')
    training={};out={}
    for study,(sites,elevations) in european_tables(folder).items():
        rows=[]
        for site in sites:
            taxa={binomial(t) for t in site['source_names']}-{None};known=taxa&values.keys()
            if len(known)<5 or len(known)/len(taxa)<.5:continue
            rows.append(dict(x=elevations[site['station']]/1000,y=float(np.mean([np.log(values[t]) for t in known])),
                coverage=len(known)/len(taxa),taxa=known))
        x=np.array([r['x'] for r in rows]);y=np.array([r['y'] for r in rows])
        training[study]=[(x,y)]
        deleted=[slope(np.delete(x,i),np.delete(y,i)) for i in range(len(x))]
        out[study]=dict(n=len(rows),coverage_range=[min(r['coverage'] for r in rows),max(r['coverage'] for r in rows)],
            slope_log_height_per_100m=slope(x,y)/10,leave_station_out_range_per_100m=[min(deleted)/10,max(deleted)/10])
        # Match species within each site when comparing trait sources, avoiding a cohort-change explanation.
        common=[]
        for r in rows:
            k=r['taxa']&leda.keys()
            if len(k)>=5:common.append((r['x'],np.mean([np.log(values[t]) for t in k]),np.mean([np.log(leda[t]) for t in k])))
        c=np.array(common);out[study]['same_species_trait_source_sensitivity']=dict(n=len(c),
            austraits_slope_per_100m=slope(c[:,0],c[:,1])/10,leda_slope_per_100m=slope(c[:,0],c[:,2])/10)
    australian_result,rows=australian_run(data)
    x=np.array([r['elevation_m']/1000 for r in rows]);y=np.array([r['log_height'] for r in rows])
    model=fit_foreign(training,minimum_sources=2)
    external={a:predict_score(x,y,model,a) for a in ['lowest','median','none']}
    return dict(definition='AusTraits adult species maximum height, contributor-balanced; mean log m',
        European_railways=out,Australian_fixed_area=australian_result['results'],
        Europe_to_Australia=external,
        limitation='Same trait definition, but no shared strict elevation range or plot area; transfer is extrapolative and exploratory.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--data',required=True);p.add_argument('--output',required=True);a=p.parse_args()
    r=run(a.data);Path(a.output).write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r['European_railways'],indent=2));print(json.dumps(r['Europe_to_Australia'],indent=2))
