"""python replication_audit.py work outputs/replication_audit.json ; tetsu 43884d11.
Required inputs and source links are listed in hypothesis_review_ja.md."""
import hashlib
import json
import re
import sys
import zipfile
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
import pymupdf
from scipy.stats import rankdata, spearmanr

data,out=map(Path,sys.argv[1:])
sys.path.insert(0,str(data/"tetsu/scripts"))
import validate_Europe_3studies_82stations_original_LEDA_2D_trait_beta_v245 as eu
import analyze_MIREN_standardized_road_edge_vegetation_v99 as road
from analyze_Slovenia_source_rail_functional_taxonomic_identity_v141 import source_sites_and_guilds
from select_two_original_rail_20_species_for_inat_v164 import original_slovenian_EPPO_codebook_binomial_mapping

from rail_frequency_sources import load_rail, names, frequencies
studies,railrich,kashmir_floras=load_rail(data)

def transfer(training,target):
    shared=sorted(target.keys() & training.keys())
    if len(shared)<3:return {"shared":len(shared)}
    a=np.array([training[t] for t in shared]);b=np.array([target[t] for t in shared])
    cuts=np.quantile(a,[.25,.75])
    lo=b[a<=cuts[0]];hi=b[a>=cuts[1]]
    # Null diagnostic only: permute target ranks within each genus.
    x=rankdata(a);y=rankdata(b)
    x=(x-x.mean())/np.linalg.norm(x-x.mean());y=(y-y.mean())/np.linalg.norm(y-y.mean())
    genera=np.array([t.split()[0] for t in shared])
    groups=[np.where(genera==g)[0] for g in sorted(set(genera)) if (genera==g).sum()>1]
    rng=np.random.default_rng(271117);null=[]
    for _ in range(999):
        perm=y.copy()
        for ids in groups:perm[ids]=rng.permutation(y[ids])
        null.append(float(x@perm))
    other=sorted(training)
    return dict(shared=len(shared),target_total=len(target),rho=float(spearmanr(a,b).statistic),
                lower_mean=lo.mean(),upper_mean=hi.mean(),quartile_count=[len(lo),len(hi)],
                movable=sum(map(len,groups)),within_genus_rank=(1+sum(v>=float(x@y)-1e-12 for v in null))/1000,
                rho_with_unrecorded_training_species_as_zero=float(spearmanr(
                    [training[t] for t in other],[target.get(t,0) for t in other]).statistic))
results={tag:transfer(studies['CZ39'],val) for tag,val in studies.items() if tag!='CZ39'}
with zipfile.ZipFile(data/'splot.zip') as z:
    head=pd.read_csv(z.open('sPlotOpen_header(2).txt'),sep='\t',low_memory=False)
    dt=pd.read_csv(z.open('sPlotOpen_DT(1).txt'),sep='\t',usecols=['PlotObservationID','Species'])
subset=head[head.Resample_1_consensus & (head.is_forest==False)][['PlotObservationID','Country']]
dt=dt.merge(subset,on='PlotObservationID')
lookup={s:eu.binomial(s) for s in dt.Species.dropna().unique()}
dt['sp']=dt.Species.map(lookup)
background={}
for tag,countries in [('SK30',['Slovak Republic']),('IT28',['Italy']),('HR_BA15',['Croatia','Bosnia-Herzegovina']),
                      ('SI89',['Slovenia']),('Kashmir16',['India'])]:
    part=dt[dt.Country.isin(countries)]
    plots=part.PlotObservationID.nunique()
    if plots<20:
        background[tag]={'plots':plots};continue
    freq=part.drop_duplicates(['PlotObservationID','sp']).sp.value_counts()/plots
    shared=sorted(studies[tag].keys() & studies['CZ39'].keys())
    x,y,b=[rankdata(vec) for vec in [[studies['CZ39'][t] for t in shared],[studies[tag][t] for t in shared],
                                     [freq.get(t,0) for t in shared]]]
    design=np.column_stack([np.ones(len(b)),b])
    def residual(vec):return vec-design@np.linalg.lstsq(design,vec,rcond=None)[0]
    background[tag]=dict(plots=plots,shared=len(shared),seen_in_background=sum(t in freq for t in shared),
        target_background_rho=float(np.corrcoef(y,b)[0,1]),partial_rail_rho=float(np.corrcoef(residual(x),residual(y))[0,1]))
    only=[t for t in shared if t in freq]
    if len(only)>=10:
        u,v,w=[rankdata(vec) for vec in [[studies['CZ39'][t] for t in only],[studies[tag][t] for t in only],[freq[t] for t in only]]]
        design=np.column_stack([np.ones(len(w)),w])
        background[tag]['partial_seen_only']=float(np.corrcoef(residual(u),residual(v))[0,1])
allpairs={a+'_'+b:transfer(studies[a],studies[b]) for i,a in enumerate(studies) for b in list(studies)[i+1:]}
# Out-of-sample countercheck: median across other European floras, requiring two contributors.
ensemble={}
for tag,test in studies.items():
    train=[v for name,v in studies.items() if name not in (tag,'Kashmir16')]
    pred={t:np.median([v[t] for v in train if t in v]) for t in set.union(*(set(v) for v in train))
          if sum(t in v for v in train)>=2}
    ensemble[tag]=transfer(pred,test)

file=data/'miren.csv'
assert hashlib.md5(file.read_bytes()).hexdigest()==road.SOURCE_MD5
clean,_=road.native_status_safe(pd.read_csv(file,low_memory=False))
_,p,_=road.matched_siteyears(road.plot_species_richness(clean))
key=list(road.SITE)
def estimates(frame):
    columns=['Alien_roadside_minus_interior','alien_fraction_roadside_minus_interior']
    means=frame.groupby(['Region','Road'])[columns].mean().groupby(level=0).mean()
    res={}
    for col in columns:
        vals=means[col].dropna()
        res[col]=dict(n_regions=len(vals),positive=int((vals>0).sum()),region_equal=vals.mean(),
                      region_values=vals.to_dict(),drop_region_range=[min((vals.sum()-vals)/(len(vals)-1)),max((vals.sum()-vals)/(len(vals)-1))])
    return dict(sites=len(frame),roads=frame.groupby(['Region','Road']).ngroups,metrics=res)
p=p.sort_values('Year')
first=p.drop_duplicates(key);assert len(first)==1005
repeat=p.groupby(key).filter(lambda g:len(g)>=2)
modern=p[p.Year>=2017].drop_duplicates(key)
largest=clean[clean.Status=='Alien'].groupby(['Region','Accepted.Name.MIREN']).size().groupby(level=0).idxmax()
rem=set(largest)
# Preserve every surveyed plot, including plots reduced to zero aliens after removal.
mask=[status=='Alien' and (reg,t) in rem for reg,t,status in clean[['Region','Accepted.Name.MIREN','Status']].itertuples(index=False,name=None)]
clean.loc[mask,'Status']='Unknown'
_,p2,_=road.matched_siteyears(road.plot_species_richness(clean))
roads={name:estimates(frame) for name,frame in [('first',first),('modern_first',modern),
      ('repeat_first',repeat.drop_duplicates(key)),('repeat_last',repeat.drop_duplicates(key,keep='last')),
      ('remove_most_recorded_alien',p2.sort_values('Year').drop_duplicates(key))]}
out.write_text(json.dumps(dict(rail_czech_transfer=results,rail_ensemble_holdout=ensemble,rail_allpairs=allpairs,background=background,
            kashmir_year_zone_transfer={tag:transfer(studies['CZ39'],v) for tag,v in kashmir_floras.items()},
            kashmir_station_embankment_minus_track=railrich,road=roads),indent=2),encoding='utf8')
