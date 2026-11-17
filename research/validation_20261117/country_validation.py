"""Country-held-out comparison. python country_validation.py work output.json"""
import hashlib
import io
import json
import re
import sys
import zipfile
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import pymupdf
from scipy.stats import rankdata

data, output = map(Path, sys.argv[1:])
sys.path.insert(0, str(data / "tetsu/scripts"))
from rail_frequency_sources import load_rail, eu, frequencies, names

studies, _, _ = load_rail(data)
target = {c: studies[k] for c, k in [('CZ','CZ39'),('IT','IT28'),('SK','SK30'),('SI','SI89')]}
adriatic = eu.publisher_hr_ba(data/'hrba.pdf')[0]
for c, rows in [('HR',adriatic[:13]),('BA',adriatic[13:])]:
    target[c] = frequencies([names(s['source_names']) for s in rows])
matrix = pd.read_excel(data/'Matrix_plot.xlsx')
def india(year):
    subset = matrix[(matrix['Year of sampling']==year) & ~matrix.iloc[:,0].str.startswith('Between')].copy()
    subset = subset[subset.Location.isin(['Track','Embankment'])]
    groups = subset.groupby(subset.iloc[:,0].str.strip())
    assert len(groups)==16 and all(set(g.Location)=={'Track','Embankment'} and len(g)==2 for _,g in groups)
    return frequencies([names(g.columns[6:][g.iloc[:,6:].gt(0).any()]) for _,g in groups])
target['IN'] = india(2014)
with zipfile.ZipFile(data/'perth.bin') as z:
    file = next(n for n in z.namelist() if n.endswith('.xlsx'))
    flora = pd.read_excel(io.BytesIO(z.read(file)), sheet_name=1, header=None)
def corridor(col):
    sites = []
    for v in flora[col].dropna():
        if str(v).startswith('Site '):
            sites.append(set())
        elif sites and (key := eu.binomial(v)):
            sites[-1].add(key)
    assert len(sites)==3
    return frequencies(sites)
target['AU'], target['CN'] = corridor(0), corridor(7)

control = {}
metadata = {}
def add(c, kind, freq, **meta):
    control[c,kind] = freq
    metadata[c+'_'+kind] = meta

m = pd.read_csv(data/'miren.csv', low_memory=False)
m = m[(m.Plot=='1') & m.Protocol.isin(['v.2017','v.2022'])].copy()
m = m[m.Year==m.groupby(['Region','Road','Site']).Year.transform('min')]
for c, regions in {'CZ':['CZE'],'IT':['ITA'],'IN':['IND'],'AU':['AUN','AUV'],'CN':['CNN']}.items():
    sub = m[m.Region.isin(regions)]
    roads = [frequencies([names(g['Accepted.Name.MIREN']) for _,g in r.groupby('Site')])
             for _,r in sub.groupby(['Region','Road'])]
    add(c,'road', {t:np.mean([r.get(t,0) for r in roads]) for t in set().union(*roads)},
        roads=len(roads), plots=sub.groupby(['Region','Road','Site']).ngroups,
        years=sorted(map(int, sub.Year.unique())), lat=[sub.Latitude.min(),sub.Latitude.max()],
        lon=[sub.Longitude.min(),sub.Longitude.max()])

trees = pd.read_excel(data/'tree1.xlsx',header=None).iloc[4:].copy()
for c,cols in {'HR':[2,7],'BA':[3],'SI':[1]}.items():
    vals = defaultdict(list)
    for _,r in trees.iterrows():
        key = eu.binomial(r[0])
        if key:
            key = {'Ochlopoa annua':'Poa annua'}.get(key,key)
            vals[key].append(float(np.mean([r[k] for k in cols]))/100)
    add(c,'urban',{k:v[0] for k,v in vals.items() if len(v)==1 and v[0]>0},
        plots=len(cols)*50, cities=len(cols), ambiguous=sum(len(v)>1 for v in vals.values()))

# Table 1, 17 association-A columns then its printed % and one association-B column.
with pymupdf.open(data/'sk_urban.bin') as doc:
    passage = doc[35].get_text().split('Tab. 1. Relev')[1]
lines = [l.strip() for l in passage.splitlines() if l.strip()]
plots = [set() for _ in range(18)]
parsed = 0
for j,l in enumerate(lines[:-19]):
    cells = lines[j+1:j+20]
    if not (re.fullmatch(r'\d+', cells[17]) and all(re.fullmatch(r'[.+rab12345]',v) for v in cells[:17]+cells[18:])):
        continue
    assert abs(100*sum(v!='.' for v in cells[:17])/17-int(cells[17]))<=0.51
    parsed += 1
    if key := eu.binomial(l):
        for plot,v in zip(plots,cells[:17]+cells[18:]):
            if v!='.': plot.add(key)
single = passage.split('Taxa found only in one relev')[1].split('E1:')[1]
for taxon,num in re.findall(r'([A-Z][a-z]+\s+[a-z]+)\s*\((\d+)\)', ' '.join(single.split())):
    plots[int(num)-1].add(taxon)
assert parsed==24 and len(re.findall(r'\(\d+\)',single))==8, (parsed, single)
add('SK','urban',frequencies([p for i,p in enumerate(plots,1) if i!=15]),plots=17,source='Malacky Table 1; rail plot 15 excluded')

h = pd.read_csv(data/'headers.csv',usecols=['COUNTRY','HABITAT','LOCALITY','RELEVE_NR'])
def fold(s):
    return s.fillna('').str.normalize('NFKD').str.encode('ascii',errors='ignore').str.decode('ascii').str.lower()
hab,loc = fold(h.HABITAT),fold(h.LOCALITY)
urban = (hab+' '+loc).str.contains(r'\b(?:praha|brno|ostrava|plzen|olomouc)\b|mestsk|sidlist|intravil')
disturbed = hab.str.contains('ruder|rumist|navaz|staveni|skladk|bouran|seslap')
rail = (hab+' '+loc).str.contains('zeleznic|nadraz|kolej|railway')
ids = set(h.loc[(h.COUNTRY=='CZ') & urban & disturbed & ~rail,'RELEVE_NR'])
rows=[]
for chunk in pd.read_csv(data/'species.csv',chunksize=200000,usecols=['RELEVE_NR','PLANT_GROUP','KAPLAN_PLUS_NAME_SIMPLE']):
    rows.append(chunk[chunk.RELEVE_NR.isin(ids) & (chunk.PLANT_GROUP=='vascular')])
veg = pd.concat(rows)
vplots=[names(g.KAPLAN_PLUS_NAME_SIMPLE) for _,g in veg.groupby('RELEVE_NR')]
add('CZ','urban',frequencies(vplots), plots=len(vplots), header_matches=len(ids))

station_groups={'CZ':['CZ'],'IT':['IT'],'SK':['SK'],'AD':['HR','BA'],'IN':['IN']}
def score(a,b):
    x,y=rankdata(a),rankdata(b)
    if x.std()==0 or y.std()==0: return None
    return float(np.corrcoef(x,y)[0,1])

def evaluate(yr, rail_pred, domestic, panel, bootstrap):
    taxa=sorted(panel)
    a,b,y=np.array([[v.get(t,0) for t in taxa] for v in (rail_pred, domestic, yr)])
    r,s=score(a,y),score(b,y)
    out=dict(n=len(taxa), recorded=int(sum(y>0)), foreign=r, domestic=s,
             delta=None if r is None or s is None else r-s,
             mae_foreign=float(np.mean(abs(a-y))),mae_domestic=float(np.mean(abs(b-y))))
    if bootstrap and out['delta'] is not None:
        groups=defaultdict(list)
        for i,t in enumerate(taxa): groups[t.split()[0]].append(i)
        groups=list(groups.values())
        draws=[]
        rng=np.random.default_rng(20261117)
        for _ in range(999):
            idx=np.concatenate([groups[i] for i in rng.integers(len(groups),size=len(groups))])
            rr,ss=score(a[idx],y[idx]),score(b[idx],y[idx])
            if rr is not None and ss is not None: draws.append(rr-ss)
        if draws: out['genus_interval']=np.quantile(draws,[.025,.975]).tolist()
    return out

results={}
candidates=defaultdict(dict)
rank_advantage=defaultdict(dict)
for c,yr in target.items():
    train={g:(studies['HR_BA15'] if g=='AD' else target[cc[0]]) for g,cc in station_groups.items() if c not in cc}
    support=Counter(t for v in train.values() for t in v)
    panel={t for t,n in support.items() if n>=2}
    pred={t:float(np.mean([v.get(t,0) for v in train.values()])) for t in panel}
    for (country,kind),bg in control.items():
        if c!=country: continue
        common=panel & bg.keys()
        result=dict(training=list(train), common=evaluate(yr,pred,bg,common,True),
                    extended=evaluate(yr,pred,bg,panel,False))
        if c=='CZ':
            both=panel & control['CZ','urban'].keys() & control['CZ','road'].keys()
            result['same_species_both_domestic_controls']=evaluate(yr,pred,bg,both,True)
        reduced=[]
        for omitted in train:
            smaller=[v for g,v in train.items() if g!=omitted]
            reduced_pred={t:float(np.mean([v.get(t,0) for v in smaller])) for t in common}
            reduced.append(evaluate(yr,reduced_pred,bg,common,False)['delta'])
        result['delta_excluding_training_study_fixed_panel']=dict(zip(train,reduced))
        aligned=sorted(common)
        # Percentile differences put whole stations and small plots on a relative scale.
        ranks=(rankdata([yr.get(t,0) for t in aligned])-rankdata([bg[t] for t in aligned]))/len(aligned)
        rank_lookup=dict(zip(aligned,ranks))
        if c in ['IN','CN','AU']:
            repeat=india(2017) if c=='IN' else corridor(2 if c=='AU' else 9)
            result['other_year']=evaluate(repeat,pred,bg,common,False)
        results[c+'_'+kind]=result
        for t in common:
            if pred[t]>=.2 and support[t]>=3:
                rank_advantage[t].setdefault(c,{})[kind]=float(rank_lookup[t])
                candidates[t].setdefault(c,{})[kind]=[yr.get(t,0),bg[t]]
nsite={'CZ':39,'IT':28,'SK':30,'HR':13,'BA':2,'IN':16,'SI':89,'CN':3,'AU':3}
country_candidates={}
for t,row in sorted(candidates.items()):
    eligible=[c for c,v in row.items() if c in ['CZ','IT','SK','HR','IN'] and
              len(v)==sum(cc==c for cc,_ in control)]
    hits=[c for c in eligible if all(a>=.2 and a-b>=.1-1e-10 for a,b in row[c].values())]
    country_candidates[t]=dict(eligible=eligible,hits=hits,frequencies=row,
        relative_hits=[c for c in hits if all(d>=.10 for d in rank_advantage[t][c].values())],
        rank_differences=rank_advantage[t])
files=['miren.csv','perth.bin','tree1.xlsx','headers.csv','species.csv','sk_urban.bin','Matrix_plot.xlsx']
def digest(path):
    sha = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024*1024), b''):
            sha.update(block)
    return sha.hexdigest()
hashes={f:digest(data/f) for f in files}
output.write_text(json.dumps(dict(sites=nsite,controls=metadata,comparisons=results,species=country_candidates,
                                sha256=hashes),ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8')
