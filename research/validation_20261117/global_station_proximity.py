"""python global_station_proximity.py SOURCE_DIR RESULT_JSON

Official inputs: https://idata.idiv.de/ddm/Data/DownloadZip/3474?version=5047 (splot.zip)
and https://download.geonames.org/export/dump/allCountries.zip (geonames.zip).
Software: numpy, pandas, scipy. GeoNames CC BY 4.0; cite sPlotOpen and its contributing datasets.
"""
import hashlib
import json
import math
import sys
import zipfile
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.spatial.distance import pdist

data, out = map(Path, sys.argv[1:])
sys.path.insert(0,str(data/"tetsu"/"scripts"))
from audit_SI_CZ_original_continuous_trait_coverage_v220 import binomial,read_leda
h,_=read_leda(data/"canopy_height.txt","height_m")
s,_=read_leda(data/"seed_mass.txt","seed_mass_mg")
leda={t:np.array([(math.log(h[t])+.5802365754696401)/1.3801000545693773,
                  (math.log(s[t])-.2936113186962433)/2.1657660559995824]) for t in h.keys() & s.keys()}
def xyz(lat, lon):
    a,b=np.radians(lat),np.radians(lon)
    return np.column_stack((np.cos(a)*np.cos(b),np.cos(a)*np.sin(b),np.sin(a)))

def km(chord):
    return 12742*np.arcsin(np.minimum(1,chord/2))

def digest(path):
    with path.open("rb") as f:
        h=hashlib.sha256()
        for block in iter(lambda:f.read(1024*1024),b""):
            h.update(block)
        return h.hexdigest()

stations=[]
with zipfile.ZipFile(data/"geonames.zip") as z, z.open("allCountries.txt") as f:
    for line in f:
        col=line.split(b"\t",9)
        if col[7]==b"RSTN" and col[6]==b"S":
            stations.append((float(col[4]),float(col[5])))
station_xyz=np.unique(xyz(*np.array(stations).T),axis=0)
tree=cKDTree(station_xyz)
with zipfile.ZipFile(data/"splot.zip") as z:
    def table(name, **kw):
        return pd.read_csv(z.open(next(p for p in z.namelist() if p.startswith(name))),
                           sep="\t",low_memory=False,**kw)
    header=table("sPlotOpen_header")
    traits=table("sPlotOpen_CWM_CWV",usecols=["PlotObservationID","TraitCoverage_cover",
               "PlantHeight_CWM","SeedMass_CWM"])
    dt=table("sPlotOpen_DT",usecols=["PlotObservationID","Species"])
all_species=dt.groupby("PlotObservationID").Species.agg(lambda x:set(x.dropna()))
df=header.merge(traits,on="PlotObservationID",validate="one_to_one")
axes=["PlantHeight_CWM","SeedMass_CWM"]
df[axes]=(df[axes]-df[axes].mean())/df[axes].std(ddof=0)
df["year"]=pd.to_numeric(df.Date_of_recording.str[:4],errors="coerce")
df=df.dropna(subset=["Latitude","Longitude","year","Releve_area","is_forest","Biome","Dataset","Country"]+axes)
df=df.loc[(df.Releve_area>0)&(df.Location_uncertainty.between(0,250))&
          (df.TraitCoverage_cover>=.5)].reset_index(drop=True)
v=xyz(df.Latitude,df.Longitude)
chord,ids=tree.query(v)
df["nearest"]=km(chord)
df["station"]=ids
df["err"]=df.Location_uncertainty/1000

def match(band,strict):
    nearby=df[df.nearest+df.err<=band].sort_values(["nearest","PlotObservationID"])
    control=df[(df.nearest-df.err>=2*band)&(df.nearest+df.err<=10)]
    cols=["Country","Dataset","Biome","Plant_recorded","is_forest"]
    partitions=control.groupby(cols).indices
    pos=control.index.to_numpy()
    edges=[]
    for i in nearby.index:
        row=df.loc[i]
        ix=partitions.get(tuple(row[c] for c in cols))
        if ix is None: continue
        j=pos[ix]
        dist=km(np.linalg.norm(v[j]-station_xyz[int(row.station)],axis=1))
        area=df.Releve_area.values[j]/row.Releve_area
        ok=(dist+df.err.values[j]<=10)&(np.abs(df.year.values[j]-row.year)<=10)&(area>=.5)&(area<=2)
        if strict:
            ok &= (np.abs(df.Elevation.values[j]-row.Elevation)<=100)&(np.abs(
                df.year.values[j]-row.year)<=5)&(area>=.8)&(area<=1.25)&(row.year>=1990)&(df.year.values[j]>=1990)
        for c,d in zip(j[ok],dist[ok]):
            edges.append((row.nearest,d,int(i),int(c),int(row.station)))
    used=set(); used_stations_seen=set(); selected=[]; anchors=[]
    for _,__,i,c,st in sorted(edges):
        if c in used or st in used_stations_seen: continue
        if anchors and np.min(km(np.linalg.norm(np.array(anchors)-station_xyz[st],axis=1)))<2: continue
        used.add(c);used_stations_seen.add(st);anchors.append(station_xyz[st]);selected.append((i,c))
    return len(nearby),selected

def compare(pairs,rng):
    a,b=np.array(pairs).T
    f=df[axes].to_numpy()
    sets=[[all_species[df.PlotObservationID.iloc[j]] for j in side] for side in (a,b)]
    ids={s:k for k,s in enumerate(sorted(set.union(*(sets[0]+sets[1]))))}
    arrays=[[np.array(sorted(ids[s] for s in p)) for p in side] for side in sets]
    widths=np.minimum(*[[len(x) for x in side] for side in arrays])
    def tax(rare):
        means=[]
        for side in arrays:
            mat=np.zeros((len(side),len(ids)),bool)
            for j,sp in enumerate(side):
                mat[j,rng.choice(sp,widths[j],replace=False) if rare else sp]=True
            means.append(pdist(mat,"jaccard").mean())
        return means[0]-means[1]
    joint=[]
    for left,right in zip(*sets):
        raw=[{binomial(t) for t in side}-{None} for side in (left,right)]
        match=[sorted(t for t in side if t in leda) for side in raw]
        if all(len(m)>=5 and len(m)>=.25*len(v) for m,v in zip(match,raw)):
            joint.append([np.array([leda[t] for t in side]) for side in match])
    extra={"n":len(joint)}
    if len(joint)>=4:
        def leda_distance(rare):
            points=[[],[]]
            for sides in joint:
                n=min(map(len,sides))
                for k,side in enumerate(sides):
                    vals=side[rng.choice(len(side),n,replace=False)] if rare else side
                    points[k].append(vals.mean(axis=0))
            return pdist(points[0]).mean()-pdist(points[1]).mean()
        extra.update(raw=leda_distance(False),equal=np.mean([leda_distance(True) for _ in range(99)]))
    return {"stations":len(a),"leda":extra,"cwm":pdist(f[a]).mean()-pdist(f[b]).mean(),
            "jaccard":tax(False),"equal_jaccard":np.mean([tax(True) for _ in range(99)])}

result={}
for band,strict in [(.25,False),(1.0,False),(2.0,False),(1.0,True)]:
    n, pairs=match(band,strict)
    strata=defaultdict(list)
    for i,c in pairs:
        r=df.iloc[i]
        strata[r.Country,r.Continent,r.Dataset].append((i,c))
    rng=np.random.default_rng(171126)
    supported=[dict(country=co,continent=ct,dataset=ds,**compare(p,rng))
               for (co,ct,ds),p in sorted(strata.items()) if len(p)>=4]
    result[str(band)+("_strict" if strict else "")]=dict(candidate_plots=n,matched_stations=len(pairs),
                           matched_countries=len({key[0] for key in strata}),supported=supported)
out.write_text(json.dumps({"inputs":{p:digest(data/p) for p in ["splot.zip","geonames.zip"]},
    "source_plots":len(header),"source_station_records":len(stations),"distinct_station_coords":len(station_xyz),
    "eligible_plots":len(df),"bands_km":result},indent=2),encoding="utf8")
