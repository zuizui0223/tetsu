"""python global_corridors.py TETSU_REPO SOURCE_DIR OUTPUT_JSON

Use tetsu commit 43884d11, Zenodo MIREN v3 record 15749536 as miren.csv,
and v24.5 university LEDA files canopy_height.txt and seed_mass.txt.
"""
import hashlib
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.spatial.distance import pdist

repo, data, out = map(Path, sys.argv[1:])
sys.path.insert(0, str(repo / "scripts"))
import analyze_global_MIREN_road_interior_LEDA_2D_functional_beta_v249 as model

records, source = model.read_source_miren(data / "miren.csv")
for file, trait in [("canopy_height.txt","height"),("seed_mass.txt","seed")]:
    assert hashlib.sha256((data/file).read_bytes()).hexdigest() == model.LEDASHA[trait]
h, _ = model.read_leda(data/"canopy_height.txt","height_m")
s, _ = model.read_leda(data/"seed_mass.txt","seed_mass_mg")
traits, coverage = model.traits_from_source(records,h,s)
surveys = defaultdict(list)
for reg,road,site,year,z in sorted(records):
    if z != 1:
        continue
    triple = [records.get((reg,road,site,year,i)) for i in (1,2,3)]
    if all(triple):
        surveys[reg,road,year].append(triple)

def select(mode):
    groups = defaultdict(list)
    for (reg,road,year), rows in surveys.items():
        keep = []
        for zones in rows:
            if mode == "all_binomials":
                keep.append(zones)
            else:
                known = [z.intersection(traits) for z in zones]
                n, frac = (10,.5) if mode == "strict_traits" else (5,.25)
                if all(len(k)>=n and len(k)>=frac*len(z) for k,z in zip(known,zones)):
                    keep.append(zones if mode=="full_taxa_on_trait_sites" else known)
        if len(keep)>=4:
            groups[reg,road].append((year,keep))
    regions = defaultdict(list)
    for (reg,road), items in sorted(groups.items()):
        year, rows = min(items,key=lambda item:(-len(item[1]),item[0]))
        regions[reg].append((road,year,rows))
    return {r:g for r,g in regions.items() if len(g)>=2}

def road_metrics(rows,functional,rng):
    species = sorted(set.union(*(z for row in rows for z in row)))
    index = {t:i for i,t in enumerate(species)}
    axes = np.array([traits[t] for t in species]) if functional else None
    cells = [[np.array(sorted(index[t] for t in z)) for z in row] for row in rows]
    minima = [min(map(len,row)) for row in cells]
    def distances(rarefy):
        zone_metrics = []
        for zone in range(3):
            matrix = np.zeros((len(rows),len(species)),dtype=bool)
            for i,row in enumerate(cells):
                ids = row[zone]
                if rarefy and len(ids)>minima[i]:
                    ids=rng.choice(ids,minima[i],replace=False)
                matrix[i,ids]=True
            metrics=[pdist(matrix,"jaccard").mean()]
            if functional:
                metrics.append(pdist((matrix @ axes)/matrix.sum(axis=1)[:,None]).mean())
            zone_metrics.append(metrics)
        return np.array(zone_metrics)
    raw=distances(False)
    equal=np.mean([distances(True) for _ in range(99)],axis=0)
    return dict(raw=raw.tolist(),equal=equal.tolist())

result={}
for mode in ("all_binomials","traits","strict_traits","full_taxa_on_trait_sites"):
    rng=np.random.default_rng(170126)
    regions=select(mode)
    groups={}
    for reg, roads in regions.items():
        groups[reg]=[dict(road=road,year=y,sites=len(rows),**road_metrics(rows,mode in ("traits","strict_traits"),rng))
                     for road,y,rows in roads]
    summary={}
    for key in ("raw","equal"):
        regional=np.array([np.mean([road[key] for road in roads],axis=0) for roads in groups.values()])
        for k,metric in enumerate(["jaccard","cwm"][:regional.shape[-1]]):
            v=regional[:,:,k]
            for a,b in [(0,2),(0,1),(1,2)]:
                d=v[:,a]-v[:,b]
                summary[f"{key}_{metric}_{a+1}minus{b+1}"]={
                    "mean":d.mean(),"negative":int((d<0).sum()),"n":len(d),
                    "median":np.median(d),"range":[min(d),max(d)],
                    "delete_one_region":[min((sum(d)-d)/(len(d)-1)),max((sum(d)-d)/(len(d)-1))] if len(d)>1 else None}
            summary[f"{key}_{metric}_strictly_ordered_regions"]=int(((v[:,0]<v[:,1])&(v[:,1]<v[:,2])).sum())
    result[mode]=dict(regions=groups,summary=summary)
out.write_text(json.dumps(dict(source=source,trait_coverage=coverage,results=result),indent=2),encoding="utf8")
