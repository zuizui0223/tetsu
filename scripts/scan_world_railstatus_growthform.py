#!/usr/bin/env python3
"""Pinned public growthform lookup -> exact iNat 2025 source replay.
Woodiness is observed source W/H, never inferred from missing values.
Output aggregate cell counts and within-family external baseline only.
"""
import argparse,csv,hashlib,io,json,math
from collections import Counter,defaultdict
from pathlib import Path
from audit_world_NHM_rail_trait_match import norm,request,active_names,scan
from scan_fullworld_2025_iNat_plant_macro_grid_v369 import source_client,taxonomy_source
SHA="6a245a72ee48d658c255fce6a25f662156ad07e3"
URL="https://raw.githubusercontent.com/traitecoevo/growthform/"+SHA+"/growth_form.csv"
W=100000
SUPPORT={"F","C","E","P","M","A","H"}
WOOD={"W","H"}

def source_rows(raw):
    rr=csv.DictReader(io.StringIO(raw.decode("utf-8-sig"),newline=""))
    if tuple(rr.fieldnames or ())!=("sp","support","woodiness","nfixing"):raise ValueError("Publisher schema changed")
    rows=list(rr)
    if len(rows)<140000 or any(not norm(r["sp"]) for r in rows):raise ValueError("Incomplete published growthform table")
    return rows

def match(rows,names):
    grouped=defaultdict(list)
    for r in rows:grouped[norm(r["sp"])].append(r)
    out={};qa=Counter()
    for name,rr in grouped.items():
        if len(rr)!=1:qa["ambiguous_source_names"]+=1;continue
        ids=names.get(name,[])
        if len(ids)!=1:qa["unmatched_or_ambiguous_iNat"]+=1;continue
        r=rr[0];wo=r["woodiness"].strip();su=r["support"].strip()
        if wo not in WOOD and wo not in ("","NA"):qa["other_woodiness_values"]+=1
        if su not in SUPPORT and su not in ("","NA"):qa["other_support_values"]+=1
        out[ids[0]]={"woodiness":wo if wo in WOOD else None,"support":su if su in SUPPORT else None}
    qa.update(source_rows=len(rows),distinct_source_names=len(grouped),matched_species=len(out),
              known_woodiness=sum(r["woodiness"] is not None for r in out.values()),
              woody=sum(r["woodiness"]=="W" for r in out.values()),
              herbaceous=sum(r["woodiness"]=="H" for r in out.values()),
              known_support=sum(r["support"] is not None for r in out.values()))
    return out,dict(qa)

def family_reference(joined,family):
    # External source taxonomy-matched species pool (not occurrence selected).
    # No imputation of unknown woodiness, no model trained on railway effects.
    counts=defaultdict(lambda:[0,0])
    for tid,r in joined.items():
        f=family.get(tid)
        if f is None or r["woodiness"] is None:continue
        a=counts[f];a[0]+=1;a[1]+=r["woodiness"]=="W"
    return {f:tuple(v) for f,v in counts.items()}

def aggregate(species,joined,family,fref):
    known=[t for t in species if t in joined and joined[t]["woodiness"] is not None]
    woody=sum(joined[t]["woodiness"]=="W" for t in known)
    ns=Counter(joined[t]["support"] for t in species if t in joined and joined[t]["support"] is not None)
    # For each focal species, remove itself from its external family pool.
    # Restrict both observed and expected to same known species with >=10
    # other source-known family members.
    obs=exp=0.;n=0
    for t in known:
        f=family.get(t);N,K=fref.get(f,(0,0))
        if N<11:continue
        w=joined[t]["woodiness"]=="W"
        n+=1;obs+=w;exp+=(K-w)/(N-1)
    return {"species":len(species),"matched":sum(t in joined for t in species),
            "known_woodiness":len(known),"woody":woody,"herbaceous":len(known)-woody,
            "known_support":sum(ns.values()),"support_counts":dict(ns),
            "family_reference_known":n,"family_observed_woody":obs,
            "family_expected_woody":exp}

def main():
    p=argparse.ArgumentParser()
    for n in ("plant","out"):p.add_argument("--"+n,type=Path,required=True)
    a=p.parse_args();plant=json.loads(a.plant.read_text())
    raw=request(URL,cap=25_000_000);rows=source_rows(raw)
    client=source_client();canonical,_,family,_,sha=taxonomy_source(client)
    if sha!=plant["source_summary"]["source_iNat_original_2026_10_27_taxa_gzip_SHA256_full"]:raise ValueError("Frozen taxonomy mismatch")
    names=active_names(client,canonical,sha)
    joined,qa=match(rows,names)
    if qa["known_woodiness"]<10000:raise ValueError("Insufficient source-known woodiness matches")
    source={"repository":"traitecoevo/growthform","commit":SHA,"url":URL,"source_bytes":len(raw),"source_SHA256":hashlib.sha256(raw).hexdigest(),"taxonomy_SHA256":sha,"qa":qa}
    print(json.dumps({"growthform_source":source}),flush=True)
    fref=family_reference(joined,family)
    del rows,raw,names
    data,replay=scan(client,canonical,plant)
    allspecies=set().union(*(v[1] for v in data.values()))
    total=aggregate(allspecies,joined,family,fref)
    cells=[]
    for (i,j),(n,s,o) in sorted(data.items()):
        if n<50 or len(o)<5 or len(s)<20:continue
        cells.append({"key":[i,j],"observations":n,"observers":len(o),**aggregate(s,joined,family,fref)})
    gates={}
    for key in ("known_woodiness","family_reference_known","known_support"):
        gates[key]={str(g):sum(r[key]>=20 and r[key]/r["species"]>=g for r in cells) for g in (.25,.5,.8)}
    out={"version":"railstatus-growthform-2025-coarse-v1","source":source,"source_replay":replay,
         "overall":total,"n_cells":len(cells),"coverage":gates,"cells":cells,
         "limits":["exact names only; no synonym resolution","source W/H annotations, unknown not herbaceous","family reference global external database, not regional available flora or phylogenetic null","known species support conditions on sampling","whole landscape, not rail microhabitat"]}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,allow_nan=False))
    print(json.dumps({k:v for k,v in out.items() if k!="cells"}),flush=True)
if __name__=="__main__":main()
