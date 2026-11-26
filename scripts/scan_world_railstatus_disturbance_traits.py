#!/usr/bin/env python3
"""Species-equal railway landscape trait summaries from authentic sources.
Only aggregated cells; exact original full and observer-half counts replay.
All source-known denominators explicit, unknown never assigned a trait.
"""
import argparse,csv,io,json,hashlib,math,shlex,zipfile
from collections import Counter,defaultdict
from pathlib import Path
from audit_railstatus_disturbance_trait_sources import (
    norm,request,figshare_metadata,get_figshare,csv_table,LEP_ROOT,LEP_PATH,LEP_BLOB,git_blob_sha,DERIVATIVE,SOURCE_SHA)
from scan_world_railstatus_cross_taxa import (
    taxonomy,source_client,open_original_source,source_coarse_cell,OBS_ETAG,OBS_BYTES,
    TAXA_ETAG,TAXA_BYTES,MAX_UNCERT_M,PROGRESS_EVERY_SOURCE_ROWS)
from scan_fullworld_2025_iNat_plant_macro_grid_v369 import source_ancestry_ids
from scan_world_railstatus_observer_split_beta import bucket,add,full_counts,split_counts,SALT,MIN_FULL,MIN_SPLIT
from analyze_world_railstatus_plant_bird_disjoint import extract
GROUPS=("plants","insects","birds")
WIDTHS=(100000,200000)
LIFE_SHA="83b23e7874d73387a0c24aeae3e2931c50ec7c78cad52ff26cdfc99d1b6623b1"
LEP_SHA="6ec35b8a31e96c971aeaa228a48aae9f107c40c33695f0d470aa4382ca6d635b"
SIZE_SHA="e1684d5d97dfc99980a7fb29391ef32467f093e80fc03dcb7a95363dcba2d3f0"
TRAITS={"plants":("annual","annual_nonwoody","woody","height_nonwoody","seed_mass_nonwoody","stem_density_nonwoody"),
        "insects":("disturbance_affinity","open_canopy","edge_affinity","multivoltine","log_wingspan_upper"),
        "birds":("grassland_vs_forest","HWI","log_wing_length","log_mass")}
MISSING={"","na","n/a","nan","null"}
def number(x):
    try:v=float(x)
    except (TypeError,ValueError):return None
    return v if math.isfinite(v) else None
def binary(v,pos,neg):
    z=norm(v)
    return 1. if z in pos else 0. if z in neg else None

def taxonomy_metadata(client,canonical,sha):
    names={g:defaultdict(list) for g in GROUPS};families={};ancestry={};roots=defaultdict(list)
    rr,raw,stream=open_original_source(client,"taxa.csv.gz",TAXA_ETAG,TAXA_BYTES)
    try:
        for r in rr:
            t=int(r["taxon_id"]);rank=r["rank"].strip().lower();n=norm(r["name"])
            active=r["active"].strip().lower() in ("t","true","1")
            if active and n in ("papilionoidea","tracheophyta"):roots[n].append(t)
            if rank=="family" and active:families[t]=n
            hit=canonical.get(t)
            if hit and hit[0] in names and hit[1]==t and rank=="species":
                names[hit[0]][n].append(t);ancestry[t]=(hit[0],tuple(map(int,source_ancestry_ids(r["ancestry"]))))
    finally:stream.close()
    if raw.bytes!=TAXA_BYTES or raw.hash.hexdigest()!=sha:raise ValueError("Exact taxonomy replay failed")
    if any(len(roots[k])!=1 for k in ("papilionoidea","tracheophyta")):raise ValueError("Ambiguous taxonomic target root")
    allowed={"plants":set(),"insects":set(),"birds":set()};fam={}
    for t,(g,anc) in ancestry.items():
        if g=="birds" or roots["tracheophyta" if g=="plants" else "papilionoidea"][0] in anc:
            allowed[g].add(t)
            f=next((x for x in reversed(anc) if x in families),None)
            if f is not None:fam[t]=f
    if len(allowed["plants"])<100000 or len(allowed["insects"])<10000 or len(allowed["birds"])<8000:raise ValueError("Taxonomic target pool incomplete")
    return names,allowed,fam,{k:v[0] for k,v in roots.items()}

def exact(rows,col,names,allowed):
    by=defaultdict(list)
    for r in rows:by[norm(r[col])].append(r)
    out={}
    for n,rr in by.items():
        ids=names.get(n,[])
        if len(rr)==1 and len(ids)==1 and ids[0] in allowed:out[ids[0]]=rr[0]
    return out

def fetch_fig(article,file,expected):
    q=figshare_metadata(article);fs=[f for f in q["files"] if f["id"]==file]
    if len(fs)!=1:raise ValueError("Missing verified file")
    raw=get_figshare(fs[0],30_000_000)
    if hashlib.sha256(raw).hexdigest()!=expected:raise ValueError("Frozen Figshare bytes changed")
    return raw

def load_life(names,allowed):
    raw=fetch_fig(20764777,37005529,LIFE_SHA);cols,rows=csv_table(raw)
    if cols!=["Family","Genus","Accepted Name","LifeCycle","Woodiness"] or len(rows)!=257401:raise ValueError("LifeForm schema changed")
    if {norm(r["LifeCycle"]) for r in rows}-MISSING-{"annual","perennial"} or {norm(r["Woodiness"]) for r in rows}-MISSING-{"non-woody","woody"}:raise ValueError("Unknown life form code")
    matched=exact(rows,"Accepted Name",names,allowed);out={}
    for t,r in matched.items():
        a=binary(r["LifeCycle"],{"annual"},{"perennial"});w=binary(r["Woodiness"],{"woody"},{"non-woody"})
        if a==1. and w==1.:continue # biologically contradictory source row; no forced label
        out[t]={"annual":a,"woody":w,"annual_nonwoody":a if w==0 else None}
    return out

def load_size(names,allowed,life):
    raw=fetch_fig(25600650,45643398,SIZE_SHA)
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        txt=z.read("Data/AllTraits.txt").decode("utf-8-sig")
    lines=txt.splitlines();cols=shlex.split(lines[0]);expected=["la","ln","ph","sla","ssd","sm","SRL","D","RTD","N","Dr","Lr"]
    if cols!=expected:raise ValueError("Size trait source columns changed")
    rows=[]
    for line in lines[1:]:
        v=shlex.split(line)
        if len(v)!=len(cols)+1:raise ValueError("Size trait row or rowname mismatch")
        rows.append(dict(zip(["species"]+cols,v)))
    if len(rows)!=39334 or len({norm(r["species"].replace("_"," ")) for r in rows})!=len(rows):raise ValueError("Original size species count mismatch")
    # Published source completeness, not a filled trait table.
    N={k:sum(number(r[k]) is not None for r in rows) for k in ("ph","sm","ssd")}
    joint=sum(number(r["ph"]) is not None and number(r["sm"]) is not None for r in rows)
    if N["ph"]!=13438 or N["sm"]!=23003 or joint!=7217:raise ValueError("Measured height/seed completeness did not reproduce publication")
    for k in ("ph","sm","ssd"):
        vv=[number(r[k]) for r in rows if number(r[k]) is not None]
        if min(vv)<=0:raise ValueError("Size file may already be transformed; review native units before fitting")
    if max(number(r["ph"]) for r in rows if number(r["ph"]) is not None)>150:raise ValueError("Height not plausible in meters")
    for r in rows:r["species"]=r["species"].replace("_"," ")
    matched=exact(rows,"species",names,allowed)
    for t,r in matched.items():
        if t not in life or life[t]["woody"]!=0:continue
        for k,f in (("ph","height_nonwoody"),("sm","seed_mass_nonwoody"),("ssd","stem_density_nonwoody")):
            v=number(r[k])
            if v is not None:life[t][f]=math.log(v)
    return {"status":"verified measured-value support","source_n":len(rows),"known":N,"height_seed_joint":joint,"exact_matches":len(matched)}

def load_lep(names,allowed):
    raw=request(LEP_ROOT+LEP_PATH,cap=8_000_000)
    if hashlib.sha256(raw).hexdigest()!=LEP_SHA or git_blob_sha(raw)!=LEP_BLOB:raise ValueError("Frozen LepTraits changed")
    cols,rows=csv_table(raw)
    if len(rows)!=12448 or not {"Species","DisturbanceAffinity","CanopyAffinity","EdgeAffinity","Voltinism","WS_U"}<=set(cols):raise ValueError("LepTraits schema changed")
    matched=exact(rows,"Species",names,allowed);out={}
    for t,r in matched.items():
        d=binary(r["DisturbanceAffinity"],{"disturbance-associated (strong)","disturbance-associated (weak)"},{"disturbance-avoidant (strong)","disturbance-avoidant (weak)"})
        c=binary(r["CanopyAffinity"],{"open canopy"},{"closed canopy","closed canopy (+edge)"})
        e=binary(r["EdgeAffinity"],{"edge-associated (strong)","edge-associated (weak)"},{"edge-avoidant (strong)","edge-avoidant (weak)"})
        v=binary(r["Voltinism"],{"m"},{"u","b"})
        wing=number(r["WS_U"])
        out[t]={"disturbance_affinity":d,"open_canopy":c,"edge_affinity":e,"multivoltine":v,
                "log_wingspan_upper":math.log(wing) if wing is not None and 0<wing<=31 else None}
    return out

def load_birds(names,allowed):
    raw=request(DERIVATIVE,cap=10_000_000)
    if hashlib.sha256(raw).hexdigest()!=SOURCE_SHA:raise ValueError("Frozen AVONET changed")
    cols,rows=csv_table(raw)
    if len(rows)!=11009 or not {"Species1","Hand-Wing.Index","Wing.Length","Mass","Habitat"}<=set(cols):raise ValueError("AVONET schema changed")
    matched=exact(rows,"Species1",names,allowed);out={}
    if len(matched)!=9872:raise ValueError("Prior bird exact match not reproduced")
    for t,r in matched.items():
        h,w,m=map(number,(r["Hand-Wing.Index"],r["Wing.Length"],r["Mass"]))
        out[t]={"grassland_vs_forest":binary(r["Habitat"],{"grassland"},{"forest"}),
                "HWI":h if h is not None and 0<=h<=100 else None,
                "log_wing_length":math.log(w) if w is not None and w>=10 else None,
                "log_mass":math.log(m) if m is not None and m>0 else None}
    return out

def family_refs(traits,family,group):
    d=defaultdict(lambda:[0,0.])
    for t,r in traits.items():
        if t not in family:continue
        for f in TRAITS[group]:
            v=r.get(f)
            if v is not None:
                a=d[f,family[t]];a[0]+=1;a[1]+=v
    return d

def summarize(species,traits,family,ref,group):
    known_nonwoody={t for t in species if traits.get(t,{}).get("woody")==0.}
    out={}
    for f in TRAITS[group]:
        den=len(known_nonwoody) if f.endswith("_nonwoody") else len(species)
        vv=[(t,traits[t][f]) for t in species if t in traits and traits[t].get(f) is not None]
        n=len(vv)
        if n>den:raise ValueError("Trait numerator outside eligible denominator")
        s=sum(v for t,v in vv);sq=sum(v*v for t,v in vv)
        obs=ex=0.;fn=0
        for t,v in vv:
            N,S=ref.get((f,family.get(t)),(0,0.))
            if N>=11:fn+=1;obs+=v;ex+=(S-v)/(N-1)
        out[f]={"denominator":den,"known":n,"mean":s/n if n else None,
                 "variance":max(0.,sq/n-(s/n)**2) if n else None,
                 "family_known":fn,"family_observed":obs/fn if fn else None,
                 "family_expected":ex/fn if fn else None,"family_excess":(obs-ex)/fn if fn else None}
    return out

def scan(client,canonical,allowed,candidates):
    from pyproj import Transformer
    tf=Transformer.from_crs(4326,6933,always_xy=True)
    rr,raw,stream=open_original_source(client,"observations.csv.gz",OBS_ETAG,OBS_BYTES)
    if not {"observer_id","latitude","longitude","positional_accuracy","taxon_id","quality_grade","observed_on"}<=set(rr.fieldnames or []):raise ValueError("Observation schema changed")
    data={w:{g:defaultdict(bucket) for g in GROUPS} for w in WIDTHS}
    target={w:{g:defaultdict(bucket) for g in GROUPS} for w in WIDTHS};qa=Counter()
    try:
        for r in rr:
            qa["rows"]+=1
            if qa["rows"]%PROGRESS_EVERY_SOURCE_ROWS==0:print(json.dumps({"trait_scan_progress":dict(qa)}),flush=True)
            if not r["observed_on"].startswith("2025-") or r["quality_grade"].strip().lower().replace(" ","_").replace("-","_")!="research":continue
            try:hit=canonical.get(int(r["taxon_id"]));u=float(r["positional_accuracy"]);lon=float(r["longitude"]);lat=float(r["latitude"])
            except (TypeError,ValueError):continue
            if hit is None or hit[0] not in GROUPS or not math.isfinite(u) or not 0<=u<=MAX_UNCERT_M:continue
            g,t=hit;k1=source_coarse_cell(lon,lat,tf,100000);k2=source_coarse_cell(lon,lat,tf,200000)
            if k1 is None or k2 is None:continue
            qa[g]+=1
            o=int(r["observer_id"]) if r["observer_id"] else None
            for w,k in ((100000,k1),(200000,k2)):
                if k not in candidates[w][g]:continue
                add(data[w][g][k],t,o)
                if t in allowed[g]:add(target[w][g][k],t,o)
    finally:stream.close()
    if raw.bytes!=OBS_BYTES or qa["rows"]!=288721938:raise ValueError("Incomplete original source")
    return data,target,dict(qa),raw.hash.hexdigest()

def main():
    p=argparse.ArgumentParser()
    for k in ("taxa","split","audit","out"):p.add_argument("--"+k,type=Path,required=True)
    a=p.parse_args();prior=json.loads(a.taxa.read_text());split=json.loads(a.split.read_text())
    if split["version"]!="railstatus-four-taxon-observer-split-beta-v1" or not split["exact_grid_replay"] or not split["exact_plant_beta_replay"] or split["split_seed"]!=SALT.decode() or split["split_gate_observations_observers_species"]!=list(MIN_SPLIT) or split["full_gate_observations_observers_species"]!=list(MIN_FULL):raise ValueError("Wrong split source")
    if any(prior[k]!=split[k] for k in ("obs_sha256","taxa_sha256")):raise ValueError("Wrong four-taxon source")
    client=source_client();canonical,_,sha=taxonomy(client)
    if sha!=split["taxa_sha256"]:raise ValueError("Taxonomy changed")
    names,allowed,fam,roots=taxonomy_metadata(client,canonical,sha)
    traits={"plants":load_life(names["plants"],allowed["plants"]),"insects":load_lep(names["insects"],allowed["insects"]),"birds":load_birds(names["birds"],allowed["birds"])}
    try:size=load_size(names["plants"],allowed["plants"],traits["plants"])
    except Exception as e:size={"status":"blocked","type":type(e).__name__,"message":str(e)[:250]}
    audit={"source_SHA256":{"LifeForm":LIFE_SHA,"LepTraits":LEP_SHA,"AVONET":SOURCE_SHA,"size_archive":SIZE_SHA},"target_roots":roots,
           "taxonomy_sha256":sha,"size_source":size,
           "pool_counts":{g:{"eligible_taxa":len(allowed[g]),"matched_taxa":len(traits[g]),"known":{f:sum(r.get(f) is not None for r in traits[g].values()) for f in TRAITS[g]}} for g in GROUPS}}
    if audit["pool_counts"]["plants"]["known"]["annual"]<100000 or audit["pool_counts"]["insects"]["known"]["disturbance_affinity"]<1000:raise ValueError("Inadequate exact life/disruption matching")
    a.audit.parent.mkdir(parents=True,exist_ok=True);a.audit.write_text(json.dumps(audit,allow_nan=False));print(json.dumps({"trait_scan_audit":audit},allow_nan=False),flush=True)
    candidates={w:{g:extract(split["pairs"][str(w)][g]) for g in GROUPS} for w in WIDTHS}
    for w in WIDTHS:
        for g in GROUPS:
            old={(r["source_equal_area_grid_ix"],r["source_equal_area_grid_iy"]):(r["observations"],r["observers"],r["species"]) for r in prior["grids"][str(w)][g]}
            if any(old.get(k)!=v[0] for k,v in candidates[w][g].items()):raise ValueError("Frozen full grid differs")
    refs={g:family_refs(traits[g],fam,g) for g in GROUPS}
    del names
    data,targets,qa,osha=scan(client,canonical,allowed,candidates)
    if osha!=split["obs_sha256"] or any(qa[g]!=prior["qa"][g] for g in GROUPS):raise ValueError("Original census not reproduced")
    grids={}
    for w in WIDTHS:
        grids[str(w)]={}
        for g in GROUPS:
            if set(data[w][g])!=set(candidates[w][g]):raise ValueError("Exact cell support mismatch")
            rr=[]
            for k,v in sorted(data[w][g].items()):
                f=full_counts(v);h=tuple(tuple(x) for x in split_counts(v))
                if (f,h)!=candidates[w][g][k]:raise ValueError("Exact full and split count mismatch")
                t=targets[w][g].get(k,bucket())
                sets=[t[1],*t[4]];cnt=[full_counts(t),*split_counts(t)]
                if any(any(a>b for a,b in zip(q,p)) for q,p in zip(cnt,(f,*h))):raise ValueError("Target counts exceed full")
                rows=[{"counts":c,"traits":summarize(s,traits[g],fam,refs[g],g)} for s,c in zip(sets,cnt)]
                rr.append({"key":k,"full":f,"halves":h,"target_full":rows[0],"target_halves":rows[1:]})
            grids[str(w)][g]=rr
            print(json.dumps({"trait_cell_panel":{"width":w,"group":g,"exact_cells":len(rr),
                   "coverage50_full":{f:sum(r["target_full"]["traits"][f]["known"]>=20 and r["target_full"]["traits"][f]["known"]>=.5*r["target_full"]["traits"][f]["denominator"] for r in rr) for f in TRAITS[g]}}}),flush=True)
    out={"version":"railstatus-disturbance-trait-cells-v1","source":audit,"obs_sha256":osha,"taxa_sha256":sha,"split_seed":SALT.decode(),
         "exact_full_and_split_replay":True,"qa":qa,"grids":grids}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,allow_nan=False))
if __name__=="__main__":main()
