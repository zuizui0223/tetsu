#!/usr/bin/env python3
"""Use author-documented log-transformed, scaled original size traits.
Do not log a second time, impute missing values or invert unknown scaling.
Replay frozen vascular-plant/life-form summaries and original split counts.
"""
import argparse,json,hashlib,math,io,zipfile,shlex
from pathlib import Path
from xml.etree import ElementTree as ET
import numpy as np
from scan_world_railstatus_disturbance_traits import (
    taxonomy,source_client,taxonomy_metadata,load_life,fetch_fig,exact,number,
    LIFE_SHA,SIZE_SHA,WIDTHS,GROUPS,TRAITS,extract,scan,bucket,full_counts,split_counts,
    summarize,family_refs,SALT)
AXES=("height_nonwoody","seed_mass_nonwoody","stem_density_nonwoody")
OLD=("annual","annual_nonwoody","woody")
EXPECTED={"ph":13438,"sm":23003,"ssd":11148}
def parse_scaled(raw):
    if hashlib.sha256(raw).hexdigest()!=SIZE_SHA:raise ValueError("Frozen source changed")
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        txt=z.read("Data/AllTraits.txt").decode("utf-8-sig")
        with zipfile.ZipFile(io.BytesIO(z.read("READ_ME.docx"))) as dz:
            root=ET.fromstring(dz.read("word/document.xml"))
            ns={"w":"http://schemas.openxmlformats.org/wordprocessingml/2006/main"}
            text=" ".join("".join(t.text or "" for t in p.findall(".//w:t",ns)) for p in root.findall(".//w:p",ns))
    if "AllTraits = Trait matrix used in the analysis" not in text or "Data are log-transformed and scale." not in text:raise ValueError("Author transformation documentation absent")
    lines=txt.splitlines();cols=shlex.split(lines[0])
    if cols!=["la","ln","ph","sla","ssd","sm","SRL","D","RTD","N","Dr","Lr"]:raise ValueError("Source schema changed")
    rows=[]
    for line in lines[1:]:
        v=shlex.split(line)
        if len(v)!=len(cols)+1:raise ValueError("Source row changed")
        r=dict(zip(["species"]+cols,v));r["species"]=r["species"].replace("_"," ");rows.append(r)
    if len(rows)!=39334:raise ValueError("Source row count changed")
    audit={}
    for k,n in EXPECTED.items():
        x=np.asarray([v for r in rows if (v:=number(r[k])) is not None])
        if len(x)!=n or np.any(np.abs(x)>20) or abs(x.mean())>.001 or abs(x.std(ddof=1)-1)>.001:raise ValueError("Author-scaled source distribution differs "+k)
        audit[k]={"known":len(x),"mean":float(x.mean()),"sample_sd":float(x.std(ddof=1))}
    if sum(number(r["ph"]) is not None and number(r["sm"]) is not None for r in rows)!=7217:raise ValueError("Observed joint support differs")
    return rows,{"n_rows":len(rows),"traits":audit,"height_seed_joint":7217,"transformation":"author log transformed and scaled, retained directly","source_sha256":SIZE_SHA}
def add_size(life,rows,names,allowed):
    matched=exact(rows,"species",names,allowed)
    for t,r in matched.items():
        if t not in life or life[t].get("woody")!=0.:continue
        for k,f in zip(("ph","sm","ssd"),AXES):
            v=number(r[k])
            if v is not None:life[t][f]=v
    return {"exact_matches":len(matched),"known_nonwoody":{f:sum(r.get(f) is not None for r in life.values()) for f in AXES}}
def check_old(new,old):
    for f in OLD:
        for k,v in new[f].items():
            q=old[f][k]
            if v is None or q is None:
                if v!=q:raise ValueError("Old missingness differs")
            elif abs(v-q)>1e-8:raise ValueError("Prior plant trait moment differs")
def main():
    p=argparse.ArgumentParser()
    for k in ("taxa","split","previous_cells","audit","out"):p.add_argument("--"+k,type=Path,required=True)
    a=p.parse_args();prior,split,old=map(lambda p:json.loads(p.read_text()),(a.taxa,a.split,a.previous_cells))
    if old["version"]!="railstatus-disturbance-trait-cells-v1" or not old["exact_full_and_split_replay"] or any(len({x[k] for x in (prior,split,old)})!=1 for k in ("obs_sha256","taxa_sha256")) or old["split_seed"]!=SALT.decode() or old["source"]["source_SHA256"]["LifeForm"]!=LIFE_SHA or old["source"]["source_SHA256"]["size_archive"]!=SIZE_SHA:raise ValueError("Wrong prior source")
    client=source_client();canonical,_,sha=taxonomy(client)
    if sha!=split["taxa_sha256"]:raise ValueError("Taxonomy changed")
    names,allowed,fam,roots=taxonomy_metadata(client,canonical,sha)
    life=load_life(names["plants"],allowed["plants"])
    rows,audit=parse_scaled(fetch_fig(25600650,45643398,SIZE_SHA))
    audit.update(add_size(life,rows,names["plants"],allowed["plants"]))
    audit.update(taxonomy_sha256=sha,LifeForm_sha256=LIFE_SHA)
    a.audit.parent.mkdir(parents=True,exist_ok=True);a.audit.write_text(json.dumps(audit,allow_nan=False))
    print(json.dumps({"verified_plant_size_source":audit}),flush=True)
    if min(audit["known_nonwoody"].values())<100:raise ValueError("Insufficient exact nonwoody size matches")
    cand={w:{g:extract(split["pairs"][str(w)][g]) if g=="plants" else {} for g in GROUPS} for w in WIDTHS}
    oldmaps={w:{tuple(r["key"]):r for r in old["grids"][str(w)]["plants"]} for w in WIDTHS}
    for w in WIDTHS:
        if {k:(tuple(v["full"]),tuple(tuple(x) for x in v["halves"])) for k,v in oldmaps[w].items()}!=cand[w]["plants"]:raise ValueError("Frozen plant support changed")
    ref=family_refs(life,fam,"plants")
    del rows,names
    data,targets,qa,osha=scan(client,canonical,allowed,cand)
    if osha!=prior["obs_sha256"] or any(qa[g]!=prior["qa"][g] for g in GROUPS):raise ValueError("Original source not replayed")
    grids={}
    for w in WIDTHS:
        rr=[]
        if set(data[w]["plants"])!=set(cand[w]["plants"]):raise ValueError("Cell geography changed")
        for k,v in sorted(data[w]["plants"].items()):
            f=full_counts(v);h=tuple(tuple(x) for x in split_counts(v))
            if (f,h)!=cand[w]["plants"][k]:raise ValueError("Original split replay failed")
            t=targets[w]["plants"].get(k,bucket())
            cnt=[full_counts(t),*split_counts(t)];ss=[t[1],*t[4]]
            o=oldmaps[w][k];previous=[o["target_full"],*o["target_halves"]]
            qq=[]
            for s,c,pr in zip(ss,cnt,previous):
                if tuple(pr["counts"])!=tuple(c):raise ValueError("Prior target counts differ")
                z=summarize(s,life,fam,ref,"plants");check_old(z,pr["traits"])
                qq.append({"counts":c,"traits":z})
            rr.append({"key":k,"full":f,"halves":h,"target_full":qq[0],"target_halves":qq[1:]})
        grids[str(w)]={"plants":rr}
        print(json.dumps({"verified_size_cell_support":{"width":w,"exact_cells":len(rr),"coverage":{f:{str(g):sum(z["target_full"]["traits"][f]["known"]>=20 and z["target_full"]["traits"][f]["known"]>=g*z["target_full"]["traits"][f]["denominator"] for z in rr) for g in (.25,.5,.8)} for f in AXES}}}),flush=True)
    out={"version":"railstatus-verified-scaled-plant-size-cells-v1","source":audit,"source_SHA256":{"LifeForm":LIFE_SHA,"size_archive":SIZE_SHA},
        "obs_sha256":osha,"taxa_sha256":sha,"split_seed":SALT.decode(),"exact_full_split_target_and_life_replay":True,
        "prior_cells_sha256":hashlib.sha256(a.previous_cells.read_bytes()).hexdigest(),"grids":grids}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,allow_nan=False))
if __name__=="__main__":main()
