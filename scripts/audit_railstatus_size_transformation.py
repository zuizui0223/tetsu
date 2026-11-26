#!/usr/bin/env python3
"""Check the pinned plant size archive before accepting transformed traits.
No trait association fitting and no individual species values exported.
"""
import hashlib,io,json,math,shlex,zipfile,re
from pathlib import Path
from xml.etree import ElementTree as ET
from audit_railstatus_disturbance_trait_sources import figshare_metadata,get_figshare
from scan_world_railstatus_disturbance_traits import SIZE_SHA,number
def main():
    q=figshare_metadata(25600650)
    fs=[f for f in q["files"] if f["id"]==45643398]
    if len(fs)!=1:raise ValueError("Archive identity missing")
    raw=get_figshare(fs[0],20000000)
    if hashlib.sha256(raw).hexdigest()!=SIZE_SHA:raise ValueError("Archive differs")
    out={"version":"railstatus-plant-size-transform-audit-v1","archive_sha256":SIZE_SHA,"tables":{},"code_excerpts":{}}
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        for name in ("Data/AllTraits.txt","Data/Herb_traits.txt","Data/Woody_traits.txt","Data/Traits_Woodyness.txt"):
            lines=z.read(name).decode("utf-8-sig").splitlines();cols=shlex.split(lines[0])
            rr=[]
            for line in lines[1:]:
                v=shlex.split(line)
                if len(v)==len(cols)+1:v=v[1:]
                if len(v)!=len(cols):raise ValueError("Unexpected row schema "+name)
                rr.append(dict(zip(cols,v)))
            a={"n":len(rr),"columns":cols,"traits":{}}
            for k in ("ph","sm","ssd"):
                if k not in cols:continue
                vv=sorted(x for r in rr if (x:=number(r[k])) is not None)
                a["traits"][k]={"n":len(vv),"nonpositive":sum(v<=0 for v in vv),"quantiles":[vv[round((len(vv)-1)*p)] for p in (0,.01,.25,.5,.75,.99,1)] if vv else []}
            a["height_seed_joint"]=sum(number(r.get("ph")) is not None and number(r.get("sm")) is not None for r in rr)
            out["tables"][name]=a
        rx=re.compile(r"AllTraits|Herb_traits|Woody_traits|Traits_Woodyness|log[ ]*[(]|log10[ ]*[(]|scale[ ]*[(]|transform|imput|standardiz|standardis|read[.]table|write[.]table|height|seed mass|stem specific|logarithm|normaliz|normalis",re.I)
        for name in sorted(x for x in z.namelist() if x.endswith(".R")):
            lines=z.read(name).decode("utf-8-sig",errors="replace").splitlines()
            ix={j for i,l in enumerate(lines) if rx.search(l) for j in range(max(0,i-1),min(len(lines),i+2))}
            out["code_excerpts"][name]=[{"line":i+1,"text":lines[i][:600]} for i in sorted(ix)][:100]
        with zipfile.ZipFile(io.BytesIO(z.read("READ_ME.docx"))) as dz:
            root=ET.fromstring(dz.read("word/document.xml"));ns={"w":"http://schemas.openxmlformats.org/wordprocessingml/2006/main"}
            paras=["".join(t.text or "" for t in p.findall(".//w:t",ns)) for p in root.findall(".//w:p",ns)]
            out["readme_relevant"]=[v for v in paras if rx.search(v)][:40]
    path=Path("results/railstatus/plant_size_transform_audit.json");path.parent.mkdir(parents=True,exist_ok=True);path.write_text(json.dumps(out,allow_nan=False))
    print(json.dumps({"size_tables":out["tables"]}),flush=True)
    print(json.dumps({"size_readme":out["readme_relevant"]}),flush=True)
    for k,v in out["code_excerpts"].items():
        print(json.dumps({"size_code":k,"excerpts":v}),flush=True)
if __name__=="__main__":main()
