#!/usr/bin/env python3
"""Generate minimal synthetic resistor networks and verify ACTUAL Julia
Circuitscape.jl NETWORK pairwise output against independent SciPy sparse
Laplacian values. No Overture, GBIF or biological outcomes in this job.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path
import sys
import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
from tetsu.circuit_graph import add_conductance,pairwise_effective_resistance,write_circuitscape_network

CASES={
    "triangle":{
        "edges":[("a","b",1.),("b","c",1.),("a","c",1.)],
        "focals":["a","b","c"]},
    "diamond":{
        "edges":[("a","b",1.),("b","d",1.),("a","c",1.),("c","d",1.)],
        "focals":["a","d"]},
    "series":{
        "edges":[("a","b",1.),("b","c",2.)],
        "focals":["a","c"]},
}


def prepare(root):
    root=Path(root)
    root.mkdir(parents=True,exist_ok=True)
    manifest={"cases":{},"status":"synthetic source only, no biological results"}
    for name,spec in CASES.items():
        g={}
        for u,v,r in spec["edges"]:
            add_conductance(g,u,v,1./r)
        values,_=pairwise_effective_resistance(g,spec["focals"])
        files=write_circuitscape_network(g,spec["focals"],root/name,stem=name)
        integer_id={node:i+1 for i,node in enumerate(sorted(g))}
        expected=[{"a":integer_id[a],"b":integer_id[b],"resistance":float(r)}
                  for (a,b),r in values.items()]
        manifest["cases"][name]={
            "config":str(Path(files["config"]).resolve()),
            "output_matrix":str((root/name/(name+"_julia_matrix.tsv")).resolve()),
            "expected":expected,
            "focals":[integer_id[k] for k in spec["focals"]],
        }
    target=root/"manifest.json"
    target.write_text(json.dumps(manifest,indent=2),encoding="utf-8")
    print(json.dumps({"prepared":list(CASES),"manifest":str(target)}),flush=True)


def validate(root):
    root=Path(root)
    inputs=json.loads((root/"manifest.json").read_text(encoding="utf-8"))
    outcome={}
    for name,record in inputs["cases"].items():
        data=np.loadtxt(record["output_matrix"])
        n=len(record["focals"])
        if data.shape!=(n+1,n+1):
            raise ValueError(f"Julia Circuitscape returned unexpected matrix {name}: {data.shape}")
        col=[int(round(x)) for x in data[0,1:]]
        rows=[int(round(x)) for x in data[1:,0]]
        if rows!=col or set(col)!=set(record["focals"]):
            raise ValueError(f"Julia focal IDs differ from exported source: {name} {col}/{rows}")
        maxerror=0.
        for v in record["expected"]:
            a,b=col.index(v["a"])+1,col.index(v["b"])+1
            measured=float(data[a,b])
            error=abs(measured-v["resistance"])
            maxerror=max(maxerror,error)
            if not np.isfinite(measured) or error>1e-5*(1+abs(v["resistance"])):
                raise ValueError(f"Julia Circuitscape resistance does not match SciPy {name}: {measured}!={v['resistance']}")
        outcome[name]={"n_checked_pairs":len(record["expected"]),
                       "largest_abs_difference":round(maxerror,10)}
    report={"actual_Julia_Circuitscape_network_executed":True,
            "synthetic_mathematical_parity":True,
            "cases":outcome,
            "NO_replicated_plant_seed_dispersal_result":True}
    (root/"parity_result.json").write_text(json.dumps(report,indent=2),encoding="utf-8")
    print(json.dumps(report),flush=True)
    return report


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--outdir",default="results/audit/actual_julia_circuitscape_v41")
    p.add_argument("--validate",action="store_true")
    a=p.parse_args()
    if a.validate:validate(a.outdir)
    else:prepare(a.outdir)


if __name__=="__main__":
    main()
