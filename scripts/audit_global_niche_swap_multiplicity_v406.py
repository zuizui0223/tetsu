#!/usr/bin/env python3
"""Reproduce v406 arithmetic from the verified v404 aggregate; not calibrated inference."""
import argparse,json
from pathlib import Path
def audit(src):
    if src["version"]!="v404-original-shared-rarefied-conditional-rail-swaps":
        raise ValueError("Wrong original v404 artifact")
    q=src["source_QA"]
    if q["all_source_rows"]!=288721938 or q["source_observation_SHA256"]!="1c48a9646f3bf6d56ffdfdc91a1a176c173c98014830ad48fa162b2803c3b279":
        raise ValueError("Wrong underlying original source")
    rows=[]
    for p in src["panels"]:
        for f in p["fits"]:
            for stat,entries in f["rank"].items():
                for target,d in entries.items():
                    if d["valid"]!=199:raise ValueError("Incomplete technical rank")
                    lo,hi=d["lower_tail_plus_one"],d["upper_tail_plus_one"]
                    if not 0<lo<=1 or not 0<hi<=1:raise ValueError("Invalid tail rank")
                    rows.append(dict(land=p["land"],support=p["support"],reference=f["reference"],stat=stat,target=target,actual=d["actual"],
                        lower_tail=lo,upper_tail=hi,two_tail_technical=min(1,2*min(lo,hi)),valid=d["valid"]))
    if len(rows)!=112 or len({(r["land"],r["support"],r["reference"],r["stat"],r["target"]) for r in rows})!=112:
        raise ValueError("Incomplete original rank family")
    rows.sort(key=lambda r:r["two_tail_technical"])
    running=0.
    for i,r in enumerate(rows):
        running=max(running,min(1.,r["two_tail_technical"]*(len(rows)-i)))
        r["holm_style_numeric_diagnostic"]=running
    return dict(n_rank_entries=len(rows),n_either_one_tail_at_most_05=sum(min(r["lower_tail"],r["upper_tail"])<=.05 for r in rows),
      n_two_tail_at_most_05=sum(r["two_tail_technical"]<=.05 for r in rows),
      n_holm_style_at_most_05=sum(r["holm_style_numeric_diagnostic"]<=.05 for r in rows),
      minimum_two_tail_technical=rows[0]["two_tail_technical"],rows=rows,
      warning="Numerical multiplicity sensitivity only: technical ranks are not valid causal or calibrated p-values.")
def main():
    p=argparse.ArgumentParser();p.add_argument("--source",type=Path,required=True);p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    d=audit(json.loads(a.source.read_text()))
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(d,allow_nan=False))
    print(json.dumps({k:v for k,v in d.items() if k!="rows"}))
if __name__=="__main__":main()
