#!/usr/bin/env python3
"""Reproduce v417 numerical sensitivity from compact verified v416 temporal log record."""
import argparse,json,math
from pathlib import Path
def holm(rows):
    rows=sorted((dict(r) for r in rows),key=lambda r:r["p_for_diagnostic"])
    m=0.
    for i,r in enumerate(rows):
        m=max(m,min(1.,r["p_for_diagnostic"]*(len(rows)-i)))
        r["holm_style_diagnostic"]=m
    return rows
def audit(source):
    if source["source_run"]!=45583858584 or source["source_job"]!=139100335398 or len(source["panels"])!=2:
        raise ValueError("Wrong verified source")
    rows=[]
    for p in source["panels"]:
        for axis,a in p["axes"].items():
            for method in ("block_FE","endpoint_FE"):
                for metric in ("difference_jaccard","difference_simpson"):
                    f=a["block"][metric] if method=="block_FE" else a["endpoint"]
                    for term in ("active","inactive","active_minus_inactive"):
                        for period in ("validation_2023_2024","descriptive_all_years"):
                            r=dict(land=p["land"],axis=axis,method=method,metric=metric,term=term,period=period,n=a["n"],status=f["status"])
                            if f["status"]!="estimated":
                                r["p_for_diagnostic"]=1.
                            else:
                                if method=="endpoint_FE" and metric!="difference_jaccard":
                                    raise ValueError("Compact source lacks supported endpoint Simpson")
                                v=f[period][term];lo,hi=v["interval95_approx"];n=len(v["years"]);cov=v["annual_covariance"]
                                se=math.sqrt(math.fsum(x for row in cov for x in row))/n
                                if abs((hi-lo)/2-1.96*se)>1e-10:raise ValueError("Original CI mismatch")
                                z=v["estimate"]/se
                                r.update(estimate=v["estimate"],ci=[lo,hi],annual_covariance=cov,se=se,z=z,p_for_diagnostic=math.erfc(abs(z)/math.sqrt(2)))
                            rows.append(r)
    full=holm(rows);focus=holm([r for r in rows if r["method"]=="block_FE" and r["term"]!="inactive"])
    if len(full)!=144 or sum(r["status"]=="estimated" for r in full)!=60 or len(focus)!=48:
        raise ValueError("Original comparison family mismatch")
    return dict(expanded_family=full,focused_family=focus,warning="Posthoc normal-tail sensitivity only, not independent confirmation or study-wide multiplicity control.")
def main():
    p=argparse.ArgumentParser();p.add_argument("--source",type=Path,required=True);p.add_argument("--out",type=Path,required=True)
    a=p.parse_args();d=audit(json.loads(a.source.read_text()));a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(d,allow_nan=False))
if __name__=="__main__":main()
