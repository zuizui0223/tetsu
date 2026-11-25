#!/usr/bin/env python3
"""Same-support richness/trait contrast and conservative missing-trait bounds.
These are conditional deterministic identification bounds, not confidence
intervals, and relax cross-cell species consistency. No missing trait imputation.
"""
import math,heapq
from collections import defaultdict
import numpy as np
from analyze_world_railstatus_growthform import FULL,POP,ACTIVE,INACTIVE,local_support,block_weights,fit,flexible,EXTRA
FIELDS=FULL+(POP,)
CLASSES=("F","C","E","P","M","A","H")

def gini_bounds(counts,total):
    if type(total)!=int or total<2 or not set(counts)<=set(CLASSES) or any(type(x)!=int or x<0 for x in counts.values()):
        raise ValueError("Invalid support counts")
    a=[counts.get(k,0) for k in CLASSES]
    u=total-sum(a)
    if u<0:raise ValueError("Known counts exceed total")
    lo=a.copy();lo[lo.index(max(lo))]+=u
    hi=a.copy();heapq.heapify(hi)
    # Each added count increases sum of squares by 2*k+1; greedy
    # allocation to a least occupied category is exact for integer counts.
    for _ in range(u):heapq.heapreplace(hi,hi[0]+1)
    def g(c):return (total*total-sum(x*x for x in c))/(total*(total-1))
    l,h=g(lo),g(hi)
    if not 0<=l<=h<=1:raise ValueError("Invalid diversity bounds")
    return l,h

def woodiness_bounds(woody,known,total):
    if any(type(x)!=int for x in (woody,known,total)) or not 0<=woody<=known<=total or total<=0:
        raise ValueError("Invalid woodiness counts")
    return woody/total,(woody+total-known)/total

def linear_weights(rows,fields):
    # Exact coefficient weights for the same block-weighted WLS design
    # as the known-outcome model; constants dropped before scaling.
    groups=defaultdict(list)
    for i,r in enumerate(rows):groups[r["continent"],r["block"]].append(i)
    if len(rows)<100 or len(groups)<20 or len(rows)-len(groups)-len(fields)<30:return None
    w=block_weights(rows)
    X=np.array([[r[k] for k in fields] for r in rows],float)
    if not np.isfinite(X).all():raise ValueError("Invalid design")
    for ix in groups.values():X[ix]-=np.average(X[ix],axis=0,weights=w[ix])
    A=X*np.sqrt(w)[:,None]
    scale=np.sqrt(np.sum(A*A,axis=0));keep=scale>1e-10
    names=[f for f,k in zip(fields,keep) if k]
    if not {ACTIVE,INACTIVE}<=set(names):return None
    Z=A[:,keep]/scale[keep];s=np.linalg.svd(Z,compute_uv=False)
    if s[-1]<1e-8*s[0]:return None
    inv=np.linalg.inv(Z.T@Z)
    out={}
    for label,terms in (("active",{ACTIVE:1.}),("inactive",{INACTIVE:1.}),("active_minus_inactive",{ACTIVE:1.,INACTIVE:-1.})):
        c=np.array([terms.get(k,0.) for k in names])
        a=np.sqrt(w)*(Z@(inv@(c/scale[keep])))
        # y is demeaned with the same weights; enforce the exact linear
        # mapping to raw cell outcomes (removes numerical rounding only).
        for ix in groups.values():a[ix]-=w[ix]/w[ix].sum()*a[ix].sum()
        if abs(float(a.sum()))>1e-8:raise ValueError("Coefficient weights do not annihilate constant")
        out[label]=a
    return out

def bound(a,l,h):
    a=np.asarray(a,float);l=np.asarray(l,float);h=np.asarray(h,float)
    if not (a.shape==l.shape==h.shape) or np.any(l>h):raise ValueError("Invalid bound inputs")
    return [float(np.where(a>=0,a*l,a*h).sum()),float(np.where(a>=0,a*h,a*l).sum())]

def conjunction(model,target,contrast):
    if model["status"]!="estimated":return "unsupported"
    z=model["effects"][contrast]
    # Joint directional conclusion requires both separate intervals in the
    # hypothesized directions; no claim of a combined p-value.
    a=z["log_species"]["interval95_approx"];b=z[target]["interval95_approx"]
    return "both directions nominally supported" if a[0]>0 and b[1]<0 else "not both directions supported"

def evaluate(rows,traits):
    rr=local_support(rows)
    r=[]
    for x in rr:
        t=traits[x["key"]];N=x["species"];n=t["known_woodiness"]
        if n<20 or t["known_support"]<20:raise ValueError("Common known support required")
        lo,hi=gini_bounds(t["support_counts"],N);wl,wh=woodiness_bounds(t["woody"],n,N)
        r.append(dict(x,woody_fraction=t["woody"]/n,wood_coverage=n/N,
                      support_lower=lo,support_upper=hi,wood_lower=wl,wood_upper=wh))
    out={"n_common_local":len(r),"models":{}}
    for label,design,fields in (
       ("common_design",r,FIELDS),
       ("flexible_effort_and_two_coverages",[flexible(x) for x in r],FIELDS+EXTRA+("coverage","wood_coverage"))):
        z={"richness_support_diversity":fit(design,fields,("log_species","support_diversity")),
           "richness_woody_fraction":fit(design,fields,("log_species","woody_fraction")),
           "coverage_diagnostics":fit(design,tuple(f for f in fields if f not in ("coverage","wood_coverage")),("coverage","wood_coverage"))}
        z["joint_direction"]={k:{t:conjunction(z[m],t,k) for t,m in
             (("support_diversity","richness_support_diversity"),("woody_fraction","richness_woody_fraction"))}
             for k in ("active","active_minus_inactive")}
        a=linear_weights(design,fields)
        z["all_species_missingness_outer_bounds"]={}
        if a is not None:
            for k,v in a.items():
                for target,model in (("support_diversity","richness_support_diversity"),("woody_fraction","richness_woody_fraction")):
                    m=z[model]
                    if m["status"]!="estimated":raise ValueError("Bound design differs from fitted design")
                    pred=float(v@np.array([x[target] for x in design]))
                    if abs(pred-m["effects"][k][target]["estimate"])>1e-8:raise ValueError("Coefficient weights fail exact fitted replay")
                z["all_species_missingness_outer_bounds"][k]={
                    "support_diversity":bound(v,[x["support_lower"] for x in design],[x["support_upper"] for x in design]),
                    "woody_fraction":bound(v,[x["wood_lower"] for x in design],[x["wood_upper"] for x in design])}
        out["models"][label]=z
    # Richness is an outcome above, never used as its own predictor.
    # A separate conditional trait model answers a different question.
    out["trait_only_conditional_on_total_richness"]=fit(r,FIELDS+("log_species",),("support_diversity","woody_fraction"))
    return out
