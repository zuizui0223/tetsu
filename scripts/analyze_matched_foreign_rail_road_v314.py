"""Fixed donor-country rail-vs-road frequency ranking in independent surveys.

No recipient observations train either score. Export aggregate statistics only.
See the effect-frozen design fd24bfe77fba3004af568552953ce24b626a07b4.
"""
from pathlib import Path
import argparse
import hashlib
import json
import re
import numpy as np
from scipy.stats import rankdata

PRIMARY_STUDIES=("SI","HRBA","SK","IN")
DONOR_COUNTRIES=("CZ","IT")
PLAN_COMMIT="fd24bfe77fba3004af568552953ce24b626a07b4"
FREQUENCY_DECIMALS=12


def normalize_name(raw):
    """Existing two-word rule plus explicit ambiguous-concept exclusion."""
    if not isinstance(raw,str): return None
    raw=" ".join(raw.split())
    if re.search(r"\b(?:agg|aggregate|group|complex|sect|section|cf|aff)\b",raw,re.I):
        return None
    if re.search(r"\bs\.\s*(?:l|lat)\.",raw,re.I) or "sensu lato" in raw.lower():
        return None
    if "×" in raw or re.search(r"(?:^|\s)x(?:\s|$)",raw):
        return None
    match=re.match(r"^([A-Z][A-Za-z-]+)\s+([a-z][a-z-]+)(?=\s|$)",raw)
    if not match or match[2] in {"sp","spp","subsp","var","indet","spec","species","s"}:
        return None
    return match[1]+" "+match[2]


def donor_scores(rail,road,unavailable):
    blocked=set().union(*(set(unavailable.get(c,())) | set(unavailable.get("road:"+c,()))
                          for c in DONOR_COUNTRIES))
    by_habitat=[]
    for vectors in (rail,road):
        keys=set().union(*(set(t for t,v in vectors[c].items() if v is not None and v>0)
                           for c in DONOR_COUNTRIES))
        by_habitat.append(keys-blocked)
    union=sorted(by_habitat[0]|by_habitat[1])
    shared=sorted(by_habitat[0]&by_habitat[1])
    result={"shared":shared,"union":union,"excluded_uncertain_donor_keys":len(blocked)}
    for label,vectors in (("rail",rail),("road",road)):
        result[label]={t:float(np.mean([vectors[c].get(t,0.) for c in DONOR_COUNTRIES])) for t in union}
    result["country"]={c:{label:{t:float(vectors[c].get(t,0.)) for t in union}
                            for label,vectors in (("rail",rail),("road",road))}
                       for c in DONOR_COUNTRIES}
    return result


def canonical_frequency(values):
    """Preserve rational ties after floating arithmetic, before all rank gates.

    Every source or weighted score is on a rational grid with spacing at least
    1/3000. Twelve decimals remove arithmetic noise without joining distinct
    observed frequencies; no recipient outcomes determine this precision.
    """
    return np.round(np.asarray(values,dtype=float),FREQUENCY_DECIMALS)


def rho(x,y):
    a,b=rankdata(canonical_frequency(x)),rankdata(canonical_frequency(y))
    if np.ptp(a)==0 or np.ptp(b)==0:
        return None
    return float(np.corrcoef(a,b)[0,1])


def evaluate(f,r,target,taxa,unavailable=(),recorded_only=False):
    candidates=sorted(set(taxa))
    bad=set(unavailable) | {t for t in candidates if t in target and target[t] is None}
    available=[t for t in candidates if t not in bad]
    if recorded_only:
        selected=[t for t in available if target.get(t,0.)>0]
    else:
        selected=available
    y=np.array([target.get(t,0.) for t in selected],float)
    x=np.array([f[t] for t in selected],float)
    z=np.array([r[t] for t in selected],float)
    if not all(np.isfinite(v).all() and np.all((v>=0)&(v<=1)) for v in (x,y,z)):
        raise ValueError("Nonfinite or out-of-range frequency")
    x,y,z=map(canonical_frequency,(x,y,z))
    positive=y[y>0]
    all_values=[v for v in target.values() if v is not None and v>0]
    total_mass=sum(all_values)
    errors=[]
    if len(selected)<(10 if recorded_only else 20): errors.append("insufficient_evaluated_taxa")
    if len(positive)<10: errors.append("fewer_than_10_recorded_taxa")
    if len(np.unique(positive))<3: errors.append("fewer_than_3_positive_frequency_levels")
    for label,v in (("rail_score",x),("road_score",z),("target",y)):
        if len(v)==0 or np.ptp(v)==0: errors.append("constant_"+label)
    a,b=(rho(x,y),rho(z,y)) if not errors else (None,None)
    return {
        "estimable":not errors,"support_failures":errors,
        "n_donor_candidates":len(candidates),"excluded_unavailable":len(set(candidates)&bad),
        "n_evaluated":len(selected),"n_recorded":len(positive),
        "positive_frequency_levels":len(np.unique(positive)),
        "all_target_named_taxa":len(all_values),
        "named_taxon_coverage":len(positive)/len(all_values) if all_values else None,
        "named_frequency_mass_coverage":float(positive.sum()/total_mass) if total_mass else None,
        "rho_rail":a,"rho_road":b,"delta_rho":a-b if a is not None else None,
        "taxon_pool_sha256":hashlib.sha256(("\n".join(selected)+"\n").encode()).hexdigest(),
    }


def study_summary(panels):
    good=[p for p in panels if p["estimable"]]
    values=[p["delta_rho"] for p in good]
    leave=[float(np.mean(values[:i]+values[i+1:])) for i in range(len(values))] if len(values)>1 else []
    return {
        "included_studies":[p["study"] for p in good],
        "excluded_studies":[p["study"] for p in panels if not p["estimable"]],
        "n_studies":len(good),"n_positive":sum(v>0 for v in values),
        "mean_delta_rho":float(np.mean(values)) if values else None,
        "leave_one_study_mean_range":[min(leave),max(leave)] if leave else None,
    }


def habitat_contrast(f,r,rail,background,taxa,rail_unavailable=(),background_unavailable=()):
    common_bad=set(rail_unavailable)|set(background_unavailable)
    common_bad|={t for t in taxa if (t in rail and rail[t] is None) or
                                    (t in background and background[t] is None)}
    a=evaluate(f,r,rail,taxa,common_bad)
    b=evaluate(f,r,background,taxa,common_bad)
    if a["taxon_pool_sha256"]!=b["taxon_pool_sha256"]:
        raise AssertionError("Habitat comparison taxon sets differ")
    return {
        "estimable":a["estimable"] and b["estimable"],
        "rail":a,"background":b,
        "difference_in_advantage":a["delta_rho"]-b["delta_rho"] if a["estimable"] and b["estimable"] else None,
    }


def analyze(source):
    scores=donor_scores(source["rail"],source["road"],source["unavailable"])
    output={
        "design_commit":PLAN_COMMIT,
        "frequency_tie_decimals":FREQUENCY_DECIMALS,
        "interpretation":"Retrospective ranking of identifiable named records, including survey nondetections; not regional absence or causal habitat effects.",
        "source_metadata":source["metadata"],
        "donor_pool":{
            "shared_taxa":len(scores["shared"]),"union_taxa":len(scores["union"]),
            "uncertain_donor_keys_excluded":scores["excluded_uncertain_donor_keys"],
            "shared_sha256":hashlib.sha256(("\n".join(scores["shared"])+"\n").encode()).hexdigest(),
            "union_sha256":hashlib.sha256(("\n".join(scores["union"])+"\n").encode()).hexdigest(),
        },
        "panels":[],"background_contrasts":[],
    }
    scenarios=[
        ("shared","both",False),("shared","both",True),("union","both",False),
        ("shared","CZ",False),("shared","IT",False),
    ]
    for study in (*PRIMARY_STUDIES,"AU_WA"):
        for pool,mode,recorded in scenarios:
            donor=scores if mode=="both" else scores["country"][mode]
            result=evaluate(donor["rail"],donor["road"],source["rail"][study],scores[pool],
                            source["unavailable"].get(study,()),recorded_only=recorded)
            output["panels"].append(dict(study=study,role="primary" if study in PRIMARY_STUDIES else "reserve_extrapolation",
                                         pool=pool,donor_mode=mode,recorded_only=recorded,**result))
    output["study_summaries"]=[]
    for pool,mode,recorded in scenarios:
        rows=[p for p in output["panels"] if p["role"]=="primary" and p["pool"]==pool and
              p["donor_mode"]==mode and p["recorded_only"]==recorded]
        output["study_summaries"].append(dict(pool=pool,donor_mode=mode,recorded_only=recorded,**study_summary(rows)))
    for study,groups in source["background"].items():
        for label,background in groups.items():
            if study not in PRIMARY_STUDIES: continue
            result=habitat_contrast(scores["rail"],scores["road"],source["rail"][study],background,scores["shared"],
                source["unavailable"].get(study,()),
                source["unavailable"].get("background:"+study+":"+label,()))
            output["background_contrasts"].append(dict(study=study,reference=label,**result))
    output["name_mapping_audit"]={
        dataset:{"raw_labels":len(mapping),"mapped_labels":sum(v is not None for v in mapping.values()),
                 "mapped_keys":len({v for v in mapping.values() if v is not None}),
                 "rejected_labels":sum(v is None for v in mapping.values())}
        for dataset,mapping in source["raw_name_audit"].items()}
    return output


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--data",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    from source_matched_rail_road_v314 import load
    source=load(args.data.resolve(),normalize_name)
    result=analyze(source)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    print(json.dumps({"donor_pool":result["donor_pool"],"primary_summary":result["study_summaries"][0]},indent=2))


if __name__=="__main__":
    main()

