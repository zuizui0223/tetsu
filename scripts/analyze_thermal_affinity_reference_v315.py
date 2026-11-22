"""Aggregate-only thermal-affinity comparison under a pre-effect frozen design."""
from pathlib import Path
from collections import Counter,defaultdict
import argparse
import hashlib
import json
import math
import statistics

PLAN_COMMIT="12688cc91aab2a3844c47d512f3ce2b90e9086b8"
EUROPEAN_STUDIES=("CZ","IT","SI","HRBA","SK")


def _hash_names(names):
    return hashlib.sha256(("\n".join(sorted(names))+"\n").encode()).hexdigest()


def _frequencies(freq,exclude=()):
    bad=set(exclude)|{key for key,value in freq.items() if value is None}
    result={}
    for key,value in sorted(freq.items()):
        if key in bad: continue
        if isinstance(value,bool) or not math.isfinite(value) or not 0<=value<=1:
            raise ValueError("Frequency must be a finite fraction")
        if value>0: result[key]=float(value)
    return result


def _validate_traits(trait):
    low,high=trait["limits"]
    if not math.isfinite(low) or not math.isfinite(high) or low>=high:
        raise ValueError("Invalid common trait scale")
    for key,value in trait["values"].items():
        if isinstance(value,bool) or not math.isfinite(value) or not low<=value<=high:
            raise ValueError("Trait outside its declared scale")
        if trait["statuses"].get(key)!="numeric":
            raise ValueError("Numeric trait has a conflicting status")
    if any(status=="numeric" and key not in trait["values"]
           for key,status in trait["statuses"].items()):
        raise ValueError("Numeric status lacks a trait value")


def summarize(freq,trait,exclude=(),min_species=10,species_equal=False):
    _validate_traits(trait)
    if min_species<1: raise ValueError("Minimum support must be positive")
    f=_frequencies(freq,exclude)
    numeric={key:value for key,value in f.items() if key in trait["values"]}
    mass=math.fsum(f.values())
    numeric_mass=math.fsum(numeric.values())
    status_counts=Counter(trait["statuses"].get(key,"unmatched") for key in f)
    status_values=defaultdict(list)
    for key,value in f.items():
        status_values[trait["statuses"].get(key,"unmatched")].append(value)
    status_mass={status:math.fsum(values) for status,values in sorted(status_values.items())}
    x_mass=status_mass.get("indifferent",0.)
    known_status=numeric_mass+x_mass
    weights={key:1. if species_equal else value for key,value in numeric.items()}
    denominator=math.fsum(weights.values())
    estimable=len(numeric)>=min_species and denominator>0
    mean=math.fsum(weight*trait["values"][key] for key,weight in weights.items())/denominator if estimable else None
    return {
        "estimable":estimable,"support_failures":[] if estimable else ["insufficient_numeric_taxa"],
        "weighting":"species_equal" if species_equal else "recorded_frequency",
        "n_named":len(f),"n_numeric":len(numeric),"mean":mean,
        "total_frequency_mass":mass,"numeric_frequency_mass":numeric_mass,
        "numeric_weight_sum":denominator,
        "numeric_taxon_coverage":len(numeric)/len(f) if f else None,
        "numeric_frequency_coverage":numeric_mass/mass if mass else None,
        "status_taxon_counts":dict(sorted(status_counts.items())),
        "status_frequency_fraction":{key:value/mass for key,value in status_mass.items()} if mass else {},
        "indifferent_fraction_among_known_status":x_mass/known_status if known_status else None,
        "named_taxon_sha256":_hash_names(f),"numeric_taxon_sha256":_hash_names(numeric),
    }


def contrast(rail,reference,trait,exclude=(),**kwargs):
    bad=set(exclude)|{key for values in (rail,reference) for key,value in values.items() if value is None}
    a=summarize(rail,trait,bad,**kwargs)
    b=summarize(reference,trait,bad,**kwargs)
    estimable=a["estimable"] and b["estimable"]
    return {"estimable":estimable,"rail":a,"reference":b,
            "difference":a["mean"]-b["mean"] if estimable else None,
            "unavailable_mask_keys":len(bad)}


def shared_unknown_envelope(rail,reference,values,limits,exclude=()):
    low,high=limits
    _validate_traits({"values":values,"statuses":{key:"numeric" for key in values},"limits":limits})
    bad=set(exclude)|{key for f in (rail,reference) for key,value in f.items() if value is None}
    r,b=_frequencies(rail,bad),_frequencies(reference,bad)
    rm,bm=math.fsum(r.values()),math.fsum(b.values())
    if not rm or not bm:
        return {"estimable":False,"lower":None,"upper":None,"reason":"empty_named_assemblage"}
    coefficients={key:r.get(key,0.)/rm-b.get(key,0.)/bm for key in sorted(r.keys()|b.keys())}
    known=math.fsum(a*values[key] for key,a in coefficients.items() if key in values)
    unknown=[a for key,a in coefficients.items() if key not in values]
    lower=known+math.fsum(min(a*low,a*high) for a in unknown)
    upper=known+math.fsum(max(a*low,a*high) for a in unknown)
    return {"estimable":True,"known_component":known,"lower":lower,"upper":upper,
            "unknown_coefficient_l1":math.fsum(abs(a) for a in unknown),
            "unknown_taxa":len(unknown),
            "sign":"positive" if lower>1e-12 else "negative" if upper< -1e-12 else "not_fixed",
            "estimand":"Hypothetical all-named-taxa EIVE contrast, not numeric-only mean or confidence interval"}


def common_numeric_traits(first,second):
    common=first["values"].keys()&second["values"].keys()
    selected=[]
    for trait in (first,second):
        selected.append({**trait,"values":{key:trait["values"][key] for key in sorted(common)},
            "statuses":{key:("outside_common_numeric" if status=="numeric" and key not in common else status)
                        for key,status in trait["statuses"].items()}})
    return tuple(selected)


def paired_summary(pairs,trait,min_species=3,min_pairs=8,exclude=()):
    rows=[]
    for rail,reference in pairs:
        rows.append(contrast({key:1. for key in rail},{key:1. for key in reference},trait,
                             exclude=exclude,min_species=min_species))
    good=[row for row in rows if row["estimable"]]
    differences=[row["difference"] for row in good]
    estimable=len(good)>=min_pairs
    return {
        "estimable":estimable,"n_candidates":len(rows),"n_pairs":len(good),
        "mean_difference":statistics.mean(differences) if estimable else None,
        "median_difference":statistics.median(differences) if estimable else None,
        "n_positive":sum(d>1e-12 for d in differences) if estimable else None,
        "n_negative":sum(d< -1e-12 for d in differences) if estimable else None,
        "n_zero":sum(abs(d)<=1e-12 for d in differences) if estimable else None,
        "difference_range":[min(differences),max(differences)] if estimable else None,
        "mean_rail_numeric_taxon_coverage":statistics.mean(row["rail"]["numeric_taxon_coverage"] for row in good) if good else None,
        "mean_reference_numeric_taxon_coverage":statistics.mean(row["reference"]["numeric_taxon_coverage"] for row in good) if good else None,
        "retained_pair_indices_sha256":_hash_names(str(i) for i,row in enumerate(rows) if row["estimable"]),
        "interpretation":"Equal station weighting within one external region; not independent climatic replication",
    }


def analyze(source,traits,pairs):
    datasets={"rail:"+study:value for study,value in source["rail"].items()}
    datasets.update({"road:"+study:value for study,value in source["road"].items()})
    datasets.update({"background:"+study+":"+label:value for study,groups in source["background"].items()
                     for label,value in groups.items() if label!="miren_road"})
    def mask(key):
        return source["unavailable"].get(key[5:] if key.startswith("rail:") else key,())
    t,e=traits["tichy"],traits["eive"]
    tc,ec=common_numeric_traits(t,e)
    scenarios=[("tichy_numeric",t,False),("eive_all",e,False),
               ("tichy_common_numeric",tc,False),("eive_common_numeric",ec,False),
               ("tichy_species_equal",t,True)]
    references=[("CZ","road:CZ","mountain_road_reference"),
                ("IT","road:IT","mountain_road_reference"),
                ("SI","background:SI:urban","urban_tree_pit_reference"),
                ("HRBA","background:HRBA:urban","urban_tree_pit_reference"),
                ("SK","background:SK:urban","urban_reference"),
                ("IN","background:IN:interior","same_station_same_year_same_area"),
                ("IN","road:IN","external_mountain_road_reference")]
    result={
        "design_commit":PLAN_COMMIT,
        "primary_estimand":"Recorded-frequency weighted regional mean among taxa with numeric Tichy T",
        "source_metadata":source["metadata"],"trait_metadata":traits["metadata"],
        "trait_audits":{key:traits[key]["audit"] for key in ("tichy","eive")},
        "climate_hypothesis_A":{"status":"not_identified","local_paired_independent_regions":1,
            "fitted":False,"reason":"Only IN has eligible paired local rail/nonrail data; regional references are not climate-matched and no common site climate is linked"},
        "assemblages":[],"reference_contrasts":[],"eive_unknown_envelopes":[],
        "paired_india":[],"sampling_uncertainty":"Not estimated; no species-level pseudo-replication or marginal binomial reconstruction",
    }
    for scenario,table,equal in scenarios:
        for name,frequency in datasets.items():
            result["assemblages"].append({"dataset":name,"scenario":scenario,"scale":table["limits"],
                **summarize(frequency,table,exclude=mask(name),species_equal=equal)})
        for study,ref,design in references:
            rail="rail:"+study
            bad=set(mask(rail))|set(mask(ref))
            result["reference_contrasts"].append({"study":study,"reference_dataset":ref,"scenario":scenario,
                "comparison_design":design,"role":"European_descriptive" if study in EUROPEAN_STUDIES else "European_indicator_extrapolation",
                "scale":table["limits"],**contrast(datasets[rail],datasets[ref],table,exclude=bad,species_equal=equal)})
    for study,ref,design in references:
        rail="rail:"+study
        result["eive_unknown_envelopes"].append({"study":study,"reference_dataset":ref,
            **shared_unknown_envelope(datasets[rail],datasets[ref],e["values"],e["limits"],
                                      exclude=set(mask(rail))|set(mask(ref)))})
    for name,table in (("tichy_numeric",t),("eive_all",e)):
        result["paired_india"].append({"scenario":name,"scale":table["limits"],
            **paired_summary(pairs,table,exclude=set(mask("rail:IN"))|set(mask("background:IN:interior")))})
    return result


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--data",type=Path,required=True)
    parser.add_argument("--traits",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    from analyze_matched_foreign_rail_road_v314 import normalize_name
    from source_matched_rail_road_v314 import load
    from source_thermal_affinity_v315 import load_traits,kashmir_pairs
    source=load(args.data.resolve(),normalize_name)
    traits=load_traits(args.traits.resolve(),normalize_name)
    pairs=kashmir_pairs(args.data.resolve(),normalize_name)
    result=analyze(source,traits,pairs)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    print(json.dumps({"primary_contrasts":[{key:row[key] for key in ("study","reference_dataset","difference")}
        for row in result["reference_contrasts"] if row["scenario"]=="tichy_numeric"],
        "climate_hypothesis_A":result["climate_hypothesis_A"]},ensure_ascii=False,indent=2))


if __name__=="__main__":
    main()
