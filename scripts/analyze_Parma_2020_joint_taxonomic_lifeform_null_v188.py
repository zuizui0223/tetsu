#!/usr/bin/env python3
"""v18.8: independent-source Parma 28-station taxonomy × life-form diagnostics.

Replicate the Croatia-Bosnia v18.5 estimands WITHOUT pooling study n: two
distinct source original studies and different site survey protocols.
All 438 original Parma species incl 5 helophytes (He) analyzed as six
lifeform groups, plus common-five sensitivity excluding He.
Internal null fixes station richness and taxon occupancy, not the independent
non-rail *regional species pool*. 5 known extra-effort stations excluded in
sensitivity. This does NOT estimate railway habitat causal effects.
"""
from __future__ import annotations
import argparse,hashlib,json,random,statistics,sys
from collections import Counter
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
from extract_Parma_2020_28station_species_matrix_v187 import (
    parse_pdf, station_stats, EXTRA_EFFORT_STATIONS, KNOWN_BENCHMARKS)
from audit_Parma_2020_28rail_station_flora_original_v186 import STATION_CODES
SEED=20261115
CHAINS=2
DRAWS=120
BURNIN=19000
GAP=1850
PLAIN_11=set("PAR SPT FID BUS CMA CLO MEZ CGU SOR VAI CHI".split())
assert len(PLAIN_11)==11 and len(set(STATION_CODES)-PLAIN_11)==17
FORM6=("T","H","P","G","Ch","He")

def prepare(rows, exclude_extra_effort=False, five_forms=False):
    station_ids=[x for x in STATION_CODES
                 if not(exclude_extra_effort and x in EXTRA_EFFORT_STATIONS)]
    filtered=[r for r in rows if not(five_forms and r["lifeform_coarse"]=="He")]
    if len(rows)!=438:raise ValueError("Requires fully source-reconstructed 438 Parma taxa")
    if any(r["lifeform_coarse"] not in FORM6 for r in filtered):
        raise ValueError("Unexpected life-form class in original Parma source")
    sets=[]
    for code in station_ids:
        sets.append(set(i for i,r in enumerate(filtered) if code in r["station_codes"]))
    if any(not s for s in sets):raise ValueError("Observed station with zero qualified taxa")
    groups=["plain" if x in PLAIN_11 else "other" for x in station_ids]
    return sets,[r["lifeform_coarse"] for r in filtered],station_ids,groups,filtered

def margin(sets,n):
    return tuple(map(len,sets)),tuple(sum(i in site for site in sets) for i in range(n))

def group_margin(sets,n,groups):
    return {g:tuple(sum(i in site for site,group in zip(sets,groups) if group==g)
                    for i in range(n)) for g in sorted(set(groups))}

def metric(sets,life,groups):
    forms=sorted(set(life))
    profiles=[]
    for s in sets:
        n=len(s); counts=Counter(life[i] for i in s)
        profiles.append({k:counts[k]/n for k in forms})
    tax,func=[],[]
    tax_same,func_same=[],[]
    tax_cross,func_cross=[],[]
    for i in range(len(sets)):
        for j in range(i+1,len(sets)):
            sim=len(sets[i]&sets[j])/len(sets[i]|sets[j])
            overlap=sum(min(profiles[i][k],profiles[j][k]) for k in forms)
            tax.append(sim);func.append(overlap)
            if groups[i]==groups[j]:
                tax_same.append(sim);func_same.append(overlap)
            else:
                tax_cross.append(sim);func_cross.append(overlap)
    def mean(a):return statistics.mean(a) if a else None
    return dict(
        tax_mean=mean(tax), functional_mean=mean(func),
        tax_within_group_mean=mean(tax_same),
        functional_within_group_mean=mean(func_same),
        tax_between_group_mean=mean(tax_cross),
        functional_between_group_mean=mean(func_cross),
        n_site_pairs=len(tax),n_within_group_pairs=len(tax_same),
        n_between_group_pairs=len(tax_cross)
    )

MEASURES=("tax_mean","functional_mean","tax_within_group_mean",
          "functional_within_group_mean","tax_between_group_mean","functional_between_group_mean")
def swaps(sites,rng,n,groups,constrained):
    site_groups=[list(range(len(sites)))] if not constrained else [
        [i for i,v in enumerate(groups) if v==g] for g in sorted(set(groups))]
    acc=0
    for _ in range(n):
        pool=rng.choice(site_groups)
        i,j=rng.sample(pool,2)
        a=sites[i]-sites[j];b=sites[j]-sites[i]
        if not a or not b:continue
        x=rng.choice(tuple(a));y=rng.choice(tuple(b))
        sites[i].remove(x);sites[i].add(y)
        sites[j].remove(y);sites[j].add(x)
        acc+=1
    return acc

def evaluate(sets,life,groups,constrained,chains=CHAINS,draws=DRAWS,
             burnin=BURNIN,gap=GAP):
    if chains<2 or draws<2:raise ValueError("Need two technical chains")
    observed=metric(sets,life,groups)
    orig_margin=margin(sets,len(life))
    regional=group_margin(sets,len(life),groups)
    draws_rows=[];chains_rows=[];rates=[]
    for chain in range(chains):
        rng=random.Random(SEED+chain*7297+int(constrained)*900001+
                          len(life)*17+len(sets)*3)
        s=[set(x) for x in sets]
        ac=swaps(s,rng,burnin,groups,constrained)
        chain_stats=[]
        for _ in range(draws):
            ac+=swaps(s,rng,gap,groups,constrained)
            if margin(s,len(life))!=orig_margin:raise AssertionError("Site richness or taxon occupancy changed")
            if constrained and group_margin(s,len(life),groups)!=regional:
                raise AssertionError("Group-specific taxon occupancy changed")
            measure=metric(s,life,groups)
            draws_rows.append(measure)
            chain_stats.append(measure)
        chains_rows.append(chain_stats);rates.append(ac/(burnin+gap*draws))
    result={}
    for field in MEASURES:
        vals=sorted(d[field] for d in draws_rows)
        avg=statistics.mean(vals)
        cmeans=[statistics.mean(x[field] for x in chain) for chain in chains_rows]
        result[field]={
            "source_observed":round(observed[field],7),
            "source_null_mean":round(avg,7),
            "observed_minus_internal_null":round(observed[field]-avg,7),
            "null_quantile025":round(vals[int(0.025*len(vals))],7),
            "null_quantile975":round(vals[int(0.975*len(vals))],7),
            "technical_lower_tail_not_field_p":round((1+sum(v<=observed[field] for v in vals))/(1+len(vals)),7),
            "technical_upper_tail_not_field_p":round((1+sum(v>=observed[field] for v in vals))/(1+len(vals)),7),
            "chain_mean_gap":round(max(cmeans)-min(cmeans),7)
        }
    return {"null_model":"same_publisher_plain_vs_other_fixed_margins" if constrained else "whole_publisher_railway_fixed_margins",
            "study_n":1,"site_n":len(sets),
            "lifeform_count":len(set(life)),
            "n_correlated_technical_draws":chains*draws,
            "n_chains":chains,"draws_per_chain":draws,
            "accepted_proposal_rates":rates,
            "number_within_group_station_pairs":observed["n_within_group_pairs"],
            "number_between_group_station_pairs":observed["n_between_group_pairs"],
            "not_a_regional_non_rail_flora_pool":True,
            "same_original_taxon_occurrences_and_station_richness_fixed":True,
            "source_group_occurrences_fixed_if_stratified":constrained,
            "stats":result}

def baseline_station_group_descriptives(sets,life,stations,groups,rows):
    info={}
    for g in sorted(set(groups)):
        idx=[i for i,gg in enumerate(groups) if g==gg]
        rich=[len(sets[i]) for i in idx]
        tprop=[sum(life[k]=="T" for k in sets[i])/len(sets[i]) for i in idx]
        nonnative=[sum(rows[k]["native_status_original"]!="Autoct." for k in sets[i])/len(sets[i]) for i in idx]
        info[g]={"site_count":len(idx),"original_codes":[stations[i] for i in idx],
                 "mean_rail_station_taxa":round(statistics.mean(rich),5),
                 "mean_T_fraction":round(statistics.mean(tprop),6),
                 "mean_non_native_fraction":round(statistics.mean(nonnative),6)}
    return info

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--original-pdf",required=True);p.add_argument("--output",required=True)
    a=p.parse_args()
    rows,err,audit=parse_pdf(a.original_pdf)
    richness,_,lfs,_=station_stats(rows)
    if len(rows)!=438 or err or lfs!=Counter({"P":71,"T":166,"H":157,"He":5,"G":27,"Ch":12}):
        raise ValueError("Require publisher-original 438 and all 6 life forms; no source fallback")
    if any(richness[k]!=v for k,v in KNOWN_BENCHMARKS.items()):
        raise ValueError("Not the v18.7 source-verified full Parma station flora")
    configs={}
    for five in [False,True]:
        for exclude in [False,True]:
            label=("five" if five else "six")+("_remove_5_extra_effort" if exclude else "_all28")
            sets,life,stations,groups,r=prepare(rows,exclude,five)
            configs[label]={
                "study":"Adorni_Ghillani_2020","site_richness":dict(zip(stations,map(len,sets))),
                "site_group_original_source":dict(zip(stations,groups)),
                "descriptives":baseline_station_group_descriptives(sets,life,stations,groups,r),
                "all_station_null":evaluate(sets,life,groups,False),
                "publisher_plain_vs_other_null":evaluate(sets,life,groups,True)
            }
    out={"version":"v18.8",
         "original_bibliography":"Adorni and Ghillani 2020 Annali Museo Civico di Rovereto 36:149-200",
         "publisher_original_pdf_sha256":audit["source_pdf_sha256"],
         "full_438_by_28_original_source_gate_passed":True,
         "original_species_number_mismatches_from_publisher_printed_n":audit["publisher_frequency_disagreements_count"],
         "station_type_not_independent_microhabitat":"mixed_station_scope",
         "station_publisher_plains_source":"Adorni and Ghillani 2024 uses subset of 11 plain stations from same 2020 source, NOT a new field study",
         "five_extra_effort_codes":sorted(EXTRA_EFFORT_STATIONS),
         "source_THS_category_comparison_to_2017_Croatia_is_DESCRIPTIVE_only":True,
         "study_n_this_Italy_cohort":1,
         "source_cases":configs,
         "warn":"Whole-station mixed habitat, not independent regional flora; this source alone cannot test universal artificial habitat override of biogeographic pools"}
    f=Path(a.output);f.parent.mkdir(parents=True,exist_ok=True)
    f.write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding="utf8")
    summary={k:{"descriptives":v["descriptives"],
                "constrained_tax":v["publisher_plain_vs_other_null"]["stats"]["tax_mean"],
                "constrained_func":v["publisher_plain_vs_other_null"]["stats"]["functional_mean"],
                "unconstrained_tax":v["all_station_null"]["stats"]["tax_mean"],
                "unconstrained_func":v["all_station_null"]["stats"]["functional_mean"]}
             for k,v in configs.items()}
    print(json.dumps({"version":"v18.8","source_sha256":audit["source_pdf_sha256"],
                      "diagnostics":summary},ensure_ascii=False,indent=2),flush=True)

if __name__=="__main__":main()
