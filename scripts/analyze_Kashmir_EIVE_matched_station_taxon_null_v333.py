#!/usr/bin/env python3
"""v333 post-hoc robustness of authentic Kashmir geographically matched stations.

Replicate v332's site/year/zone, exact EIVE trait gate and source-geo matching
before any trait-shuffle, one-taxon omission or equal-known-species technical
rarefaction. No source field site GPS or plot-by-species list is redistributed.
All numerical diagnostics are exploratory; source taxon traits remain globally
consistent within each technical permutation.
"""
import argparse
import json
import math
import random
from collections import defaultdict,Counter
from pathlib import Path
from statistics import mean

from analyze_Kashmir_same_sites_reference_year_v254 import parse_xlsx_matrix
from analyze_Kashmir_rail_interior_EIVE_thermal_v327 import species_names_from_publisher
from analyze_Slovenia_EIVE_temperature_width_v321 import eive_table
from analyze_Kashmir_geographic_station_between_EIVE_v332 import (
    site_metadata,site_panel,source_match,TRAIT_GATES
)

PANELS=(
 (2014,"strict_3systems_5species_50pct","G1_geographic_altitude_primary"),
 (2014,"relaxed_1system_3species_30pct","G1_geographic_altitude_primary"),
 (2014,"relaxed_1system_3species_30pct","G2_tight_geographic_altitude"),
 (2014,"relaxed_1system_3species_30pct","G3_geo_altitude_richness_coverage"),
 (2014,"sparse_3systems_1species","G1_geographic_altitude_primary"),
 (2017,"sparse_3systems_1species","G1_geographic_altitude_primary")
)
REFERENCE="Interior Plot 2"
METRIC="position"
DRAWS=1999
RARE_DRAWS=999
SEED=333

def aggregate(values):
    if not values:return {"n":0,"mean":None,"min":None,"max":None}
    return {"n":len(values),"mean":mean(values),
            "min":min(values),"max":max(values)}

def source_plots(plots,names,eive,selected,paired,year,metric,system_count):
    reverse={id(x):site for site,x in selected.items()}
    out=[]
    for pair in paired:
        between=reverse[id(pair["between"])]
        station=reverse[id(pair["station"])]
        item={}
        for label,site in (("B",between),("S",station)):
            for zone,suffix in (("Track","T"),(REFERENCE,"I")):
                original=plots.get((site,year,zone))
                if original is None:
                    raise ValueError("Eligible source field site unexpectedly missing")
                taxa={names[i] for i in original
                    if names[i] in eive and eive[names[i]]["systems"]>=system_count}
                if not taxa:
                    raise ValueError("Matched source has no accepted EIVE trait species")
                item[label+suffix]=taxa
        out.append(item)
    return out

def coefficients(pairs):
    if not pairs:raise ValueError("Empty matched set cannot be analyzed")
    source=defaultdict(float)
    n=len(pairs)
    for p in pairs:
        for field,sign in (("BT",1),("BI",-1),("ST",-1),("SI",1)):
            taxa=p[field]
            if not taxa:raise ValueError("Missing numeric trait species in paired zone")
            w=sign/(n*len(taxa))
            for sp in taxa:source[sp]+=w
    return dict(source)

def paired_mean(pairs,trait):
    diffs=[]
    for p in pairs:
        av=lambda label:mean(trait[s] for s in p[label])
        diffs.append(av("BT")-av("BI")-av("ST")+av("SI"))
    return mean(diffs)

def numeric_observed(pairs,trait):
    co=coefficients(pairs)
    total=sum(co.get(sp,0)*v for sp,v in trait.items())
    if not math.isclose(total,paired_mean(pairs,trait),rel_tol=0,abs_tol=1e-10):
        raise AssertionError("EIVE source trait linear coefficients do not preserve matched CWM difference")
    return total

def shuffled_value_distribution(coeff,traits,kind,draws,seed):
    taxa=sorted(coeff)
    groups=defaultdict(list)
    for sp in taxa:groups[sp.split()[0]].append(sp)
    source={s:traits[s] for s in taxa}
    obs=sum(coeff[s]*source[s] for s in taxa)
    rng=random.Random(seed)
    samples=[]
    for _ in range(draws):
        permuted={}
        if kind=="global":
            vals=[source[s] for s in taxa]
            rng.shuffle(vals)
            permuted=dict(zip(taxa,vals))
        elif kind=="within_genus":
            for g,members in sorted(groups.items()):
                vals=[source[s] for s in members]
                rng.shuffle(vals)
                permuted.update(zip(members,vals))
        else:raise ValueError("Unknown trait label shuffle mode")
        samples.append(sum(coeff[s]*permuted[s] for s in taxa))
    samples.sort()
    extreme=sum(abs(x)>=abs(obs)-1e-12 for x in samples)
    return {"label_shuffle":kind,
        "n_original_source_mapped_taxa_in_matched_comparison":len(taxa),
        "n_source_genera":len(groups),
        "n_genera_with_at_least_two_species":sum(len(a)>1 for a in groups.values()),
        "n_taxa_in_permutable_genera":sum(len(a) for a in groups.values() if len(a)>1),
        "n_draws":draws,"source_trait_value_observed":obs,
        "technical_null_mean":mean(samples),
        "technical_null_95percentile":[samples[int(.025*draws)],samples[int(.975*draws)]],
        "technical_two_sided_fraction_plus_one":(extreme+1)/(draws+1),
        "NOT_field_experiment_or_multiplicity_corrected":True}

def leave_one_taxon_out(pairs,traits):
    baseline=paired_mean(pairs,traits)
    pool=sorted(set().union(*(set().union(*p.values()) for p in pairs)))
    estimates=[];failures=[]
    for species in pool:
        reduced=[]
        for p in pairs:
            q={label:taxa-{species} for label,taxa in p.items()}
            if any(not vals for vals in q.values()):
                continue
            reduced.append(q)
        if not reduced:
            failures.append(species);continue
        new=paired_mean(reduced,traits)
        estimates.append({"taxon":species,"delta_after_removal":new,
            "original_minus_removed_effect":baseline-new,
            "n_retained_pairs":len(reduced)})
    full=[x for x in estimates if x["n_retained_pairs"]==len(pairs)]
    changes=sorted(full,key=lambda x:-x["original_minus_removed_effect"])
    changes_negative=sorted(full,key=lambda x:x["original_minus_removed_effect"])
    return {"n_taxa_seen_in_matched_source_plots":len(pool),
        "taxa_whose_removal_has_estimable_CWM":len(estimates),
        "taxa_preserving_same_all_original_pairs_after_removal":len(full),
        "source_species_omission_delta_with_same_all_pairs":aggregate(
            [x["delta_after_removal"] for x in full]),
        "number_of_full_pair_taxon_deletions_making_difference_nonpositive":
            sum(x["delta_after_removal"]<=0 for x in full),
        "min_matched_pairs_after_any_evaluable_deletion":
            min((x["n_retained_pairs"] for x in estimates),default=None),
        "top_five_taxa_supporting_positive_delta_when_present":[
            {"taxon":x["taxon"],
             "change_in_matched_mean_if_taxon_excluded":-x["original_minus_removed_effect"],
             "remaining_matched_pairs":x["n_retained_pairs"]}
            for x in changes[:5]],
        "top_five_taxa_opposing_positive_delta_when_present":[
            {"taxon":x["taxon"],
             "change_in_matched_mean_if_taxon_excluded":-x["original_minus_removed_effect"],
             "remaining_matched_pairs":x["n_retained_pairs"]}
            for x in changes_negative[:5]],
        "source_missing_or_less_supported_species_not_imputed_or_selected":True}

def rarefy_known_taxa(pairs,trait,draws=RARE_DRAWS,seed=SEED):
    if not pairs:raise ValueError("No source-matched pairs")
    rng=random.Random(seed)
    results=[]
    n_same_pairs=len(pairs)
    for _ in range(draws):
        ds=[]
        for p in pairs:
            means={}
            for suffix in ("T","I"):
                n=min(len(p["B"+suffix]),len(p["S"+suffix]))
                if n<1:raise AssertionError("Source originally eligible matching pair had no taxa")
                for group in ("B","S"):
                    source=sorted(p[group+suffix])
                    vals=(rng.sample(source,n) if len(source)>n else source)
                    means[group+suffix]=mean(trait[s] for s in vals)
            ds.append(means["BT"]-means["BI"]-means["ST"]+means["SI"])
        results.append(mean(ds))
    results.sort()
    return {"n_source_original_pairs":n_same_pairs,
            "n_species_equalized_diagnostic_draws":draws,
            "mean_across_technical_draws":mean(results),
            "technical_95percentile":[results[int(.025*draws)],results[int(.975*draws)]],
            "technical_draws_positive":sum(x>0 for x in results),
            "original_exact_species_site_presence_fixed_except_technical_subsampling":True,
            "technical_percentiles_are_NOT_independent_source_field_confidence_intervals":True}

def analyze(kashmir_original,eive_original):
    plots,originalqa=parse_xlsx_matrix(Path(kashmir_original))
    names=species_names_from_publisher(kashmir_original)
    eive,eiveqa=eive_table(eive_original)
    geo,geoqa=site_metadata(kashmir_original)
    output={"version":"v333",
       "source_kashmir_sha256":originalqa["original_file_sha256"],
       "source_eive_md5":eiveqa["EIVE_original_md5"],
       "source_field_sites":geoqa["n_field_sites"],
       "independent_original_studies":1,
       "status":"SOURCE_VERIFIED_AFTER_V332_POST_HOC_DIAGNOSTIC",
       "panels":[]}
    for i,(year,gate,matching) in enumerate(PANELS):
        candidates,eligibility=site_panel(plots,names,eive,geo,year,
                                        REFERENCE,METRIC,gate)
        matched,matching_qc=source_match(candidates,matching)
        item={"year":year,"zone_pair":["Track",REFERENCE],"metric":METRIC,
              "trait_source_gate":gate,"geographical_matching_gate":matching,
              "source_eligibility":eligibility,"source_matching_qa":matching_qc}
        if not matched:
            item["estimable"]=False
            output["panels"].append(item);continue
        systems=TRAIT_GATES[gate][0]
        pairs=source_plots(plots,names,eive,candidates,matched,year,METRIC,systems)
        trait={sp:float(eive[sp][METRIC])
               for sp in set().union(*(set().union(*p.values()) for p in pairs))}
        obs=numeric_observed(pairs,trait)
        expected=mean(x["between"]["temp_total"]-x["station"]["temp_total"]
                      for x in matched)
        if not math.isclose(obs,expected,abs_tol=1e-10,rel_tol=0):
            raise AssertionError("v332 original-matched geographic baseline drift")
        coeff=coefficients(pairs)
        item.update({
          "estimable":True,"n_unique_original_geographical_pairs":len(matched),
          "baseline_v332_same_source_matched_thermal_mean":obs,
          "global_trait_label_conditional_null":shuffled_value_distribution(
              coeff,trait,"global",DRAWS,SEED+100*i),
          "within_genus_trait_label_conditional_null":shuffled_value_distribution(
              coeff,trait,"within_genus",DRAWS,SEED+100*i+1),
          "source_one_taxon_leaveout":leave_one_taxon_out(pairs,trait),
          "technical_numeric_species_equalization":rarefy_known_taxa(
              pairs,trait,RARE_DRAWS,SEED+100*i+2)})
        output["panels"].append(item)
    output["limitations"]=[
      "v333 is post hoc after v332 exploratory positive difference, so these are diagnostics NOT confirmatory tests",
      "All source plots are from one Kashmir corridor, temporally repeated site IDs do not create independent studies",
      "Global species EIVE label shuffle destroys genus and phylogenetic niche association, not a genuine habitat no-effect null",
      "Within-genus shuffle can have few movable taxa and low power; null nonrejection does not prove no ecological selection",
      "Technical rarefaction equalizes species *already present* and by linearity does not correct selective taxon composition bias",
      "Taxon-name leave-one-out is descriptive source influence, not direct proof of adaptation or propagule delivery",
      "Station versus between is a nonrandom contrast and may differ in land use, light, vegetation and maintenance",
      "Original site-specific coordinates, site names and plot-species source incidence not redistributed",
      "Microclimate temperatures, physiological niche limits, train seed flux, genotypes or recruitment are not measured"]
    return output

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--kashmir-original",required=True)
    p.add_argument("--eive-original",required=True)
    p.add_argument("--output",required=True)
    args=p.parse_args()
    report=analyze(args.kashmir_original,args.eive_original)
    dst=Path(args.output);dst.parent.mkdir(parents=True,exist_ok=True)
    dst.write_text(json.dumps(report,indent=2,ensure_ascii=False),encoding="utf-8")
    compact=[]
    for p in report["panels"]:
        compact.append({"year":p["year"],"gate":p["trait_source_gate"],
             "match":p["geographical_matching_gate"],
             "n_pairs":p["source_matching_qa"]["matched_count"],
             **({k:p[k] for k in ("baseline_v332_same_source_matched_thermal_mean",
                  "global_trait_label_conditional_null",
                  "within_genus_trait_label_conditional_null",
                  "source_one_taxon_leaveout",
                  "technical_numeric_species_equalization")} if p.get("estimable") else {})})
    print(json.dumps({"version":"v333","panels":compact},
                     ensure_ascii=False),flush=True)

if __name__=="__main__":main()
