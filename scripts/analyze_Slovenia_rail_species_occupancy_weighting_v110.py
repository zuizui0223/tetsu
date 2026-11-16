#!/usr/bin/env python3
"""v11.0 real Slovenian railway species rarity x wind trait occupancy sensitivity.

Original Šabić et al. 2026 Slovenian railway field plots, DOI
10.1007/s44473-026-00155-z, CC BY 4.0 original supplement. 89 sampled
railway field locations nested in 15 railway sections, genuine species
labels, exact same-location BS/BT source matrices. No railway causal claims.

Quantifies estimand dependence: equal species weighting versus taxon
occurrence-frequency weighting, keeping the SAME source species set.
Measures conditional wind x log-frequency interaction descriptively
and uses source lifespan x occurrence-band restricted trait permutations.
No claim of a priori discovery: occupancy bins were motivated by v10.9.
"""
from __future__ import annotations
import argparse,json,math,random
from collections import defaultdict
from pathlib import Path
from statistics import mean

from analyze_Slovenia_rail_microhabitats_v103 import (
    SOURCE_DOI,SOURCE_SHA256,extract_original_xlsx
)
from analyze_Slovenia_rail_species_occupancy_permutation_v109 import (
    STRATA,SEED,original_location_matched_sets,taxa_preference,
)
from analyze_Czech_2024_rail_station_area_aliens_v108 import solve_linear

BINS=("1-2","3-5","6-9","10+")
PERM_DRAWS=3999


def frequency_bin(n):
    if n<1:raise ValueError("Cannot use unobserved species")
    if n<=2:return "1-2"
    if n<=5:return "3-5"
    if n<=9:return "6-9"
    return "10+"


def group_estimator(taxa,min_occ=1,weight_mode="species"):
    """Same species set and phenotypes, change only within-guild weights."""
    if weight_mode not in ("species","opportunities"):
        raise ValueError("Unknown source species weighting definition")
    rows=sorted([r for r in taxa if
                 r["n_original_locations_present_either_zone"]>=min_occ],
                key=lambda r:r["taxon_code"])
    groups={s:{"wind":[],"other":[]} for s in STRATA}
    for r in rows:
        w=1 if weight_mode=="species" else r["n_original_locations_present_either_zone"]
        groups[r["stratum"]]["wind" if r["wind"] else "other"].append(
            (r["occupancy_preference_BS_minus_BT"],w))
    if any(min(len(groups[s]["wind"]),len(groups[s]["other"]))<3 for s in STRATA):
        return {"identified":False,"min_original_site_occupancy":min_occ}
    totalw=sum(w for s in STRATA for v in groups[s].values() for _,w in v)
    contrast=0
    by={}
    for s in STRATA:
        d=groups[s]
        subgroup={}
        for name, vals in d.items():
            sw=sum(w for _,w in vals)
            mu=sum(x*w for x,w in vals)/sw
            subgroup[name]={"n_original_taxa":len(vals),
                            "total_weighted_occupancy_sites":sw,
                            "mean_BS_minus_BT_preference":mu}
        wts=sum(v["total_weighted_occupancy_sites"] for v in subgroup.values())/totalw
        diff=subgroup["wind"]["mean_BS_minus_BT_preference"]-subgroup["other"]["mean_BS_minus_BT_preference"]
        contrast+=wts*diff
        by[s]={"weight":wts,"wind_minus_other_source_preference":diff,
               "guilds":subgroup}
    return {"identified":True,"weight_mode":weight_mode,
            "n_source_taxa":len(rows),"minimum_site_occupancy":min_occ,
            "wind_minus_other_preference_lifespan_adjusted":contrast,
            "by_stratum":by,
            "not_proof_of_seed_transport_or_selection":True}


def rank_logfreq_ols(source):
    """Species-unit (not site-unit) conditional descriptive linear contrast."""
    rows=sorted(source,key=lambda r:r["taxon_code"])
    freqcenter=mean(math.log(r["n_original_locations_present_either_zone"]) for r in rows)
    x=[]
    for r in rows:
        life=int(r["stratum"]=="strict_annual")
        wind=int(r["wind"])
        z=math.log(r["n_original_locations_present_either_zone"])-freqcenter
        x.append([1.0,life,float(wind),z,wind*z])
    y=[r["occupancy_preference_BS_minus_BT"] for r in rows]
    k=len(x[0])
    gram=[[sum(v[i]*v[j] for v in x) for j in range(k)] for i in range(k)]
    rhs=[sum(v[j]*yi for v,yi in zip(x,y)) for j in range(k)]
    b=solve_linear(gram,rhs)
    return {"n_original_source_plant_taxa":len(rows),
            "mean_source_log_occupancy":freqcenter,
            "wind_minus_other_at_mean_log_occupancy":b[2],
            "logfreq_slope_nonwind":b[3],
            "additional_logfreq_slope_for_wind":b[4],
            "coefficient_vector":[round(v,9) for v in b]}


def permute_source_traits_within_lifespan_x_occurrence_band(source,draws=PERM_DRAWS,seed=SEED):
    """Restricted label randomization, preserve lifespan and commonness band.

    Contrast is diagnostic; species are not exchangeable phylogenetically.
    """
    observed=rank_logfreq_ols(source)["additional_logfreq_slope_for_wind"]
    rows=sorted(source,key=lambda r:r["taxon_code"])
    groups=defaultdict(list)
    for i,r in enumerate(rows):
        group=(r["stratum"],frequency_bin(r["n_original_locations_present_either_zone"]))
        groups[group].append(i)
    counts={}
    for group,indices in sorted(groups.items()):
        nw=sum(rows[i]["wind"] for i in indices)
        counts[group[0]+"__"+group[1]]={"n_source_species":len(indices),
                                       "n_source_wind_species":nw}
    rng=random.Random(seed)
    extremes=0
    vals=[]
    for _ in range(draws):
        perm=[dict(r) for r in rows]
        for group,indices in sorted(groups.items()):
            labels=[rows[i]["wind"] for i in indices]
            rng.shuffle(labels)
            for idx,label in zip(indices,labels):
                perm[idx]["wind"]=label
        slope=rank_logfreq_ols(perm)["additional_logfreq_slope_for_wind"]
        vals.append(slope)
    nullcenter=mean(vals)
    extremes=sum(abs(v-nullcenter)+1e-14>=abs(observed-nullcenter) for v in vals)
    return {"n_source_species":len(rows),
            "n_source_lifespan_x_frequency_band_cells":len(counts),
            "original_source_guild_cell_counts":counts,
            "observed_source_wind_logfreq_interaction":observed,
            "restricted_permutation_null_mean":nullcenter,
            "restricted_permutation_two_sided_centered_p":(extremes+1)/(draws+1),
            "n_permutations":draws,
            "frequency_groups_preselected_after_v109_exploration_not_confirmatory":True}


def descriptive_freq_bands(rows):
    groups=defaultdict(list)
    for r in rows:
        groups[(r["stratum"],frequency_bin(r["n_original_locations_present_either_zone"]),
                "wind" if r["wind"] else "other")].append(r["occupancy_preference_BS_minus_BT"])
    out={}
    for s in STRATA:
        out[s]={}
        for band in BINS:
            out[s][band]={}
            for category in ("wind","other"):
                arr=groups[(s,band,category)]
                out[s][band][category]={"n_species":len(arr),
                    "mean_species_BS_minus_BT_preference":mean(arr) if arr else None}
    return out


def run(original,verified109):
    if (verified109.get("complete") is not True or
        verified109.get("version")!="v10.9" or
        verified109.get("original_publisher_field_XLSX_sha256")!=SOURCE_SHA256):
        raise ValueError("Original field source v10.9 verified receipt required")
    sites,traits,provenance=original_location_matched_sets(original)
    data,n=taxa_preference(sites,traits)
    if n!=89 or len(data)!=verified109["n_source_taxa_with_both_trait_codes_present"]:
        raise ValueError("Source original 89-site source taxa differ from v10.9")
    by_occ={}
    for threshold in (1,3,5,10):
        equal=group_estimator(data,threshold,"species")
        weighted=group_estimator(data,threshold,"opportunities")
        if equal["identified"]!=weighted["identified"]:
            raise ValueError("Identical taxa became unidentified solely from weights")
        by_occ[str(threshold)]={"equal_species":equal,"site_occupancy_weighted":weighted}
        old=verified109["species_equal_diagnostics_by_min_site_occupancy"][str(threshold)]
        if equal["identified"] and not math.isclose(
            equal["wind_minus_other_preference_lifespan_adjusted"],
            old["equal_species_within_lifespan_wind_minus_other_preference"],abs_tol=1e-12):
            raise ValueError("v10.9 equal species estimate not reproduced")
    sections=sorted({x["section_number"] for x in sites.values()})
    loo=[]
    for sect in sections:
        subset,_=taxa_preference(sites,traits,exclude_section=sect)
        w=group_estimator(subset,1,"opportunities")
        e=group_estimator(subset,1,"species")
        loo.append({"omitted_source_railway_section":sect,
                    "wind_minus_other_pref_species":e.get(
                        "wind_minus_other_preference_lifespan_adjusted"),
                    "wind_minus_other_pref_occupancy":w.get(
                        "wind_minus_other_preference_lifespan_adjusted")})
    contrast_w=by_occ["1"]["site_occupancy_weighted"][
        "wind_minus_other_preference_lifespan_adjusted"]
    contrast_e=by_occ["1"]["equal_species"][
        "wind_minus_other_preference_lifespan_adjusted"]
    return {
        "complete":True,"version":"v11.0",
        "published_original_source_DOI":SOURCE_DOI,
        "original_source_XLSX_sha256":SOURCE_SHA256,
        "n_true_source_independent_railway_field_locations":n,
        "n_true_railway_sections":len(sections),
        "n_original_taxa_known_both_lifespan_diaspore":len(data),
        "full_1plus_species_equal_wind_minus_other_preference":contrast_e,
        "full_1plus_occupancy_weighted_wind_minus_other_preference":contrast_w,
        "same_taxa_weighting_reversal_occurred":contrast_e*contrast_w<0,
        "equal_species_versus_occupancy_weighted_on_same_source_taxa":by_occ,
        "source_frequency_band_descriptive_2by4":descriptive_freq_bands(data),
        "species_unit_original_logfreq_interaction":rank_logfreq_ols(data),
        "posthoc_frequency_band_restricted_wind_label_permutation":
            permute_source_traits_within_lifespan_x_occurrence_band(data),
        "original_railway_section_leave_one_out_weighting":loo,
        "number_leave_one_section_negative_occ_weight":sum(
            (x["wind_minus_other_pref_occupancy"] is not None and
             x["wind_minus_other_pref_occupancy"]<0) for x in loo),
        "number_leave_one_section_negative_species_equal":sum(
            (x["wind_minus_other_pref_species"] is not None and
             x["wind_minus_other_pref_species"]<0) for x in loo),
        "status_noncausal_exploratory":(
            "The selected frequency bands (1-2,3-5,6-9,10+) follow v10.9; "
            "do not treat p as prospective validation. Species composition "
            "depends on original field site selection, unequal plot area, "
            "source taxon labelling, species phylogenetic affiliation and rarity. "
            "No observation of actual train cargo seeds, viable propagules, "
            "gust flow, maintenance cutting, flowering phenology or genetic adaptation."
        ),
    }


def main():
    p=argparse.ArgumentParser()
    for n in ("original-xlsx","source-audit102","verified109","output"):
        p.add_argument("--"+n,required=True)
    a=p.parse_args()
    data=extract_original_xlsx(a.original_xlsx,a.source_audit102)
    orig109=json.loads(Path(a.verified109).read_text(encoding="utf-8"))
    result=run(data,orig109)
    dest=Path(a.output);dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding="utf-8")
    print(json.dumps(result,ensure_ascii=False),flush=True)


if __name__=="__main__":main()
