#!/usr/bin/env python3
"""tetsu v16.2 — source-authenticated 2-country × 2×2 functional-group
rail-vs-road species identity prediction. POSTHOC pilot; NOT a rail-seed meta.

Independent ORIGINAL author botanical observations:
  SI: Sabic et al 2026 (DOI 10.1007/s44473-026-00155-z),
      89 real field 100-m rail sites / 15 original line sections, 2022–23.
      Original author STRICT annual lifespan × source diaspore strategy.
      TRUE original field transect GPS and 2024-01 Geofabrik OSM graph.
  CZ: Kutlvasr et al 2024 (DOI 10.23855/preslia.2024.247),
      original 39 unequal-area stations, 35 strict-name matched to
      2026 SR70 station CENTROID PROXIES, not field GPS; 5 line stretches.
      Original Pladias Raunkiaer THEROPHYTE proxy × Sádlo strategy.
      2023-01 Geofabrik Czech physical rail+drivable road graph.

4 biological guilds are independently coded from ORIGINAL source values:
annual/proxy × wind-STRATEGY ARCHETYPE, annual×other,
nonannual/proxy × wind, nonannual×other. Wind strategy is NOT
observed seed carried by train wind; annual is NOT observed mowing
tolerance or evolved plasticity.

Main prior v15.1 and v15.6 full-trait-known Jaccard road+rail CV MSEs
are recomputed from original sources and exactly checked before
any posthoc subgroup output.

Original site-pair rail/road geography/section masks and baseline source
covariates are inherited exactly, never tuned to any ecological response.
Each subgroup Jaccard is UNDEFINED if BOTH sites have no original
species of that guild: drop such pairs, do NOT set similarity=0.
A group observed at one endpoint only yields Jaccard=0 legitimately.
Each source original full-assemblage Jaccard metric is a separate
known-taxon baseline, not a species-independent seed-distance effect.

For all subgroups, two richness controls:
  group_log1p: group-species richness of each original site; zeros allowed;
  whole_flora_source: the EXACT ORIGINAL v15.1/15.6 site richness controls.
All other geography, station-area (Czech), source functional profile
and real road/rail path predictors unchanged, hence comparisons can be
made WITHIN original study; do not pool country betas/mean squared errors.

Minimum support gate is identical to parent source study's base CV
(>=180 pairs/50 sites/8 original sections SI; >=250 pairs/30 sites/
4 original stretches CZ). Report not-estimable if group too sparse.
Original heldout folds remain dependent; there is NO independent
randomized seed-flow, mowing, train traffic OD or common-garden evidence.
"""
from __future__ import annotations
import argparse,hashlib,json,math
from collections import Counter
from pathlib import Path

from analyze_Slovenia_rail_microhabitats_v103 import extract_original_xlsx
from analyze_Slovenia_true_OSM_rail_road_paths_v150 import run_network as si_network
from analyze_Slovenia_rail_road_heldout_species_identity_v151 import (
    model_records as si_records,section_heldout_validation as si_cv
)
from analyze_Czech_2023_true_rail_road_osm_routes_v155 import run as cz_network
from analyze_Czech_2023_true_rail_road_heldout_taxon_identity_v156 import (
    pair_records as cz_records,cv_section_heldout as cz_cv
)
from analyze_Czech_2024_39station_functional_identity_v143 import jaccard

SI_DOI="10.1007/s44473-026-00155-z"
CZ_DOI="10.23855/preslia.2024.247"
SI_PBF_SHA="b572258fe000d6dfc76dc63afe5456121c5f463995e132a5dcd79957d819ae39"
CZ_PBF_SHA="26ff09ded5fa72dc9556b543b532d7df011f2e3c08af461ea2f537a9b0b42747"
FROZEN_REFERENCE={
  "SI":{"road":0.0037366620160146245,"road_rail":0.0037567314602607142,
        "n_original_sites":89,"n_original_segments":15,"n_supported_site_pairs":3896},
  "CZ":{"road":0.0027606838643546936,"road_rail":0.0028333736388869778,
        "n_original_sites":35,"n_original_segments":5,"n_supported_site_pairs":594},
}
MAIN_GROUPS={
  "annual_wind":{
      "SI":"strict_annual_wind","CZ":"th_wind"},
  "annual_nonwind":{
      "SI":"strict_annual_other","CZ":"th_other"},
  "nonannual_wind":{
      "SI":"other_known_wind","CZ":"nonth_wind"},
  "nonannual_nonwind":{
      "SI":"other_known_other","CZ":"nonth_other"},
}
RICHNESS_CONTROLS=("guild_specific_log1p","whole_original_known_flora")

def source_references(study,topo,base,sites):
    if study not in FROZEN_REFERENCE:
        raise ValueError("No support for unknown original railway study")
    known=FROZEN_REFERENCE[study]
    if len(sites)!=known["n_original_sites"] or len(base)!=known["n_supported_site_pairs"]:
        raise ValueError("Source original site and pair original sample cohort drifted")
    if len({s["section"] for s in sites})!=known["n_original_segments"]:
        raise ValueError("Original author study line section units drifted")
    if study=="SI":
        if topo["original_vegetation_DOI"]!=SI_DOI or topo["original_OSM_pbf_SHA256"]!=SI_PBF_SHA:
            raise ValueError("Slovenian original published DOI and genuine fixed 2024 OSM SHA drifted")
        cv=si_cv(base)
        if not cv["identified"] or cv["n_heldout_sections_estimable"]!=15:
            raise ValueError("Source 15 original Slovenian field line sections did not reproduce")
        m=cv["source_section_equal_weight_MSE_by_model"]
        a,b=m["plus_road"],m["plus_road_and_rail"]
    else:
        if topo["original_botanical_source_DOI"]!=CZ_DOI or topo["actual_download_source_SHA256"]!=CZ_PBF_SHA:
            raise ValueError("Czech original published DOI and genuine fixed 2023 OSM SHA drifted")
        cv=cz_cv(base)
        if not cv["identified"] or cv["n_holdout_original_stretches"]!=5:
            raise ValueError("Source five Czech station route stretches did not reproduce")
        m=cv["equal_original_stretch_weighted_mse"]
        a,b=m["M1_plus_actual_2023_road_route"],m["M2_plus_actual_2023_road_and_rail_ratio"]
    if not all(math.isclose(float(x),known[y],rel_tol=1e-9,abs_tol=1e-12)
               for x,y in ((a,"road"),(b,"road_rail"))):
        raise ValueError("Original author-source v15.1/v15.6 full-flora railway-vs-road HOLDOUT result not exactly reproducible")
    return {"identified":True,"author_original_DOI":SI_DOI if study=="SI" else CZ_DOI,
            "n_original_study_sites":len(sites),
            "n_original_line_sections":known["n_original_segments"],
            "n_original_comparable_rail_and_road_site_pairs":len(base),
            "reference_road_only_holdout_MSE":float(a),
            "reference_road_plus_rail_holdout_MSE":float(b),
            "original_field_source_2023_or_2024_OSM_SHA":SI_PBF_SHA if study=="SI" else CZ_PBF_SHA,
            "no_source_assemblage_reference_change":True}


def guild_sets(sites,guild,study):
    if study not in ("SI","CZ"):
        raise ValueError("Unknown source study and original trait ontology")
    source_labels={g[study] for g in MAIN_GROUPS.values()}
    if set(guild.values())!=source_labels:
        raise ValueError("Original field author four-guild source code ontology changed")
    known=set(guild)
    if any(not s["known"]<=known for s in sites):
        raise ValueError("Original known site species code absent from author guild mapping")
    classes={}
    for group,source_name in MAIN_GROUPS.items():
        label=source_name[study]
        taxa_set={name for name,tag in guild.items() if tag==label}
        original=[s["known"]&taxa_set for s in sites]
        classes[group]={
            "source_label":label,
            "observed_taxa":len(set().union(*original)),
            "source_site_guild_incidences":sum(len(x) for x in original),
            "n_original_site_with_guild_taxa":sum(bool(x) for x in original),
            "n_source_site_with_zero_guild_species":sum(not x for x in original),
            "site_sets":original
        }
    original_union=set().union(*(s["known"] for s in sites))
    reconstructed=set().union(*(set().union(*c["site_sets"]) for c in classes.values()))
    if original_union!=reconstructed or sum(c["source_site_guild_incidences"] for c in classes.values())!=sum(len(s["known"]) for s in sites):
        raise ValueError("Four source guild partitions do not exactly cover original known plant incidence")
    return classes


def replace_response(study,base,sites,group_sites,richness="guild_specific_log1p"):
    """Keep original source physical site pair/sampling frame and all covariates.
    Drop double-zero guild pairs ONLY, not presumed sampled absences.
    """
    if study not in ("SI","CZ") or richness not in RICHNESS_CONTROLS:
        raise ValueError("No supported original source/species richness mode")
    if len(sites)!=len(group_sites) or len(sites)<2:
        raise ValueError("Original site source arrays not aligned")
    mean_key,abs_key=(("site_log_richness_mean","site_log_richness_absdiff")
                      if study=="SI" else
                      ("site_ln_species_count_mean","site_ln_species_count_absdiff"))
    excluded_double_zero=0
    rows=[]
    for r in base:
        i,j=r["i"],r["j"]
        x,y=group_sites[i],group_sites[j]
        if not x and not y:
            excluded_double_zero+=1
            continue
        val=jaccard(x,y)
        if val is None or not(0<=val<=1):
            raise ValueError("Source guild Jaccard missing/invalid; cannot fabricate zero-zero similarity")
        new=dict(r)
        f=dict(r["features"])
        if richness=="guild_specific_log1p":
            ni,nj=math.log1p(len(x)),math.log1p(len(y))
            f[mean_key]=(ni+nj)/2
            f[abs_key]=abs(ni-nj)
        new["features"]=f
        new["jaccard"]=float(val)
        new["source_original_4guild_posthoc_ontology"]=study
        rows.append(new)
    return rows,{
       "source_full_original_comparable_site_pairs":len(base),
       "source_pairs_excluded_because_both_sites_have_zero_original_guild_species":excluded_double_zero,
       "source_guild_pair_comparisons_with_at_least_one_original_species":len(rows),
       "unobserved_original_species_NOT_inserted_as_absences":True,
       "other_source_geography_road_rail_functional_covariates_unchanged":True
    }


def group_cv(study,rows):
    if study=="SI":
        cv=si_cv(rows)
        if not cv["identified"]:
            return {
              "identified":False,
              "source_support":cv.get("support") or cv.get("source_pair_support"),
              "reason":cv.get("reason","Original study section-fold/source site eligibility insufficient")
            }
        vals=cv["source_section_equal_weight_MSE_by_model"]
        m0,m1,m2=vals["geo_environment_section_control"],vals["plus_road"],vals["plus_road_and_rail"]
        folds=cv["details_per_original_section"]
        npos=cv["n_source_sections_rail_increment_positive"]
        total=cv["n_heldout_sections_estimable"]
        fold_info=[{
            "original_author_stretch":d["source_section"],
            "source_n_test_dependent_site_pairs":d["n_holdout"],
            "rail_increment_after_road_source_fold_MSE_reduction":d["rail_increment_beyond_road_mse_reduction"]
          } for d in folds if d["estimated"]]
        support=cv["source_pair_support"]
    elif study=="CZ":
        cv=cz_cv(rows)
        if not cv["identified"]:
            return {"identified":False,"source_support":cv.get("source_overlap_support"),
                    "reason":cv.get("reason","Original source station section-fold eligibility insufficient")}
        vals=cv["equal_original_stretch_weighted_mse"]
        m0,m1,m2=vals["M0_geography_functional_area_station_size"],vals["M1_plus_actual_2023_road_route"],vals["M2_plus_actual_2023_road_and_rail_ratio"]
        folds=cv["source_stretch_holdout"]
        npos=cv["rail_improves_how_many_stretches"]
        total=cv["n_holdout_original_stretches"]
        fold_info=[{
            "original_author_stretch":d["heldout_original_source_stretch"],
            "source_n_test_dependent_site_pairs":d["n_test_pairs"],
            "rail_increment_after_road_source_fold_MSE_reduction":d["rail_after_road_increment_MSE_gain"]
          } for d in folds if d["identified"]]
        support=cv["source_overlap_support"]
    else:
        raise ValueError("Unknown original study code")
    return {"identified":True,
       "source_original_section_holdout_MSE_M0":m0,
       "source_original_section_holdout_MSE_M1_road":m1,
       "source_original_section_holdout_MSE_M2_road_and_rail":m2,
       "source_original_section_equal_weight_rail_after_road_MSE_gain_pct":100*(m1-m2)/m1,
       "n_original_source_sections_rail_addition_improved":npos,
       "n_original_source_sections_eligible":total,
       "original_source_section_heldout_guild_fold_score":fold_info,
       "source_support":support,
       "not_assumed_independent_site_pair_tests":True}


def run(study,si_xlsx=None,si_audit=None,cz_s1=None,cz_s2=None,
        cz_disp=None,cz_life=None,cz_sr70=None,pbf=None):
    if study=="SI":
        if not(si_xlsx and si_audit and pbf):
            raise ValueError("Slovenia true original XLSX, audit and original OSM PBF required")
        original=extract_original_xlsx(si_xlsx,si_audit)
        topo,pairs,sites,guild=si_network(original,pbf,True)
        base,omitted=si_records(pairs,sites,True)
    elif study=="CZ":
        if not all((cz_s1,cz_s2,cz_disp,cz_life,cz_sr70,pbf)):
            raise ValueError("Czech Preslia S1/S2, Pladias source, SR70 and historic OSM required")
        topo,pairs,sites,guild=cz_network(cz_s1,cz_s2,cz_disp,cz_life,cz_sr70,pbf,True)
        base,omitted=cz_records(pairs,sites,guild,"both_Pladias_known")
    else:raise ValueError("Unrecognized original field study")
    reference=source_references(study,topo,base,sites)
    classes=guild_sets(sites,guild,study)
    result={}
    for name,grp in classes.items():
        values={}
        for mode in RICHNESS_CONTROLS:
            rows,support=replace_response(study,base,sites,grp["site_sets"],mode)
            values[mode]={"source_comparison_eligibility":support,
                          "original_source_heldout_rail_vs_road_prediction":group_cv(study,rows)}
        result[name]={
            "source_original_lifeform_diaspore_code":grp["source_label"],
            "source_distinct_observed_taxa":grp["observed_taxa"],
            "source_guild_site_taxon_incidences":grp["source_site_guild_incidences"],
            "source_sites_with_at_least_one_group_taxon":grp["n_original_site_with_guild_taxa"],
            "source_sites_with_none_of_guild_taxa":grp["n_source_site_with_zero_guild_species"],
            "richness_source_sensitivity":values,
        }
    return {
      "complete":True,
      "version":"v16.2_TWO_ORIGINAL_COUNTRIES_DISTINCT_2x2_FUNCTIONAL_GUILD_RAILROAD",
      "source_original_country":study,
      "source_authentic_original_botanical_and_geofabrik_provenance":topo,
      "source_reference_assemblage_road_vs_rail_CV_reproduced":reference,
      "original_source_fully_known_4guild_partition":result,
      "n_biologically_interpretable_classes":4,
      "two_original_studies_NOT_identical_lifespan_traits_or_survey_grain":True,
      "posthoc_novel_4guild_hypothesis_NOT_independent_validation":True,
      "strict_limits":[
        "SI original STRICT annual lifespan is not CZ Pladias THEROPHYTE lifeform proxy; never compare as identical effect size.",
        "Source diaspore wind strategy label is a trait guild archetype, not observed strong wind/turbulence exposure, measured true flight, train-borne seed capture or adaptation.",
        "Original site-by-taxon incidence is not actual same cohort arrival/germination/recruitment or individual phenotype; no mowing history or genetic common-garden.",
        "CZ original 2020-22 whole-station field data matched to 2026 derived SR70 STATION CENTROID, not author field plot GPS; SI uses original 2022-23 real field GPS.",
        "Source 2023/2024 historic Geofabrik OSM maps physical shortest connectivity, not operational freight OD/actual train traffic, and national border routes may be clipped.",
        "Double-zero guild pairs dropped, including differently selected regional species pools: posthoc supported pair populations vary across groups and must be shown.",
        "Only original 15 SI sections/5 CZ stretches are within-country route clusters; pairs and heldout fold errors not independent field studies.",
        "The two original study cohorts cannot be combined as a causal or standardized global rail-seed transfer meta effect.",
      ],
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--study",required=True,choices=("SI","CZ"))
    for name in ("si-original-xlsx","si-source-audit102","cz-s1","cz-s2",
                 "cz-diaspore","cz-lifeform","cz-sr70","osm-pbf","output"):
        p.add_argument("--"+name,required=(name in ("osm-pbf","output")))
    a=p.parse_args()
    result=run(a.study,si_xlsx=a.si_original_xlsx,si_audit=a.si_source_audit102,
               cz_s1=a.cz_s1,cz_s2=a.cz_s2,cz_disp=a.cz_diaspore,
               cz_life=a.cz_lifeform,cz_sr70=a.cz_sr70,pbf=a.osm_pbf)
    dst=Path(a.output)
    dst.parent.mkdir(parents=True,exist_ok=True)
    dst.write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding="utf-8")
    print(json.dumps({
      "version":result["version"],
      "study":result["source_original_country"],
      "original_source_v15_whole_guild_control_reproduction":result["source_reference_assemblage_road_vs_rail_CV_reproduced"],
      "four_source_functional_guild_posthoc_heldout":{
         group:{
            "original_source_label":data["source_original_lifeform_diaspore_code"],
            "observed_taxa":data["source_distinct_observed_taxa"],
            "source_site_taxon_incidences":data["source_guild_site_taxon_incidences"],
            "n_actual_source_sites_with_guild":data["source_sites_with_at_least_one_group_taxon"],
            "mode_results":{
                mode:{
                  "n_original_guild_evaluable_site_pairs":z["source_comparison_eligibility"]["source_guild_pair_comparisons_with_at_least_one_original_species"],
                  "double_zero_source_pairs_excluded":z["source_comparison_eligibility"]["source_pairs_excluded_because_both_sites_have_zero_original_guild_species"],
                  "rail_after_road_CV":z["original_source_heldout_rail_vs_road_prediction"]
                } for mode,z in data["richness_source_sensitivity"].items()
             }
        } for group,data in result["original_source_fully_known_4guild_partition"].items()
      },
      "never_pool_distinct_sources_or_claim_wind_selection":True
    },ensure_ascii=False))


if __name__=="__main__":main()
