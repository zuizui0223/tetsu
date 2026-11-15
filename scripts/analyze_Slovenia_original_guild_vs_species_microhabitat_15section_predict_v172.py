#!/usr/bin/env python3
"""v17.2 ORIGINAL Šabić et al. 2026 Slovenia 89 genuine railway 100m
field sites in 15 author line sections: test OUT-OF-SECTION transportability
of taxon identity vs life-history x wind-diaspore guild for where original
plants occur inside an actual railway microhabitat.

Actual same-location original BS=ballast shoulder and BT=ballast toe plots,
left/right source sides UNIONED WITHIN each genuine one of 89 field sites
(non-independent site sides not treated as separate cross-validation folds).
A species at a site is observed in exactly one state:
  BS_only / BT_only / both BS_and_BT.
Never manufacture an absent species as an observed BT/BS absence: every
analyzed site×species row has at least one real original observation.

AUTHENTIC publisher source SHA-locked (401,279 original bytes, source
DOI 10.1007/s44473-026-00155-z), 312 source taxa/traits; source-trait
known both lifespan and diaspore strategy for 261 genuine observed source
taxa across 89 sites. Exact original taxonomy, no imputed unknown trait.
Original adjacent microhabitat zone widths/areas NOT proven equal, no
matched nonrail control, no freight cargo/ontrain seed or actual mowing.

Does the original source taxon and its 4-guild source-assigned trait
archetype predict a *different original rail route section*?
  M0: global state frequencies on 14 training author line sections
  M1: M0 plus guild conditional state probability
  M2: M1 plus source species identity shrinkage toward guild.
All formulae fixed BEFORE inspecting this new holdout result:
  M0 Dirichlet(1,1,1)
  M1 counts guild + 10 × M0 pseudo-events
  M2 counts source species + 10 × M1 pseudo-events.
For source guild, strict_annual differs from other KNOWN lifespan, and
wind source class is published Phragmites/Epilobium diaspore STRATEGY
archetype, NOT a measured railway slipstream or true wind seed capture.

Primary heldout endpoint is 3 categorical states above, using log-loss
and multiclass Brier. This retains shared-BS/BT plant observations rather
than discarding them. Secondary is conditional only among observed
BS-only vs BT-only at a site; "both" is explicitly excluded and cannot be
a fake label. Both have 15 original railway section leaveout folds.
Report per original section (equal section weight), per site (equal site
weight within section), and occurrence-event-weighted score separately.
Never treat 2,729 repeated site×species incidences as 2,729 independently
randomized rail corridors, or the same original taxa across folds as
independent studies.

Species-level crossfold per-taxon score improvements and rank listings
are EXPLORATORY only: derived from same data, not unbiased tests of an
independent list of railway-selected species or a priori theta optimal.
Pre-adaptation by means of mean source TRAIT labels is assessed only as a
conditional predictive model, NOT shown to be a genetic cause or evolved
response. This study contains NO within-species phenotype measurements.
"""
from __future__ import annotations
import argparse,json,math
from collections import Counter,defaultdict
from pathlib import Path
from statistics import mean

from analyze_Slovenia_rail_microhabitats_v103 import (
    extract_original_xlsx, SOURCE_DOI,SOURCE_SHA256,validate_original_keys
)
from analyze_Slovenia_rail_species_occupancy_permutation_v109 import (
    original_location_matched_sets,taxa_preference
)
from analyze_Slovenia_source_rail_functional_taxonomic_identity_v141 import (
    source_sites_and_guilds
)
from select_two_original_rail_20_species_for_inat_v164 import (
    original_slovenian_EPPO_codebook_binomial_mapping
)

N_SITES=89
N_SECTIONS=15
N_JOINT_TRAIT_OBSERVED_TAXA=261
N_SITE_X_TAXON_INCIDENCES=2729
STATES=("BS_only","BT_only","both")
MODES=("all_observed_three_state","exclusive_BS_or_BT_only")
MODELS=("M0_global","M1_source_lifespan_by_diaspore_guild","M2_original_source_species_identity_shrunk_to_guild")
GUILD_PSEUDO=10.
SPECIES_PSEUDO=10.
GLOBAL_CATEGORY_PSEUDO=1.
N_BOOTSTRAP=0  # Not a clustered bootstrap; 15 source folds reported openly.

def author_source_rows(original):
    source_sites,trait_codes,original_audit=source_sites_and_guilds(original)
    actual,publisher_traits,integrity=original_location_matched_sets(original)
    code_to_name,cb_quality=original_slovenian_EPPO_codebook_binomial_mapping(original)
    if len(source_sites)!=N_SITES or len(actual)!=N_SITES or len(trait_codes)<N_JOINT_TRAIT_OBSERVED_TAXA:
        raise ValueError("Source original 89 field sites or traits missing")
    # original site indices are private to this in-memory computation.
    source_by_location={x["site"]:x for x in source_sites}
    if set(source_by_location)!={str(x) for x in actual}:
        raise ValueError("Original published source site identifiers changed")
    records=[]
    for site_id,obj in actual.items():
        site=source_by_location[str(site_id)]
        sBS=obj["BS"]; sBT=obj["BT"]
        if site["section"]!=str(obj["section_number"]):
            raise ValueError("One site falsely placed into two source rail sections")
        known=(sBS|sBT)&set(trait_codes)
        if known!=site["known"]:
            raise ValueError("Original source taxon incidence of location not identical to source v14.1")
        for taxon in sorted(known):
            state=("both" if taxon in sBS and taxon in sBT else
                   "BS_only" if taxon in sBS else "BT_only")
            records.append({
                "internal_original_site":str(site_id),
                "author_original_section":site["section"],
                "source_taxon_EPPO":taxon,
                "original_taxon_name_binomial_if_codebook_unambiguous":code_to_name.get(taxon),
                "source_four_lifespan_diaspore_guild":trait_codes[taxon],
                "observed_original_microhabitat_state":state,
            })
    if (len(records)!=N_SITE_X_TAXON_INCIDENCES or
        len({x["internal_original_site"] for x in records})!=N_SITES or
        len({x["author_original_section"] for x in records})!=N_SECTIONS or
        len({x["source_taxon_EPPO"] for x in records})!=N_JOINT_TRAIT_OBSERVED_TAXA or
        any(x["observed_original_microhabitat_state"] not in STATES for x in records)):
        raise ValueError("Original SI BS/BT 89×261 source 2729 observed site×taxon incidence drifted")
    if set(x["source_four_lifespan_diaspore_guild"] for x in records)!={
       "strict_annual_wind","strict_annual_other","other_known_wind","other_known_other"
    }:
        raise ValueError("Original four guild source strategy/lifespan ontology changed")
    prior_taxa,_=taxa_preference(actual,publisher_traits)
    if {x["taxon_code"] for x in prior_taxa}!={x["source_taxon_EPPO"] for x in records}:
        raise ValueError("Original v10.9 fully source-known source species taxonomy not reproduced")
    return records,{
       "source_DOI":SOURCE_DOI,"source_original_XLSX_SHA256":SOURCE_SHA256,
       "n_original_field_sites":N_SITES,"n_original_railway_line_sections":N_SECTIONS,
       "n_actual_taxa_with_both_author_lifespan_and_diaspore_source_codes":N_JOINT_TRAIT_OBSERVED_TAXA,
       "n_actual_original_site_taxon_presence_records":len(records),
       "n_source_taxon_name_mapped_to_explicit_binomial":sum(
          x["source_taxon_EPPO"] in code_to_name for x in prior_taxa),
       "original_author_species_codebook_QA":cb_quality,
       "n_original_BS_only_taxon_site_events":sum(x["observed_original_microhabitat_state"]=="BS_only" for x in records),
       "n_original_BT_only_taxon_site_events":sum(x["observed_original_microhabitat_state"]=="BT_only" for x in records),
       "n_original_taxon_site_BOTH_BS_BT_events":sum(x["observed_original_microhabitat_state"]=="both" for x in records),
       "original_taxon_site_events_cannot_infer_tillage_or_railway_true_land_area":True,
    }

def mode_data(rows,mode):
    if mode not in MODES:
        raise ValueError("Unrecognized biological source selection contrast")
    selected=[r for r in rows
              if mode=="all_observed_three_state" or r["observed_original_microhabitat_state"]!="both"]
    return selected, STATES if mode=="all_observed_three_state" else STATES[:2]

def source_smoothed_model(train,states):
    n_states=len(states)
    count=Counter(r["observed_original_microhabitat_state"] for r in train)
    assert set(count)<=set(states)
    n=sum(count.values())
    if n<=0:
        raise ValueError("No true original BS/BT source training observations")
    p0={s:(count[s]+GLOBAL_CATEGORY_PSEUDO)/(n+GLOBAL_CATEGORY_PSEUDO*n_states)
        for s in states}
    gc=defaultdict(Counter)
    sc=defaultdict(Counter)
    species_guild={}
    for r in train:
        state=r["observed_original_microhabitat_state"]
        guild=r["source_four_lifespan_diaspore_guild"]
        taxon=r["source_taxon_EPPO"]
        if taxon in species_guild and species_guild[taxon]!=guild:
            raise ValueError("A source species cannot change fixed published trait guild between sites")
        species_guild[taxon]=guild
        gc[guild][state]+=1
        sc[taxon][state]+=1
    pg={}
    for g,counts in gc.items():
        total=sum(counts.values())
        pg[g]={s:(counts[s]+GUILD_PSEUDO*p0[s])/(total+GUILD_PSEUDO)
               for s in states}
    return {"global":p0,"guild":pg,"taxon":sc,
            "source_taxon_guild":species_guild,"states":states}

def predict(row,fit):
    guild=row["source_four_lifespan_diaspore_guild"]
    taxon=row["source_taxon_EPPO"]
    p0=fit["global"]
    pg=fit["guild"].get(guild,p0)
    if taxon in fit["source_taxon_guild"] and fit["source_taxon_guild"][taxon]!=guild:
        raise ValueError("Original heldout author's source trait code differs from training")
    sc=fit["taxon"].get(taxon,Counter())
    total=sum(sc.values())
    ps={state:(sc[state]+SPECIES_PSEUDO*pg[state])/(total+SPECIES_PSEUDO)
        for state in fit["states"]}
    for p in (p0,pg,ps):
        if not math.isclose(sum(p.values()),1.,abs_tol=1e-10) or min(p.values())<=0:
            raise ValueError("Invalid genuine source multinomial predictive probability")
    return {
      MODELS[0]:p0,MODELS[1]:pg,MODELS[2]:ps,
      "original_heldout_species_has_prior_original_train_observations":total>0,
      "n_training_source_taxon_zone_events":total
    }

def scores(label,p,states):
    if label not in states:
        raise ValueError("Source observation category absent from original stratified response")
    logloss=-math.log(p[label])
    brier=sum((p[s]-int(s==label))**2 for s in states)/len(states)
    return logloss,brier

def score_one_holdout(train,test,states):
    fit=source_smoothed_model(train,states)
    events=[]; unseen=0
    for r in test:
        pred=predict(r,fit)
        if not pred["original_heldout_species_has_prior_original_train_observations"]:
            unseen+=1
        event={
          "internal_original_site":r["internal_original_site"],
          "source_taxon_EPPO":r["source_taxon_EPPO"],
          "original_taxon_name_binomial_if_codebook_unambiguous":
             r["original_taxon_name_binomial_if_codebook_unambiguous"],
          "source_guild":r["source_four_lifespan_diaspore_guild"],
          "state":r["observed_original_microhabitat_state"],
          "previously_observed_same_taxon_in_OTHER_source_sections":
             pred["original_heldout_species_has_prior_original_train_observations"],
          "out_of_original_source_section_log_loss":{},
          "out_of_original_source_section_multiclass_Brier":{},
        }
        for name in MODELS:
            ll,b=scores(event["state"],pred[name],states)
            event["out_of_original_source_section_log_loss"][name]=ll
            event["out_of_original_source_section_multiclass_Brier"][name]=b
        events.append(event)
    if not events:
        return None,[],0
    grouped=defaultdict(list)
    for e in events:grouped[e["internal_original_site"]].append(e)
    site_loss={
       metric:{
          m:mean(mean(e[metric][m] for e in siteevents)
                 for siteevents in grouped.values())
          for m in MODELS
       }
       for metric in ("out_of_original_source_section_log_loss",
                      "out_of_original_source_section_multiclass_Brier")
    }
    event_loss={
       metric:{m:mean(e[metric][m] for e in events) for m in MODELS}
       for metric in ("out_of_original_source_section_log_loss",
                      "out_of_original_source_section_multiclass_Brier")
    }
    return {
       "n_source_original_test_sites":len(grouped),
       "n_source_original_test_site_x_taxon_events":len(events),
       "n_test_events_taxon_not_seen_in_other_sections":unseen,
       "site_mean_equal_weight_scores":site_loss,
       "event_weight_scores":event_loss,
       "heldout_source_state_frequencies":dict(Counter(e["state"] for e in events))
    },events,unseen

def interval_from_fold_omission(heldout,metric="out_of_original_source_section_log_loss"):
    if len(heldout)!=N_SECTIONS:
        raise ValueError("Not 15 original site section folds")
    contributions=[(r["original_source_rail_section"],
                    r["site_mean_equal_weight_scores"][metric][MODELS[1]]-
                    r["site_mean_equal_weight_scores"][metric][MODELS[2]])
                   for r in heldout]
    n=len(contributions)
    total=sum(x for _,x in contributions)
    remaining=[{"omitted_source_section":s,
                "M2_source_species_shrink_better_than_guild_only_remaining_nats":
                (total-x)/(n-1)}
               for s,x in contributions]
    return {
      "original_equal_section_source_score_gain_M1_minus_M2_nats":
        total/n,
      "n_original_scored_fold_omissions_maintaining_positive_species_added_gain":
        sum(x["M2_source_species_shrink_better_than_guild_only_remaining_nats"]>0 for x in remaining),
      "min_after_removing_one_scored_section":
        min(x["M2_source_species_shrink_better_than_guild_only_remaining_nats"] for x in remaining),
      "max_after_removing_one_scored_section":
        max(x["M2_source_species_shrink_better_than_guild_only_remaining_nats"] for x in remaining),
      "original_section_scores_excluded_without_MODEL_RETRAINING":True,
      "site_species_pairs_are_NOT_statistically_independent":True,
      "source_remaining_section_details":remaining
    }

def evaluate_study(rows,mode):
    train_rows,states=mode_data(rows,mode)
    sections=sorted({r["author_original_section"] for r in rows})
    if len(sections)!=N_SECTIONS:
        raise ValueError("Source real author 15 line sections unavailable for holdout")
    holdouts=[];events=[]
    for section in sections:
        train=[r for r in train_rows if r["author_original_section"]!=section]
        test=[r for r in train_rows if r["author_original_section"]==section]
        original_score,scored,unseen=score_one_holdout(train,test,states)
        if original_score is None:
            raise ValueError("Original section has ZERO source eligible events in "+mode)
        holdouts.append({"original_source_rail_section":section,**original_score})
        events.extend(scored)
    if len(events)!=len(train_rows):
        raise ValueError("Some original source site-taxonomy events lost from heldout scoring")
    metric_names=("out_of_original_source_section_log_loss",
                  "out_of_original_source_section_multiclass_Brier")
    scoring={
      met:{
        "original_source_section_equal_site_equal_section":{m:mean(
            sec["site_mean_equal_weight_scores"][met][m] for sec in holdouts
        ) for m in MODELS},
        "original_source_event_weighted":{m:mean(e[met][m] for e in events)
                                         for m in MODELS},
        "original_source_site_equal_weight":{m:mean(
           mean(e[met][m] for e in rows_site) for rows_site in
             [ [ev for ev in events if ev["internal_original_site"]==sid]
               for sid in sorted({ev["internal_original_site"] for ev in events}) ]
        ) for m in MODELS},
      } for met in metric_names
    }
    # True original out-of-section scores by species; species-level summaries
    # may have n=1 and are exploratory NOT independently validated species tests.
    by_taxon=defaultdict(list)
    for e in events:by_taxon[e["source_taxon_EPPO"]].append(e)
    taxon_diag=[]
    for t,trs in by_taxon.items():
        mean_ll={model:mean(x["out_of_original_source_section_log_loss"][model]
                            for x in trs) for model in MODELS}
        taxon_diag.append({
          "source_taxon_name":trs[0]["original_taxon_name_binomial_if_codebook_unambiguous"],
          "source_taxon_original_EPPO_code":t,
          "author_original_lifespan_by_diaspore_guild":trs[0]["source_guild"],
          "n_original_heldout_test_site_x_taxon_events":len(trs),
          "n_original_source_rail_sections_with_real_taxon_observation":len({
               x["author_original_section"] for x in train_rows
               if x["source_taxon_EPPO"]==t
          }),
          "n_original_source_test_events_with_species_seen_in_other_training_sections":
             sum(x["previously_observed_same_taxon_in_OTHER_source_sections"] for x in trs),
          "source_mean_out_of_section_negative_logloss_by_fixed_model":mean_ll,
          "M2_minus_M1_mean_heldout_logloss":mean_ll[MODELS[2]]-mean_ll[MODELS[1]],
          "only_an_EXPLORATORY_species_priority_score_NOT_adaptation":True
        })
    taxon_diag.sort(key=lambda x:(
        -x["n_original_heldout_test_site_x_taxon_events"],
        x["source_taxon_original_EPPO_code"]))
    output_species=[
      x for x in taxon_diag
      if x["n_original_heldout_test_site_x_taxon_events"]>=10
      and x["source_taxon_name"]
    ][:25]
    main=scoring["out_of_original_source_section_log_loss"]["original_source_section_equal_site_equal_section"]
    base=main[MODELS[0]]
    return {
       "source_response_mode":mode,
       "n_real_original_observed_source_site_x_species_model_events":len(train_rows),
       "n_original_unique_sites_source":len({r["internal_original_site"] for r in train_rows}),
       "n_original_source_sections":len(holdouts),
       "n_original_authors_distinct_species_original_observed":len(by_taxon),
       "n_heldout_events_for_taxa_unseen_in_other_training_sections":sum(
           sec["n_test_events_taxon_not_seen_in_other_sections"] for sec in holdouts),
       "n_original_BS_only_BT_only_BOTH_events":dict(Counter(r["observed_original_microhabitat_state"]
                                                          for r in train_rows)),
       "n_proper_source_observed_count_site_events_NOT_261_independent_species_studies":True,
       "model_smoothing":{"global_dirichlet_per_state":GLOBAL_CATEGORY_PSEUDO,
                          "guild_source_effective_prior_n":GUILD_PSEUDO,
                          "species_identity_effective_prior_n":SPECIES_PSEUDO},
       "original_source_study_holdout_scores":scoring,
       "source_equal_section_improvement_percents_LOGLOSS":{
          "guild_M1_vs_global_M0":100*(main[MODELS[0]]-main[MODELS[1]])/base,
          "species_M2_vs_guild_M1":100*(main[MODELS[1]]-main[MODELS[2]])/main[MODELS[1]],
          "species_M2_vs_global_M0":100*(main[MODELS[0]]-main[MODELS[2]])/base,
       },
       "source_n_sections_M1_guild_beats_global_M0":sum(
          sec["site_mean_equal_weight_scores"]["out_of_original_source_section_log_loss"][MODELS[1]]<
          sec["site_mean_equal_weight_scores"]["out_of_original_source_section_log_loss"][MODELS[0]]
          for sec in holdouts),
       "source_n_sections_M2_species_beats_guild_M1":sum(
          sec["site_mean_equal_weight_scores"]["out_of_original_source_section_log_loss"][MODELS[2]]<
          sec["site_mean_equal_weight_scores"]["out_of_original_source_section_log_loss"][MODELS[1]]
          for sec in holdouts),
       "original_heldout_section_scores":holdouts,
       "source_scored_fold_omission_only_M2_beyond_guild_influence":
           interval_from_fold_omission(holdouts),
       "source_ranked_high_frequency_species_outofsample_descriptive_NOT_global_rail_preference":
          output_species,
       "not_a_measurement_of_within_species_genetic_or_plastic_adaptation":True,
       "not_an_original_nonrail_community_control":True,
       "no_independent_seed_traffic_mowing_wind_or_recruitment":True
    }

def run(original):
    rows,source=author_source_rows(original)
    modes={mode:evaluate_study(rows,mode) for mode in MODES}
    return {
      "complete":True,
      "version":"v17.2_SI_original_microhabitat_species_vs_guild_15route_heldout",
      "authentic_original_botanist_source":source,
      "source_true_observed_89site_261taxon_15section_data":True,
      "holdout_species_vs_guild_two_distinct_response_definitions":modes,
      "study_limitations":[
        "Only source 2022–23 Slovenia railway BS vs BT within same sampled 100m rail field location, not genuine RAIL-vs-ROAD or RAIL-vs-nonrail control.",
        "Original source BS vs BT microhabitat sampling widths/areas not standardized and detection/sampling nonrandom; no causal selection coefficient.",
        "Published Phragmites/Epilobium diaspore.strategy archetypes are species-level labels, not measured individual dispersal, train slipstream or railway cargo seed deposition.",
        "The source 89 original sites within 15 selected railway route sections, with repeated taxa across sites and folds; no independent 2729 field replicates or p-values.",
        "Model M1 and M2 smoothing constants pre-fixed, with no tuning or posthoc null p-values; improvements are only out-of-section conditional prediction.",
        "Conditional exclusive-only BS-vs-BT endpoint excludes original plants present BOTH habitats. Main three-state retains both; divergent conclusions should remain documented.",
        "Species taxon identity shrunk model, even if predictive, does not demonstrate inherited genetic adaptation, within-species traits, plasticity, railway-only occupancy or biological trait optimum.",
        "High-frequency species per-taxon posthoc holdout priority list selected for descriptive future source study only, not independent test.",
      ],
      "NOT_estimated_species_within_trait_change_or_evolved_adaptation":True,
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--original-xlsx",required=True)
    p.add_argument("--source-audit102",required=True)
    p.add_argument("--output",required=True)
    a=p.parse_args()
    raw=extract_original_xlsx(a.original_xlsx,a.source_audit102)
    result=run(raw)
    dst=Path(a.output);dst.parent.mkdir(parents=True,exist_ok=True)
    dst.write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding="utf-8")
    summary={
      "version":result["version"],
      "true_source":result["authentic_original_botanist_source"],
      "outcomes":{name:{
         "n_real_original_source_events":part["n_real_original_observed_source_site_x_species_model_events"],
         "n_taxa":part["n_original_authors_distinct_species_original_observed"],
         "n_source_sites":part["n_original_unique_sites_source"],
         "n_source_sections":part["n_original_source_sections"],
         "n_original_BS_BT_and_BOTH":part["n_original_BS_only_BT_only_BOTH_events"],
         "score_equal_source_sections_and_sites":
            part["original_source_study_holdout_scores"]["out_of_original_source_section_log_loss"]["original_source_section_equal_site_equal_section"],
         "effect_guild_then_species_vs_baseline_pct":
            part["source_equal_section_improvement_percents_LOGLOSS"],
         "n_section_guild_improves":part["source_n_sections_M1_guild_beats_global_M0"],
         "n_section_species_improves":part["source_n_sections_M2_species_beats_guild_M1"],
         "species_over_guild_scored_heldout_deletion_sensitivity":
            {key:value for key,value in part["source_scored_fold_omission_only_M2_beyond_guild_influence"].items()
             if key!="source_remaining_section_details"},
         "top_stable_source_species_research_candidates":
            part["source_ranked_high_frequency_species_outofsample_descriptive_NOT_global_rail_preference"][:15]
       } for name,part in result["holdout_species_vs_guild_two_distinct_response_definitions"].items()},
      "no_source_genetic_adaptation_or_actual_freight_seed_arrival":True
    }
    print(json.dumps(summary,ensure_ascii=False))

if __name__=="__main__":main()
