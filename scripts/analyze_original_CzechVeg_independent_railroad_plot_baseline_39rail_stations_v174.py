#!/usr/bin/env python3
"""v17.4 independent original CzechVeg population baseline for the 39 actual
Kutlvašr 2020-22 Czech railway STATION observed plant floras.

DO NOT compute naive 39 whole-station/117,739 small-plot enrichment as a
true ecological railway-selection coefficient. Original station survey
grain includes multiple habitat types over hectares, while CzechVeg
relevés are 1-1000m2. 2020-22 rail station geographic location GPS
is NOT in author's original S1/S2, so no source rail-station-to-background
local neighborhood claim can be made. Instead:

1) Preserve as independent CONTEXT ONLY full-country CzechVeg-Open 2026
   original 117,739 survey headers / 2,510,340 site×taxon source records.
   For each original Czech botanical railway station species, record its
   frequency in all source area-valid botanical vascular relevés, and
   note why this prevalence cannot be directly divided into station
   occupancy.
2) Within CzechVeg's OWN harmonized real source PLOT observational grain,
   compare textual RAIL vs ROAD vegetation plot labels within same
   PUBLISHED STUDY TITLE, phytogeographic region, sampling decade and
   measured plot AREA CLASS. Frozen v12.3 SOURCE S2 rule requires >=2
   source RAIL and ROAD plots per block. Compute for each taxon a
   block-standardized (sum w*(nrail_sp/nrail - nroad_sp/nroad)/sum w)
   where w=min(nrail,nroad), and use only original botanical records.
   This describes railway-associated TEXT labels, NOT GPS-verified active
   railway right-of-way or causal response to train/vegetation cutting.
3) Independent station-Preslia source-2020–22 39 actual stations floras
   are joined by EXACT botanical normalized Latin name, never fuzzy or
   genus-level taxonomy, and retain actual number of the 39 stations
   occupied. Candidates selected using CzechVeg source matched q, NOT
   photos or station success; require >=5 station presences, >=5 observed
   plot incidences (rail+road), and >=2 *different author source STUDIES*
   with original taxon observed within corresponding retained matched
   blocks. If fewer than20 candidate species eligible, report fewer
   rather than fill with guesses, negative or unmatched names.
4) Report 2010+ (recent) and 2020+ original survey-year matched
   panel as PREDECLARED sensitivity without silently using 1800s
   railway/road relevés to describe modern 2020-22 stations. Same
   publication×region×decade×area-class field coverage and min2
   each original plot arm required; not all years may have support.
5) Also strict only-in-source-HABITAT label as more conservative
   railway-context sensitivity (often severely underpowered).
   Source botanical 'road' is context text label not true plot matched
   track management. Never treat no label as 'nonrail' and never
   infer variation in train slipstream gust, freight seeds or genetics.

This is not actual railway preference of individual plants and NOT a
species-level causal rail preference. It is an independently BASELINED
conditional priority index for future PHOTO response cohort selection,
deliberately keeping plot grain, region & publication bias explicit.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,math
from collections import Counter,defaultdict
from pathlib import Path

from analyze_CzechVeg_Pladias_rail_road_diaspore_v119 import (
    verify_file,plot_labels,normalize_taxon,HEADERS_NAME,SPECIES_NAME
)
from analyze_CzechVeg_rail_road_region_year_area_overlap_v121 import (
    header_source_context
)
from analyze_CzechVeg_rail_road_same_study_trait_v123 import (
    source_study_meta,common_original_study_blocks,SPECS
)
from analyze_Czech_2024_39station_functional_identity_v143 import (
    read_original_czech_publisher,PAPER
)

ZENO_DOI="10.5281/zenodo.19231152"
SOURCE_CZECHVEG_HEADERS=117739
SOURCE_CZECHVEG_PLANT_ROWS=2510340
N_SOURCE_STATIONS=39
N_SOURCE_STATION_TAXA=761
N_SOURCE_STATION_SITE_OCCURRENCES=6941
MATCHING_S2=next(x for x in SPECS if x[0]=="S2_study_region_decade_area4_min2")
FRAMES=(
   ("BROAD_original_same_STUDY_region_DECADE_area_class_all_years","broad",None),
   ("BROAD_original_same_STUDY_region_DECADE_area_class_2010_onward","broad",2010),
   ("BROAD_original_same_STUDY_region_DECADE_area_class_2020_onward","broad",2020),
   ("STRICT_HABITAT_original_same_STUDY_region_DECADE_area_class_all_years","strict",None),
)
# source selection excludes photo outcomes by design.
MIN_ORIGINAL_STATION_OCCURRENCES=5
MIN_BACKGROUND_MATCHED_PLOT_SPECIES_INCIDENCES=5
MIN_DISTINCT_ORIGINAL_BOTANICAL_STUDIES_WITH_TAXON=2

def station_source(s1,s2):
    sites,meta,proof=read_original_czech_publisher(s1,s2)
    if (len(sites)!=N_SOURCE_STATIONS or
        len({t for st in sites.values() for t in st})!=N_SOURCE_STATION_TAXA or
        sum(map(len,sites.values()))!=N_SOURCE_STATION_SITE_OCCURRENCES or
        proof["original_source_Czech_railway_Preslia_DOI"]!=PAPER):
        raise ValueError("Original real botanist 39 railway stations / 761 taxa / 6941 incidence source not authentic")
    occ=Counter(t for st in sites.values() for t in st)
    if max(occ.values())>N_SOURCE_STATIONS:
        raise ValueError("An original author taxon was counted more than once in one railway station")
    return occ,{
      "original_source_DOI":PAPER,
      "n_genuine_original_botanist_station_survey_sites":len(sites),
      "n_distinct_original_current_station_source_taxa":len(occ),
      "n_source_station_x_taxon_presences":sum(occ.values()),
      "original_Preslia_S1_S2_source_verification":proof,
      "station_observations_2020_2022_are_NOT_matched_by_exact_railway_subplot_with_CzechVeg":True,
    }

def source_vegetation_scan(original_species_csv,valid_source_plot_meta,
                           matched_text_candidates,original_station_taxa):
    # Keep *all* source candidate plot vascular species for validating true
    # CzechVeg plot survey support. Store only requested 761 station taxon
    # background station comparison counts, not full 2.5m author matrix.
    name_counts=Counter()
    site_valid_vascular=set()
    original_text_plot_taxa=defaultdict(set)
    context_seen=set()
    context_hits=0
    scanned=0
    with open(original_species_csv,encoding="utf-8-sig",newline="") as f:
        reader=csv.DictReader(f)
        expected={"RELEVE_NR","KAPLAN_PLUS_NAME_SIMPLE","PLANT_GROUP"}
        if expected-set(reader.fieldnames or ()):
            raise ValueError("CzechVeg original botanical source plot/species header schema changed")
        for row in reader:
            scanned+=1
            if str(row.get("PLANT_GROUP") or "").strip().lower()!="vascular":
                continue
            id=str(row.get("RELEVE_NR") or "").strip()
            if id not in valid_source_plot_meta:
                continue
            taxon=normalize_taxon(row.get("KAPLAN_PLUS_NAME_SIMPLE"))
            if not taxon or taxon.lower() in ("nan","na","null"):
                continue
            site_valid_vascular.add(id)
            if id in matched_text_candidates:
                original_text_plot_taxa[id].add(taxon)
            if taxon in original_station_taxa:
                # Deduplicate exact ORIGINAL plot×taxon; do not count
                # multiple vegetation layers as repeated independent plots.
                key=(id,taxon)
                if key not in context_seen:
                    context_seen.add(key)
                    name_counts[taxon]+=1
                    context_hits+=1
    if scanned!=SOURCE_CZECHVEG_PLANT_ROWS:
        raise ValueError("Full CzechVeg original 2.51034 million plot-taxonomic source records not read")
    if not site_valid_vascular or len(site_valid_vascular)>len(valid_source_plot_meta):
        raise ValueError("Official source CzechVeg vascular site survey denominator invalid")
    return original_text_plot_taxa,name_counts,{
        "source_original_CzechVeg_site_x_taxon_CSV_rows_scanned":scanned,
        "source_original_plot_headers_with_recorded_plot_area_1_to_1000m2":
             len(valid_source_plot_meta),
        "source_original_valid_area_plots_with_actual_vascular_species_records":
             len(site_valid_vascular),
        "source_station_taxa_matched_original_CzechVeg_eligible_area_plot_by_species_incidence":
             context_hits,
        "source_original_published_study_text_labeled_valid_area_plots_with_vascular_records":
             len(original_text_plot_taxa),
        "original_full_country_background_site_x_species_duplicated_life_layers_removed":True,
        "original_full_country_context_CAN_NOT_be_ratioed_to_39_multi_zone_station_survey":True,
    },len(site_valid_vascular)

def retain_blocks(meta,plot_species,label,min_year):
    selected={key:value for key,value in meta.items()
              if min_year is None or value["source_year_int"]>=min_year}
    blocks,quality=common_original_study_blocks(
        selected,plot_species,label,MATCHING_S2)
    source_studies={key[0] for key in blocks}
    # no source outcome selected before the matching scheme
    return blocks,{
       **quality,
       "original_source_year_cutoff":min_year,
       "original_source_taxon_agnostic_matching_studies_with_both_real_labelled_arms":
           len(source_studies),
       "ORIGINAL_source_railway_and_road_labels_are_text_NOT_verified_right_of_way":True
    }

def standardized_taxon_panel(blocks,source_plot_taxa,station_occs,background,denom):
    # Pre-generate plot×names for each eligible original matched block.
    # For all 761 station source taxa, missing from a surveyed relevé
    # is a recorded non-detection, NOT confirmed biological absence.
    original_station_names=set(station_occs)
    block_profiles=[]
    source_taxon_total_rail=Counter()
    source_taxon_total_road=Counter()
    by_study=defaultdict(set)
    by_block=defaultdict(set)
    total_w=0
    all_rail=set();all_road=set()
    for k,gr in sorted(blocks.items()):
        rail=gr["rail"];road=gr["road"]
        nr=len(rail);nd=len(road)
        if nr<2 or nd<2:raise ValueError("Original source block did not preserve min2 railway and roadside plots")
        w=min(nr,nd);total_w+=w
        rc=Counter();dc=Counter()
        for plot in rail:
            if plot in all_rail or plot in all_road:raise ValueError("Original source vegetation plot reused in more than one rail-road match block")
            all_rail.add(plot)
            rc.update(source_plot_taxa[plot]&original_station_names)
        for plot in road:
            if plot in all_rail or plot in all_road:raise ValueError("Original source vegetation plot reused in more than one rail-road match block")
            all_road.add(plot)
            dc.update(source_plot_taxa[plot]&original_station_names)
        source_taxon_total_rail.update(rc)
        source_taxon_total_road.update(dc)
        for t in rc.keys()|dc.keys():
            by_study[t].add(k[0])
            by_block[t].add(k)
        block_profiles.append((k,nr,nd,w,rc,dc))
    if blocks and total_w<=0:
        raise ValueError("Cannot have original real rail-road matching blocks without source plot opportunity")
    weighted_num=Counter()
    study_weight=Counter()
    original_study_species_numerators=defaultdict(Counter)
    for k,nr,nd,w,rc,dc in block_profiles:
        study_weight[k[0]]+=w
        for t in rc.keys()|dc.keys():
            increment=w*(rc[t]/nr-dc[t]/nd)
            weighted_num[t]+=increment
            original_study_species_numerators[t][k[0]]+=increment
    species={}
    for name in station_occs:
        rail_n=source_taxon_total_rail[name]
        road_n=source_taxon_total_road[name]
        n_both=rail_n+road_n
        delta=(weighted_num[name]/total_w if total_w>0 and n_both else None)
        if delta is not None and not(-1-1e-10<=delta<=1+1e-10):
            raise ValueError("Original block-standardized taxon rail/road plot prevalence difference outside [-1,1]")
        # Posthoc influence, not independent validation: remove each whole
        # original CzechVeg author study's already matched plot blocks, and
        # recompute the exact weighted prevalence contrast. DO NOT infer
        # independent station-scale railway adaptation from the sign.
        omit=[]
        for pub,total_study_weight in sorted(study_weight.items()):
            remaining=total_w-total_study_weight
            if remaining<=0:continue
            score=(weighted_num[name]-original_study_species_numerators[name][pub])/remaining if n_both else None
            if score is not None:omit.append(score)
        positive_loo=sum(x>0 for x in omit if x is not None)
        species[name]={
          "station_original_2020_22_n_of_39":station_occs[name],
          "station_original_2020_22_occupancy_fraction":station_occs[name]/39,
          "CzechVeg_all_eligible_plot_background_species_n":background[name],
          "CzechVeg_all_eligible_plot_background_species_prevalence_CONTEXT_NOT_COMPARABLE_TO_STATION":
               background[name]/denom if denom else None,
          "CzechVeg_same_study_region_decade_area_source_N_RAIL_text_plot_with_this_taxon":rail_n,
          "CzechVeg_same_study_region_decade_area_source_N_ROAD_text_plot_with_this_taxon":road_n,
          "CzechVeg_same_study_source_N_distinct_publications_with_taxon":len(by_study[name]),
          "CzechVeg_same_study_source_N_matched_strata_with_taxon":len(by_block[name]),
          "CzechVeg_all_original_source_AUTHOR_STUDY_scored_loo_full_source_studies":len(omit),
          "CzechVeg_original_AUTHOR_STUDY_scored_loo_N_with_positive_rail_minus_road":positive_loo,
          "CzechVeg_original_AUTHOR_STUDY_scored_loo_min_delta":min(omit) if omit else None,
          "CzechVeg_original_AUTHOR_STUDY_scored_loo_max_delta":max(omit) if omit else None,
          "CzechVeg_original_AUTHOR_STUDY_scored_loo_POSITIVE_for_ALL_eligible_omissions":
              bool(omit and len(omit)==positive_loo and delta is not None and delta>0),
          "CzechVeg_source_loo_deletes_entire_original_publication_relevés_not_station_cohort":True,
          "CzechVeg_study_region_decade_area_STANDARDIZED_rail_minus_road_plot_prevalence":
               delta,
          "RAIL_text_labels_NOT_externally_validated_right_of_way":True,
          "station_occurrence_ratio_NOT_a_genuine_rail_affinity_effect":True
        }
    return species,{
       "n_authentic_same_study_plot_balance_blocks":len(blocks),
       "n_original_rail_text_plots":len(all_rail),
       "n_original_road_text_plots":len(all_road),
       "n_original_source_field_studies":len({k[0] for k in blocks}),
       "n_original_plot_pair_MIN_support_weight_mass":total_w,
       "n_original_station_source_taxa_with_CzechVeg_matched_rail_or_road_identified":sum(
           x["CzechVeg_study_region_decade_area_STANDARDIZED_rail_minus_road_plot_prevalence"] is not None
           for x in species.values()),
       "n_original_39station_source_taxa":len(species),
       "n_original_station_taxa_present_somewhere_in_CzechVeg_plot_source":
           sum(background[x]>0 for x in station_occs)
    }

def rank_independent_rail_context(main,modern=None,strict=None):
    eligible=[]
    for name,x in main.items():
        score=x["CzechVeg_study_region_decade_area_STANDARDIZED_rail_minus_road_plot_prevalence"]
        if (score is None or score<=0 or
            x["station_original_2020_22_n_of_39"]<MIN_ORIGINAL_STATION_OCCURRENCES or
            x["CzechVeg_same_study_region_decade_area_source_N_RAIL_text_plot_with_this_taxon"]+
            x["CzechVeg_same_study_region_decade_area_source_N_ROAD_text_plot_with_this_taxon"]<
                  MIN_BACKGROUND_MATCHED_PLOT_SPECIES_INCIDENCES or
            x["CzechVeg_same_study_source_N_distinct_publications_with_taxon"]<
                MIN_DISTINCT_ORIGINAL_BOTANICAL_STUDIES_WITH_TAXON):
            continue
        z={
          "scientific_name_original_exact_normalized":name,
          "station_author_present_n_of39":x["station_original_2020_22_n_of_39"],
          "independent_CzechVeg_same_publication_region_decade_area_plot_STANDARDIZED_rail_minus_road_text_label_prevalence":score,
          "n_original_studies_with_taxon_source":x["CzechVeg_same_study_source_N_distinct_publications_with_taxon"],
          "n_matched_original_plots_with_taxon_RAIL_text":x["CzechVeg_same_study_region_decade_area_source_N_RAIL_text_plot_with_this_taxon"],
          "n_matched_original_plots_with_taxon_ROAD_text":x["CzechVeg_same_study_region_decade_area_source_N_ROAD_text_plot_with_this_taxon"],
          "source_original_author_study_LOO_sign_always_positive":
             x["CzechVeg_original_AUTHOR_STUDY_scored_loo_POSITIVE_for_ALL_eligible_omissions"],
          "source_original_author_study_LOO_min_rail_minus_road":
             x["CzechVeg_original_AUTHOR_STUDY_scored_loo_min_delta"],
          "source_original_author_study_LOO_N_positive":
             x["CzechVeg_original_AUTHOR_STUDY_scored_loo_N_with_positive_rail_minus_road"],
          "source_original_author_study_LOO_N_evaluable":
             x["CzechVeg_all_original_source_AUTHOR_STUDY_scored_loo_full_source_studies"],
          "modern_2010plus_source_same_plot_type_score":(
             modern.get(name,{}).get(
                "CzechVeg_study_region_decade_area_STANDARDIZED_rail_minus_road_plot_prevalence")
             if modern else None),
          "strict_source_HABITAT_only_plot_score":(
             strict.get(name,{}).get(
                "CzechVeg_study_region_decade_area_STANDARDIZED_rail_minus_road_plot_prevalence")
             if strict else None),
          "CzechVeg_observational_proxy_NOT_independently_confirmed_active_railway_affinity":True,
        }
        eligible.append(z)
    eligible.sort(key=lambda r:(
        -r["independent_CzechVeg_same_publication_region_decade_area_plot_STANDARDIZED_rail_minus_road_text_label_prevalence"],
        -r["n_original_studies_with_taxon_source"],
        r["scientific_name_original_exact_normalized"]))
    return eligible[:20],{
      "n_eligible_station_taxa_with_positive_independent_CzechVeg_railroad_proxy_and_source_support":len(eligible),
      "source_selected_top20_max_not_forced_to_fill":min(20,len(eligible)),
      "new_selection_does_NOT_depend_on_iNaturalist_photograph_counts_or_flower_traits":True,
      "original_same_study_full_country_observational_rail_TEXT_NOT_railway_environment_causality":True,
    }

def run(zenodo,header_csv,species_csv,s1,s2):
    original_head=verify_file(zenodo,HEADERS_NAME,header_csv)
    original_sp=verify_file(zenodo,SPECIES_NAME,species_csv)
    station,station_proof=station_source(s1,s2)
    meta,plotqc=plot_labels(header_csv)
    withyear,yearqc=header_source_context(header_csv,meta)
    withpaper,studyqc=source_study_meta(header_csv,withyear)
    candidate={id for id,m in withpaper.items()
               if m["broad"] in ("rail","road") or m["strict"] in ("rail","road")}
    plot_taxa,background,scan,denom=source_vegetation_scan(
        species_csv,meta,candidate,set(station))
    result={};source_plots={}
    for name,label,yearcut in FRAMES:
        blocks,summary=retain_blocks(withpaper,plot_taxa,label,yearcut)
        species,orig=standardized_taxon_panel(blocks,plot_taxa,station,background,denom)
        result[name]=species
        source_plots[name]={**summary,**orig}
    original=FRAMES[0][0]
    panel=result[original]
    ranking,elig=rank_independent_rail_context(
        panel,
        modern=result[FRAMES[1][0]],
        strict=result[FRAMES[3][0]]
    )
    base_q=source_plots[original]
    # Source v12.3 frozen common study exact area4/decade min2
    # observed real 44 blocks 229 RAIL and 462 ROAD for broad labels;
    # preserve the field cohort and do not silently use all 117k plots.
    if (base_q["n_original_rail_text_plots"]!=229 or
        base_q["n_original_road_text_plots"]!=462 or
        base_q["n_authentic_same_study_plot_balance_blocks"]!=44):
        raise ValueError("Independent original CzechVeg v12.3 common study/decade/area panel did not reproduce 229 rail /462 road /44 blocks")
    if source_plots[FRAMES[3][0]]["n_original_rail_text_plots"]!=20 or source_plots[FRAMES[3][0]]["n_original_road_text_plots"]!=25:
        raise ValueError("Original CzechVeg STRICT only HABITAT current same-study matched 20 rail /25 road did not reproduce")
    return {
      "complete":True,
      "version":"v17.4_CzechVeg_original_plot_baselined_39Preslia_station_taxa_railroad_source_proxy",
      "independent_original_botanist_sources":[
         "10.5281/zenodo.19231152","10.23855/preslia.2024.247"],
      "CzechVeg_full_original_source_HEADER_MD5_verified":original_head,
      "CzechVeg_full_original_source_2510340_species_MD5_verified":original_sp,
      "original_39station_publisher_Czech_flora_verified":station_proof,
      "original_CzechVeg_full_country_plot_and_selected_taxa_source_counts":scan,
      "original_CzechVeg_plot_metadata_and_study_checks":{
         "source_original_plot_habitat_label_quality":plotqc,
         "source_original_header_area_year_phyto_region":yearqc,
         "source_original_published_author_study_reference":studyqc
      },
      "source_original_observational_rail_vs_road_comparisons":source_plots,
      "n_station_source_taxa":len(station),
      "n_unique_publisher_stations":39,
      "source_original_station_taxa_independent_CzechVeg_reference_context_per_species":result,
      "independent_SOURCE_TAXON_RAIL_TEXT_CONTEXT_ranked_TOP20_NOT_true_railroad_affinity":ranking,
      "source_eligibility_criteria_AND_no_forced_rank20":elig,
      "NOT_station_frequency_divided_by_incompatible_background_plot_prevalence":True,
      "NOT_station_plots_geo_matched_to_CzechVeg_vegetation_plots":True,
      "original_CzechVeg_2026_vegetation_plots_and_2020_22_39whole_station_flora_differ_in_sampling_grain":True,
      "predefined_2010_2020_recent_and_strict_HABITAT_source_sensitivity_retained_even_if_negative":True,
      "limits":[
         "The 39 original Preslia station survey units aggregate many habitats over hectares, incomparable to single 1–1000m2 CzechVeg vascular vegetation relevés; no station/background prevalence ratio or causal rail preference inferred.",
         "Only CzechVeg's OWN plot-grain source railway/roadside text-labeled relevés matched within original published study×phytogeographic-region×survey decade×area class are used to rank observational RAIL context. Text not actual validated right-of-way.",
         "Same original published study reference does not guarantee identical sampling observer/effort, or two plots assigned railway exposure. Urbanization, maintenance and reporting bias remain.",
         "CzechVeg modern 2010+ and 2020+ are separately reported; older plots cannot be assumed contemporaneous with Preslia 2020–22 actual flora.",
         "A positive adjusted plot prevalence is observational presence enrichment in labelled railway *context*, not within-species genetic adaptation, mowing, survival or train seed dispersal.",
         "The original Preslia taxa are used as a station filter and original CzechVeg as independently measured plot context; each originates from different independent field authors, NOT original railway-vs-nonrail paired station plots.",
         "Species source list ranked by CzechVeg conditional observational rail text plot prevalence, station occupancy only minimum support; cannot infer optimum theta rail or inherited phenotypic traits.",
         "The CzechVeg survey species CSV contains 2,510,340 original rows but possible multiple taxonomic and layer records, collapsed to unique original plot×species; 2,510,340 is not number of independent plant sampling sites.",
         "Strict habitat text-only matched source study panel has only 20 rail and 25 road plots, likely underpowered to rank 761 station taxa; missing coverage is not negative preference."
      ]
    }

def main():
    p=argparse.ArgumentParser()
    for arg in ("zenodo-record","source-CzechVeg-headers","source-CzechVeg-species",
                "Preslia-original-S1","Preslia-original-S2","output"):
        p.add_argument("--"+arg,required=True)
    a=p.parse_args()
    obj=run(json.loads(Path(a.zenodo_record).read_text(encoding="utf-8")),
            a.source_CzechVeg_headers,a.source_CzechVeg_species,
            a.Preslia_original_S1,a.Preslia_original_S2)
    out=Path(a.output);out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(obj,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps({
      "version":obj["version"],
      "original39_rail_stations":obj["n_unique_publisher_stations"],
      "source_station_taxa":obj["n_station_source_taxa"],
      "source_full_CzechVeg_plot_and_species":obj["original_CzechVeg_full_country_plot_and_selected_taxa_source_counts"],
      "source_4_original_study_area_year_frame_support":obj["source_original_observational_rail_vs_road_comparisons"],
      "source_eligible_taxa_count":obj["source_eligibility_criteria_AND_no_forced_rank20"],
      "source_ranked_candidate_taxa":obj["independent_SOURCE_TAXON_RAIL_TEXT_CONTEXT_ranked_TOP20_NOT_true_railroad_affinity"],
      "no_biologically_measured_rail_vs_nonnear_rail_adaptation":True
    },ensure_ascii=False))

if __name__=="__main__":main()
