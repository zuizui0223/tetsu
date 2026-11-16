#!/usr/bin/env python3
"""v23.8 CzechVeg source-original functional convergence vs independent text-nontransport controls.

SOURCE AUTHENTICITY:
 CzechVeg-Open v1.1 true 2026 original 117739 relevé metadata (MD5) and
 2,510,340 ORIGINAL species x plot incidences (MD5); unchanged source.
 LEDA Oldenburg ORIGINAL canopy height and seed mass (SHA256 locked).

Use SAME rail-text plot and SAME Czech source phytogeographic region,
survey decade, original sampled plot-area band; choose separate-original-
publication nontransport-text controls. Physical matching: measured
ALTITUDE within 100 m, independent of vegetation composition and SYNTAXON.
Second reference: all original nontransport-text plots on the SAME rail
and region support. The latter is not confirmed natural or rail-free;
the broad original "rail" label is from original HABITAT or LOCALITY
text, NOT verified operational trackbed. Reference sites are independent
field plots and studies *within a compiled multi-study vegetation DB*,
NOT independent-country ecological replication. No rail causes proven.

Source opportunity from existing v19.5 exact source code:
 - get_meta() verifies 76,386,498-byte 117,739-row CzechVeg headers MD5.
 - source_arms() pre-selects controls via original region, decade, area,
   other publication, and measured altitude, BEFORE species content.
 - source_species() verifies 382,939,947-byte 2,510,340-row species CSV MD5.
 - LEDA exact species-level medians; no genus-level filling or synthetic
   mean values. Reject aggregate/affinity taxa for 2D traits.

Output:
 Per source region sample one RAIL plot and SAME RAIL-eligible paired
 reference (elevation matched vs all nontransport) 399 technical draws.
 Compare pairwise original-region CWM Euclidean functional distance,
 same regional species Jaccard, and mean within-plot 2D trait variance.
 Each region is one draw unit, pairwise region distances are dependent.
 All comparisons use exactly SAME actual rail IDs and source regions.

Panels:
  source_5plus_constrained (>=5 jointly trait known species per plot)
  source_richness_matched (original species-count 2/3..1.5 for controls)
  source_recent_survey2010_era (eligible >2010 original year, only if support)
  source_HABITAT_only_rail_label (strict label only, only if support).
 Hard stop on >=3 original source phyto regions, min 5 original source
 reference opportunities for each eligible RAIL plot after trait QC.
 Do not create a fake zero-effect result if a panel lacks support.

No direct contemporaneous train counts, substrate/mowing/wind exposure,
direct seed carriage, identical station-level botanical survey plots,
or rail adaptive genotype is observed.
"""
from __future__ import annotations
import argparse,hashlib,json,math,random,re
from collections import Counter,defaultdict
from pathlib import Path
from statistics import mean,median

from audit_CzechVeg_original_nonrail_analog_reference_pool_v190 import (
    DOI,get_meta,matched_rail_source
)
from audit_CzechVeg_abiotic_regional_reference_opportunity_v194 import abiotic_rows
from analyze_CzechVeg_abiotic_matched_reference_taxonomic_v195 import (
    source_arms
)
from analyze_CzechVeg_original_independent_reference_taxonomic_v191 import (
    source_species,jaccard
)
from audit_SI_CZ_original_continuous_trait_coverage_v220 import read_leda,binomial
from analyze_Slovenia_original_LEDA_height_seed_breadth_v224 import (
    components,quart
)

LEDA_SHA={
 "height_m":"cbf5e761390d78a843512dfc1139acfc1154d46d6c1995a827354840a9f78a92",
 "seed_mass_mg":"34aa12c333fc1ba6de2e560a495ac5d693aa0722a09ecf68418ad94ac1c19cef",
}
MIN_REF=5
MIN_KNOWN_TAXA=5
MIN_TRAIT_COVERAGE=0.55
N_DRAWS=399
N_LOO_DRAWS=79
SEED=20261117+238
PANELS=("all_source_rail_text","original_richness_matched",
        "strict_original_HABITAT_rail_only","actual_survey2010plus")
METRICS=("cwm_distance","taxon_Jaccard","mean_alpha_trace")
# Exclude source morphospecies/groups which are NOT the same as a
# single LEDA species median; avoid any silent "agg." -> exact species.
GROUP_PATTERN=re.compile(r"(?i)(?:\bagg\.?\b|\bgroup\b|\bspp?\.?\b|\bcf\.?\b|\baff\.?\b|\bs\.?\s*l\.?\b)")

def exact_source_name(name):
    if not isinstance(name,str) or GROUP_PATTERN.search(name):return None
    return binomial(name)

def trait_db(height,seed,names):
    # The same fixed source trait axes standardize ALL arms and panels.
    matched={}
    for raw in names:
        tax=exact_source_name(raw)
        if tax and tax in height and tax in seed:
            matched[raw]=(math.log(height[tax]),math.log(seed[tax]))
    if len(matched)<100:
        raise ValueError("Original CzechVeg study species lack enough exact LEDA overlap")
    mu=[mean(v[i] for v in matched.values()) for i in (0,1)]
    sd=[math.sqrt(mean((v[i]-mu[i])**2 for v in matched.values())) for i in (0,1)]
    if min(sd)<.01:raise ValueError("Source species 2D LEDA standardization invalid")
    return {t:tuple((v[i]-mu[i])/sd[i] for i in (0,1)) for t,v in matched.items()},{
        "n_unique_unambiguous_CzechVeg_taxon_names_both_LEDA_traits":len(matched),
        "reference_zlog_mean_unweighted_unique_source_species":mu,
        "reference_zlog_sd_unweighted_unique_source_species":sd,
        "excluded_agg_group_cf_aff_source_names_not_substituted":True}

def profile_plots(plants,values):
    profiles={};audit=Counter()
    for id,src in plants.items():
        audit["source_vegetation_plots_with_vascular_records"]+=1
        valid=src&set(values)
        if len(valid)<MIN_KNOWN_TAXA:
            audit["rejected_fewer_than_5_exact_joint_LEDA_species"]+=1
            continue
        cover=len(valid)/len(src)
        if cover<MIN_TRAIT_COVERAGE:
            audit["rejected_below_55pct_true_species_trait_coverage"]+=1
            continue
        p=components(valid,values)
        profiles[id]={
            "n_original_species":len(src),
            "n_both_exact_trait_species":len(valid),
            "source_dual_trait_coverage":cover,
            "CWM_height":p["height_mean"],
            "CWM_seed_mass":p["seed_mean"],
            "joint_trace":p["joint_trace"],
            "exact_species":src
        }
    audit["source_plots_qualified_true_original_continuous_traits"]=len(profiles)
    return profiles,dict(audit)

def candidate_arms(meta,preselected,profiles,panel):
    if panel not in PANELS:raise ValueError("Unknown preregistered ecological source panel")
    allarm=preselected["all"];abiotic=preselected["elev100"]
    selected=defaultdict(list);qa=Counter()
    for rid in sorted(set(allarm)&set(abiotic)):
        if rid not in profiles:
            qa["rail_text_source_plot_no_sufficient_trait_join"]+=1;continue
        study=meta[rid]["study"]
        year=meta[rid]["year"]
        if panel=="actual_survey2010plus" and year<2010:continue
        if panel=="strict_original_HABITAT_rail_only" and meta[rid]["strict"]!="rail":
            continue
        controls={}
        for arm,source in (("altitude",abiotic),("all",allarm)):
            ids=[]
            for ref in source[rid]:
                if ref not in profiles:continue
                if meta[ref]["study"]==study:raise ValueError("Original control publication same as rail author")
                if panel=="actual_survey2010plus" and meta[ref]["year"]<2010:continue
                if panel=="original_richness_matched":
                    ratio=profiles[ref]["n_original_species"]/profiles[rid]["n_original_species"]
                    if not 2/3<=ratio<=1.5:continue
                ids.append(ref)
            controls[arm]=ids
        if min(len(controls["altitude"]),len(controls["all"]))<MIN_REF:
            qa["rail_text_plot_fewer_than5_trait_qualified_independent_reference_plots"]+=1
            continue
        selected[meta[rid]["region"]].append((rid,controls))
        qa["rail_text_plot_has_both_true_control_types"]+=1
    return dict(selected),{
        "panel":panel,"source_rail_text_plots_matched_to_two_types":qa["rail_text_plot_has_both_true_control_types"],
        "source_phyto_regions_supported":len(selected),
        "source_original_rail_text_plots_qualifying_before_characteristic_filter":len(set(allarm)&set(abiotic)),
        "original_phyto_region_source_plots_no_ID":{
             str(reg):len(v) for reg,v in sorted(selected.items())},
        "source_eligible_original_rail_plots_excluded_due_to_insufficient_reference":qa[
             "rail_text_plot_fewer_than5_trait_qualified_independent_reference_plots"],
        "all_original_reference_plots_preselected_before_species":True,
        "original_true_railway_habitat_NOT_confirmed":True}

def pairwise_mean(v,kind):
    if len(v)<3:raise ValueError("Minimum three source independent phytoregions")
    d=[]
    for i in range(len(v)):
        for j in range(i+1,len(v)):
            if kind=="CWM":
                a,b=v[i],v[j]
                d.append(math.hypot(a["CWM_height"]-b["CWM_height"],
                                    a["CWM_seed_mass"]-b["CWM_seed_mass"]))
            elif kind=="species":
                d.append(jaccard(v[i]["exact_species"],v[j]["exact_species"]))
            else:raise ValueError("Unknown pairwise metric")
    return mean(d)

def one_draw(frame,prof,regions,rng):
    selections={"rail":[],"altitude":[],"all":[]}
    for reg in regions:
        rail,control=rng.choice(frame[reg])
        ids={"rail":rail,"altitude":rng.choice(control["altitude"]),
             "all":rng.choice(control["all"])}
        for k,id in ids.items():selections[k].append(prof[id])
    out={}
    for label,groups in selections.items():
        out[f"{label}_cwm_distance"]=pairwise_mean(groups,"CWM")
        out[f"{label}_taxon_Jaccard"]=pairwise_mean(groups,"species")
        out[f"{label}_mean_alpha_trace"]=mean(p["joint_trace"] for p in groups)
        out[f"{label}_source_species_n"]=mean(p["n_original_species"] for p in groups)
        out[f"{label}_actual_trait_species_n"]=mean(p["n_both_exact_trait_species"] for p in groups)
        out[f"{label}_actual_trait_coverage"]=mean(p["source_dual_trait_coverage"] for p in groups)
    for control in ("altitude","all"):
        for typ in METRICS:
            out[f"rail_minus_{control}_{typ}"]=out[f"rail_{typ}"]-out[f"{control}_{typ}"]
    return out

def summarize(v):
    if not v:raise ValueError("No source technical draws")
    return {"mean":mean(v),"q025_technical_draws_NOT_field_CI":quart(v,.025),
      "q975_technical_draws_NOT_field_CI":quart(v,.975),
      "n_source_technical_draws":len(v),
      "negative_draw_fraction":sum(x<0 for x in v)/len(v)}

def eval_panel(frame,profiles,seed,draws=N_DRAWS):
    regions=sorted(frame)
    if len(regions)<3:
        return {"identified":False,"n_phyto_regions":len(regions),
                "reason":"Fewer than three actual original source regions after true trait and control gates",
                "underpowered_rail_vs_regional_control_not_zero_effect":True}
    rng=random.Random(seed)
    samples=[one_draw(frame,profiles,regions,rng) for _ in range(draws)]
    keys=sorted(samples[0])
    summarized={k:summarize([s[k] for s in samples]) for k in keys}
    important=tuple(k for k in keys if k.startswith("rail_minus_"))
    loo={}
    if len(regions)>=4:
        for i,drop in enumerate(regions):
            remain=[x for x in regions if x!=drop]
            rr=random.Random(seed+10000+i*113)
            res=[one_draw(frame,profiles,remain,rr) for _ in range(N_LOO_DRAWS)]
            loo[str(drop)]={k:mean(row[k] for row in res) for k in important}
    return {"identified":True,
        "source_geographical_regions":len(regions),
        "n_region_pairs_per_draw_DEPENDENT":len(regions)*(len(regions)-1)//2,
        "number_source_rail_plots":sum(len(x) for x in frame.values()),
        "source_technical_resampling_not_independent_field_replicates":draws,
        "source_results":summarized,
        "leave_one_original_phyto_region_out":loo,
        "loo_signs_for_source_rail_minus_control":{
            k:{"n_loo":len(loo),"positive":sum(p[k]>0 for p in loo.values()),
               "negative":sum(p[k]<0 for p in loo.values())}
            for k in important},
        "site_and_region_reuse_not_a_study_level_field_CI":True}

def run(headers,species,height,seed,draws=N_DRAWS):
    source_meta,original_header=get_meta(headers)
    matching_blocks,blockqa=matched_rail_source(source_meta)
    abiotic,abioticqa=abiotic_rows(headers,source_meta)
    arms,read_ids,armqa=source_arms(source_meta,matching_blocks,abiotic)
    original_plants,srcqa=source_species(species,read_ids)
    allnames=set().union(*original_plants.values())
    x,zqa=trait_db(height,seed,allnames)
    profiles,profileqa=profile_plots(original_plants,x)
    panels={}
    for i,panel in enumerate(PANELS):
        frame,q=candidate_arms(source_meta,arms,profiles,panel)
        panels[panel]={
            "original_source_plot_support":q,
            "statistical_functional_result":eval_panel(frame,profiles,SEED+137*i,draws=draws)}
    return {
      "version":"v23.8",
      "status":"SOURCE_VERIFIED_CZECHVEG_RAIL_TEXT_VS_INDEPENDENT_NONTRANSPORT_FUNCTIONAL_REFS",
      "authentic_data_sources":{
        "original_CzechVeg_DOI":DOI,
        "CzechVeg_2026_v1_1_plot_headers_md5":original_header["publisher_original_headers_md5"],
        "CzechVeg_2026_v1_1_full_species_csv_md5":srcqa["publisher_species_csv_md5"],
        "CzechVeg_original_header_n":117739,
        "CzechVeg_original_plant_species_records_n":srcqa["all_original_species_rows"],
        "source_full_selected_study_reference_plot_IDs_n":len(read_ids),
        "original_selected_original_vascular_plant_plot_species_counts":srcqa,
        "source_Pladias_CzechVeg_abiotic_original_header_QA":abioticqa,
        "source_dual_LEDA_trait_species_QA":zqa,
        "source_real_plant_2D_plot_profile_QA":profileqa,
        "all_reference_candidates_original_metadata_preselected_before_plant_assemblages":armqa},
      "panels":panels,
      "interpretation_limitations":[
        "CzechVeg railway category is HABITAT or LOCALITY text (not physically confirmed operating ballast), so it is not direct controlled railway field treatment",
        "Absence of transportation terms is not proof that reference plots are natural/nonrail; these are unlabeled nontransport-text original vegetation plots",
        "Reference plots are from a different original study than each rail source plot, but within the same compiled country dataset",
        "Plot geographic regions, decade, sampled-area band matched by original metadata; elevation within ±100m, not exact rock substrate, management, climate or human density",
        "CzechVeg RAIL-text plots have far fewer plots than all nontransport vegetation; original sampling is nonrandom and survey study contexts differ",
        "Czech station-wide Preslia 39 study and CzechVeg plot-site grain are NOT identical survey units or paired physical railway station vs nonrail measurements",
        "LEDA canopy height and seed mass are species medians, not individual phenotypes or genotype evolution",
        "Resampling draws and phyto-region pairs are not independent biological studies; output quantiles are not formal 95% confidence intervals",
        "SYNTAXON plant-composition-derived analogous classes are deliberately NOT used to select main reference groups",
        "Reported functional distances and alpha breadth different estimands; mean annual fraction from Slovenia source not directly pooled",
        "No actual freight train OD and 2022/2023 same-year train traffic linked to this CzechVeg plot context",
        "No causal universal railroad homogenization/divergence, plant dispersal or selection inference can be made"]}

def main():
    p=argparse.ArgumentParser()
    for k in ("headers","species","height","seed","output"):
        p.add_argument("--"+k,required=True)
    p.add_argument("--draws",type=int,default=N_DRAWS)
    a=p.parse_args()
    for file,dig in ((a.height,LEDA_SHA["height_m"]),(a.seed,LEDA_SHA["seed_mass_mg"])):
        if hashlib.sha256(Path(file).read_bytes()).hexdigest()!=dig:
            raise ValueError("Original University Oldenburg LEDA table source checksum changed")
    h,hq=read_leda(a.height,"height_m")
    s,sq=read_leda(a.seed,"seed_mass_mg")
    ans=run(a.headers,a.species,h,s,draws=a.draws)
    target=Path(a.output);target.parent.mkdir(parents=True,exist_ok=True)
    target.write_text(json.dumps(ans,ensure_ascii=False,indent=2),encoding="utf8")
    quick={
      "status":ans["status"],
      "authentic_data_sources":{k:v for k,v in ans["authentic_data_sources"].items()
                               if k not in ("source_Pladias_CzechVeg_abiotic_original_header_QA",)},
      "panels":{panel:{
          "source_original_rail_plots":p["original_source_plot_support"]["source_rail_text_plots_matched_to_two_types"],
          "regions":p["original_source_plot_support"]["source_phyto_regions_supported"],
          "identified":p["statistical_functional_result"]["identified"],
          "reason":p["statistical_functional_result"].get("reason"),
          "results":{k:v for k,v in p["statistical_functional_result"].get("source_results",{}).items()
                     if k.startswith("rail_minus_") or k in
                      ("rail_cwm_distance","altitude_cwm_distance","all_cwm_distance",
                       "rail_mean_alpha_trace","altitude_mean_alpha_trace","all_mean_alpha_trace",
                       "rail_actual_trait_species_n","altitude_actual_trait_species_n","all_actual_trait_species_n")},
          "LOO_signs":p["statistical_functional_result"].get("loo_signs_for_source_rail_minus_control",{}),
      } for panel,p in ans["panels"].items()}
    }
    print("V238_CZECHVEG_REAL_FUNCTIONAL_REGIONAL_REFERENCE="+json.dumps(quick,ensure_ascii=False,separators=(",",":")),flush=True)
if __name__=="__main__":main()
