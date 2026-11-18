#!/usr/bin/env python3
"""v27.4 official G/main vs R/regional rail section local plant-coupling audit.

Authentic Šabić et al. (2026) railway plant XLSX, 89 sites (15 original
national line sections), 173 true BS/BT side pairs. Same original single-site
union of species on left/right sides as source v27.3 (not doubled as samples).
Original government G(main) / R(regional) classification from Slovenian
government plan (2024, page 7) via verified original v22.5 function:
official_class(original source national line). G/R is NOT train traffic.

Within EACH original line, keep all BS source species sets frozen and use
2x2 checkerboard swaps on BT species incidence: exactly preserve each BT
site's plant richness and per-taxon frequency on that national source line.
Outcome local-coupling excess = shuffled-null Jaccard beta - observed
Jaccard beta (positive => more same-site taxa overlap than null).
Aggregate site -> original section -> official G/R (equal section weights).
Track layout SINGLE versus DOUBLE is strongly confounded with G/R;
separate STRICT SAME single-track-only source panels (within each
original section) and never claim main/branch effects from unequal tracks.

This is not independent rail traffic, no causal departure frequency,
network connectivity, freight transfer, train-induced wind, mowing,
seed transport, individual plant fitness or adaptation measured.
99 correlated technical Markov draws are not randomized field p-values.
Only aggregate original line-class effects and counts exported.
"""
from __future__ import annotations
import argparse,json,math
from collections import Counter,defaultdict
from pathlib import Path
from statistics import mean

from analyze_Slovenia_rail_microhabitats_v103 import (
    extract_original_xlsx,validate_original_keys
)
from analyze_rail_local_site_coupling_v257 import slovenia_location_sets
from analyze_Slovenia_gov_official_main_regional_traits_v225 import (
    official_class,KNOWN
)
from audit_RAIL_ROAD_source_frequency_fixed_reference_null_v270 import (
    sampler,finish_null,METRICS
)

SEED=274118
DRAWS=99
BURNIN=2500
SPACING=350
KNOWN_COUNTS={"main_G":(6,34),"regional_R":(9,55)}


def original_grouped(original:dict)->tuple[dict,dict]:
    source,qc=slovenia_location_sets(original)
    loc,zones,sp,tr,cb,qa=validate_original_keys(original)
    if len(source)!=89 or len({v["section"] for v in source.values()})!=15:
        raise ValueError("Genuine publisher 89 sites / 15 line sections drift")
    by=defaultdict(list)
    for lid,q in source.items():
        section=q["section"]
        typ=str(loc[lid]["track_type"])
        if typ not in ("SINGLE","DOUBLE"):
            raise ValueError(f"Original track type not verified: {typ}")
        cls=official_class(section)
        if cls is None:
            raise ValueError(f"Original government G/R section not mapped: {section}")
        by[section].append({"id":lid,"BS":q["BS"],"BT":q["BT"],
                            "layout":typ,"class":cls})
    if set(by)!=set(KNOWN):
        raise ValueError("Source 15 sections not equal true national government G/R list")
    counts={typ:(len({sec for sec,rs in by.items() if rs[0]["class"]==typ}),
                 sum(len(rs) for sec,rs in by.items() if rs[0]["class"]==typ))
            for typ in KNOWN_COUNTS}
    if counts!=KNOWN_COUNTS:
        raise ValueError(f"Original 6/9 lines, 34/55 field-site true G/R gate changed: {counts}")
    qout={
        "original_source_89_field_sites_173_same_side_pairs":True,
        "national_government_2024_classification_2013_legal_source_not_actual_2020_train_count":True,
        "G_R_original_section_field_site_counts":{k:{"sections":v[0],"sites":v[1]} for k,v in counts.items()},
        "original_track_layout_G_R":{
            k:dict(Counter(s["layout"] for rs in by.values() for s in rs if s["class"]==k))
            for k in KNOWN_COUNTS
        },
        "all_sites_original_vegetation_source_track_type_verified":True,
        "unknown_one_source_section_local_null_not_included":True
    }
    return by,qout


def per_source_section(source:dict,min_site:int=2,
                       track_only:str|None=None,seed:int=SEED,
                       draws:int=DRAWS,burnin:int=BURNIN,spacing:int=SPACING)->dict:
    output={}
    for k,(section,sites) in enumerate(sorted(source.items())):
        if track_only is not None:
            sites=[s for s in sites if s["layout"]==track_only]
        if len(sites)<min_site:
            continue
        if len({s["class"] for s in sites})!=1:
            raise ValueError("Original one true national line cannot have G/R cross-labels")
        res=sampler([s["BS"] for s in sites],[s["BT"] for s in sites],
                    draws=draws,burnin=burnin,spacing=spacing,seed=seed+1957*k)
        stat=finish_null(res["observed"],res["null_draws"])
        output[section]={
            "official_class":sites[0]["class"],
            "n_field_locations":len(sites),
            "source_track_layouts":dict(Counter(s["layout"] for s in sites)),
            "local_extra_similarity_nullminusobserved_beta":
                -stat["beta_jac"]["observed_minus_frequency_null"],
            "local_extra_shared_taxa_observedminusnull":
                stat["shared_taxa"]["observed_minus_frequency_null"],
            "observed_jaccard_beta":stat["beta_jac"]["observed"],
            "within_source_section_occupancy_fixed_null_jaccard_beta":
                stat["beta_jac"]["frequency_null_mean"],
            "source_acceptance_swaps":res["accepted_swaps_total"]
        }
    return output


def contrasts(data:dict)->dict:
    groups={tag:[r for r in data.values() if r["official_class"]==tag]
            for tag in KNOWN_COUNTS}
    n={tag:len(vals) for tag,vals in groups.items()}
    support={"n_source_sections":len(data),
             "source_locations":sum(v["n_field_locations"] for v in data.values()),
             "source_main_G_sections":n["main_G"],
             "source_regional_R_sections":n["regional_R"],
             "all_source_section_species_occupancy_reference_margins_fixed":True,
             "source_sections_with_zero_accepted_swaps":sum(
                 v["source_acceptance_swaps"]==0 for v in data.values())}
    if not all(n.values()):
        return {"identified":False,"source_support":support,
                "main_regional_contrast_not_identifiable":True}
    results={}
    for key in ("local_extra_similarity_nullminusobserved_beta",
                "local_extra_shared_taxa_observedminusnull",
                "observed_jaccard_beta"):
        g=[x[key] for x in groups["main_G"]]
        r=[x[key] for x in groups["regional_R"]]
        actual=mean(g)-mean(r)
        leave=[]
        for section,record in data.items():
            keep={k:v for k,v in data.items() if k!=section}
            gg=[x[key] for x in keep.values() if x["official_class"]=="main_G"]
            rr=[x[key] for x in keep.values() if x["official_class"]=="regional_R"]
            if gg and rr:
                leave.append(mean(gg)-mean(rr))
        results[key]={
            "main_G_equal_original_line_mean":mean(g),
            "regional_R_equal_original_line_mean":mean(r),
            "G_minus_R_descriptive":actual,
            "leave_one_original_line_out_range":[min(leave),max(leave)] if leave else None,
            "leave_one_original_line_out_sign_consistent":all(v>0 for v in leave)
               if actual>0 else all(v<0 for v in leave) if actual<0 else False,
            "not_a_field_class_randomization_or_railway_traffic_effect":True,
        }
    return {
        "identified":True,
        "source_support":support,
        "source_section_equal_main_regional_comparisons":results
    }


def run(xlsx:Path,audit:Path,draws:int=DRAWS)->dict:
    original=extract_original_xlsx(xlsx,audit)
    by,qa=original_grouped(original)
    full=per_source_section(by,2,draws=draws)
    if len(full)!=14 or sum(x["n_field_locations"] for x in full.values())!=88:
        raise ValueError("Original v27.3 14 source section 88 site panel changed")
    all_=contrasts(full)
    if abs(mean(v["observed_jaccard_beta"] for v in full.values())-
           .6367636451153651)>1e-10:
        raise ValueError("Prior v27.3 original same-site βjac must reproduce")
    single=per_source_section(by,2,track_only="SINGLE",seed=SEED+100007,draws=draws)
    double=per_source_section(by,2,track_only="DOUBLE",seed=SEED+200019,draws=draws)
    strict_single=contrasts(single)
    strict_double=contrasts(double)
    paired={}
    for sec in sorted(set(single)&set(double)):
        paired[sec]={
            "official_class":single[sec]["official_class"],
            "single_sites":single[sec]["n_field_locations"],
            "double_sites":double[sec]["n_field_locations"],
            "single_minus_double_local_coupling":
                single[sec]["local_extra_similarity_nullminusobserved_beta"]-
                double[sec]["local_extra_similarity_nullminusobserved_beta"],
        }
    return {
        "version":"v27.4",
        "genuine_original_slovenian_author_source_DOI":"10.1007/s44473-026-00155-z",
        "government_official_G_R_designation_URL":"https://www.gov.si/assets/organi-v-sestavi/URSZR/Datoteke/Drzavni-nacrti/dn_zelezniska_azuriran_2024.pdf",
        "original_site_section_layout_source_QA":qa,
        "primary_14_original_rail_sections_88_site_source_fixed_frequency":all_,
        "strict_original_single_track_source_fixed_frequency":strict_single,
        "strict_original_double_track_source_fixed_frequency":strict_double,
        "original_sections_with_both_single_and_double_source_support":{
            "source_n_both_layout_sections_with_2plus_sites_each":len(paired),
            "descriptive_single_minus_double_local_coupling_mean":
                mean(x["single_minus_double_local_coupling"] for x in paired.values()) if paired else None
        },
        "source_section_class_counts_and_local_coupling_internal_no_plant_names":{
            sec:{"official_class":x["official_class"],
                 "n_source_sites":x["n_field_locations"],
                 "original_layout":x["source_track_layouts"],
                 "observed_Jaccard":x["observed_jaccard_beta"],
                 "fixed_margin_excess_similarity":x["local_extra_similarity_nullminusobserved_beta"]}
            for sec,x in full.items()
        },
        "source_MCMC_draws_not_field_experimental_replication":True,
        "no_actual_passenger_or_freight_frequency_in_official_G_R_class":True,
        "no_causal_rail_transport_connectivity_mowing_or_wind_adaptation":True
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--original-xlsx",required=True,type=Path)
    p.add_argument("--source-audit",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    p.add_argument("--draws",type=int,default=DRAWS)
    a=p.parse_args()
    result=run(a.original_xlsx,a.source_audit,a.draws)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf8")
    for key in ("primary_14_original_rail_sections_88_site_source_fixed_frequency",
                "strict_original_single_track_source_fixed_frequency",
                "strict_original_double_track_source_fixed_frequency"):
        x=result[key]
        print(json.dumps({"version":"v27.4","panel":key,
            "support":x["source_support"],
            "identified":x["identified"],
            "main_vs_regional":x.get("source_section_equal_main_regional_comparisons")
        },ensure_ascii=False))
if __name__=="__main__":main()
