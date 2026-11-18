#!/usr/bin/env python3
"""v27.5 Slovenia source original local plant-coupling x ACTUAL 2020 trains.

Original published author vegetation: 89 genuine 2022/23 BS/BT source sites,
173 same-side pairs, 15 national line sections; exact author XLSX SHA-gated.
Original infrastructure manager 2023 station kilometre positions from
operator Annex 2B (PDF SHA-gated). Original Eurostat 2020 2-directional goods
and passenger TRAINS per physical route interval, matched by BOTH official
station chainages and national railway line, never by rail code alone.

Only 27/89 source railway vegetation sites map unambiguously to exactly
9 true physical 2020 train intervals nested within 5 original national
lines. Unmapped sites' train counts are UNKNOWN, not zero.

To measure biology at a matched true site, compare original BS vs BT
Jaccard dissimilarity to a reference matrix null constructed over ALL
original plant sites ON THE SAME national rail line, NOT only the 27
traffic-mapped sites. Each BT plot source richness and each source taxon's
BT incidence frequency on that true line are preserved by checkerboard
swaps; BS plant lists frozen. For site i:
    local_excess_i = E_null_beta(BS_i, BT_i_reassigned) - observed_beta_i
positive means its real same-site rail microhabitats overlap more than
can be explained by original line-specific plant frequencies.
Technical MCMC draws are correlated, NOT randomized ecological p-values.

Aggregate each mapped source site within its actual physical traffic
interval (9 independent traffic labels, sites are not nine independent
plant studies), compare SOURCE segment-mean local_excess to historic
2020 goods & passenger log(1 + trains per year). Summarize ranks and
leave-one-segment/line, then EXHAUSTIVELY reassign goods labels within
actual national lines (48 possible permutations including the original)
while preserving source line, plant groups and passenger values.
Assess a within-line centered freight slope, and a passenger-adjusted
within-line slope if identifiable. Exact finite 48 arrangements are a
technical negative-control source-label falsification, not experiment p.

No actual 2022/23 traffic, origin-destination, wind, mowing or seed
transport data is claimed. The official main G vs regional R category,
SINGLE/DOUBLE layout, source route, original site year, geographic
gradients, COVID 2020 and original field selection are confounded.
No raw taxon-by-site, original latitude/longitude or publisher source
full records are republished.
"""
from __future__ import annotations
import argparse,json,math,random
from collections import Counter,defaultdict
from pathlib import Path
from statistics import mean

from analyze_Slovenia_rail_microhabitats_v103 import extract_original_xlsx
from analyze_rail_local_site_coupling_v257 import slovenia_location_sets
from audit_RAIL_ROAD_source_frequency_fixed_reference_null_v270 import (
    swaps,beta,source_freq
)
from audit_SI_2023_official_station_chainage_2020_eurostat_v231 import (
    extract_authoritative_stations,source_vegetation_sites,
    eurostat_2020_segments,segment_route_station_join,
    collapse_exact_reverse_direction_segments,assign_segment
)
from analyze_SI_2020_traffic_vs_2022_2023_flora_exploratory_v232 import (
    SOURCE_PDF_SHA,segment_traffic_sites,rank
)
from validate_SI_actual_2020_freight_vegetation_beta_holdout_v233 import (
    exact_within_line_labels
)
from analyze_Slovenia_gov_official_main_regional_traits_v225 import official_class

DRAWS=99
BURNIN=3000
SPACING=350
SEED=275118
METRICS=("local_excess_beta","source_same_site_jaccard_beta","source_taxa_BS_richness",
         "source_taxa_BT_richness")
EXPECTED={"publisher_true_rail_locations":89,
          "traffic_original_flora_locations":27,
          "traffic_original_physical_segments":9,
          "traffic_national_rail_lines":5}
TOL=1e-10


def site_specific_null(source:dict,draws:int=DRAWS,
                       burnin:int=BURNIN,spacing:int=SPACING)->tuple[dict,dict]:
    if draws<30:raise ValueError("MCMC technical draws too small")
    by=defaultdict(list)
    for sid,q in source.items():
        if not q["BS"] or not q["BT"]:
            raise ValueError("Original true side-matched BS/BT taxon sets incomplete")
        by[q["section"]].append((str(sid),q))
    if len(by)!=15 or sum(len(v) for v in by.values())!=89:
        raise ValueError("Original source full 15 line /89 site roster drift")
    site={}
    no_null=[]
    section_summary={}
    for line_ix,(line,items) in enumerate(sorted(by.items())):
        items=sorted(items,key=lambda x:x[0])
        if len(items)<2:
            no_null.extend(sid for sid,q in items)
            continue
        bs=[set(q["BS"]) for sid,q in items]
        bt=[set(q["BT"]) for sid,q in items]
        hold_size=[len(s) for s in bt]
        hold_freq=source_freq(bt)
        observed=[beta(a,b)["beta_jac"] for a,b in zip(bs,bt)]
        matrix=[set(x) for x in bt]
        rng=random.Random(SEED+1799*line_ix)
        accepted=swaps(matrix,rng,burnin)
        accum=[0.]*len(items)
        for _ in range(draws):
            accepted+=swaps(matrix,rng,spacing)
            for j in range(len(items)):
                accum[j]+=beta(bs[j],matrix[j])["beta_jac"]
        if [len(s) for s in matrix]!=hold_size or source_freq(matrix)!=hold_freq:
            raise ValueError("Source local site richness or original section taxon incidence changed")
        excess=[]
        for j,(sid,q) in enumerate(items):
            null_beta=accum[j]/draws
            observed_beta=observed[j]
            site[sid]={
                "section":line,
                "local_excess_beta":null_beta-observed_beta,
                "source_same_site_jaccard_beta":observed_beta,
                "source_taxa_BS_richness":len(bs[j]),
                "source_taxa_BT_richness":len(bt[j])
            }
            excess.append(null_beta-observed_beta)
        section_summary[line]={
            "original_independent_sites":len(items),
            "mean_source_local_coupling_excess":mean(excess),
            "source_accepted_checkerboard_swaps":accepted,
            "source_original_BT_taxa_incidence_frequencies_fixed_exact":True
        }
    if len(site)!=88 or len(no_null)!=1 or len(section_summary)!=14:
        raise ValueError("Strict v27.3 original 14 sections 88 sites + 1 no-support section changed")
    if abs(mean(v["source_same_site_jaccard_beta"] for v in
                [x for x in site.values() if x["section"] in section_summary])-.6367636451153651)>0.20:
        # Site-equal is intentionally DIFFERENT from source section-equal.
        raise ValueError("Actual original source local beta unexpectedly outside prior source range")
    return site,{
        "full_original_source_site_roster":sum(len(items) for items in by.values()),
        "original_sections_with_two_plus_field_sites":len(section_summary),
        "source_sites_with_valid_rail_section_null":len(site),
        "original_section_single_site_no_frequency_null_count":len(no_null),
        "no_unmatched_source_sites_given_fake_zero_or_null":True,
        "original_line_source_support":section_summary,
        "technical_MCMC_draws_are_not_independent_field_experiments":True
    }


def mapped_segment_data(original:dict,pdf:Path)->tuple[list[dict],dict]:
    stations,source_station_qa=extract_authoritative_stations(pdf)
    if source_station_qa["sha256_original_pdf_source"]!=SOURCE_PDF_SHA:
        raise ValueError("Original official 2023 network statement PDF SHA changed")
    plants,floraqa=source_vegetation_sites(original)
    routes,trafqa=eurostat_2020_segments()
    joined,euroqa=segment_route_station_join(stations,routes)
    complete,physicalqa=collapse_exact_reverse_direction_segments(joined)
    mapped,crossqa=assign_segment(plants,complete)
    actual,guild=segment_traffic_sites(original,mapped)
    qc={
        "official_2023_station_PDF_SHA256":source_station_qa["sha256_original_pdf_source"],
        "original_eurostat_2020_JSON_SHA256":trafqa["source"]["sha256_original_json_response"],
        "source_true_sites":len(plants),
        "source_2020_mapped_sites":len(actual),
        "source_2020_mapped_undirected_physical_segments":len({x["source_segment"] for x in actual}),
        "source_true_national_lines":len({x["section"] for x in actual}),
        "source_genuine_physical_reverse_direction_intervals":len(complete),
        "unmapped_real_plant_sites":len(plants)-len(actual),
        "source_mapped_line_site_counts":dict(Counter(x["section"] for x in actual)),
        "source_2020_covid_year_not_field_2022_2023":True,
    }
    if {k:qc[k] for k in (
        "source_2020_mapped_sites","source_2020_mapped_undirected_physical_segments",
        "source_true_national_lines") } != {
            "source_2020_mapped_sites":27,
            "source_2020_mapped_undirected_physical_segments":9,
            "source_true_national_lines":5
        }:
        raise ValueError("2020 source original 27site 9segment 5national-line crosswalk drift")
    return actual,qc


def join_and_group(traffic:list[dict],sites:dict)->tuple[list[dict],dict]:
    grouped=defaultdict(list)
    unmatched=[]
    for row in traffic:
        sid=str(row["sid"])
        site=sites.get(sid)
        if site is None:
            unmatched.append(sid)
            continue
        if site["section"]!=row["section"]:
            raise ValueError("Rail source original location and government section mismatch")
        if not isinstance(row["freight"],(int,float)) or row["freight"]<0:
            raise ValueError("Not genuine historical freight count")
        if not isinstance(row["passenger"],(int,float)) or row["passenger"]<0:
            raise ValueError("Not genuine historical passenger count")
        grouped[row["source_segment"]].append({
            "line":row["section"],"freight":row["freight"],
            "passenger":row["passenger"],"year":row["year"],
            "site_beta_excess":site["local_excess_beta"],
            "site_beta":site["source_same_site_jaccard_beta"]
        })
    segment=[]
    for segment_key,rows in sorted(grouped.items()):
        if len({x["line"] for x in rows})!=1:
            raise ValueError("One undirected physical Eurostat segment cross national lines")
        if len({(x["freight"],x["passenger"]) for x in rows})!=1:
            raise ValueError("One source segment has contradictory original historical traffic")
        segment.append({
            "original_line":rows[0]["line"],
            "class":official_class(rows[0]["line"]),
            "source_plant_site_count":len(rows),
            "goods_2020":rows[0]["freight"],
            "passenger_2020":rows[0]["passenger"],
            "mean_local_excess":mean(x["site_beta_excess"] for x in rows),
            "mean_observed_Jaccard_beta":mean(x["site_beta"] for x in rows),
            "source_survey_years":dict(Counter(str(x["year"]) for x in rows))
        })
    if len(segment)<7:
        raise ValueError("Too few 2020 physical segments after source frequency null join")
    return segment,{
        "original_27_2020_mapped_sites":len(traffic),
        "original_mapped_sites_with_valid_frequency_fixed_null":sum(len(x) for x in grouped.values()),
        "matched_original_physical_segments":len(segment),
        "matched_national_lines":len({x["original_line"] for x in segment}),
        "mapped_sites_on_original_singleton_line_without_local_reference_null":len(unmatched),
        "historical_traffic_missing_not_zero_imputed":True
    }


def centred_slope(segment:list[dict],goods_override:dict[int,float]|None=None,
                  passenger_adjusted:bool=False)->dict:
    """Only within original national line, each true physical segment one vote."""
    n=len(segment)
    goods=[(math.log1p(row["goods_2020"]) if goods_override is None else
            goods_override[i]) for i,row in enumerate(segment)]
    pax=[math.log1p(x["passenger_2020"]) for x in segment]
    y=[x["mean_local_excess"] for x in segment]
    grp=defaultdict(list)
    for i,r in enumerate(segment):grp[r["original_line"]].append(i)
    def center(x):
        v=list(x)
        for ids in grp.values():
            a=mean(x[i] for i in ids)
            for i in ids:v[i]=x[i]-a
        return v
    xc=center(goods);yc=center(y);pc=center(pax)
    xx=sum(x*x for x in xc)
    if xx<1e-12:return {"identified":False,"source_within_line_goods_variance_too_small":True}
    if passenger_adjusted:
        pp=sum(p*p for p in pc)
        if pp<1e-12:return {"identified":False,"source_within_line_passenger_variance_too_small":True}
        xp=sum(x*p for x,p in zip(xc,pc))/pp
        yp=sum(v*p for v,p in zip(yc,pc))/pp
        xres=[x-xp*p for x,p in zip(xc,pc)]
        yres=[v-yp*p for v,p in zip(yc,pc)]
    else:
        xres=xc;yres=yc
    residual_sxx=sum(x*x for x in xres)
    if residual_sxx<1e-10:
        return {"identified":False,"source_collinear_passenger_goods_within_line":True}
    return {"identified":True,
            "within_line_goods_log_count_local_excess_slope":
                sum(x*z for x,z in zip(xres,yres))/residual_sxx,
            "original_source_segment_n":n,
            "original_national_lines":len(grp),
            "lines_with_2plus_physical_2020_segments":sum(len(v)>1 for v in grp.values()),
            "source_2020_freight_intensity_not_network_origin_destination":True}


def traffic_diagnostic(segment:list[dict])->dict:
    if len(segment)!=9 or len({x["original_line"] for x in segment})!=5:
        raise ValueError("Only original 9 source physical segments in five actual national lines supported")
    goods=[math.log1p(x["goods_2020"]) for x in segment]
    pax=[math.log1p(x["passenger_2020"]) for x in segment]
    y=[x["mean_local_excess"] for x in segment]
    meta={i:{"national_line":s["original_line"],"freight_2020":s["goods_2020"]}
          for i,s in enumerate(segment)}
    base=centred_slope(segment,passenger_adjusted=False)
    adjusted=centred_slope(segment,passenger_adjusted=True)
    counts=[]
    remapped_results={False:[],True:[]}
    for n,swap in exact_within_line_labels(meta):
        if swap is None:
            counts.append(n);continue
        for adjusted_flag in (False,True):
            v=centred_slope(segment,swap,passenger_adjusted=adjusted_flag)
            if v["identified"]:
                remapped_results[adjusted_flag].append(v["within_line_goods_log_count_local_excess_slope"])
    if counts!=[48]:
        raise ValueError("Original v23.3 source within-national-line only 48 freight assignments not reproduced")
    for mode,ref in ((False,base),(True,adjusted)):
        if ref["identified"] and len(remapped_results[mode])!=48:
            raise ValueError("Some true physical segment freight label assignments not evaluable")
    def permutation_summary(obs,data):
        if not obs["identified"]:
            return {"identified":False,"source_permutation_attempts":48,
                    "reason":"Within national line baseline or passenger-adjusted rank deficient"}
        effect=obs["within_line_goods_log_count_local_excess_slope"]
        return {
            **obs,
            "source_exact_2020_freight_permutations_within_original_national_lines":48,
            "n_same_or_more_extreme_absolute_source_slopes":sum(abs(a)>=abs(effect)-1e-12 for a in data),
            "original_observed_source_slope_rank_fraction_two_sided_not_field_p":
                sum(abs(a)>=abs(effect)-1e-12 for a in data)/48,
            "source_reassignment_min_max":[min(data),max(data)],
            "source_traffic_reassignment_NOT_randomized_field_treatment":True
        }
    loo_segment=[]
    loo_line=[]
    for i in range(9):
        x=[item for j,item in enumerate(segment) if j!=i]
        rho=rank([math.log1p(q["goods_2020"]) for q in x],
                 [q["mean_local_excess"] for q in x])
        if rho is not None:loo_segment.append(rho)
    for line in sorted({x["original_line"] for x in segment}):
        x=[item for item in segment if item["original_line"]!=line]
        rho=rank([math.log1p(q["goods_2020"]) for q in x],
                 [q["mean_local_excess"] for q in x])
        if rho is not None:loo_line.append(rho)
    counts_gr=dict(Counter(x["class"] for x in segment))
    return {
        "source_2020_historical_goods_trains_year_9segment_spearman_local_excess":
            rank(goods,y),
        "source_2020_historical_passenger_trains_year_9segment_spearman_local_excess":
            rank(pax,y),
        "source_9_true_segments":len(segment),
        "source_5_original_national_lines":5,
        "source_main_G_vs_regional_R_segment_counts":counts_gr,
        "historical_2020_goods_both_directions_minmax":[min(x["goods_2020"] for x in segment),
                                                          max(x["goods_2020"] for x in segment)],
        "leave_one_true_traffic_segment_out_spearman_goods_minmax":
            [min(loo_segment),max(loo_segment)] if loo_segment else None,
        "leave_one_true_national_line_out_spearman_goods_minmax":
            [min(loo_line),max(loo_line)] if loo_line else None,
        "within_line_goods_slope_with_source_exact_48_label_redistributions":
            permutation_summary(base,remapped_results[False]),
        "within_line_goods_adjusted_for_passenger_slope_with_exact_48":
            permutation_summary(adjusted,remapped_results[True]),
        "no_source_site_or_segment_taxon_matrix_republished":True
    }


def run(xlsx:Path,audit:Path,pdf:Path,draws:int=DRAWS)->dict:
    original=extract_original_xlsx(xlsx,audit)
    locations,source=slovenia_location_sets(original)
    sites,localqa=site_specific_null(locations,draws=draws)
    actual,trafficqa=mapped_segment_data(original,pdf)
    phys,joinqa=join_and_group(actual,sites)
    if len(phys)!=9:
        raise ValueError("Not all nine genuine source physical segments have local source test support")
    association=traffic_diagnostic(phys)
    return {
        "version":"v27.5",
        "source_author_original_plant_DOI":"10.1007/s44473-026-00155-z",
        "source_historical_actual_train_traffic":"Eurostat official rail_tf_ns20_si, true 2020 physical directions paired",
        "official_2023_station_PDF":"NS_2023_1.2_Annex_2B_Stations.pdf",
        "source_source_site_species_frequency_exact_null_QA":localqa,
        "source_traffic_station_route_geo_and_actual_train_2020_QA":trafficqa,
        "genuine_physical_traffic_segments_and_source_coupling_QA":joinqa,
        "true_original_2020_segment_local_coupling_descriptive_association":association,
        "no_true_2022_2023_same_field_year_traffic":True,
        "nonmapped_original_plant_site_freight_is_unknown_not_zero":True,
        "historical_freight_intensity_NOT_direct_train_origin_destination_or_seed_arrival":True,
        "no_rail_traffic_seed_transport_causal_or_evolved_adaptation_claim":True,
        "no_author_source_plant_or_GPS_site_records_exported":True
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--original-xlsx",type=Path,required=True)
    p.add_argument("--source-audit",type=Path,required=True)
    p.add_argument("--official-station-pdf",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    p.add_argument("--draws",type=int,default=DRAWS)
    a=p.parse_args()
    out=run(a.original_xlsx,a.source_audit,a.official_station_pdf,a.draws)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(out,ensure_ascii=False,indent=2)+"\n",encoding="utf8")
    print(json.dumps({
        "version":"v27.5",
        "original_source_QA":out["genuine_physical_traffic_segments_and_source_coupling_QA"],
        "local_source_support":{k:v for k,v in out["source_source_site_species_frequency_exact_null_QA"].items()
          if k!="original_line_source_support"},
        "true_historical_goods_and_passenger_local_excess_association":
            out["true_original_2020_segment_local_coupling_descriptive_association"]
    },ensure_ascii=False))
if __name__=="__main__":main()
