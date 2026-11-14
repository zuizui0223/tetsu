#!/usr/bin/env python3
"""v12.8: ORIGINAL published Austrian 60-site × 2014/2015 field sample
count stratified comparison with strict real-source outcome labeling.

SOURCE: Pascher K, Hainz-Renetzeder C & Schneeweiss GM (2016).
Risk of seed spillage of imported oilseed rape along transport routes.
Austrian Federal Ministry of Health report, 107 pp., Table 14 pp. 53-54.
Official original report:
https://www.verbrauchergesundheit.gv.at/dam/jcr:e7da8178-d364-4ca4-8e1d-5fa87f4cd1dc/rapsstudie_risk_of_speed_spillage_201707071.pdf
Original 2017 peer-reviewed PERSPECTIVE summarizing source original field study:
Pascher et al. (2017), DOI 10.3389/fevo.2017.00104.

60 original field sites independently transcribed from Table 14 (ONLY
site code and published NUMBERS OF SAMPLED FERAl PLANTS, NOT numbers of
seeds), and verified to match the paper's authentic categorical sample
sizes BHR 10, BHK 10, GBH 6, FB 2, DH 6, SR 11, SK 11, OEM 4; original
1001 sampled feral plants in 2014, 1112 in 2015, 2113 total, and
44 out of 60 source sites with feral sampled OSR in at least one year.

Crucial scientific distinctions:
- Original source collected PLANT LEAVES for DNA; 2014/2015 numbers
  are sampled established feral Brassica napus PLANTS. They are NOT
  directly trapped seeds, viable-cargo propagules, germination rate,
  total living plant abundance, establishment probability conditional
  on arriving seed, or railway cargo flow.
- Station vs roadside source survey units/length/observer constraints
  differ; within/without oilseed rape cultivating areas is an observed
  source stratum, NOT randomized confounder control.
- Source rail *stations* BHR/BHK vs *road segments* SR/SK; target
  rail sections and motorways were excluded by original authors
  for access/safety. Cargo traffic is not measured.
- No pooling original targeted hotspots (border station/switchyard/
  port/oil mill) with randomly selected ordinary rail stations/roads
  for a freight causal claim.
- Site positive if and only if at least one source feral plant SAMPLE
  was collected in either year. A sample count of zero is not a
  known zero entire vegetation community nor proof of absent seeds.
- Exact tests are conditional description under a hypothetical
  exchangeability assumption that original source design DOES NOT
  establish. Report exploratory diagnostics, not treatment causal p.

Original table is manually transcribed and total-verified. It is NOT
claimed as a byte-identical original 107-page PDF or a government CSV.
"""
from __future__ import annotations
import argparse,csv,json,math
from collections import Counter,defaultdict
from pathlib import Path
from statistics import mean
from math import comb

SOURCE_PDF="https://www.verbrauchergesundheit.gv.at/dam/jcr:e7da8178-d364-4ca4-8e1d-5fa87f4cd1dc/rapsstudie_risk_of_speed_spillage_201707071.pdf"
PAPER_DOI="10.3389/fevo.2017.00104"
REPORTED_SITE_CATEGORIES={"BHR":10,"BHK":10,"GBH":6,"FB":2,
                          "DH":6,"SR":11,"SK":11,"OEM":4}
ACTUAL_BORDER_CODES=("GBH1","GBH2","GBH3","GBH4","GBH5","GBH7")
ORIGINAL_YEAR_TOTALS=(1001,1112)
ACTUAL_SOURCE_SITES_WITH_SAMPLED_FERAL_OSR=44


def canonical_source_codes():
    out=set()
    for group,n in REPORTED_SITE_CATEGORIES.items():
        if group=="GBH":out.update(ACTUAL_BORDER_CODES)
        elif group in ("SR","SK"):
            out.update(group+str(i).zfill(2) for i in range(1,n+1))
        else:out.update(group+str(i) for i in range(1,n+1))
    return out


def source_category(code):
    if code.startswith(("BHR","BHK","GBH","OEM")):
        return code[:3]
    if code.startswith(("FB","DH","SR","SK")):
        return code[:2]
    raise ValueError("Not an original author-coded Austrian 2014-2015 field site")


def validate_original_table(transcribed_csv):
    with open(transcribed_csv,encoding="utf-8-sig",newline="") as file:
        reader=csv.DictReader(file)
        if reader.fieldnames!=["source_site_id","sampled_OSR_2014","sampled_OSR_2015"]:
            raise ValueError("Original published Table14 transcribed columns changed")
        rows=[]
        ids=set()
        for x in reader:
            code=x["source_site_id"].strip()
            if code in ids or code not in canonical_source_codes():
                raise ValueError("Original 60 field survey site codes missing, duplicate or invented")
            ids.add(code)
            a=int(x["sampled_OSR_2014"])
            b=int(x["sampled_OSR_2015"])
            if a<0 or b<0 or a>10000 or b>10000:
                raise ValueError("Original site sample counts invalid")
            category=source_category(code)
            if code not in canonical_source_codes():
                raise ValueError("Unrecognized source category/site key")
            rows.append({"source_site_code":code,"source_category":category,
                         "sampled_2014":a,"sampled_2015":b,
                         "feral_sampled_2014":a>0,"feral_sampled_2015":b>0,
                         "feral_sampled_2014_or_2015":a+b>0})
    if len(rows)!=60 or ids!=canonical_source_codes():
        raise ValueError("Published source 60 site IDs have not been recovered")
    counts=Counter(x["source_category"] for x in rows)
    if dict(counts)!=REPORTED_SITE_CATEGORIES:
        raise ValueError("Source category sample design differs from original 60-site research")
    original_totals=(sum(x["sampled_2014"] for x in rows),
                     sum(x["sampled_2015"] for x in rows))
    if original_totals!=ORIGINAL_YEAR_TOTALS:
        raise ValueError("Source transcribed original Table14 sample counts do not reproduce published totals")
    if sum(x["feral_sampled_2014_or_2015"] for x in rows)!=ACTUAL_SOURCE_SITES_WITH_SAMPLED_FERAL_OSR:
        raise ValueError("Original 44/60 established-plant positive source sites do not reconcile")
    return rows


def observed_count_summary(rows):
    return {
        "n_original_author_field_sites":len(rows),
        "n_site_sources_with_feral_OSR_sampled_in_either_year":sum(
            x["feral_sampled_2014_or_2015"] for x in rows),
        "n_site_sources_with_feral_OSR_sampled_2014":sum(
            x["feral_sampled_2014"] for x in rows),
        "n_site_sources_with_feral_OSR_sampled_2015":sum(
            x["feral_sampled_2015"] for x in rows),
        "n_original_sites_with_plant_sampling_in_BOTH_years":sum(
            x["feral_sampled_2014"] and x["feral_sampled_2015"] for x in rows),
        "n_original_sites_with_plant_sampled_ONLY_2014":sum(
            x["feral_sampled_2014"] and not x["feral_sampled_2015"] for x in rows),
        "n_original_sites_with_plant_sampled_ONLY_2015":sum(
            not x["feral_sampled_2014"] and x["feral_sampled_2015"] for x in rows),
        "n_original_sites_with_no_sampled_feral_OSR_in_EITHER_year":sum(
            not x["feral_sampled_2014_or_2015"] for x in rows),
        "original_2014_sampled_plant_individuals_NOT_plant_census":sum(
            x["sampled_2014"] for x in rows),
        "original_2015_sampled_plant_individuals_NOT_plant_census":sum(
            x["sampled_2015"] for x in rows),
        "original_author_total_sampled_plant_individuals_NOT_seeds":sum(
            x["sampled_2014"]+x["sampled_2015"] for x in rows),
    }


def category_report(rows):
    out={}
    for code in REPORTED_SITE_CATEGORIES:
        sub=[x for x in rows if x["source_category"]==code]
        d=observed_count_summary(sub)
        d["feral_plant_sampling_site_positive_fraction_either_year"]=(
            d["n_site_sources_with_feral_OSR_sampled_in_either_year"]/len(sub))
        out[code]=d
    return out


def source_presence(rows,category):
    group=[x for x in rows if x["source_category"]==category]
    return sum(x["feral_sampled_2014_or_2015"] for x in group),len(group)


def hypergeom_prob(N,K,n,a):
    if not (max(0,n-(N-K))<=a<=min(n,K)):return 0.0
    return comb(K,a)*comb(N-K,n-a)/comb(N,n)


def fisher_exact_two_sided(a,b,c,d):
    """Pure-standard-library Fisher probability-ordering exact test."""
    if min(a,b,c,d)<0:
        raise ValueError("Nonnegative observed original site counts required")
    N=a+b+c+d;K=a+c;n=a+b
    p_obs=hypergeom_prob(N,K,n,a)
    p=sum(hypergeom_prob(N,K,n,x) for x in range(
        max(0,n-(N-K)),min(n,K)+1)
        if hypergeom_prob(N,K,n,x)<=p_obs+1e-12)
    return min(1.0,p)


def conditional_stratified_exact(specs):
    """Exact convolution under hypothetical within-growing-area rail
    station vs road site exchangeability, NOT real randomized treatment.
    Test statistic is total rail positive sites across source strata.
    """
    dist={0:1.0}
    observed=0
    expectation=0.
    checks=[]
    for rail_positive,rail_n,road_positive,road_n in specs:
        if min(rail_positive,road_positive)<0 or rail_positive>rail_n or road_positive>road_n:
            raise ValueError("Original stratified source positive site count impossible")
        N=rail_n+road_n;K=rail_positive+road_positive
        nxt=defaultdict(float)
        for old,val in dist.items():
            for x in range(max(0,rail_n-(N-K)),min(rail_n,K)+1):
                nxt[old+x]+=val*hypergeom_prob(N,K,rail_n,x)
        dist=dict(nxt)
        observed+=rail_positive
        expectation+=rail_n*K/N
        checks.append({"n_source_rail_station_sites":rail_n,
                       "n_source_road_section_sites":road_n,
                       "n_source_sites_with_feral_plants_sampled_total":K,
                       "rail_source_site_positive":rail_positive,
                       "road_source_site_positive":road_positive})
    p=sum(prob for x,prob in dist.items()
          if abs(x-expectation)+1e-12>=abs(observed-expectation))
    if abs(sum(dist.values())-1)>1e-10:
        raise ValueError("Exact conditional null probabilities failed normalization")
    return {
        "original_category_source_strata":checks,
        "source_station_successes_observed":observed,
        "source_station_successes_exchangeability_expected":expectation,
        "exact_two_sided_centered_stratified_tail_probability":p,
        "exact_one_sided_station_higher_tail_probability":sum(
            prob for x,prob in dist.items() if x>=observed),
        "hypothetical_null":"exchange original railway STATION versus ROAD SEGMENT labels WITHIN original OSR-growing-area source strata",
        "NOT_randomized_railway_seed_carriage_treatment_test":True,
    }


def risk_and_source_sensitivity(rows):
    combos=(("within_OSR_cultivation_area","BHR","SR"),
            ("outside_OSR_cultivation_area","BHK","SK"))
    stats={}
    source_specs=[]
    for name,rail,road in combos:
        a,ra=source_presence(rows,rail)
        c,rd=source_presence(rows,road)
        source_specs.append((a,ra,c,rd))
        stats[name]={
            "railway_station_feral_sampled_sites":a,
            "n_source_railway_station_sites":ra,
            "road_segment_feral_sampled_sites":c,
            "n_source_road_segments":rd,
            "railway_station_site_detection_fraction":a/ra,
            "road_site_detection_fraction":c/rd,
            "observational_railway_station_minus_road_site_detection_fraction":a/ra-c/rd,
            "exploratory_two_sided_Fisher_exact_conditional_site_label_p":
                fisher_exact_two_sided(a,ra-a,c,rd-c),
            "source_site_type_different_station_vs_road_segment_area_and_safety_access":True,
        }
    a=sum(s[0] for s in source_specs)
    ra=sum(s[1] for s in source_specs)
    c=sum(s[2] for s in source_specs)
    rd=sum(s[3] for s in source_specs)
    source_rd=(a/ra-c/rd)
    bal_rd=mean((s[0]/s[1]-s[2]/s[3]) for s in source_specs)
    leaveout=[]
    for source in sorted(rows,key=lambda z:z["source_site_code"]):
        if source["source_category"] not in ("BHR","SR","BHK","SK"):continue
        reduced=[r for r in rows if r["source_site_code"]!=source["source_site_code"]]
        p=[]
        for _,rail,road in combos:
            rr,nr=source_presence(reduced,rail);dr,nd=source_presence(reduced,road)
            p.append(rr/nr-dr/nd)
        leaveout.append({"excluded_source_site_id":source["source_site_code"],
                         "stratum_balanced_rail_minus_road_detection":mean(p)})
    return {
      "observed_original_source_paired_sampling_frame":"10 source railway stations + 11 road sections in growing areas, 10 source railway stations + 11 road sections outside growing areas",
      "original_source_station_vs_road_comparison_by_cultivation_stratum":stats,
      "n_original_non_hotspot_station_source_sites":ra,
      "n_original_non_hotspot_road_source_sections":rd,
      "n_source_station_sites_with_OSR_samples":a,
      "n_source_road_sections_with_OSR_samples":c,
      "pooled_original_station_sampled_fraction":a/ra,
      "pooled_original_road_sampled_fraction":c/rd,
      "original_station_minus_road_observed_detection_fraction":source_rd,
      "source_cultivation_stratum_equal_station_minus_road_detection_fraction":bal_rd,
      "exploratory_pooled_Fisher_exact_2sided_p":
          fisher_exact_two_sided(a,ra-a,c,rd-c),
      "source_stratified_exact_probability_NOT_field_randomization":
          conditional_stratified_exact(source_specs),
      "original_42_site_leave_one_out_n":len(leaveout),
      "original_42_site_leave_one_out_min_balanced_difference":min(x[
          "stratum_balanced_rail_minus_road_detection"] for x in leaveout),
      "original_42_site_leave_one_out_max_balanced_difference":max(x[
          "stratum_balanced_rail_minus_road_detection"] for x in leaveout),
      "original_42_site_leave_one_out_opposite_sign_vs_original":
          sum(x["stratum_balanced_rail_minus_road_detection"]*bal_rd<0 for x in leaveout),
      "source_42_sites_rail_and_road_categories_ONLY_other_hotspots_excluded":True,
    }


def run(csv_path):
    rows=validate_original_table(csv_path)
    counts=observed_count_summary(rows)
    by=category_report(rows)
    stat=risk_and_source_sensitivity(rows)
    return {
        "complete":True,"version":"v12.8",
        "original_source_publication":"Pascher et al. 2016 Austrian Federal Ministry of Health report pp107 TABLE 14 pp53–54",
        "original_author_report_official_URL":SOURCE_PDF,
        "published_peer_reviewed_perspective_DOI":PAPER_DOI,
        "source_table_status":"manual exact numeric transcription of original author Table 14, NOT source-original PDF byte export, independently reconciled 60 site codes and official 1001+1112 samples",
        "original_author_60_site_2yr_sample_validation":counts,
        "original_author_plant_sampling_by_true_source_site_type":by,
        "original_stratified_regular_railway_station_vs_road_segment_source_site_comparison":stat,
        "original_site_count_not_total_species_flora_or_cargo_seeds":True,
        "critical_causal_limitations":[
            "The source counts are SAMPLED feral OSR plants for leaf DNA, not complete population abundance, propagule seed catches or a count of seeds that germinated.",
            "Sampling area/length and railway station versus road SEGMENT field access differ; source sites were observationally sampled, so exact test exchangeability assumption is not established.",
            "Railway line sections and motorways were NOT directly sampled because of safety and permits; source rail sites are STATIONS, not random rail-route kilometer sections.",
            "The 2014 and 2015 samples do not track individual plant survival; two-year repeat detection is not demographic persistence or seed bank survival.",
            "Source cultivating area is only one coarse ecological supply factor; source cargo rates, local farming, neighboring landscape, herbicide and plant search effort are not controlled.",
            "Original source 6 border stations, 2 switchyards, 6 ports and 4 oil-mill/processing facilities are targeted high-risk infrastructure, not exchangeable with ordinary stations and road sections; excluded from stratified station-road contrast.",
            "Austrian observation among Brassica napus cannot be generalized as trait-wide seed dispersal or biodiversity response to railway networks.",
            "The source is a transcribed PUBLISHED numeric table, not verified digit-for-digit original PDF bytes, the original 107-page report is cited for verification.",
        ],
        "integration_with_previous_tetsu_research":{
           "arrival_measured_directly_by_Pascher2016":False,
           "established_feral_plant_sample_detected":True,
           "actual_CzechVeg_and_Slovenia_community_surveys_are_NOT_same_units":True,
           "actual_train_slipstream_mowing_seed_viability_adaptation_not_measured":True,
        }
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--original-table14-transcription",required=True)
    p.add_argument("--output",required=True)
    a=p.parse_args()
    res=run(a.original_table14_transcription)
    out=Path(a.output);out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(res,ensure_ascii=False,indent=2),encoding="utf-8")
    core=res["original_stratified_regular_railway_station_vs_road_segment_source_site_comparison"]
    print(json.dumps({
        "complete":True,"version":"v12.8",
        "source":res["original_author_60_site_2yr_sample_validation"],
        "by_category":{
            k:{"sites":v["n_original_author_field_sites"],
               "positive":v["n_site_sources_with_feral_OSR_sampled_in_either_year"],
               "both_years":v["n_original_sites_with_plant_sampling_in_BOTH_years"],
               "sampled_2014":v["original_2014_sampled_plant_individuals_NOT_plant_census"],
               "sampled_2015":v["original_2015_sampled_plant_individuals_NOT_plant_census"]}
            for k,v in res["original_author_plant_sampling_by_true_source_site_type"].items()},
        "comparison":{
           k:v for k,v in core.items()
           if k not in ("original_source_station_vs_road_comparison_by_cultivation_stratum",)
        },
        "source_strata":core["original_source_station_vs_road_comparison_by_cultivation_stratum"],
        "NO_direct_measured_train_seed_transport":True,
    },ensure_ascii=False),flush=True)


if __name__=="__main__":main()
