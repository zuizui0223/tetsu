#!/usr/bin/env python3
"""v24.5: independent 3-study / 4-country railway station 2D trait beta audit.

True original station plant source cohorts, from THREE separate studies:
 CZ: Kutlvašr et al. (2024), Preslia: 39 stations sampled 2020-22,
     761 current rail-station original species, 5 author rail stretches,
     publisher S1 and corrected S2, original SHA verified.
 IT: Adorni & Ghillani (2020), Ann Museo Rovereto: 28 actual railway
     stations, 438 author taxa; original field 2012-13, 11 plain and
     17 other publisher geography groups, original PDF SHA verified.
 HR/BA: Jasprica et al. (2017), Natura Croatica: 15 railway stations,
     359 author taxa, original TWS/TSS/MED climate classes, original
     33p PDF and SHA verified. HR+BA = ONE independent publication,
     not two studies.

All station cohorts are floristic presence from ENTIRE station/rail
areas, not area-equal local field relevés and not a matched nonrail
reference. Slovenia 89 source railway BALLAST plots remain a separate
grain and are not misreported as additional 89 railway stations.

Continuous traits: identical University Oldenburg LEDA genuine
canopy height (m) and seed mass (mg), two original SHA-pinned source
files. Strict Latin binomial matching; collapse author subspecies to
species within a station rather than treating same species' traits as
different genotypes. No genus/aggregate imputation.
Log transform and z-standardize ONCE across pooled matched species
from ALL THREE authentic independent railway studies, so numerics
are on the SAME 2D trait axes, unlike different study-specific scaling.

For EACH source study separately:
 (1) observed CWM Euclidean distance over original station pairs,
     plus within publisher geography class and across class;
 (2) observed species-equal within-station 2D variance trace;
 (3) observed original matched-species Jaccard (different scale);
 (4) conditional checkerboard fixed-margin null swaps within EACH
     publisher original geography group only. Preserve all original
     trait-known site richness AND each exact matched taxon's
     occupancy frequency within geography group.
 2 independent technical swap chains, 100 draws each with burnin
 18k + thinning 1700 attempts (MCMC diagnostic NOT independent
 country-wide p-values).

No independent nonrail vegetation sampled at stations here, no
contemporaneous train frequency, freight OD, seed capture, wind,
mowing or plasticity/genetic adaptation measured.
No pooling source study station pairs as independent field n.
Primary question is whether v23.6 Czech internal fixed-margin
station functional divergence also appears in independently
published IT and HR/BA, NOT whether all rail globally override
natural regional flora.
"""
from __future__ import annotations
import argparse,hashlib,json,math,random
from collections import Counter,defaultdict
from itertools import combinations
from pathlib import Path
from statistics import mean,median

from analyze_Czech_2024_39station_functional_identity_v143 import read_original_czech_publisher
from extract_Parma_2020_28station_species_matrix_v187 import parse_pdf
from audit_Parma_2020_28rail_station_flora_original_v186 import STATION_CODES
from analyze_Parma_2020_joint_taxonomic_lifeform_null_v188 import PLAIN_11
from extract_Croatia_Bosnia_2017_original_rotated_table_matrix_v183 import rows_from_original
from analyze_Croatia_Bosnia_15station_lifeform_fixedmargin_null_v184 import STATION_MACROCLIMATE
from audit_SI_CZ_original_continuous_trait_coverage_v220 import binomial,read_leda
from analyze_Slovenia_original_LEDA_height_seed_breadth_v224 import quart

LEDASHA={
    "height":"cbf5e761390d78a843512dfc1139acfc1154d46d6c1995a827354840a9f78a92",
    "seed":"34aa12c333fc1ba6de2e560a495ac5d693aa0722a09ecf68418ad94ac1c19cef"}
PUBLISHER_HASH={
    "italy_original_pdf":"c09bd96ff7307a936030195cc4510bf2591576a713de54019625c04aa9333787",
    "hr_ba_original_pdf":"90d489ca8cbc29464d3e6a2b1d8b5d9568c9940e953b1ce720964991366ab956"}
CZ_EXPECT=39
IT_EXPECT=28
HR_EXPECT=15
SOURCE_STUDY_N=3
OBSERVED_N=82
SEED=20261117+245
CHAINS=2
DRAWS=100
BURNIN=18000
GAP=1700
METRICS=("functional_CWM_centroid_distance_all_stationpairs",
         "functional_CWM_centroid_distance_within_geography",
         "functional_CWM_centroid_distance_between_geography",
         "source_station_2D_trait_trace_mean",
         "source_matched_species_Jaccard_all_stationpairs",
         "source_matched_species_Jaccard_within_geography",
         "source_matched_species_Jaccard_between_geography")

def publisher_italy(pdf):
    rows,errors,source=parse_pdf(pdf)
    if len(rows)!=438 or errors or len(STATION_CODES)!=28 or not(
         source["source_pdf_sha256"]==PUBLISHER_HASH["italy_original_pdf"]):
        raise ValueError("Original Parma 2020 publisher PDF/438 taxa source integrity drift")
    stations=[]
    for station in STATION_CODES:
        raw={row['taxon_original'] for row in rows if station in row['station_codes']}
        stations.append({"station":station,"source_names":raw,
          "region_group":"plain" if station in PLAIN_11 else "other"})
    original_counts={station:len({i for i,r in enumerate(rows) if station in r["station_codes"]})
                     for station in STATION_CODES}
    if original_counts.get("SAL")!=157 or original_counts.get("BVT")!=151 or original_counts.get("CHI")!=32:
        raise ValueError("Authentic original Parma station richness crosschecks failed")
    return stations,{"original_n_taxa":len(rows),"study":"IT_Parma_Adorni_Ghillani_2020",
       "original_stations":len(stations),"source_geography_groups":["plain","other"],
       "field_survey_original_years":"2012–2013 (five stations with additional observations)",
       "author_original_pdf_sha256":source["source_pdf_sha256"],
       "original_table_438_rows_station_frequency_discrepancies":source[
           "publisher_frequency_disagreements_count"]}

def publisher_hr_ba(pdf):
    rows,audit=rows_from_original(pdf)
    if len(rows)!=359 or not audit.get("full_matrix_verified") or audit[
        "publisher_original_pdf_sha256"]!=PUBLISHER_HASH["hr_ba_original_pdf"]:
        raise ValueError("Original Jasprica2017 15station rotated table2 not source-verifiable")
    stations=[]
    for id in range(1,16):
        stations.append({"station":f"source_station_{id}",
           "source_names":{row["taxon_original"] for row in rows
                           if id in row["source_station_numbers"]},
           "region_group":STATION_MACROCLIMATE[id]})
    if set(s['region_group'] for s in stations)!={"MED","TWS","TSS"}:
        raise ValueError("Original macroclimate source group assignment changed")
    return stations,{"original_n_taxa":len(rows),
       "study":"HR_BA_Jasprica_2017_SINGLE_STUDY",
       "original_stations":len(stations),"original_2017_publication_2_countries_1_study":True,
       "field_original_survey_years":"2013–2014 (unequal season and subplot area)",
       "original_pdf_sha256":audit["publisher_original_pdf_sha256"],
       "source_Canna_indica_frequency_contradiction_retained":True}

def publisher_czech(s1,s2):
    floras,meta,source=read_original_czech_publisher(s1,s2)
    taxa=set().union(*floras.values())
    if len(floras)!=39 or len(taxa)!=761 or sum(len(x) for x in floras.values())!=6941:
        raise ValueError("Original Preslia 39 railway station 761 taxa source changed")
    st=[{"station":k,"source_names":floras[k],
         "region_group":meta[k]["stretch"]} for k in sorted(floras)]
    if len({x["region_group"] for x in st})!=5:
        raise ValueError("Original Czech railway source 5 stretches missing")
    return st,{"original_n_taxa":len(taxa),"study":"CZ_Kutlvasr_2024",
       "original_stations":39,"n_source_study_rail_stretches":5,
       "field_original_survey_years":"2020–2022",
       "source_original_author_publisher_QA":source}

def traits_join(studies,h,s):
    species={binomial(raw) for sites in studies.values() for plot in sites
             for raw in plot["source_names"]}
    species.discard(None)
    joined={t:(math.log(h[t]),math.log(s[t])) for t in species if t in h and t in s}
    if len(joined)<300:
        raise ValueError("Three real original railway study sources do not contain enough authentic LEDA 2D taxa")
    mu=[mean(x[k] for x in joined.values()) for k in (0,1)]
    sd=[math.sqrt(mean((x[k]-mu[k])**2 for x in joined.values())) for k in (0,1)]
    if min(sd)<.05:raise ValueError("LEDA 2D axes have no genuine between taxon variance")
    return {t:((x[0]-mu[0])/sd[0],(x[1]-mu[1])/sd[1]) for t,x in joined.items()},{
         "source_unique_original_species_binomial_names_3_studies":len(species),
         "source_original_LEDA_both_canopy_height_seed_mass_species_n":len(joined),
         "trait_scale":"z log canopy_height[m] and z log seed_mass[mg], same global pooled original matched species mu/SD",
         "species_equal_global_log_mean":mu,"species_equal_global_log_SD":sd,
         "taxon_concepts_collapsed_to_binomial_to_prevent_duplicate_subspecies_trait_pseudoreplication":True}

def source_cohort(studies,traits):
    cohorts={};qa={}
    for study,src in studies.items():
        source_raw_taxa=set().union(*(p["source_names"] for p in src))
        included=[]
        site_info=[]
        group=defaultdict(list)
        valid_all=set()
        n_total_inc=sum(len(p["source_names"]) for p in src)
        n_matched_inc=0
        for i,plot in enumerate(src):
            binned={binomial(n) for n in plot["source_names"]}
            binned.discard(None)
            valid=binned&set(traits)
            if len(valid)<10:
                raise ValueError("Actual source original station had fewer than ten source-authentic LEDA joint-trait taxa: "+str((study,i,len(valid))))
            n_matched_inc+=len(valid)
            valid_all.update(valid)
            included.append({"station":plot["station"],"study":study,
                "group":plot["region_group"],"taxa":set(valid),
                "n_original_noncollapsed_taxon_names":len(plot["source_names"])})
            group[plot["region_group"]].append(i)
            site_info.append(len(valid))
        if len(group)<2 or min(len(v) for v in group.values())<2:
            raise ValueError("Original source geographical regional classes not large enough for conditional null")
        cohorts[study]=included
        qa[study]={
          "source_total_author_station_n":len(src),
          "source_original_unique_noncollapsed_taxon_names":len(source_raw_taxa),
          "source_exact_joint_LEDA_known_species_n":len(valid_all),
          "source_original_author_station_by_taxon_incidence":n_total_inc,
          "source_exact_LEDA_joint_study_station_taxon_incidence":n_matched_inc,
          "source_incidence_joint_trait_coverage_fraction":n_matched_inc/n_total_inc,
          "source_station_trait_species_min_median_max":[min(site_info),median(site_info),max(site_info)],
          "source_publisher_geography_classes_sample_sizes":{k:len(v) for k,v in group.items()},
          "all_source_sites_have_at_least10_actual_dual_trait_species":True}
    return cohorts,qa

def population_trace(taxa,trait):
    xx=[trait[i][0] for i in taxa]
    yy=[trait[i][1] for i in taxa]
    ax=mean(xx);ay=mean(yy)
    return sum((x-ax)**2 for x in xx)/len(xx)+sum((y-ay)**2 for y in yy)/len(yy),ax,ay

def get_metrics(sites,traits):
    if len(sites)<10:raise ValueError("Original minimum independent-station source cohort not met")
    profiles=[]
    for s in sites:
        tr,h,sd=population_trace(s["taxa"],traits)
        profiles.append((tr,h,sd))
    paired={m:[] for m in METRICS if "trace" not in m}
    for i,j in combinations(range(len(sites)),2):
        xi,yi=sites[i]["taxa"],sites[j]["taxa"]
        cwm=math.hypot(profiles[i][1]-profiles[j][1],profiles[i][2]-profiles[j][2])
        jac=len(xi&yi)/len(xi|yi)
        for prefix,value in (("functional_CWM_centroid_distance",cwm),
                             ("source_matched_species_Jaccard",jac)):
            paired[prefix+"_all_stationpairs"].append(value)
            key=("within" if sites[i]["group"]==sites[j]["group"] else "between")
            paired[prefix+"_"+key+"_geography"].append(value)
    if any(len(a)==0 for a in paired.values()):
        raise ValueError("Published source has no within/between independent-study railway station pairs")
    return {m:(mean(x) if x else None) for m,x in paired.items()} | {
       "source_station_2D_trait_trace_mean":mean(v[0] for v in profiles)}

def regional_margin(sites):
    group=defaultdict(Counter)
    for src in sites:
        group[src["group"]].update(src["taxa"])
    return {k:dict(v) for k,v in group.items()}

def one_null(sites,traits,seed,draws=DRAWS,burnin=BURNIN,gap=GAP):
    original_station_rich=[len(s["taxa"]) for s in sites]
    original_region_margin=regional_margin(sites)
    groups=defaultdict(list)
    for i,src in enumerate(sites):groups[src["group"]].append(i)
    rng=random.Random(seed)
    state=[set(s["taxa"]) for s in sites]
    gnames=sorted(groups)
    accept=0
    def swap(n):
        nonlocal accept
        for attempt in range(n):
            g=rng.choice(gnames)
            i,j=rng.sample(groups[g],2)
            a=state[i]-state[j];b=state[j]-state[i]
            if not a or not b:continue
            ta=rng.choice(sorted(a));tb=rng.choice(sorted(b))
            state[i].remove(ta);state[i].add(tb)
            state[j].remove(tb);state[j].add(ta)
            accept+=1
    swap(burnin)
    out=[]
    for rep in range(draws):
        swap(gap)
        pert=[{**s,"taxa":state[i]} for i,s in enumerate(sites)]
        if [len(x) for x in state]!=original_station_rich or regional_margin(pert)!=original_region_margin:
            raise AssertionError("Region fixed-margin original source abundance/occupancy changed")
        out.append(get_metrics(pert,traits))
    if accept<1000:raise ValueError("Original source group marginal swap chain did not move")
    return out,{"accepted_true_source_two_by_two_incidence_swaps":accept,
                "source_null_proposal_count":burnin+draws*gap,
                "original_station_counts_and_original_geography_group_species_frequencies_held_constant":True}

def summarize(observed,draw):
    ans={}
    for key in METRICS:
        val=observed[key];vec=[v[key] for v in draw]
        avg=mean(vec)
        ans[key]={"original_observed":val,
          "internal_source_geography_fixed_margin_null_mean":avg,
          "observed_minus_null":val-avg,
          "technical_null_q025_NOT_field_CI":quart(vec,.025),
          "technical_null_q975_NOT_field_CI":quart(vec,.975),
          "technical_lower_tail_rank_plus_one_NOT_field_p":(
              1+sum(x<=val for x in vec))/(1+len(vec)),
          "direction":(
             "positive means more differentiated than rail-only conditional null"
             if "distance" in key else
             "negative means less within-station 2D functional breadth"
             if "trace" in key else
             "positive means more taxonomically similar than rail-only conditional null")}
    return ans

def run(cz_s1,cz_s2,it_pdf,hr_pdf,rawheight,rawseed,draws=DRAWS,burnin=BURNIN,gap=GAP):
    studies={};source_qc={}
    for key,fn,args in (
        ("CZ_Czechia",publisher_czech,(cz_s1,cz_s2)),
        ("IT_Italy_Parma",publisher_italy,(it_pdf,)),
        ("HR_BA_2017_one_study",publisher_hr_ba,(hr_pdf,))):
        studies[key],source_qc[key]=fn(*args)
    if sum(map(len,studies.values()))!=OBSERVED_N or len(studies)!=SOURCE_STUDY_N:
        raise ValueError("3 truly independently published railway cohort totals 82 source stations failed")
    traits,tq=traits_join(studies,rawheight,rawseed)
    valid,cohortqa=source_cohort(studies,traits)
    results={}
    for a,(study,cohort) in enumerate(valid.items()):
        observed=get_metrics(cohort,traits)
        allnull=[];chains=[]
        for ch in range(CHAINS):
            sims,mix=one_null(cohort,traits,seed=SEED+10000*a+ch*7117,draws=draws,burnin=burnin,gap=gap)
            allnull.extend(sims);chains.append(mix)
        stats=summarize(observed,allnull)
        results[study]={
            "source_pub_original_field_stations_n":len(cohort),
            "source_study_country_not_independent_station_pairs":True,
            "source_geography_group_conditional_null_N":len(allnull),
            "source_technical_chain_QA":chains,
            "source_source_trait_observed_null":stats}
    signs={"CWM_distance_positive_rail_internal_source_only_n_studies":sum(
       r["source_source_trait_observed_null"][
          "functional_CWM_centroid_distance_all_stationpairs"]["observed_minus_null"]>0
       for r in results.values()),
       "2D_within_station_breadth_negative_n_studies":sum(
       r["source_source_trait_observed_null"][
          "source_station_2D_trait_trace_mean"]["observed_minus_null"]<0
       for r in results.values()),
       "number_original_independent_publisher_studies":len(results),
       "original_station_number_across_three_studies":sum(len(x) for x in valid.values())}
    return {"version":"v24.5",
       "status":"SOURCE_VERIFIED_CZ_IT_HR_BA_82_INDEPENDENT_RAIL_STATIONS_2D_LEDA_INTERNAL_TRAIT_NULL",
       "source_true_original_publisher_3_studies":source_qc,
       "source_all_three_study_taxa_shared_trait_axis_QA":tq,
       "source_3_real_study_rail_station_matched_2D_trait_coverage":cohortqa,
       "original_three_study_separate_continuous_functional_inference":results,
       "original_three_study_descriptive_sign_coherence":signs,
       "not_an_independent_natural_regional_nonrail_pool_test":True,
       "important_limitations":[
          "82 source railway stations arise from 3 original studies in Czechia, Italy, Croatia/Bosnia; HR+BA only one original study despite two countries",
          "Czech station-wide source 2020-22, Italy original 2012-13, Croatia/Bosnia 2013-14, different site size, management and survey effort",
          "All null models preserve only railway study species pool and author region classes, not independent natural local nonrail flora",
          "Two technical checkerboard chains with 100 draws are not two field studies nor a formal p-value",
          "Cannot combine 82 author stations with 89 Slovenian 100m plot locations as independent identical units",
          "LEDA NW European botanical trait coverage may be poorer and taxonomically biased in SE Europe; no genus imputation",
          "Species-level normalization merges authored subspecific occurrences within same station, not adaptive alleles",
          "No 2020-22 field station traffic OD, 2012-13 historical freight exposure, seed rain, local wind, station mowing or adaptive genotype",
          "Do not infer external railway ecosystem functional convergence or causal homogenization from the within-rail fixed-margin null alone"]}

def main():
    p=argparse.ArgumentParser()
    for name in ("cz_s1","cz_s2","it_pdf","hr_pdf","height","seed","output"):
        p.add_argument("--"+name.replace("_","-"),required=True)
    p.add_argument("--draws",type=int,default=DRAWS)
    p.add_argument("--burnin",type=int,default=BURNIN)
    p.add_argument("--gap",type=int,default=GAP)
    a=p.parse_args()
    for file,digest in ((a.height,LEDASHA["height"]),(a.seed,LEDASHA["seed"])):
        if hashlib.sha256(Path(file).read_bytes()).hexdigest()!=digest:
            raise ValueError("Source original LEDA file SHA changed")
    h,hq=read_leda(a.height,"height_m")
    s,sq=read_leda(a.seed,"seed_mass_mg")
    data=run(a.cz_s1,a.cz_s2,a.it_pdf,a.hr_pdf,h,s,
             draws=a.draws,burnin=a.burnin,gap=a.gap)
    f=Path(a.output);f.parent.mkdir(parents=True,exist_ok=True)
    f.write_text(json.dumps(data,indent=2,ensure_ascii=False),encoding="utf8")
    compact={"status":data["status"],
        "source_all_three_study_taxa_shared_trait_axis_QA":data["source_all_three_study_taxa_shared_trait_axis_QA"],
        "source_3_real_study_rail_station_matched_2D_trait_coverage":data[
             "source_3_real_study_rail_station_matched_2D_trait_coverage"],
        "original_three_study_separate_continuous_functional_inference":data[
             "original_three_study_separate_continuous_functional_inference"],
        "original_three_study_descriptive_sign_coherence":data[
             "original_three_study_descriptive_sign_coherence"]}
    print("V245_GENUINE_THREE_RAILWAY_STUDIES_82_STATIONS="+json.dumps(compact,ensure_ascii=False,separators=(",",":")),flush=True)
if __name__=="__main__":main()
