#!/usr/bin/env python3
"""v14.3: independent REAL Czech 2024 Preslia 39 station contemporary
flora taxonomic vs Pladias functional-profile section clustering.

Original field authors Kutlvašr et al. (2024), Preslia 96:247-266,
doi 10.23855/preslia.2024.247, 2020–2022 actual field surveys
39 railway stations, 5 original railway stretches, three managed
habitat zones per STATION, 11139 species observations/764 recent taxa.
The original publisher S1 79261-byte XLSX source and corrected
erratum S2 13533-byte original XLSX are SHA-pinned and checked,
not republished as raw tables.
S1 1053 taxa include historical-only source records; ONLY actual
39 "rec" station columns with 6941 literal "x" marks are used.
S2 original 39 station codes and source stretch labels drive groups.

Exact Pladias publicly retrievable source 3147 strategy and 4599
lifeform records (original Sádlo wind archetype and Raunkiaer
"therophyte" proxy), with strict normalized exact full taxon names.
A Pladias therophyte proxy is NOT Slovenian original strict annual
lifespan, so these two sources do NOT share exact guild definitions,
and effects will NEVER be naively pooled into common beta.
Missing Pladias traits remain missing and no genus guessing.

First identify 39-station original sample taxon+lifeform coverage.
If coverage too weak to estimate meaningful original site profiles,
return not estimable (never equate missing to non-wind or native).
Estimate source within-stretch taxonomic overlap versus other
stretch and coarsened species functional profile overlap.
Use species-by-station frequency- and guild-by-station richness-
preserving checkerboard swaps, a descriptive null without
geographic distance, freight, climate or management controls.
No contemporaneous station GPS in author's S1/S2: extremely
strong spatial and unequal-area survey confounding. This is
one-country station-scale exploratory validation ONLY and is
NOT true train seed carriage/flight or local evolutionary adaptation.
"""
from __future__ import annotations
import argparse,copy,hashlib,json,math,random,zipfile
from collections import Counter,defaultdict
from itertools import combinations
from pathlib import Path
from statistics import mean,median

from analyze_Czech_2024_rail_station_area_aliens_v108 import (
    EXPECTED_S1_SHA,EXPECTED_S2_SHA
)
from analyze_CzechVeg_Pladias_rail_road_diaspore_v119 import (
    parse_traits,normalize_taxon
)

PAPER="10.23855/preslia.2024.247"
GROUPS=("th_wind","th_other","nonth_wind","nonth_other")
WIND={"Epilobium","Phragmites"}
DRAWS=300
BURNIN=24000
SWAP_BETWEEN=3000
SEED=20261115


def read_original_czech_publisher(s1,s2):
    # Validate publisher-original fingerprints BEFORE loading an optional
    # Excel parser; even the minimal built-in CI must fail closed on
    # counterfeit source workbooks without requiring openpyxl installed.
    expected=((s1,EXPECTED_S1_SHA,79261),(s2,EXPECTED_S2_SHA,13533))
    for path,sha,length in expected:
        raw=Path(path).read_bytes()
        if len(raw)!=length or hashlib.sha256(raw).hexdigest()!=sha:
            raise ValueError("Czech 2024 original publisher S1 or corrected S2 source SHA/bytes changed")
        if not zipfile.is_zipfile(path):
            raise ValueError("Czech original publisher source Excel invalid")
        with zipfile.ZipFile(path) as z:
            if z.testzip() is not None:raise ValueError("Czech publisher source XLSX CRC failed")
    from openpyxl import load_workbook
    ws1=load_workbook(s1,read_only=True,data_only=True)
    ws2=load_workbook(s2,read_only=True,data_only=True)
    try:
        rows=list(next(iter(ws1.worksheets)).values)
        meta=list(next(iter(ws2.worksheets)).values)
    finally:
        ws1.close();ws2.close()
    recent=[]
    current=None
    for col in range(2,len(rows[4])):
        if rows[4][col] is not None and str(rows[4][col]).strip():
            current=str(rows[4][col]).strip()
        mark=rows[5][col]
        if mark is not None and str(mark).strip().lower()=="rec":
            recent.append((col,current))
    if len(recent)!=39 or len({x[1] for x in recent})!=39:
        raise ValueError("Exact 39 Czech original recent railway station columns not identified")
    stats={}
    for row in meta[5:]:
        if len(row)<11 or row[1] is None or not str(row[1]).strip():
            continue
        code=str(row[1]).strip()
        if code in stats:raise ValueError("Original Czech source S2 station duplicated")
        stats[code]={
           "stretch":str(row[2]).strip(),
           "original_published_total_taxa":int(row[3]),
           "original_source_station_total_area_ha":float(row[10]) if row[10] is not None and str(row[10]).strip()!="NA" else None,
        }
    if len(stats)!=39 or {code for _,code in recent}!=set(stats):
        raise ValueError("Czech original publisher S1 versus corrected S2 station codes mismatch")
    if len({x["stretch"] for x in stats.values()})!=5:
        raise ValueError("Original 5 Czech source rail stretches not recovered")
    floras={code:set() for _,code in recent}
    source_rows=0;source_rec_markers=0;source_taxa={}
    for row in rows[7:]:
        if len(row)<2 or row[1] is None or not str(row[1]).strip():continue
        name=normalize_taxon(row[1])
        if name in source_taxa:raise ValueError("Original S1 1053 source taxon names duplicated after strict normalization")
        source_taxa[name]=str(row[0] or "").strip()
        source_rows+=1
        for col,code in recent:
            m=str(row[col] or "").strip().lower()
            if m and m!="x":raise ValueError("Czech S1 current field survey non-X marker, do not recode historical")
            if m=="x":
                floras[code].add(name)
                source_rec_markers+=1
    if source_rows!=1053 or source_rec_markers!=6941:
        raise ValueError("Original Czech 1053 taxon rows + 6941 contemporary source site incidences not conserved")
    known_recent=set().union(*floras.values())
    if len(known_recent)!=761:
        # Original SOURCE S1 has 1053 published taxonomy rows but only 761
        # appear in contemporary rec columns; paper ABSTRACT reports 764
        # taxa across a finer three-zone survey grain. Those 3 absent
        # from the released rec presence columns CANNOT be fabricated.
        raise ValueError("Original Czech released S1 recent 761 taxa + 6941 rec X source sample differs; n="+str(len(known_recent)))
    source_errors={c:len(s)-stats[c]["original_published_total_taxa"] for c,s in floras.items()
                   if len(s)!=stats[c]["original_published_total_taxa"]}
    if source_errors!={"DV":4,"SV":1,"TI":1}:
        raise ValueError("Exact known original Czech S1 vs erratum S2 3 station richness discrepancies changed")
    return floras,stats,{
      "original_source_Czech_railway_Preslia_DOI":PAPER,
      "source_S1_SHA256":EXPECTED_S1_SHA,
      "source_S2_corrected_erratum_SHA256":EXPECTED_S2_SHA,
      "source_original_recent_station_IDs":len(recent),
      "source_original_railway_stretches":5,
      "n_original_recent_real_taxa_in_public_S1_station_rec_columns":len(known_recent),
      "n_original_author_paper_reported_species_from_all_zone_records":764,
      "observed_source_vs_paper_taxon_count_difference":764-len(known_recent),
      "DO_NOT_INVENT_3_TAXA_MISSING_FROM_PUBLIC_REC_COLUMNS":True,
      "source_all_historical_plus_current_taxon_rows":source_rows,
      "original_39_station_contemporary_exact_presence_marks":source_rec_markers,
      "known_original_publisher_S1_S2_station_richness_disagreements":source_errors,
      "original_recent_and_historical_not_pooled":True,
      "original_habitat_zones_aggregated_at_station_scale_NOT_true_Slovenia_BS_BT_comparison":True,
    }


def czech_guild(name,strategies,lifes):
    d=strategies.get(name);l=lifes.get(name)
    if d is None or l is None:return None
    wind=d in WIND
    th=l.strip().lower().startswith("thero") or l.strip().lower() in ("th","t")
    return ("th_" if th else "nonth_")+("wind" if wind else "other")


def current_sites(floras,meta,strategies,lifes):
    observed=set().union(*floras.values())
    codes={k:czech_guild(k,strategies,lifes) for k in observed}
    known={k:v for k,v in codes.items() if v is not None}
    sites=[]
    for code in sorted(floras):
        taxa=floras[code]
        kept=taxa&set(known)
        sites.append({"site":code,"section":meta[code]["stretch"],
                      "known":set(kept),"all_original_taxa":set(taxa),
                      "area_ha":meta[code]["original_source_station_total_area_ha"]})
    known_n=len(known);nobs=len(observed)
    nknowninc=sum(len(s["known"]) for s in sites)
    allinc=sum(len(s["all_original_taxa"]) for s in sites)
    qc={"n_original_rail_stations":len(sites),
        "n_source_recent_species_in_original_publisher_S1_rec":nobs,
        "n_Pladias_both_strategy_therophyte_proxy_known_recent_taxa":known_n,
        "n_original_contemporary_site_taxon_incidence_all_species":allinc,
        "n_original_incidence_with_both_exact_Pladias_traits":nknowninc,
        "actual_original_species_and_site_trait_coverage":nknowninc/allinc,
        "source_field_site_min_known_taxa":min(len(s["known"]) for s in sites),
        "n_source_stretches":len({s["section"] for s in sites}),
        "source_Pladias_lifeform_therophyte_IS_PROXY_NOT_Strict_Annual_Lifespan":True}
    return sites,known,qc


def jaccard(x,y):
    return len(x&y)/len(x|y) if x|y else None


def function_sim(x,y,guild):
    if not x or not y:return None
    a=Counter(guild[k] for k in x);b=Counter(guild[k] for k in y)
    return 1-0.5*sum(abs(a[g]/len(x)-b[g]/len(y)) for g in GROUPS)


def site_pairs(sites,guild):
    wa,be,fw,fb,allw,allb=[],[],[],[],[],[]
    for i,j in combinations(range(len(sites)),2):
        a,b=sites[i],sites[j]
        same=a["section"]==b["section"]
        p=jaccard(a["known"],b["known"]) if a["known"] and b["known"] else None
        f=function_sim(a["known"],b["known"],guild)
        raw=jaccard(a["all_original_taxa"],b["all_original_taxa"])
        if p is not None:
            (wa if same else be).append(p)
            (fw if same else fb).append(f)
        (allw if same else allb).append(raw)
    def result(a,b):
        return {"n_pairs_same":len(a),"n_pairs_other":len(b),
                "mean_same":mean(a) if a else None,
                "mean_other":mean(b) if b else None,
                "gap":mean(a)-mean(b) if a and b else None}
    return {"taxonomic_pladias_known":result(wa,be),
            "functional_four_coarse_guild_overlap":result(fw,fb),
            "taxonomic_all_original_recent_S1_761_taxa":result(allw,allb)}


def swap_null(sites,guild,draws=DRAWS,burnin=BURNIN,thin=SWAP_BETWEEN,seed=SEED):
    # Global prevalence and each source station's exact functional
    # guild richness fixed, no false geographic/area matching.
    universe={g:{t for t,c in guild.items() if c==g} for g in GROUPS}
    data={g:[s["known"]&universe[g] for s in sites] for g in GROUPS}
    initial=({g:tuple(len(q) for q in data[g]) for g in GROUPS},
             Counter(t for g in GROUPS for st in data[g] for t in st))
    rng=random.Random(seed)
    def advance(proposals):
        accept=0
        for k in range(proposals):
            i,j=rng.sample(range(len(sites)),2)
            rows=data[rng.choice(GROUPS)]
            a=rows[i]-rows[j];b=rows[j]-rows[i]
            if not a or not b:continue
            u=rng.choice(tuple(a));v=rng.choice(tuple(b))
            rows[i].remove(u);rows[i].add(v);rows[j].remove(v);rows[j].add(u)
            accept+=1
        return accept
    acc=advance(burnin)
    n_pairs=len(sites)*(len(sites)-1)//2
    if n_pairs!=741:raise ValueError("39 author station source matrix cannot be relabelled as independent vegetation plots")
    observed=site_pairs(sites,guild)["taxonomic_pladias_known"]["gap"]
    if observed is None:return {"identified":False,"cause":"no original station dual-known Pladias source pairs"}
    null=[]
    for _ in range(draws):
        acc+=advance(thin)
        rebuilt=[set().union(*(data[g][i] for g in GROUPS)) for i in range(len(sites))]
        total=({g:tuple(len(q) for q in data[g]) for g in GROUPS},
               Counter(t for g in GROUPS for st in data[g] for t in st))
        if total!=initial:raise ValueError("Original Czech source station guild or overall taxon prevalence not preserved")
        within=[];between=[]
        for i,j in combinations(range(len(sites)),2):
            score=jaccard(rebuilt[i],rebuilt[j])
            if score is None:continue
            (within if sites[i]["section"]==sites[j]["section"] else between).append(score)
        null.append(mean(within)-mean(between))
    z=sorted(null)
    return {"identified":True,
            "observed_original_stretch_species_identity_gap":observed,
            "guild_and_species_frequency_fixed_null_expected_gap":mean(null),
            "excess_original_over_guild_frequency_reference":observed-mean(null),
            "null_chain_2_5":z[round((len(z)-1)*.025)],
            "null_chain_97_5":z[round((len(z)-1)*.975)],
            "plus_one_conditional_shuffle_upper_tail_NOT_causal_p":
                (sum(v>=observed for v in null)+1)/(len(null)+1),
            "draws_autocorrelated":len(null),
            "swap_proposals":burnin+thin*draws,
            "accepted_source_guild_checkerboard_swaps":acc,
            "source_original_station_size_and_species_frequency_margins_preserved":True,
            "no_GPS_road_or_active_train_freight_adjustment_possible_from_S1_S2":True}


def run(s1,s2,disp_xlsx,life_xlsx):
    floras,meta,source=read_original_czech_publisher(s1,s2)
    d,l,traitqc=parse_traits(disp_xlsx,life_xlsx)
    sites,guild,coverage=current_sites(floras,meta,d,l)
    stat=site_pairs(sites,guild)
    if coverage["actual_original_species_and_site_trait_coverage"]<.25 or coverage["source_field_site_min_known_taxa"]<3:
        null={"identified":False,"cause":"Pladias source strategy/lifeform exact join yields insufficient source site coverage"}
    else:null=swap_null(sites,guild)
    return {
      "complete":True,"version":"v14.3",
      "original_Preslia_2024_czech_station_publisher_provenance":source,
      "original_Pladias_strategy_and_lifeform_author_provenance":traitqc,
      "original_2020_22_39_stations_Pl adias_exact_guild_coverage".replace(" ",""):coverage,
      "true_author_source_recent_flora_taxonomic_and_functional_comparison":stat,
      "original_stretch_plant_identity_shuffle_diagnostic":null,
      "not_a_slovenian_study_replication":{
          "original_2024_Czech_site_unit":"39 unequal-area stations 3 zones aggregated, 5 railway stretches",
          "original_2026_Slovenian_site_unit":"89 100-m 2022/23 railway transects BS and BT union, 15 sections",
          "Pladias_lifeform_therophyte_vs_original_Slovenia_strict_annual_not_same_trait_ontology":True,
          "no_Czech_original_station_GPS_or_road_control_in_these_S1_S2_appendices":True},
      "limits":[
          "Czech source S1 is 2020–22 station-combined contemporary, not original 39 station×3 habitat×plot observations; all habitat zones pooled per station.",
          "Czech source S1 rec x marks are actual 39 publisher-station recent observations; 761 unique recent S1 taxa vs paper 764 taxa across three-zone field survey (3 absent from published recent station rec matrix), never invent missing three.",
          "Historical uneven a/b/c/d records excluded.",
          "Original 39 station area ranges widely, so source taxonomic overlap can be affected by uneven effort/area; not fully adjusted in this first pilot.",
          "Five source railway stretches not independently sampled country studies. No GPS, actual freight OD, parallel road routes or measured seed carriage in S1/S2.",
          "Pladias public plant taxa traits may only match a subset of Czech original 764 recent species; no source unknown traits imputed.",
          "Pladias therophyte is lifeform proxy NOT the Slovenian strict annual lifespan, thus source crosscountry functional categories are not literally identical.",
          "Randomized species identities inside functional groups preserve source species occurrence totals; it cannot rule out geographic source pool, landcover, microhabitat heterogeneity or urbanization.",
          "No event-specific freight vector, mowing schedule, living plant fitness, phenotypic plasticity or heritable adaptation measured.",
          "No combined country-level effect estimated from dissimilar Czech station versus Slovenian railway transect units."
      ]
    }


def main():
    p=argparse.ArgumentParser()
    for n in ("source-s1","source-s2","pladias-diaspore","pladias-lifeform","output"):
        p.add_argument("--"+n,required=True)
    a=p.parse_args()
    out=run(a.source_s1,a.source_s2,a.pladias_diaspore,a.pladias_lifeform)
    dest=Path(a.output);dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(out,indent=2,ensure_ascii=False),encoding="utf-8")
    print(json.dumps({
       "complete":True,"version":"v14.3",
       "source":out["original_Preslia_2024_czech_station_publisher_provenance"],
       "trait_coverage":out["original_2020_22_39_stations_Pladias_exact_guild_coverage"],
       "actual_study":out["true_author_source_recent_flora_taxonomic_and_functional_comparison"],
       "guild_null":out["original_stretch_plant_identity_shuffle_diagnostic"],
       "no_geo_matched_rail_freight_replication":True},ensure_ascii=False),flush=True)


if __name__=="__main__":main()
