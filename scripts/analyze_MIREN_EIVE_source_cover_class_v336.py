#!/usr/bin/env python3
"""v336: MIREN mountain ROAD EIVE thermal position/width, presence vs Cover class.

The source's Cover variable is an ORDINAL COVER CLASS 1..8, not % cover.
Original long source CSV v3 checksum is fixed. This is an independent road
benchmark (4 European region codes) and MUST NOT count as railway replication.
Original site/plot coordinates/species observations are never republished.
"""
import argparse,csv,hashlib,json,math
from collections import Counter,defaultdict
from pathlib import Path
from statistics import mean
from analyze_original_MIREN_EIVE_road_vs_interior_v323 import (
    MIREN_MD5,read_miren,select_first_complete)
from analyze_MIREN_EIVE_genus_position_width_v324 import eligible_roads
from analyze_Slovenia_EIVE_temperature_width_v321 import eive_table,normal_name
from audit_SI_CZ_original_continuous_trait_coverage_v220 import binomial
from analyze_Slovenia_thermal_generalist_section_v319 import section_test

REGION_CODES=("AUT","CHE","ITA","NOR")
METRICS=("position","width")
WEIGHTS=("presence","cover_class","sqrt_cover_class")
COVER_GATES=(0.8,1.0)
MIN_NUMERIC_SPECIES=5
MIN_ROAD_SITES=4
MIN_REGION_ROADS=2

def numeric_class(v):
    if v is None or str(v).strip().casefold() in ("","na","nan","null"):
        return None
    try:
        x=float(v)
    except ValueError:
        raise ValueError("Unexpected MIREN source ordinal Cover label") from None
    if not math.isfinite(x) or x<1 or x>8:
        raise ValueError("Original ordinal cover class outside published 1..8")
    return x

def cover_by_source_plant_plot(path,only_keys):
    if hashlib.md5(Path(path).read_bytes()).hexdigest()!=MIREN_MD5:
        raise ValueError("Not the original Zenodo MIREN source v3")
    out=defaultdict(dict);qa=Counter()
    with Path(path).open("r",encoding="utf-8-sig",newline="") as f:
        rows=csv.DictReader(f)
        for k in ("Region","Road","Site","Year","Plot","Plot_id",
                  "Accepted.Name.MIREN","Cover","Protocol"):
            if k not in (rows.fieldnames or []):
                raise ValueError("Required MIREN original source ordinal cover field absent")
        for row in rows:
            qa["source_total_species_rows"]+=1
            reg=(row["Region"] or "").strip()
            road=(row["Road"] or "").strip()
            site=(row["Site"] or "").strip()
            year=(row["Year"] or "").strip()
            plot=(row["Plot"] or "").strip()
            if reg not in REGION_CODES or plot not in ("1","3") or not year.isdigit():
                continue
            if not str(row["Plot_id"] or "").strip().endswith("."+plot):
                continue
            key=(reg,road,site,int(year),plot)
            if key not in only_keys:continue
            raw=binomial(str(row["Accepted.Name.MIREN"] or "").strip())
            if raw is None:
                qa["source_selected_rows_without_unambiguous_binomial"]+=1
                continue
            sp=normal_name(raw)
            cls=numeric_class(row["Cover"])
            qa["selected_original_species_rows"]+=1
            if cls is None:
                qa["selected_rows_without_cover_class"]+=1
            else:
                qa["selected_rows_with_cover_class"]+=1
            if sp in out[key]:
                qa["duplicate_same_taxon_same_plot_year_records"]+=1
                old=out[key][sp]
                if old is not None and cls is not None and old!=cls:
                    qa["duplicate_taxon_conflicting_cover_class"]+=1
                if cls is not None:
                    out[key][sp]=max(old,cls) if old is not None else cls
            else:
                out[key][sp]=cls
    if qa["source_total_species_rows"]!=167804:
        raise ValueError("Original MIREN source field row count changed")
    return out,dict(qa)

def source_communities(roads):
    keys=set()
    for (region,road),sites in roads.items():
        for x in sites:
            keys.add((region,road,x["site"],x["year"],"1"))
            keys.add((region,road,x["site"],x["year"],"3"))
    return keys

def cwm(plot_species,cover_by_taxon,eive,metric,scheme):
    source_species=sorted(plot_species)
    eligible={s:cover_by_taxon.get(s) for s in source_species
              if cover_by_taxon.get(s) is not None}
    support=len(eligible)
    coverage=support/len(source_species) if source_species else 0
    if support==0:return {"mean":None,"source_trait_species":len(source_species),
          "covered_trait_species":0,"cover_completeness":coverage,"neff":None}
    if scheme=="presence":
        w={s:1 for s in eligible}
    elif scheme=="cover_class":
        w=eligible
    elif scheme=="sqrt_cover_class":
        w={s:math.sqrt(v) for s,v in eligible.items()}
    else:raise ValueError("Not one of original MIREN cover-class sensitivities")
    sw=sum(w.values())
    p={s:v/sw for s,v in w.items()}
    return {"mean":sum(p[s]*eive[s][metric] for s in p),
            "source_trait_species":len(source_species),
            "covered_trait_species":support,
            "cover_completeness":coverage,
            "neff":1/sum(v*v for v in p.values())}

def collect_sitepanels(roads,cover,eive,minimum_coverage):
    result=defaultdict(list)
    srcqa=Counter()
    for (reg,road),pairs in sorted(roads.items()):
        for p in pairs:
            k1=(reg,road,p["site"],p["year"],"1")
            k3=(reg,road,p["site"],p["year"],"3")
            a=cover.get(k1,{})
            b=cover.get(k3,{})
            ca=cwm(p["roadside"],a,eive,"position","presence")
            cb=cwm(p["inland"],b,eive,"position","presence")
            srcqa["original_qualified_first_visit_site_pairs_before_cover_gates"]+=1
            if (min(ca["covered_trait_species"],cb["covered_trait_species"])<MIN_NUMERIC_SPECIES
                 or min(ca["cover_completeness"],cb["cover_completeness"])<minimum_coverage):
                srcqa["dropped_cover_missingness_or_min_species"]+=1
                continue
            result[(reg,road)].append({"cover_roadside":a,"cover_inland":b,
                "source_roadside":p["roadside"],"source_inland":p["inland"],
                "year":p["year"]})
    roads_good={k:v for k,v in result.items() if len(v)>=MIN_ROAD_SITES}
    regions=Counter(k[0] for k in roads_good)
    filtered={k:v for k,v in roads_good.items()
              if regions[k[0]]>=MIN_REGION_ROADS}
    return filtered,dict(srcqa)

def mean_by_region(roads,eive,metric,scheme):
    region=defaultdict(list)
    road_counts=defaultdict(int)
    site_counts=defaultdict(int)
    eff_n=defaultdict(list)
    for (reg,road),sites in sorted(roads.items()):
        d=[]
        for row in sites:
            a=cwm(row["source_roadside"],row["cover_roadside"],eive,metric,scheme)
            b=cwm(row["source_inland"],row["cover_inland"],eive,metric,scheme)
            if a["mean"] is None or b["mean"] is None:raise ValueError("Invalid filtered source cover plot")
            d.append(a["mean"]-b["mean"])
            eff_n[reg].append((a["neff"],b["neff"]))
        region[reg].append(mean(d))
        road_counts[reg]+=1;site_counts[reg]+=len(sites)
    stats={reg:{"n_road_ids":road_counts[reg],"n_original_first_visit_sites":site_counts[reg],
             "road_equal_cwm_roadside_minus_inland":mean(rr),
             "positive_road_ids":sum(x>0 for x in rr),
             "negative_road_ids":sum(x<0 for x in rr),
             "mean_effective_species_number_roadside":mean(x[0] for x in eff_n[reg]),
             "mean_effective_species_number_inland":mean(x[1] for x in eff_n[reg])}
          for reg,rr in sorted(region.items())}
    vals=[v["road_equal_cwm_roadside_minus_inland"] for v in stats.values()]
    return {"metric":metric,"species_weighting":scheme,
        "n_source_regions":len(stats),"n_source_roads":sum(road_counts.values()),
        "n_original_source_site_year_first_visits":sum(site_counts.values()),
        "region_equal_mean_roadside_minus_inland_EIVE":mean(vals) if vals else None,
        "positive_regions":sum(x>0 for x in vals),
        "negative_regions":sum(x<0 for x in vals),
        "source_region_signflip_diagnostic":section_test(region),
        "region_specific_results":stats}

def analyze(miren_file,eive_file):
    original,originalqa=read_miren(miren_file)
    eive,eiveqa=eive_table(eive_file)
    paired=select_first_complete(original,2010)
    road_source=eligible_roads(paired,eive)
    raw_keys=source_communities(road_source)
    coverage,coverageqa=cover_by_source_plant_plot(miren_file,raw_keys)
    results=[]
    for mincov in COVER_GATES:
        cohort,subqa=collect_sitepanels(road_source,coverage,eive,mincov)
        panel={"minimum_cover_class_completeness":mincov,
           "source_initial_road_site_counts":{"roads":len(road_source),
                "sites":sum(len(sites) for sites in road_source.values()),
                "regions":len(set(x[0] for x in road_source))},
           "source_cover_site_eligibility":subqa,
           "source_result_original_roads":len(cohort),
           "source_result_original_sites":sum(len(x) for x in cohort.values()),
           "panels":[]}
        for metric in METRICS:
            for scheme in WEIGHTS:
                panel["panels"].append(mean_by_region(cohort,eive,metric,scheme))
        results.append(panel)
    return {"version":"v336","original_MIREN_md5":originalqa["source_md5"],
        "original_EIVE_md5":eiveqa["EIVE_original_md5"],
        "source_cover_class_qa":coverageqa,
        "results":results,
        "limitations":[
         "MIREN source cover is 1-8 ordinal class, NOT actual percent cover",
         "Cover class numeric weighting is exploratory rank weighting, not calibrated live biomass",
         "Original source 2007/08 data lack cover, use first complete paired observation 2010 onward",
         "Only original AUT/CHE/ITA/NOR region codes, one MIREN global field programme, not 4 independent studies",
         "No claim of railway replication or direct causal difference between railway BS/BT and road/inland contrasts",
         "Missing cover values cause eligible source cohort changes; compare occurrence and cover on exactly same retained cover-known taxa",
         "Duplicate same-plot species cover classes are conservatively summarized by maximum; conflicts audited",
         "No original site-by-taxon data or original coordinates redistributed",
         "Source data paper recommends notifying MIREN Steering Committee for substantial publications"]}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--miren-csv",required=True)
    p.add_argument("--eive-xlsx",required=True)
    p.add_argument("--output",required=True)
    a=p.parse_args()
    out=analyze(a.miren_csv,a.eive_xlsx)
    dest=Path(a.output);dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps(out,ensure_ascii=False),flush=True)

if __name__=="__main__": main()
