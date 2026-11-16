#!/usr/bin/env python3
"""v19.0 CzechVeg source-based, non-rail regional reference-pool opportunity audit.

This is a HEADER-only audit, not a plant compositional test. True
CzechVeg-Open v1.1 official 117739 vegetation relevés with original study,
geobotanical region, year and area. Natural analogue candidates are text-
nontransport plots in taxonomically published TH (dry grasslands), TF
(pioneer sand/shallow soils), and TG (sandy steppes) vegetation classes.
These classes include secondary or managed sites; natural origin NOT assured.
Rail/road labels are original study TEXT, not verified operating track.
No RAIL/CAR causal interpretation or unobserved vegetation species inferred.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,re,sys
from collections import Counter,defaultdict
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
from analyze_CzechVeg_Pladias_rail_road_diaspore_v119 import plot_labels
from analyze_CzechVeg_rail_road_region_year_area_overlap_v121 import (
    area_key,year_from_source
)
from audit_CzechVeg_shared_railroad_reference_studies_v122 import source_title

DOI="10.5281/zenodo.19231152"
SOURCE_MD5="0bc920462d1f9eac473e7b2a231e6782"
SOURCE_N=76386498
ROWS=117739
ANALOG_CODES={"TH":"dry grassland","TF":"pioneer sand/shallow soil","TG":"sandy steppe"}
CLASSES=["TH","TF","TG"]
MIN_ARM=2

def md5_file(path):
    md=hashlib.md5()
    with open(path,"rb") as f:
        for b in iter(lambda:f.read(4*1024*1024),b""):
            md.update(b)
    return md.hexdigest()

def class_two_char(v):
    val=str(v or "").strip()
    if val in ("NA","", "N/A","NULL"):
        return None
    m=re.fullmatch(r"([A-Z]{2})(?:[A-Z]\d{0,3}[a-z]?|)",val)
    return m.group(1) if m else None

def get_meta(header):
    if Path(header).stat().st_size!=SOURCE_N or md5_file(header)!=SOURCE_MD5:
        raise ValueError("Not exact publisher original CzechVeg-Open 2026 v1.1 headers")
    source_labels,label_receipt=plot_labels(header)
    sites={}
    candidate=Counter()
    roots=Counter()
    with open(header,encoding="utf-8-sig",newline="") as fh:
        rdr=csv.DictReader(fh)
        expected={"RELEVE_NR","YEAR","SURF_AREA","SKALICKY1988","REFERENCE_TITLE",
                  "HABITAT","LOCALITY","SYNTAXON"}
        if expected-set(rdr.fieldnames or []):
            raise ValueError("Missing source Czech vegetation geography, study or syntaxon fields")
        for row in rdr:
            candidate["original_headers_scanned"]+=1
            key=str(row.get("RELEVE_NR") or "").strip()
            if key not in source_labels:
                candidate["rejected_plot_area"]+=1;continue
            year=year_from_source(row.get("YEAR"))
            region=str(row.get("SKALICKY1988") or "").strip()
            paper=source_title(row.get("REFERENCE_TITLE"))
            if year is None or not region or region=="NA":
                candidate["rejected_year_or_region"]+=1;continue
            if paper is None:
                candidate["rejected_missing_study"]+=1;continue
            area=area_key(source_labels[key]["area"],"class4")
            cls=class_two_char(row.get("SYNTAXON"))
            if cls:roots[cls]+=1
            sites[key]={"study":paper,"region":region,"decade":year//10*10,
                        "year":year,"area":area,"label":source_labels[key]["broad"],
                        "strict":source_labels[key]["strict"],"source_class":cls}
            candidate["eligible_all_sample_metadata"]+=1
            if cls in ANALOG_CODES:candidate["analogous_class_"+cls]+=1
    if candidate["original_headers_scanned"]!=ROWS:
        raise ValueError("Original full source header row count changed")
    return sites,{
        "publisher_DOI":DOI,"publisher_original_headers_md5":SOURCE_MD5,
        "source_head_count":ROWS,
        "source_columns_exact":True,
        "source_quality":dict(candidate),
        "source_vegetation_standalone_class_root_counts":{k:roots[k] for k in CLASSES},
        "plot_label_source_receipt":label_receipt["counts"]
    }

def matched_rail_source(sites):
    groups=defaultdict(lambda:{"rail":[],"road":[]})
    for sid,meta in sites.items():
        label=meta["label"]
        if label not in ("rail","road"):continue
        k=(meta["study"],meta["region"],meta["decade"],meta["area"])
        groups[k][label].append(sid)
    common={k:v for k,v in groups.items()
            if len(v["rail"])>=MIN_ARM and len(v["road"])>=MIN_ARM}
    return common,{
        "matched_source_railroad_study_region_decade_area_blocks":len(common),
        "studies":len({k[0] for k in common}),
        "regions":len({k[1] for k in common}),
        "n_rail_plots":sum(len(v["rail"]) for v in common.values()),
        "n_road_plots":sum(len(v["road"]) for v in common.values()),
        "source_publication_grain_requires_independent_plant_records":True
    }

def opportunity(sites,matched):
    bygeo=defaultdict(list)
    for pid,m in sites.items():
        # A source unlabeled relevé is a background observation *candidate*,
        # not a verified railway-free/natural habitat by negation.
        if m["label"] is None and m["strict"] is None:
            bygeo[(m["region"],m["decade"],m["area"])].append(pid)
    stats=Counter()
    region_support=defaultdict(lambda:{"rail":0,"analog_rail_with_geographic_support":0,
                                       "broad_rail_with_geographic_support":0})
    rows=[]
    for k,v in sorted(matched.items()):
        study,reg,dec,area=k
        others=[pid for pid in bygeo[(reg,dec,area)] if sites[pid]["study"]!=study]
        analog=[pid for pid in others if sites[pid]["source_class"] in CLASSES]
        analog_classes=Counter(sites[pid]["source_class"] for pid in analog)
        row={"region":reg,"decade":dec,"area_band":area,
             "rail_n":len(v["rail"]),"road_n":len(v["road"]),
             "independent_unlabeled_source_plots":len(others),
             "independent_analog_source_plots":len(analog),
             "analog_class_counts":dict(analog_classes),
             "independent_unlabeled_source_studies":len({sites[p]["study"] for p in others}),
             "independent_analog_source_studies":len({sites[p]["study"] for p in analog}),
             "surveyed_2010plus":dec>=2010}
        # Do not expose original exact publication titles, locality, plot IDs,
        # site coordinates or full copyrighted field raw data in GitHub report.
        rows.append(row)
        stats["blocks"]+=1
        stats["rail_plots"]+=len(v["rail"])
        stats["source_nontransport_geo_year_area_support"]+=bool(others)
        stats["source_analog_geo_year_area_support"]+=bool(analog)
        for threshold in (1,2,5,10):
            stats["rail_plots_unlabeled_pool_min"+str(threshold)]+=len(v["rail"]) if len(others)>=threshold else 0
            stats["rail_plots_analog_pool_min"+str(threshold)]+=len(v["rail"]) if len(analog)>=threshold else 0
        stats["rail_plots_2010plus"]+=len(v["rail"]) if dec>=2010 else 0
        if len(analog)>=5:
            stats["rail_plot_georegions_analog5_marker"]+=len(v["rail"])
        region_support[reg]["rail"]+=len(v["rail"])
        if len(analog)>=5:region_support[reg]["analog_rail_with_geographic_support"]+=len(v["rail"])
        if len(others)>=5:region_support[reg]["broad_rail_with_geographic_support"]+=len(v["rail"])
    stats["distinct_regions_with_analog_min5"] = sum(q["analog_rail_with_geographic_support"]>0 for q in region_support.values())
    stats["distinct_regions_with_unlabeled_min5"] = sum(q["broad_rail_with_geographic_support"]>0 for q in region_support.values())
    stats["distinct_rail_regions"] =len(region_support)
    return {
        "source_original_rail_vs_road_text_matched_publication_frame":stats["rail_plots"],
        "independent_reference_study_and_geo_decade_area_restrictions":True,
        "strict_natural_habitat_NOT_VERIFIED":True,
        "nontransport_label_is_missing_NOT_confirmed_natural":True,
        "analog_classes_are_published_TH_TF_TG_labels_not_field_measured_substrate":True,
        "country_level_natural_analogs_cannot_be_claimed_from_raw_plot_counts":True,
        "support":dict(stats),
        "regional_support_no_plot_ids":{str(k):v for k,v in sorted(region_support.items())},
        "block_support_without_study_IDs":rows,
    }

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--source-headers",required=True)
    ap.add_argument("--out",required=True)
    args=ap.parse_args()
    sites,source=get_meta(args.source_headers)
    b,m=matched_rail_source(sites)
    potential=opportunity(sites,b)
    ret={"version":"v19.0","original_source":source,
         "observed_railway_context_source_blocks":m,
         "independent_nonrail_regional_plot_opportunity":potential,
         "no_species_matrix_or_taxonomic_effect_yet":True,
         "go_no_go_gate_for_analog_cross_region":
             potential["support"]["distinct_regions_with_analog_min5"]>=3,
         "goal":"Next version 19.1 will use real species×plot rows only for supported source plots, preserve one-study blocks, and separate all nontransport background vs specific analogue classes."
    }
    dest=Path(args.out);dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(ret,ensure_ascii=False,indent=2),encoding="utf-8")
    digest={k:v for k,v in ret.items() if k not in ("original_source","independent_nonrail_regional_plot_opportunity")}
    digest["original_source_counts"]=source["source_quality"]
    digest["reference_opportunity"]=potential["support"]
    print(json.dumps(digest,ensure_ascii=False,indent=2),flush=True)

if __name__=="__main__":main()
