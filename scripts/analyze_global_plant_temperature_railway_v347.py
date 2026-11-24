#!/usr/bin/env python3
"""v347: world-sourced plant temperature envelopes × authentic 3-continent rail floras.

Global species climate preference = Harris et al (2022) EIDC corrected tmean
(Kelvin), derived from occurrence climatology, NOT rail-ballast temperature.
Original railway assemblage data from v314 source; actual source SHA audited
before reanalysis. All contrasts labelled by comparator design, not pooled.

Main source: https://doi.org/10.5285/ca339c86-3674-4030-b891-35326e71141e
Design: docs/tetsu_v346_global_plant_temperature_envelope_frozen_ja.md
"""
from __future__ import annotations
import argparse,csv,hashlib,json,math,re,statistics,sys,unicodedata
from collections import Counter
from pathlib import Path

MIN_ORIGINAL_SPECIES=15
COVERAGE=0.5
MIN_PAIRED_EIDC_SPECIES=3
MIN_PAIRED_STATIONS=8
SOURCE_MIN_RECORDS=(20,100)
GLOBAL_K_MIN=197.0
GLOBAL_K_MAX=320.0
C4_CANDIDATES=frozenset((
  "amaranthus caudatus","bothriochloa ischaemum","cynodon dactylon",
  "cyperus rotundus","echinochloa colona","euphorbia prostrata",
  "fimbristylis dichotoma","saccharum spontaneum","setaria viridis",
  "sorghum halepense"
))
REFS=(
  ("CZ","road:CZ","MIREN_mountain_ROAD_not_local_nonrail_pair"),
  ("IT","road:IT","MIREN_mountain_ROAD_not_local_nonrail_pair"),
  ("SI","background:SI:urban","URBAN_tree_pits_not_local_natural_nonrail_pair"),
  ("HRBA","background:HRBA:urban","URBAN_tree_pits_not_local_natural_nonrail_pair"),
  ("SK","background:SK:urban","URBAN_vegetation_not_local_natural_nonrail_pair"),
  ("IN","background:IN:interior","same_station_2014_nearby_interior_not_guaranteed_pristine"),
  ("IN","road:IN","MIREN_mountain_ROAD_not_local_nonrail_pair"),
)
PANELS=(("all_source_taxa",frozenset()),("exclude_Cynodon",frozenset(("cynodon dactylon",))),("exclude_10_C4_candidates",C4_CANDIDATES))

def simple_binomial(v):
    x=" ".join(unicodedata.normalize("NFKC",str(v or "")).strip().lower().split())
    if re.fullmatch(r"[a-z-]+ [a-z-]+",x):
        return x
    return None

def fieldmap(source):
    fields=source or []
    names=set(fields)
    expect={"species","numrecords","tmean"}
    if not expect<=names:
        raise ValueError("Published author source exact global species/numrecords/tmean fields missing")
    return {f:f for f in expect}

def load_global(path):
    file=Path(path)
    if not file.is_file():raise FileNotFoundError("Need authenticated original EIDC CSV bytes")
    h=hashlib.sha256()
    with file.open("rb") as f:
        for block in iter(lambda:f.read(2**20),b""):h.update(block)
    names={}
    blocked=set()
    qa=Counter()
    temp_min=math.inf;temp_max=-math.inf
    with file.open(encoding="utf-8-sig",newline="") as f:
        reader=csv.DictReader(f)
        fields=fieldmap(reader.fieldnames)
        for row in reader:
            qa["source_total_rows"]+=1
            name=simple_binomial(row.get("species"))
            if not name:
                qa["source_name_not_unambiguous_binomial"]+=1
                continue
            try:
                t=float(row["tmean"])
                n=float(row["numrecords"])
            except (ValueError,TypeError,KeyError):
                qa["source_empty_temperature_or_record_count"]+=1
                continue
            if not math.isfinite(t) or not math.isfinite(n) or
                t<GLOBAL_K_MIN or t>GLOBAL_K_MAX or n<0 or not n.is_integer():
                qa["source_invalid_temperature_or_records"]+=1
                continue
            temp_min=min(temp_min,t);temp_max=max(temp_max,t)
            tr=(t-273.15,int(n))
            if name in names and names[name]!=tr:
                blocked.add(name)
                qa["source_conflicting_duplicate_taxon_rows"]+=1
            elif name in names:qa["source_exact_duplicate_rows"]+=1
            else:names[name]=tr
    for name in blocked:names.pop(name,None)
    if qa["source_total_rows"]<100000:
        raise ValueError("EIDC global raw original source species row census under 100k: data incomplete")
    if len(names)<90000:
        raise ValueError("Too few exact species binomials: official global climate table not represented")
    return names,{"source_DOI":"10.5285/ca339c86-3674-4030-b891-35326e71141e",
                 "original_sha256":h.hexdigest(),"original_file_bytes":file.stat().st_size,
                 "original_source_CSV_header":reader.fieldnames,
                 "global_original_plant_source_qa":dict(qa),
                 "valid_original_exact_species_binomials":len(names),
                 "verified_global_K_min_max":[temp_min,temp_max],
                 "physically_K_to_C_conversion":273.15,
                 "field_definition":"Author source EVI/effort corrected species climate average tmean, not measured rail habitat surface heat"}

def original_source_value(v):
    if v is None or isinstance(v,bool):return None
    f=float(v)
    if not math.isfinite(f) or f<0 or f>1:
        raise ValueError("Original railway published site frequency invalid")
    return f

def summarize(freq,traits,records_min,exclude=frozenset(),source_exclude=frozenset(),equal=False):
    if records_min<0:raise ValueError("Invalid original occurrence support threshold")
    excluded={simple_binomial(s) for s in exclude|source_exclude}
    original={simple_binomial(s):original_source_value(v) for s,v in freq.items()
              if simple_binomial(s) is not None and original_source_value(v)>0
              and simple_binomial(s) not in excluded}
    known={name:(val,original[name]) for name,val in traits.items()
           if name in original and val[1]>=records_min}
    n=len(original);nvalid=len(known)
    mass=sum(original.values());nummass=sum(f for value,f in known.values())
    numeric_known=(nvalid>=MIN_ORIGINAL_SPECIES and nvalid/n>=COVERAGE
                    and nummass/mass>=COVERAGE) if n and mass else False
    total=sum(1 if equal else f for x,f in known.values())
    cwm=sum((1 if equal else f)*val[0] for val,f in known.values())/total if numeric_known and total else None
    return {"estimable":cwm is not None,
            "source_named_taxa_n":n,
            "source_global_temperature_covered_species_n":nvalid,
            "source_species_trait_coverage":nvalid/n if n else None,
            "source_frequency_covered_mass":nummass/mass if mass else None,
            "climate_tmean_C":cwm,"global_trait_min_original_GBIF_records":records_min,
            "species_equal_vs_published_frequency":"equal_species" if equal else "author_frequency",
            "source_unknown_species_not_imputed":True}

def pairs_kashmir(pairs,traits,records_min,exclude=frozenset()):
    vals=[]
    cover=[]
    for a,b in pairs:
        ra=summarize({v:1 for v in a},traits,records_min,exclude,equal=True)
        rb=summarize({v:1 for v in b},traits,records_min,exclude,equal=True)
        if (ra["source_global_temperature_covered_species_n"]>=MIN_PAIRED_EIDC_SPECIES
            and rb["source_global_temperature_covered_species_n"]>=MIN_PAIRED_EIDC_SPECIES
            and ra["source_species_trait_coverage"]>=COVERAGE
            and rb["source_species_trait_coverage"]>=COVERAGE):
            # The fixed ratio and n=3 support gate supersedes regional source n>=15
            t_a=[traits[simple_binomial(sp)][0] for sp in a
                  if simple_binomial(sp) in traits and simple_binomial(sp) not in exclude
                  and traits[simple_binomial(sp)][1]>=records_min]
            t_b=[traits[simple_binomial(sp)][0] for sp in b
                  if simple_binomial(sp) in traits and simple_binomial(sp) not in exclude
                  and traits[simple_binomial(sp)][1]>=records_min]
            vals.append(statistics.mean(t_a)-statistics.mean(t_b))
            cover.append((ra["source_species_trait_coverage"],rb["source_species_trait_coverage"]))
    return {"original_station_candidate_pairs":len(pairs),
            "original_paired_source_stations_eligible":len(vals),
            "min_original_eligible_stations":MIN_PAIRED_STATIONS,
            "estimable":len(vals)>=MIN_PAIRED_STATIONS,
            "same_station_equal_species_temperature_difference_C":
                statistics.mean(vals) if len(vals)>=MIN_PAIRED_STATIONS else None,
            "n_positive_sites":sum(v>0 for v in vals) if len(vals)>=MIN_PAIRED_STATIONS else None,
            "n_negative_sites":sum(v<0 for v in vals) if len(vals)>=MIN_PAIRED_STATIONS else None,
            "source_stations_not_independent_climate_regions":True,
            "source_2014_nearby_interior_not_confirmed_unimpacted_nature":True}

def analyze(source,traits,pairs,metadata):
    out={"version":"v347-genuine-global-species-tmean-K-rail-assemblages",
         "source_EIDC":metadata,"rails":[],"comparators":[],"original_Kashmir_paired":[],
         "the_input_railway_surveys_are_not_worldwide_balanced":"Europe_5_research_programmes;Asia_1;Oceania_1",
         "trait_global_species_distribution_is_not_a_global_railway_microhabitat_causal_effect":True,
         "original_published_railway_sources_are_sha_pinned":True}
    for n in SOURCE_MIN_RECORDS:
        for equ in (False,True):
            for panel,excluded in PANELS:
                for study,fr in source["rail"].items():
                    original_unavailable=source["unavailable"].get(study,frozenset())
                    x=summarize(fr,traits,n,excluded,original_unavailable,equal=equ)
                    out["rails"].append({"source_study":study,"trait_min_source_records":n,
                                      "taxon_panel":panel,**x})
                for study,ref,design in REFS:
                    r=summarize(source["rail"][study],traits,n,excluded,
                      source["unavailable"].get(study,frozenset()),equal=equ)
                    reference=(source["road"][ref.split(":")[1]] if ref.startswith("road:")
                               else source["background"][ref.split(":")[1]][ref.split(":")[2]])
                    b=summarize(reference,traits,n,excluded,
                       source["unavailable"].get(ref,frozenset()),equal=equ)
                    contrast=(r["climate_tmean_C"]-b["climate_tmean_C"]
                        if r["estimable"] and b["estimable"] else None)
                    out["comparators"].append({
                       "study":study,"documented_source_reference_kind":design,
                       "source_reference":ref,"trait_min_source_records":n,
                       "species_weighting":"equal" if equ else "author_frequency",
                       "taxon_panel":panel,
                       "rail":r,"reference":b,
                       "rail_minus_reference_global_species_tmean_C":contrast,
                       "different_reference_designs_are_not_one_pooled_effect":True,
                    })
        for label,excluded in PANELS:
            out["original_Kashmir_paired"].append({
                "taxon_panel":label,
                "min_source_records":n,
                **pairs_kashmir(pairs,traits,n,excluded)})
    return out

def main():
    p=argparse.ArgumentParser()
    for field in ("data","eidc-species-csv","output"):
        p.add_argument("--"+field,required=True)
    a=p.parse_args()
    from analyze_matched_foreign_rail_road_v314 import normalize_name
    from source_matched_rail_road_v314 import load
    from source_thermal_affinity_v315 import kashmir_pairs
    original_source=load(Path(a.data).resolve(),normalize_name)
    trait,traitqa=load_global(a.eidc_species_csv)
    pairs=kashmir_pairs(Path(a.data).resolve(),normalize_name)
    r=analyze(original_source,trait,pairs,traitqa)
    target=Path(a.output);target.parent.mkdir(parents=True,exist_ok=True)
    target.write_text(json.dumps(r,ensure_ascii=False,indent=2),encoding="utf8")
    print(json.dumps({"source":r["source_EIDC"],"rails":r["rails"],
                      "comparators":r["comparators"],
                      "Kashmir":r["original_Kashmir_paired"]},ensure_ascii=False),flush=True)

if __name__=="__main__":main()
