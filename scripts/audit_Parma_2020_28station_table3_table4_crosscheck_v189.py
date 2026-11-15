#!/usr/bin/env python3
"""v18.9 independent within-publisher source QA for Parma 2020.

Crossvalidate the fully reconstructed original 438x28 station-code matrix
against TWO DISTINCT summary tables in SAME publisher original PDF:
Table 3 (all 28 station taxon counts) and Table 4 (27 most frequent taxa).
This checks reconstruction but is NOT a new independent field experiment.
Publisher article source S4.1 has 68 bracket count inconsistencies, which
must never be used to fabricate absences/presences.
"""
from __future__ import annotations
import argparse,json,re,sys
from collections import Counter
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
from extract_Parma_2020_28station_species_matrix_v187 import parse_pdf,station_stats
from audit_Parma_2020_28rail_station_flora_original_v186 import STATION_CODES

# Publisher original Table 3, p.179, completely independent of the species
# list text within 4.1 (but SAME field source, hence statistical study n=1).
TABLE3={
"SAL":157,"BVT":151,"PAR":145,"SEL":131,"BER":130,"NOC":126,"SPT":124,
"OST":123,"FID":117,"BUS":113,"MED":112,"VIC":110,"VAL":107,"CMA":100,
"FEL":99,"SOL":98,"ROC":97,"BRG":94,"FOR":93,"CLO":83,"MEZ":83,
"OZZ":80,"CLL":77,"CGU":74,"SOR":74,"VAI":71,"GAI":63,"CHI":32}

# Publisher original Table 4, p.182, 27 listed taxa with ≥20 stations.
# Distinct from the 4.1 bracket frequencies, many of which differ from actual
# station-code enumerations. For ordinary cases, first two taxon words fix
# the lookup (full rank names are retained in the underlying CSV).
TABLE4={
"Galium aparine":27, "Cerastium glomeratum":26,
"Equisetum ramosissimum":26,"Geranium purpureum":26,
"Anisantha sterilis":24,"Arenaria serpyllifolia":24,
"Chenopodium album":24,"Lactuca sativa":24,
"Convolvulus arvensis":23,"Cynodon dactylon":23,
"Hordeum murinum":23,"Setaria italica":23,
"Stellaria media":23,"Veronica arvensis":23,
"Capsella bursa-pastoris":22,"Crepis sancta":22,
"Erigeron canadensis":22,"Papaver rhoeas":22,
"Senecio vulgaris":22,"Sonchus oleraceus":22,
"Poa annua":21,"Portulaca oleracea":21,
"Sonchus asper":21,"Avena barbata":20,
"Clematis vitalba":20,"Robinia pseudoacacia":20,
"Taraxacum sect. Taraxacum":20
}
assert len(TABLE3)==28 and set(TABLE3)==set(STATION_CODES)
assert len(TABLE4)==27

def taxon_key(name):
    if name.startswith("Taraxacum ") and "sect. Taraxacum" in name:
        return "Taraxacum sect. Taraxacum"
    m=re.match(r"^([A-Z][a-zA-ZÀ-ž-]+)\s+([a-z][a-zA-ZÀ-ž-]+)",name)
    return " ".join(m.groups()) if m else None

def source_check(rows):
    if len(rows)!=438 or len(set(r["taxon_original"] for r in rows))!=438:
        raise ValueError("Expected 438 actual source taxa")
    site_richness,_,_,_=station_stats(rows)
    source_by_key={}
    for r in rows:
        k=taxon_key(r["taxon_original"])
        if k in TABLE4:source_by_key.setdefault(k,[]).append(r)
    t4={}
    for key,n in TABLE4.items():
        hits=source_by_key.get(key,[])
        if len(hits)!=1:
            t4[key]={"publisher_expected":n,
                     "observed_count":None,
                     "source_matching_rows":len(hits),"pass":False}
        else:
            t4[key]={"publisher_expected":n,
                     "observed_count":len(hits[0]["station_codes"]),
                     "original_bracket_count":hits[0]["printed_station_total"],
                     "pass":len(hits[0]["station_codes"])==n,
                     "source_matching_rows":1}
    t3={k:{"publisher_expected":v,"observed_count":site_richness[k],
           "pass":site_richness[k]==v} for k,v in TABLE3.items()}
    out={
        "source_n_438x28":True,
        "source_table3_all_28_stations":t3,
        "source_table4_27_frequent_taxa":t4,
        "source_table3_passed":sum(r["pass"] for r in t3.values()),
        "source_table4_passed":sum(r["pass"] for r in t4.values()),
        "source_bracket_frequency_mismatches":sum(
             r["printed_station_total"]!=len(r["station_codes"]) for r in rows),
        "source_raw_station_codes_NOT_modified_to_match_printed_bracket":True,
        "audit_passed_without_fake_site_species":all(r["pass"] for r in t3.values())
                 and all(r["pass"] for r in t4.values()),
        "same_original_publisher_not_independent_field_replication":True}
    return out

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--original-pdf",required=True)
    p.add_argument("--output",required=True)
    a=p.parse_args()
    rows,err,diag=parse_pdf(a.original_pdf)
    if err:raise ValueError("PDF 4.1 species list parse errors")
    out=source_check(rows)
    out["source_pdf_sha256"]=diag["source_pdf_sha256"]
    f=Path(a.output);f.parent.mkdir(parents=True,exist_ok=True)
    f.write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding="utf8")
    print(json.dumps({"source_sha256":out["source_pdf_sha256"],
                      "source_table3_passed":out["source_table3_passed"],
                      "source_table4_passed":out["source_table4_passed"],
                      "source_bracket_frequency_mismatches":out["source_bracket_frequency_mismatches"],
                      "failed_table4":{k:v for k,v in out["source_table4_27_frequent_taxa"].items() if not v["pass"]},
                      "failed_table3":{k:v for k,v in out["source_table3_all_28_stations"].items() if not v["pass"]},
                      "audit_passed_without_fake_site_species":out["audit_passed_without_fake_site_species"]},
                      ensure_ascii=False,indent=2),flush=True)
    if not out["audit_passed_without_fake_site_species"]:
        raise ValueError("Source reconstruction did not meet original publisher independent summary tables")

if __name__=="__main__":main()
