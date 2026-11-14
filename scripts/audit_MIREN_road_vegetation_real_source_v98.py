#!/usr/bin/env python3
"""v9.8 authentic MIREN roadside standardized plot source audit.

Original: Seipel/Haider/MIREN consortium, Zenodo record 15749536,
v3, 2007–2023, direct source CSV (published 2025). Original file
MD5 from Zenodo 15749536 metadata: b640195cc44f252ccb6ddb28eb6fff13.

This is a completely INDEPENDENT measured-road-vegetation dataset:
25 mountain regions / 3364 vegetation plots / 6854 vascular spp
in the 2026 Haider et al. metadata, including original regional
native/non-native determinations and road-distance layout.

Audit raw schema, row conservation, source column missingness and
original metadata before assuming how road-distance/native status
fields are actually coded. NEVER invent mowing exposure or train use.
No railway seed dispersal/adaptation effects estimated here.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,requests
from collections import Counter
from datetime import datetime,timezone
from pathlib import Path

ZENODO_RECORD="https://zenodo.org/records/15749536"
URL="https://zenodo.org/records/15749536/files/MIREN_road_survey_%282007-2023%29_zen_v3_002507.csv?download=1"
EXPECTED_MD5="b640195cc44f252ccb6ddb28eb6fff13"
MIN_BYTES=15_000_000
MAX_BYTES=25_000_000
MAX_SOURCE_ROWS=400000
MAX_COLUMNS=100
SOURCE_ROW_SAMPLE=5

def acquire(url,dest):
    dest=Path(dest)
    dest.parent.mkdir(parents=True,exist_ok=True)
    with requests.get(url,timeout=(20,180),stream=True,
                      headers={"User-Agent":"tetsu-independent-MIREN-road-plot-2026-v98"}) as r:
        r.raise_for_status()
        with open(dest,"wb") as f:
            for chunk in r.iter_content(chunk_size=512*1024):
                if chunk:f.write(chunk)
    if not MIN_BYTES<dest.stat().st_size<MAX_BYTES:
        raise ValueError("Original independent MIREN public source CSV size does not fit metadata")
    md5=hashlib.md5(dest.read_bytes()).hexdigest()
    if md5!=EXPECTED_MD5:
        raise ValueError("MIREN raw 2007–2023 original source md5 changed since Zenodo v3; stop before ecological analysis")
    return {"source_original_zenodo_md5":md5,
            "source_raw_bytes":dest.stat().st_size,
            "source_original_Zenodo_DOI":"10.5281/zenodo.15749536"}

def inspect_csv(path):
    path=Path(path)
    with path.open("r",encoding="utf-8-sig",newline="") as f:
        first=f.read(15000)
        f.seek(0)
        dialect=csv.Sniffer().sniff(first,delimiters=",;\t")
        reader=csv.DictReader(f,dialect=dialect)
        if reader.fieldnames is None or len(reader.fieldnames)<7 or len(reader.fieldnames)>MAX_COLUMNS:
            raise ValueError("MIREN actual source CSV column schema invalid")
        headers=[h.strip() if h is not None else "" for h in reader.fieldnames]
        if any(not h for h in headers) or len(set(headers))!=len(headers):
            raise ValueError("Original MIREN 2025 CSV has duplicate or unnamed columns")
        miss=Counter(); non=Counter(); uniques={h:set() for h in headers}
        rows=0
        values_sample=[]
        for row in reader:
            if None in row:raise ValueError("Malformed CSV row width/source delimiter")
            rows+=1
            if rows>MAX_SOURCE_ROWS:raise ValueError("MIREN source row count above predeclared 400k")
            for key in headers:
                v=(row[key] or "").strip()
                if v.lower() in ("","na","nan","null"):
                    miss[key]+=1
                else:
                    non[key]+=1
                    if len(uniques[key])<500:uniques[key].add(v[:90])
            if len(values_sample)<SOURCE_ROW_SAMPLE:
                values_sample.append({h:(row[h] or "").strip()[:80]
                    for h in headers
                    if not any(z in h.lower() for z in ("latitude","longitude","lat","lon","xcoord","ycoord"))})
        if rows<50000:raise ValueError("Not original global vegetation source: <50000 rows")
        # Do not assert raw counts equal the 2026 article's deduplicated
        # 170k records/plot counts: observational row and plot grains differ.
        return {
            "n_original_csv_rows":rows,
            "n_original_csv_columns":len(headers),
            "raw_source_field_names":headers,
            "source_csv_delimiter":dialect.delimiter,
            "n_null_by_column":dict(miss),
            "n_nonempty_by_column":dict(non),
            "low_cardinality_field_values_sample":{
                h:sorted(uniques[h])[:65] for h in headers
                if len(uniques[h])<=60},
            "raw_source_first_five_rows_without_original_latlong":values_sample,
            "no_claim_original_rows_equal_species_or_plant_plot_count":True,
        }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source-csv",required=True)
    p.add_argument("--audit-output",required=True)
    a=p.parse_args()
    source=acquire(URL,a.source_csv)
    audit=inspect_csv(a.source_csv)
    out={
        "complete":True,"UTC":datetime.now(timezone.utc).isoformat(),
        "source_zenodo_record":ZENODO_RECORD,
        "source_original_publisher":"Seipel T, Haider S and MIREN consortium (2025), MIREN road survey v3",
        "independent_published_data_road_verge_not_train_rail_verge":True,
        "source_2026_paper_doi":"10.1111/geb.70284",
        "source_original_binary_metadata":source,
        "actual_CSV_source_schema":audit,
        "biological_inference_not_yet_fitted":True,
        "NO_railway_wind_mowing_seed_transport_or_adaptation_inferred":True,
    }
    dest=Path(a.audit_output)
    dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps(out,ensure_ascii=False),flush=True)

if __name__=="__main__":
    main()
