#!/usr/bin/env python3
"""v10.2 official 2026 Slovenian railway vegetation ORIGINAL supplementary
XLSX source fetch and ZIP integrity/schema-size audit.

Šabić A, Grašič M, Rovanšek A (2026), Plant Biosystems
'Vegetation patterns, functional traits and invasion dynamics across
railway microhabitats in Slovenia', DOI 10.1007/s44473-026-00155-z.
Official article has one CC-BY 4.0 XLSX supplement named
44473_2026_155_MOESM1_ESM.xlsx, about 392 KB, with plant traits and
sampling location/abundance/frequency.

Study: 92 authentic field sampling locations across 15 Slovenian
railway line segments, summer 2022–2023, 312 recorded plant taxa,
four actual railway microhabitats: sleeper zone SZ, ballast shoulder BS,
ballast toe BT, and intertrack area IA. File may NOT contain original
site×zone×taxon raw observations; inspect actual sheets before asserting
a species/site source response. No train seed carriage observed.

The official paper explicitly says the original sites were chosen by
visible developed vegetation, NOT a randomized national rail sample.
An XLSX published supplementary trait/abundance summary is not a
global railway plot survey or a counterfactual matched MIREN road plot.
"""
from __future__ import annotations
import argparse,hashlib,io,json,requests,zipfile,time
from datetime import datetime,timezone
from pathlib import Path

ARTICLE="https://link.springer.com/article/10.1007/s44473-026-00155-z"
SOURCE_URLS=(
  "https://media.springernature.com/original/springer-static/esm/art%3A10.1007%2Fs44473-026-00155-z/MediaObjects/44473_2026_155_MOESM1_ESM.xlsx",
  "https://static-content.springer.com/esm/art%3A10.1007%2Fs44473-026-00155-z/MediaObjects/44473_2026_155_MOESM1_ESM.xlsx",
)
MIN_BYTES=150_000
MAX_BYTES=1_000_000
MIN_ZIP_MEMBERS=7
ALLOWED_TYPE=(
    "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
    "application/octet-stream","binary/octet-stream",
)

def fetch_source(session):
    failures=[]
    for url in SOURCE_URLS:
        try:
            r=session.get(url,timeout=(20,160),
                headers={"User-Agent":"Mozilla/5.0 (compatible; tetsu-macro-source-audit/1.0)",
                         "Referer":"https://link.springer.com/"})
            r.raise_for_status()
            raw=r.content
            if not MIN_BYTES<len(raw)<MAX_BYTES or not zipfile.is_zipfile(io.BytesIO(raw)):
                raise ValueError("Official supplementary file is not the expected XLSX ZIP size/type")
            return url,raw
        except (requests.RequestException,ValueError) as e:
            failures.append({"source_host":url.split("/")[2],
                             "error_type":type(e).__name__})
            time.sleep(.3)
    raise RuntimeError("Official 2026 Springer supplementary XLSX download inaccessible: "+str(failures))

def audit_xlsx_container(raw):
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        if z.testzip() is not None:
            raise ValueError("Published Slovenian XLSX source ZIP CRC corrupted")
        names=z.namelist()
        if len(names)<MIN_ZIP_MEMBERS or "xl/workbook.xml" not in names or (
            "[Content_Types].xml" not in names):
            raise ValueError("Expected source supplementary XLSX workbook internal files absent")
        sheets=[x for x in names if x.startswith("xl/worksheets/sheet") and x.endswith(".xml")]
        if not sheets:raise ValueError("Original published XLSX source has no worksheet data")
        if len(sheets)>80:raise ValueError("Unexpected too many original supplementary sheet tables")
        return {
            "n_original_xlsx_zip_members":len(names),
            "n_original_source_worksheet_xml_files":len(sheets),
            "worksheet_XML_files":sheets,
            "original_shared_strings_worksheet_file_present":"xl/sharedStrings.xml" in names,
            "real_supplement_workbook_SOURCE_ONLY_not_automatically_92_site_quadrats":True,
        }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--xlsx-output",required=True)
    p.add_argument("--audit-output",required=True)
    a=p.parse_args()
    with requests.Session() as s:
        url,raw=fetch_source(s)
    provenance=audit_xlsx_container(raw)
    p=Path(a.xlsx_output);p.parent.mkdir(parents=True,exist_ok=True)
    p.write_bytes(raw)
    summary={
      "complete":True,
      "source_UTC":datetime.now(timezone.utc).isoformat(),
      "publication":"Sabic A, Grasic M, Rovans ek A (2026), Plant Biosystems",
      "original_article_DOI":"10.1007/s44473-026-00155-z",
      "original_article_url":ARTICLE,
      "original_supplement_url":url,
      "original_supplement_sha256":hashlib.sha256(raw).hexdigest(),
      "n_source_original_bytes":len(raw),
      "zip_integrity":provenance,
      "source_publication_field_sample_sites":92,
      "source_publication_rail_segments":15,
      "source_publication_recorded_plant_taxa":312,
      "source_publication_sample_years":"2022–2023",
      "no_assumed_original_plot_x_zone_species_observations_until_real_sheets_checked":True,
      "NO_actual_rail_seed_movement_mowing_or_evolutionary_adaptation":True,
    }
    q=Path(a.audit_output);q.parent.mkdir(parents=True,exist_ok=True)
    q.write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps(summary,ensure_ascii=False),flush=True)

if __name__=="__main__":main()
