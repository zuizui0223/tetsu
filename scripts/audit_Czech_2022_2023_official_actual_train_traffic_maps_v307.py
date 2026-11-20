#!/usr/bin/env python3
"""v307 genuine official Czech infrastructure operator 2022/2023 Q1 TRAINS traffic map source gate.

The official public map "Skutečné počty všech vlaků" is NOT equivalent to
station-specific annual passenger/freight counts. This process only verifies
what exact original PDFs contain. It never invents numerical exposure for
the 39 railway floristic stations surveyed 2020-2022.

All source PDFs stay transient in Actions. Output is audit metadata and
low-resolution first-page visual review only (never original map file).
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import re
import fitz

AUTHORITY="Správa železnic"
PUBLIC_SOURCE_PAGE="https://provoz.spravazeleznic.cz/Portal/ViewArticle.aspx?oid=2104272"
SOURCE_URLS={
    "2022": "https://provoz.spravazeleznic.cz/Portal/Show.aspx?path=%2FData%2FMapy%2Fpocty_vse_22.pdf",
    "2023": "https://provoz.spravazeleznic.cz/Portal/Show.aspx?path=%2FData%2FMapy%2Fpocty_vse_23.pdf",
}
ORIGINAL_PUBLISHED_PERIOD={
    "2022":"1 January–31 March 2022",
    "2023":"1 January–31 March 2023",
}
MAX_SIZE_MB=60

def inspect(path: Path, year: str, render_dir: Path|None=None)->dict:
    if year not in SOURCE_URLS:
        raise ValueError("Only original 2022/2023 public maps supported")
    raw=path.read_bytes()
    if len(raw)<10000 or len(raw)>MAX_SIZE_MB*1024*1024:
        raise ValueError("Official map source PDF unexpected byte size")
    if not raw.startswith(b"%PDF-"):
        raise ValueError("Official Czech traffic map did not return a PDF")
    try:
        d=fitz.open(stream=raw,filetype="pdf")
    except Exception as exc:
        raise ValueError("Official traffic map PDF broken or encrypted") from exc
    if not 1 <= len(d) <= 20 or d.needs_pass:
        raise ValueError("Official traffic map page support invalid")
    parts=[]
    widths=[]
    images=[]
    drawings=[]
    words=[]
    text_digest=[]
    chars=[]
    for page in d:
        txt=page.get_text("text",sort=True)
        chars.append(len(txt))
        words.append(len(page.get_text("words")))
        widths.append([round(float(page.rect.width),2),round(float(page.rect.height),2)])
        images.append(len(page.get_images(full=True)))
        drawings.append(len(page.get_drawings()))
        # Source words only in internal audit as hash: not relicense official map
        text_digest.append(hashlib.sha256(txt.encode()).hexdigest())
        parts.append(txt)
    full="\n".join(parts)
    if render_dir is not None:
        render_dir.mkdir(parents=True,exist_ok=True)
        page=d[0]
        pix=page.get_pixmap(matrix=fitz.Matrix(.75,.75),alpha=False)
        pix.save(render_dir/f"czech_traffic_{year}_firstpage_source_visual_inspection.png")
    d.close()
    num_tokens=re.findall(r"(?<!\w)(?:\d+(?:[.,]\d+)?)(?!\w)",full)
    # PDF text may encode map legend numeric bands in glyph streams. Numeric
    # tokens do NOT establish which station/track segment belongs to a band.
    return {
        "year":year,
        "official_original_source_url":SOURCE_URLS[year],
        "original_official_public_landing_page":PUBLIC_SOURCE_PAGE,
        "official_map_infrastructure_operator":AUTHORITY,
        "published_measure":"actual number of ALL trains, not verified freight-only or passenger-only",
        "published_months":ORIGINAL_PUBLISHED_PERIOD[year],
        "publisher_original_pdf_sha256":hashlib.sha256(raw).hexdigest(),
        "publisher_original_pdf_bytes":len(raw),
        "original_pdf_pages":len(parts),
        "pages_sizes_points":widths,
        "page_text_lengths":chars,
        "page_word_counts":words,
        "page_image_object_counts":images,
        "page_vector_drawings_counts":drawings,
        "page_text_sha256":text_digest,
        "numeric_text_tokens_unassigned_count":len(num_tokens),
        "has_extractable_glyph_text":any(n>0 for n in chars),
        "railway_vegetation_2020_2022_station_specific_passes_identified":False,
        "real_year_specific_botanical_site_linkage_performed":False,
        "freight_passenger_split_proven":False,
        "station_level_numeric_traffic_values_output":0,
        "status":"PUBLIC_ORIGINAL_PDF_SOURCE_VERIFIED_ONLY_NO_STATION_TRAIN_EXPOSURE_ESTIMATED",
    }

def verify_pair(records:dict)->dict:
    if set(records)!={"2022","2023"}:
        raise ValueError("2022 and 2023 must both be present")
    if any(records[k]["year"]!=k for k in records):
        raise ValueError("Traffic map publisher years do not match source file keys")
    if records["2022"]["publisher_original_pdf_sha256"]==records["2023"]["publisher_original_pdf_sha256"]:
        raise ValueError("Distinct public historical train map URLs unexpectedly identical PDF")
    return {
        "status":"PUBLIC_ORIGINAL_CZECH_MAPS_2022_AND_2023_SOURCE_FILES_AVAILABLE",
        "n_source_botanical_station_exposures_linked":0,
        "raw_train_maps_republished":False,
        "original_number_of_genuine_Czech_2020_2022_floristic_stations":39,
        "original_registry_geolocated_stations_previous_source_audit":35,
        "original_publication":"Kutlvašr et al 2024, Preslia DOI 10.23855/preslia.2024.247",
        "year_2022_2023_map_SAMPLES_ONLY_1JAN_TO_31MAR":True,
        "cannot_claim_annual_freight_trains":True,
        "source_maps":records
    }

def main():
    p=argparse.ArgumentParser()
    for year in SOURCE_URLS:p.add_argument("--map-"+year,required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    p.add_argument("--render-dir",type=Path)
    args=p.parse_args()
    r={k:inspect(getattr(args,"map_"+k),k,args.render_dir) for k in SOURCE_URLS}
    report=verify_pair(r)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report,ensure_ascii=False,indent=2)+"\n",encoding="utf8")
    print(json.dumps({
        "status":report["status"],
        "genuine_2022_2023_Q1_only":True,
        "map_2022_sha256":r["2022"]["publisher_original_pdf_sha256"],
        "map_2023_sha256":r["2023"]["publisher_original_pdf_sha256"],
        "pages_by_year":{k:v["original_pdf_pages"] for k,v in r.items()},
        "text_characters_by_year":{k:sum(v["page_text_lengths"]) for k,v in r.items()},
        "n_botanical_station_traffic_values_estimated":0,
    },ensure_ascii=False))

if __name__=="__main__":
    main()
