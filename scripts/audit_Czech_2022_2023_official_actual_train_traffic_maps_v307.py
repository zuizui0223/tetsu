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
import unicodedata
import math
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

# Official botanical station names per corrected 2024 publication; these
# names are only searched as GLYPH LABELS, not linked to train exposure.
ORIGINAL_CZECH_BOTANICAL_STATIONS=[
"Břeclav","Zaječí","Hrušovany u Brna","Brno dolní nádraží",
"Brno-Maloměřice","Skalice nad Svitavou","Svitavy","Česká Třebová",
"Ústí nad Orlicí","Choceň","Pardubice","Přelouč","Tišnov",
"Křižanov","Žďár nad Sázavou","Havlíčkův Brod","Světlá nad Sázavou",
"Čáslav","Kutná Hora","Kolín","Velký Osek","Nymburk",
"Lysá nad Labem","Všetaty","Praha-Horní Měcholupy",
"Praha-Malešice","Praha-Libeň","Praha-Modřany","Praha-Bubny",
"Kralupy nad Vltavou","Mělník","Roudnice nad Labem",
"Lovosice","Velké Žernoseky","Ústí nad Labem západ",
"Ústí nad Labem-Střekov","Děčín","Děčín-Prostřední Žleb",
"Děčín východ"
]

def normalize_source_map_name(value:str)->str:
    v=unicodedata.normalize("NFKD",value)
    v="".join(c for c in v if not unicodedata.combining(c))
    return " ".join(re.findall(r"[a-z0-9]+",v.casefold()))

def audit_botanical_name_glyph_coverage(pdf_pages_text:list[str])->dict:
    """Only determine exact normalized TEXT appearance in official page.

    A matched city name is not a mapped original station, track segment,
    survey year or numerical traffic measurement.
    """
    raw=normalize_source_map_name(" ".join(pdf_pages_text))
    labels={}
    for station in ORIGINAL_CZECH_BOTANICAL_STATIONS:
        tokens=normalize_source_map_name(station)
        labels[station]=bool(re.search(r"(?<![a-z0-9])"+re.escape(tokens)+r"(?![a-z0-9])",raw))
    return {
        "n_original_botanical_station_names_in_corrected_publisher_roster":len(labels),
        "n_exact_normalized_station_label_text_occurrences":sum(labels.values()),
        "matched_text_labels_NOT_exposure_or_station_position":[x for x,v in labels.items() if v],
        "unmatched_not_proof_absent_from_visual_map":[x for x,v in labels.items() if not v],
        "original_station_numeric_pass_counts_still_unidentified":True,
    }




def station_pdf_text_spatial_candidate_audit(d,botanical_station_names:list[str])->dict:
    """Text-only placement QA. Nearby number token NEVER becomes traffic.

    A station at a junction may have multiple services/directions; standalone
    numeric labels could be segment traffic, route marker or unrelated annotation.
    Report only short candidate snippets for manual map-based traceability.
    """
    result={}
    for station in botanical_station_names:
        boxes=[]
        for page_no,page in enumerate(d):
            raw_boxes=page.search_for(station,quads=False)
            num_words=[]
            for w in page.get_text("words"):
                txt=str(w[4]).strip()
                m=re.fullmatch(r"([0-9]{1,3})(?:[.,][0-9]+)?",txt)
                if m and int(m.group(1))<=600:
                    num_words.append((int(m.group(1)),fitz.Rect(w[:4])))
            for b in raw_boxes[:8]:
                cx=(b.x0+b.x1)/2;cy=(b.y0+b.y1)/2
                distances=[]
                for val,rect in num_words:
                    rpx=(rect.x0+rect.x1)/2;rpy=(rect.y0+rect.y1)/2
                    dist=math.hypot(rpx-cx,rpy-cy)
                    if dist<28:
                        distances.append({"number_text_unclassified":val,
                            "distance_pdf_points":round(dist,2),
                            "number_text_bbox_pdf_points":[round(x,2) for x in rect]})
                distances.sort(key=lambda x:x["distance_pdf_points"])
                boxes.append({"page_index":page_no,
                              "label_bbox_pdf_points":[round(x,2) for x in b],
                              "nearest_number_tokens_UNASSIGNED":distances[:8]})
        result[station]=boxes
    return {
      "source_original_botanical_station_names":len(botanical_station_names),
      "station_names_with_EXACT_PDF_search_boxes":sum(bool(x) for x in result.values()),
      "name_bbox_records_not_original_geo_or_actual_train_exposure":result,
      "declared_station_numeric_traffic_exposures":0,
      "all_nearby_numeric_tokens_AMBIGUOUS_not_actual_counts":True
    }

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
    reviewed_stations=station_pdf_text_spatial_candidate_audit(d,ORIGINAL_CZECH_BOTANICAL_STATIONS)
    if render_dir is not None:
        render_dir.mkdir(parents=True,exist_ok=True)
        page=d[0]
        pix=page.get_pixmap(matrix=fitz.Matrix(2.0,2.0),alpha=False)
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
        "published_measure":"daily count of actually running ALL trains on the ninth-busiest day (publisher: 9th decile) of 1 January–31 March, rounded to nearest five for values over 25; not annual freight/passenger counts",
        "ninth_decile_NOT_a_mean_day":True,
        "counts_above_25_rounded_to_multiples_of_5":True,
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
        "botanical_station_NAME_glyph_match_not_numeric_exposure":audit_botanical_name_glyph_coverage(parts),
        "station_text_bbox_and_nearby_NUMBER_CANDIDATES_NOT_EXPOSURES":reviewed_stations,
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
        "official_9th_decile_day_count_not_period_mean":True,
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
