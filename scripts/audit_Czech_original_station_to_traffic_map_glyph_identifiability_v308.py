#!/usr/bin/env python3
"""v308 original Czech botanical locality × official train-map glyph audit.

No station is ever assigned a train frequency by nearest number. This pipeline
tests matching *source identity* and geographic map text positions, not
ecological traffic effects. Existing 2022/2023 original operator-map metrics
are 9th decile of daily ALL-train counts for first quarter; NOT annual freight.
"""
from __future__ import annotations
import argparse
from collections import Counter
from pathlib import Path
import hashlib
import json
import math
import re
import fitz
from openpyxl import load_workbook
from audit_Czech_2022_2023_official_actual_train_traffic_maps_v307 import (
    ORIGINAL_CZECH_BOTANICAL_STATIONS, normalize_source_map_name,
    SOURCE_URLS as OPERATOR_PDF_URLS
)
from analyze_Czech_2024_39station_functional_identity_v143 import (
    read_original_czech_publisher,
)
from audit_Czech_2024_S2_to_SR70_original_station_points_v153 import (
    original_s2_names,sr70_records,match_czech_stations,
)
from audit_Czech_2024_real_railway_appendices_v107 import S1_URL,S2_URL

YEAR_SHA={
    "2022":"c4d2436a7e1950bc6c247926b12dbba16649491221c94fc7288e234c1c3830ee",
    "2023":"4cbd4349814ac636bda74d3389810d6100ecc9e12daef237a12fa8b62eb25b57",
}
MIN_XLSX_STATIONS=39

def read_original_botanical_station_support(s1:Path,s2:Path,sr70:Path)->tuple[list[dict],dict]:
    flora,meta,source_qa=read_original_czech_publisher(s1,s2)
    stations=original_s2_names(s2)
    if len(stations)!=MIN_XLSX_STATIONS or len(flora)!=MIN_XLSX_STATIONS:
        raise ValueError("Original botanical study has no authentic 39 stations")
    pub_names={normalize_source_map_name(x['original_source_s2_col1_label']) for x in stations}
    previous={normalize_source_map_name(x) for x in ORIGINAL_CZECH_BOTANICAL_STATIONS}
    if pub_names!=previous:
        raise ValueError("The earlier v307 hard-coded botanical station list is NOT an exact set of original corrected S2 names")
    pointrows,pointqa=sr70_records(sr70)
    matched=match_czech_stations(stations,pointrows)
    support=Counter(x['match_status'] for x in matched)
    wb=load_workbook(s2,read_only=True,data_only=True)
    try:
        ws=wb.worksheets[0]
        head=[['' if z is None else str(z)[:120] for z in row[:20]]
              for row in list(ws.values)[:5]]
        source_cols=ws.max_column
    finally:wb.close()
    # No raw geographic coordinates, station-by-taxon matrix, or exact
    # individual lines are serialized; source names already public in Fig.1.
    return matched,{
        "original_Preslia_DOI":"10.23855/preslia.2024.247",
        "original_Preslia_S1_url":S1_URL,
        "original_Preslia_S2_corrected_url":S2_URL,
        "source_S1_S2_verified_biological_incidence_markers":source_qa['original_39_station_contemporary_exact_presence_marks'],
        "n_original_source_current_station_taxa":source_qa['n_original_recent_real_taxa_in_public_S1_station_rec_columns'],
        "n_author_station_rows":len(stations),
        "n_independently_matched_2026_sr70_centroids":sum(x['matched_czech_2026_source_station_centroid_WGS84'] is not None for x in matched),
        "station_name_list_matches_original_2024_corrected_S2_exact_normalized_set":True,
        "n_original_railway_stretches":len({x['source_section'] for x in stations}),
        "source_s2_header_preview_not_original_station_individual_data":head,
        "source_S2_max_columns":source_cols,
        "source_original_stations_selected_partly_by_freight_passenger_volume_not_random":True,
        "station_registry_derived_2026_centroid_not_original_field_plot_coordinates":True,
        "station_match_status":dict(support),
        "third_party_registry_source_qa":pointqa,
    }

def check_pdf_original(p:Path,year:str)->fitz.Document:
    if year not in YEAR_SHA:raise ValueError("year outside original locked 2022/2023")
    if hashlib.sha256(p.read_bytes()).hexdigest()!=YEAR_SHA[year]:
        raise ValueError("Original operator real-train PDF SHA changed; no unverified map allowed")
    d=fitz.open(str(p))
    if len(d)!=1:raise ValueError("Expected exactly one operator diagram page")
    return d

def number_glyphs(page)->list[dict]:
    rows=[]
    # Text-word bounding boxes and map-vector geometry use an unrotated
    # coordinate system. Transform each to the page's rendered rotation to
    # inspect crops consistently; do not match to station traffic.
    for word in page.get_text("words"):
        label=str(word[4]).strip()
        if not re.fullmatch(r"\d{1,3}",label):
            continue
        n=int(label)
        if not 1<=n<=600:continue
        bbox=fitz.Rect(word[:4])
        rotated=bbox*page.rotation_matrix
        rows.append({"number":n,
                     "bbox_unrotated":bbox,
                     "bbox_rotated":rotated})
    return rows

def _center(rect):
    return ((rect.x0+rect.x1)/2,(rect.y0+rect.y1)/2)

def _distance(a,b):
    x,y=_center(a);xx,yy=_center(b)
    return math.hypot(xx-x,yy-y)

def glyph_pairs_for_station(page,station_name:str,numbers:list[dict],radius_pt=10.0):
    labels=page.search_for(station_name)
    candidate_groups=[]
    for b in labels:
        rot=b*page.rotation_matrix
        # Reject impossible rotated coords, otherwise overlaying onto
        # pixmap would create a false line attribution.
        if (rot.x0<-1 or rot.y0<-1 or rot.x1>page.rect.width+1 or
            rot.y1>page.rect.height+1):
            raise ValueError("PDF rotation transform inconsistent with page.rect")
        candidates=[{"n":w["number"],"distance_points":round(_distance(b,w["bbox_unrotated"]),3)}
                    for w in numbers if _distance(b,w["bbox_unrotated"])<=radius_pt]
        candidate_groups.append({"label_bbox_rotated_points":[round(z,2) for z in rot],
                                 "label_bbox_unrotated_points":[round(z,2) for z in b],
                                 "nearby_numeric_glyphs_UNCLASSIFIED":sorted(candidates,key=lambda z:z["distance_points"])})
    return candidate_groups

def inspect_map(original_pdf:Path,year:str,matched:list[dict])->dict:
    d=check_pdf_original(original_pdf,year)
    p=d[0]
    words=number_glyphs(p)
    station_details={}
    for source in matched:
        name=source["original_source_s2_col1_label"]
        station_details[source["original_locality_code"]]={
            "original_publisher_station_name":name,
            "original_study_stretch":source["source_section"],
            "SR70_place_match_status":source["match_status"],
            "num_station_word_search_boxes":len((q:=glyph_pairs_for_station(p,name,words))),
            "PDF_word_search_boxes_rotated_correctly_and_nearby_numeric_GLYPHS_ONLY":q,
        }
    counts=Counter()
    single_unambiguous_proximity=[]
    for site,desc in station_details.items():
        if desc["num_station_word_search_boxes"]>0:counts["station_label_detected"]+=1
        if desc["num_station_word_search_boxes"]==1:
            c=desc["PDF_word_search_boxes_rotated_correctly_and_nearby_numeric_GLYPHS_ONLY"][0]["nearby_numeric_glyphs_UNCLASSIFIED"]
            if len(c)==1:
                single_unambiguous_proximity.append({
                    "code":site,"name":desc["original_publisher_station_name"],
                    "single_numeric_GLYPH_not_verified_train_count":c[0]["n"]})
    result={
        "year":year,"original_pdf_sha256":YEAR_SHA[year],
        "operator_source_url":OPERATOR_PDF_URLS[year],
        "operator_traffic_statistic":"ninth decile (~90th percentile) of ALL trains per day during Jan-Mar; counts over 25 rounded to 5; excludes passenger/freight split",
        "pdf_page_rotation_deg":p.rotation,
        "pdf_rendered_page_rect_points":[round(p.rect.width,2),round(p.rect.height,2)],
        "pdf_unrotated_page_rect_points":[round(p.cropbox.width,2),round(p.cropbox.height,2)],
        "n_unclassified_numeric_text_glyphs":len(words),
        "n_station_labels_found_by_pdf_glyph_bbox":counts["station_label_detected"],
        "one_box_one_number_within_10pt_ONLY_not_rail_segment_assignment":single_unambiguous_proximity,
        "n_one_box_one_number_within_10pt_only":len(single_unambiguous_proximity),
        "station_line_and_number_assignment_not_proven":True,
        "n_approved_botanical_station_traffic_exposures":0,
        "original_study_station_geo_box_caution":"Station text may recur in detailed city insets and describe towns or junctions; no single route taken as exposure",
        "source_station_text_and_nearby_unclassified_map_glyphs":station_details,
    }
    d.close()
    return result

def main():
    p=argparse.ArgumentParser()
    for k in ("s1","s2","sr70","map-2022","map-2023","output"):
        p.add_argument("--"+k,required=True,type=Path)
    a=p.parse_args()
    stations,original=read_original_botanical_station_support(a.s1,a.s2,a.sr70)
    years={
        "2022":inspect_map(getattr(a,"map_2022"),"2022",stations),
        "2023":inspect_map(getattr(a,"map_2023"),"2023",stations)}
    x=years["2022"]["one_box_one_number_within_10pt_ONLY_not_rail_segment_assignment"]
    y=years["2023"]["one_box_one_number_within_10pt_ONLY_not_rail_segment_assignment"]
    px={z['name']:z['single_numeric_GLYPH_not_verified_train_count'] for z in x}
    py={z['name']:z['single_numeric_GLYPH_not_verified_train_count'] for z in y}
    both={z:{"map_2022_numeric_glyph":px[z],"map_2023_numeric_glyph":py[z]}
          for z in sorted(px.keys()&py.keys())}
    # The "Lysá nad Labem = 1/1" counterexample highlights how a
    # naive closest-number rule could capture map metadata/branch details
    # rather than through-running traffic at a busy junction.
    out={
        "version":"v308","status":"ORIGINAL_CZECH_PUBLISHER_STUDY_AND_OPERATOR_MAPS_MATCHED_AT_TEXT_ONLY",
        "one_botanical_source_publication_not_39_independent_studies":True,
        "original_botanical_publisher_and_station_geometry":original,
        "operator_maps":years,
        "cross_year_uniquely_near_num_glyphs_NOT_verifiable_exposures":both,
        "n_original_floristic_stations":39,
        "n_station_level_train_counts_approved":0,
        "no_ecological_rail_traffic_slope_fitted":True,
        "do_NOT_use_a_station_numeric_GLYPH_as_train_exposure_without_same_track_physical_line_trace":True,
    }
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(out,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    print(json.dumps({
        "status":out["status"],
        "original_botanical_station_names_proved":39,
        "original_2026_centroid_matches":original["n_independently_matched_2026_sr70_centroids"],
        "2022_station_labels":years["2022"]["n_station_labels_found_by_pdf_glyph_bbox"],
        "2023_station_labels":years["2023"]["n_station_labels_found_by_pdf_glyph_bbox"],
        "single_number_within_10pt_both_years":len(both),
        "actual_station_traffic_verified":0,
        "example_naive_numeric_GLYPH_nearest_might_be_wrong":both.get("Lysá nad Labem"),
    },ensure_ascii=False))
if __name__=="__main__":main()
