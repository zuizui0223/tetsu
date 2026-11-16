#!/usr/bin/env python3
"""v15.3 verify exact 39 Czech original botanical-station identities against
2026 third-party extract of the official Czech SR70 railway point registry.

Primary plant source: Kutlvašr et al (2024), DOI 10.23855/preslia.2024.247;
Preslia author-publisher S1 & corrected S2 SHA-checked in v14.3.
Independent spatial candidate source: the Czech railway official SR70 codes as
distributed by Sap1k/jrunify-ext-geodata /rail/SR70.csv, 2026-08-15 effective
official sheet. This is a *third-party derivative*, NOT primary author GPS.
Git blob SHA-1 pinned to known source; payload SHA256 recorded separately.

Conservative matching:
- exact accent-insensitive, punctuation-insensitive FULL names only;
- distinct registry rows at nearly the SAME place collapse as one station;
- non-identical positions for a name are ambiguous, not guessed;
- no fuzzy string match, suffix deletion, code inference or town centre fill.
- source S2 first column may carry names; if not, leave unmatched. Source
  Fig.1 explicitly lists 39 original station names, but an ordinal->code
  join may only be used after checking original S2 row semantics.
No fitting of rail vs road effects before source geography is authenticated.
"""
from __future__ import annotations
import argparse,csv,hashlib,io,json,math,re,unicodedata
from pathlib import Path
from collections import Counter
from analyze_Czech_2024_39station_functional_identity_v143 import read_original_czech_publisher
from analyze_Slovenia_rail_spatial_matched_taxon_identity_v142 import haversine_km

SR70_REPO="https://github.com/Sap1k/jrunify-ext-geodata"
SR70_URL="https://raw.githubusercontent.com/Sap1k/jrunify-ext-geodata/main/rail/SR70.csv"
SR70_GIT_BLOB_SHA1="df1f4b41a6292071b58c4a9a862bd535914706c3"
SR70_SOURCE_DATE="2026-08-15"
DOI="10.23855/preslia.2024.247"
MAX_SAME_PLACE_M=150.0


def normalized_name(v):
    if v is None:return ""
    s=unicodedata.normalize("NFKD",str(v).casefold())
    s="".join(x for x in s if not unicodedata.combining(x))
    # Czech letters ěčďňřšťůž -> e c d n r s t u z via Unicode.
    s=s.replace("–","-").replace("—","-")
    s=re.sub(r"[^a-z0-9]+"," ",s).strip()
    return re.sub(r"\s+"," ",s)


def original_s2_names(path):
    # The independent source-check function validates BOTH authentic XLSX
    # sizes and their SHA before this reader accesses source sheets.
    from openpyxl import load_workbook
    wb=load_workbook(path,read_only=True,data_only=True)
    try:
        rows=list(next(iter(wb.worksheets)).values)
    finally:wb.close()
    originals=[]
    for index,row in enumerate(rows[5:],start=6):
        if len(row)<3 or row[1] is None or not str(row[1]).strip():continue
        code=str(row[1]).strip()
        title=str(row[0]).strip() if row[0] is not None else ""
        section=str(row[2]).strip()
        originals.append({"source_s2_sheet_row":index,
                          "original_locality_code":code,
                          "original_source_s2_col1_label":title,
                          "source_section":section})
    if len(originals)!=39 or len(set(x["original_locality_code"] for x in originals))!=39:
        raise ValueError("Original publisher S2 must have 39 unique station codes")
    if len(set(x["source_section"] for x in originals))!=5:
        raise ValueError("Original publisher S2 5 original stretches drifted")
    return originals


def sr70_records(path):
    raw=Path(path).read_bytes()
    hash_git=hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()
    if hash_git!=SR70_GIT_BLOB_SHA1:
        raise ValueError("Third-party published SR70 CSV Git blob SHA drift (not original 2026 official file)")
    try: txt=raw.decode("utf-8-sig")
    except UnicodeDecodeError as e:raise ValueError("Third-party SR70 UTF8 source invalid") from e
    original=[]
    n_source_missing_coordinate_0_0=0
    for line,row in enumerate(csv.reader(io.StringIO(txt)),start=1):
        if not any(v.strip() for v in row):continue
        if len(row)!=4:raise ValueError("SR70 original derived rows not id,name,lat,lon")
        code,name,lat,lon=(x.strip() for x in row)
        if not (code.isdigit() and len(code)>=5 and name):raise ValueError("Malformed SR70 code/name")
        lat=float(lat);lon=float(lon)
        if lat==0 and lon==0:
            # Genuine SR70 derivative marks missing-source geometry with 0,0.
            # Exclude it from matching, preserving the number of omissions.
            n_source_missing_coordinate_0_0+=1
            continue
        if not (48<=lat<=52 and 12<=lon<=20 and math.isfinite(lat) and math.isfinite(lon)):
            raise ValueError("Source SR70 non-missing point outside Czech bounds")
        original.append({"sr70_code":code,"name":name,"lat":lat,"lon":lon,
                         "match_key":normalized_name(name)})
    return original,{
      "original_third_party_regenerated_rail_registry_url":SR70_URL,
      "third_party_2026_effective_version_according_to_upstream_README":SR70_SOURCE_DATE,
      "source_blob_git_sha1":hash_git,
      "downloaded_csv_sha256":hashlib.sha256(raw).hexdigest(),
      "download_bytes":len(raw),
      "n_sr70_valid_coordinate_rows":len(original),
      "n_sr70_source_rows_with_explicit_zero_zero_missing_coordinate_excluded":n_source_missing_coordinate_0_0,
      "not_original_station_botanical_sampling_coordinates":True,
      "third_party_derived_FROM_official_SR70_NOT_official_2022_original":True,
    }


def match_czech_stations(originals,srrows,max_same_place_m=MAX_SAME_PLACE_M):
    bykey={}
    for p in srrows:bykey.setdefault(p["match_key"],[]).append(p)
    result=[]
    for r in originals:
        key=normalized_name(r["original_source_s2_col1_label"])
        candidates=bykey.get(key,[]) if key else []
        groups=[]
        for c in candidates:
            assigned=False
            for group in groups:
                # A candidate can enter a group only if near EVERY member.
                if all(haversine_km(c["lat"],c["lon"],d["lat"],d["lon"])*1000
                       <=max_same_place_m for d in group):
                    group.append(c);assigned=True;break
            if not assigned:groups.append([c])
        if not key or not candidates:status="no_exact_full_name_match"
        elif len(groups)>1:status="ambiguous_different_station_locations"
        else:status="unique_place_by_identical_full_station_name"
        loc=None
        if status=="unique_place_by_identical_full_station_name":
            group=groups[0]
            loc={"lat":sum(x["lat"] for x in group)/len(group),
                 "lon":sum(x["lon"] for x in group)/len(group),
                 "n_sr70_registry_codes_collapsed_same_place":len(group),
                 "source_sr70_registry_codes":[x["sr70_code"] for x in group],
                 "source_sr70_full_names":[x["name"] for x in group]}
        result.append({**r,
            "original_source_full_name_normalized":key,
            "source_registry_full_name_candidates":len(candidates),
            "distinct_150m_place_clusters":len(groups),
            "match_status":status,
            "matched_czech_2026_source_station_centroid_WGS84":loc,
        })
    return result


def run(s1,s2,sr70_path):
    _,meta,source=read_original_czech_publisher(s1,s2)
    originals=original_s2_names(s2)
    if {p["original_locality_code"] for p in originals}!=set(meta):
        raise ValueError("Actual Czech source S2 original code mismatch with verified S1")
    points,provenance=sr70_records(sr70_path)
    matches=match_czech_stations(originals,points)
    statuses=Counter(x["match_status"] for x in matches)
    exact=[x for x in matches if x["matched_czech_2026_source_station_centroid_WGS84"] is not None]
    return {
       "version":"v15.3_original_Czech_2024_botanical_station_vs_2026_SR70_location_audit",
       "original_Czech_study_DOI":DOI,
       "publisher_S1_and_corrected_S2_authentic_provenance":source,
       "original_publisher_S2_station_count":len(originals),
       "original_publisher_source_stretches":5,
       "third_party_official_rail_registry_derivative":provenance,
       "exact_name_conservative_station_centroid_status_counts":dict(statuses),
       "n_unambiguously_matched_station_centroids":len(exact),
       "n_unmatched_or_ambiguous_station_centroids":len(originals)-len(exact),
       "all_39_station_original_source_identity_and_candidate_provenance":matches,
       "do_NOT_treat_registry_station_centroid_as_source_field_plot_GPS":True,
       "station_name_matches_may_be_spatially_off_in_wide_station_yards":True,
       "source_2026_vs_botanical_2020_2022_temporal_mismatch":True,
       "strictly_NO_movement_network_or_seed_causality_claimed":True,
    }


def main():
    p=argparse.ArgumentParser()
    for n in ("s1","s2","sr70","output"):p.add_argument("--"+n,required=True)
    a=p.parse_args()
    out=run(a.s1,a.s2,a.sr70)
    dest=Path(a.output);dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(out,indent=2,ensure_ascii=False),encoding="utf-8")
    print(json.dumps({
        "version":out["version"],
        "n_publisher_stations":out["original_publisher_S2_station_count"],
        "registry_SHA":out["third_party_official_rail_registry_derivative"]["downloaded_csv_sha256"],
        "status_counts":out["exact_name_conservative_station_centroid_status_counts"],
        "all_names_station_codes_match_flags":[{
           "station":s["original_source_s2_col1_label"],
           "code":s["original_locality_code"],"section":s["source_section"],
           "exact_candidates":s["source_registry_full_name_candidates"],
           "status":s["match_status"]
        } for s in out["all_39_station_original_source_identity_and_candidate_provenance"]],
        "not_field_GPS":True
    },ensure_ascii=False))


if __name__=="__main__":main()
