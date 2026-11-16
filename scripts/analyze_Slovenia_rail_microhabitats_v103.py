#!/usr/bin/env python3
"""v10.3: REAL original 2026 Slovenia railway vegetation by microhabitat.

Read *unchanged* official 401,279-byte Springer CC BY 4.0 Supplementary 1
XLSX (92 railway site surveys, 482 microhabitat zone records, 6331
species-by-zone observations, 312 trait-coded taxa; 15 rail segments),
verify its source SHA256/ZIP CRC and key joins, then directly compare
Ballast Shoulder BS vs Ballast Toe BT on SAME site and SAME track side.

We use a tiny standard-library read-only Office Open XML cell decoder
so that GitHub Actions does not need an editable workbook or an
unsupported spreadsheet writer; source workbook itself remains intact.
NO implicit treatment of the supplementary "cover" field as a
percentage: one original source value is 200, while most are <=100.
Do not assign physical zone area equal unless measured! BS/BT both are
~1–2m transition zones surveyed over 100m railway transects, but
actual widths vary; emphasize within-source trait composition shares,
not absolute species richness per square metre.

Annual fraction (Lifespan.categories exactly 'annual') and therophyte
fraction (Life.form exactly 'Th') use denominators of trait-labelled
taxa in the zone. Mixed 'annual/biennial' and 'Th/He' stay in source
and are included as explicit broad-coded sensitivities. Other
microhabitat zones SZ/IA are not coerced into the same-side BS/BT
pairing because the number of replicate sides differs by track type.
Nativeness.status blank is UNKNOWN: only explicitly 'alien' tagged
species are a fraction among all recorded species, NOT a true native/
alien compositional fraction. 2022–23 cross-sectional observations
follow a 2021 herbicide ban but cannot estimate time-treatment effects.

Results independent per 89 railway survey locations, then 15 source
line segments; do not treat 173 sides as independent railway networks.
This is a secondary reanalysis of original work by Šabić, Grašič &
Rovanšek (2026), DOI 10.1007/s44473-026-00155-z, NOT tetsu field data.
Actual freight vector, seed carriage, mowing/gust selection, fitness,
plasticity or evolved adaptation are NOT observed.
"""
from __future__ import annotations
import argparse,hashlib,io,json,math,re,zipfile
from collections import defaultdict,Counter
from datetime import datetime,timedelta,timezone
from pathlib import Path
from statistics import mean,median
from xml.etree import ElementTree as ET

SOURCE_SHA256="26be7ca9513dd54136d4ec499d8413c1cb288045f87fe3cf8814731f2edac99f"
SOURCE_DOI="10.1007/s44473-026-00155-z"
SOURCE_BYTES=401279
SHEETS=("README","locations","zones","species","traits","species_codebook")
SHEET_ROWS={"locations":92,"zones":482,"species":6331,"traits":312,"species_codebook":312}
AREAS=("BS","BT")
MIN_SOURCE_LOCATIONS_PER_SEGMENT=3
NS={"s":"http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
SURI=NS["s"]
NOT_ACTUAL_MOWING_WIND_OR_CARGO_SEEDS=True


def extract_original_xlsx(source_path,source_audit):
    raw=Path(source_path).read_bytes()
    h=hashlib.sha256(raw).hexdigest()
    audit=json.loads(Path(source_audit).read_text(encoding="utf-8"))
    if (h!=SOURCE_SHA256 or len(raw)!=SOURCE_BYTES or audit.get("complete") is not True
        or audit.get("original_supplement_sha256")!=h or
        audit.get("original_article_DOI")!=SOURCE_DOI):
        raise ValueError("Source 2026 Slovenian railway Springer supplement SHA256/provenance mismatch")
    return decode_original_xlsx_bytes(raw)


def decode_original_xlsx_bytes(raw):
    """Read-only source XML extraction; never write or edit original workbook.
    The SOURCE SHA is independently checked by extract_original_xlsx.
    """
    if not zipfile.is_zipfile(io.BytesIO(raw)):
        raise ValueError("Source Supplementary 1 is not a valid Office Open XML ZIP")
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        if z.testzip() is not None:
            raise ValueError("Original publisher XLSX ZIP CRC failed")
        root=ET.fromstring(z.read("xl/workbook.xml"))
        found=[e.attrib.get("name") for e in root.findall(".//s:sheets/s:sheet",NS)]
        if found!=list(SHEETS):
            raise ValueError("Publisher source workbook original sheet order/names changed")
        shared=ET.fromstring(z.read("xl/sharedStrings.xml"))
        strings=["".join((t.text or "") for t in node.iter(f"{{{SURI}}}t"))
                 for node in shared.findall("s:si",NS)]
        data={}
        for ix,name in enumerate(SHEETS,1):
            sheet=ET.fromstring(z.read(f"xl/worksheets/sheet{ix}.xml"))
            cells=[]
            for row in sheet.findall(".//s:sheetData/s:row",NS):
                x={}
                for cell in row.findall("s:c",NS):
                    address=cell.attrib.get("r","")
                    match=re.match(r"^([A-Z]+)[0-9]+$",address)
                    if not match:raise ValueError("Malformed original XLSX source address")
                    col=0
                    for char in match.group(1):
                        col=col*26+(ord(char)-64)
                    col-=1
                    value=cell.find("s:v",NS)
                    inline=cell.find("s:is",NS)
                    if value is None and inline is None:
                        v=None
                    elif inline is not None:
                        v="".join(n.text or "" for n in inline.iter(f"{{{SURI}}}t"))
                    elif cell.attrib.get("t")=="s":
                        v=strings[int(value.text)]
                    elif cell.attrib.get("t")=="b":
                        v=value.text=="1"
                    elif cell.attrib.get("t")=="e":
                        raise ValueError("Original published XLSX contains spreadsheet error cell")
                    else:
                        try:
                            num=float(value.text)
                        except (ValueError,TypeError):
                            v=value.text
                        else:
                            v=int(num) if num.is_integer() else num
                    x[col]=v
                if x:cells.append(x)
            if name=="README":
                continue
            if not cells:raise ValueError("Original railway source XLSX sheet is empty")
            width=max(cells[0])+1
            headers=[str(cells[0].get(i,"")) for i in range(width)]
            if any(not h for h in headers) or len(set(headers))!=len(headers):
                raise ValueError("Original railway data sheet header malformed")
            entries=[{h:row.get(i) for i,h in enumerate(headers)}
                     for row in cells[1:] if any(row.get(i) is not None for i in range(width))]
            if len(entries)!=SHEET_ROWS[name]:
                raise ValueError("Original railway sample zone/taxon/trait source row count changed for "+name)
            data[name]=entries
        return data


def float_field(value,name,allow200=False):
    if type(value) not in (int,float) or not math.isfinite(value) or value<0:
        raise ValueError("Original Slovenian rail zone "+name+" numeric source invalid")
    if not allow200 and value>100000:
        raise ValueError("Unexpected original source value of "+name)
    return float(value)


def validate_original_keys(original):
    locs=original["locations"]
    zones=original["zones"]
    taxa=original["species"]
    traits=original["traits"]
    names=original["species_codebook"]
    loc={r["location_id"]:r for r in locs}
    tr={r["EPPO_code"]:r for r in traits}
    cb={r["EPPO_code"]:r for r in names}
    if len(loc)!=92 or len(tr)!=312 or len(cb)!=312:
        raise ValueError("Original Slovenian 92-source-location or 312-species source row counts inconsistent")
    # Original published workbook has EPPO_code LOTHE in traits and
    # species observations, but corresponding Lotus herbaceus is DOYPH
    # in the codebook. Preserve source disagreement rather than silently
    # matching codes by name or deleting the 4 actual occurrence rows.
    trait_not_codebook=sorted(set(tr)-set(cb))
    codebook_not_trait=sorted(set(cb)-set(tr))
    keys={}
    for row in zones:
        k=(row["location_id"],row["zone_id"])
        if k in keys or k[0] not in loc or row.get("area") not in ("SZ","BS","BT","IA"):
            raise ValueError("Duplicated or unknown railway source location × zone")
        keys[k]=row
        float_field(row["richness"],"richness")
        if row.get("cover") is not None:float_field(row["cover"],"cover",allow200=True)
    species=defaultdict(set)
    n_original_repeated_EPPO_within_zone=0
    for r in taxa:
        k=(r["location_id"],r["zone_id"])
        t=r["EPPO_code"]
        if k not in keys or t not in tr:
            raise ValueError("Railway taxon original source zone or trait foreign key invalid")
        if t in species[k]:
            # Original Supplementary source contains repeated EPPO codes
            # for the SAME field zone. Biological richness is the UNIQUE
            # taxon set; keep raw row provenance and report duplicate count.
            n_original_repeated_EPPO_within_zone+=1
        species[k].add(t)
    for k,row in keys.items():
        if int(row["richness"])!=len(species.get(k,set())):
            raise ValueError("Original zone richness does not conserve actual observed species EPPO source")
    segments={r["track_section_number"] for r in locs}
    if len(segments)!=15:
        raise ValueError("Original source line segment universe changed from published 15")
    overs=[(row["location_id"],row["zone_id"]) for row in zones
           if row.get("cover") is not None and row["cover"]>100]
    return loc,keys,species,tr,cb,{
        "n_original_field_locations":len(loc),
        "n_original_published_zone_rows":len(keys),
        "n_original_species_zone_occurrence_rows":len(taxa),
        "n_original_duplicate_species_EPPO_rows_within_field_zone_deduplicated":n_original_repeated_EPPO_within_zone,
        "n_unique_original_EPPO_species_times_zone_source_incidence":sum(len(x) for x in species.values()),
        "n_original_312_taxa_complete_trait_codebook":len(tr),
        "trait_EPPO_codes_absent_from_species_codebook":trait_not_codebook,
        "codebook_EPPO_codes_absent_from_trait_source":codebook_not_trait,
        "n_original_occurrence_rows_with_trait_EPPO_not_in_codebook":sum(1 for r in taxa if r["EPPO_code"] in trait_not_codebook),
        "do_not_silently_reconcile_original_LOTHE_vs_DOYPH_codebook_mismatch":True,
        "n_original_railway_line_segments":len(segments),
        "n_source_cover_values_above_100_not_assumed_percent":len(overs),
        "all_482_original_zone_richness_equals_unique_species_codes":True,
        "original_supplement_had_no_verified_mowing_or_wind_exposure":True,
    }


def trait_shares(codes,traits):
    if not codes:return {"no_taxa":True}
    valid_life=[traits[k]["Life.form"] for k in codes
                if traits[k].get("Life.form") is not None]
    valid_year=[traits[k]["Lifespan.categories"] for k in codes
                if traits[k].get("Lifespan.categories") is not None]
    alien=[traits[k].get("Nativeness.status") for k in codes]
    # Blank nativeness remains unknown, NOT native.
    return {
        "species_richness_per_original_zone_not_area_standardized":len(codes),
        "strict_therophyte_share_of_known_lifeforms":
            sum(v=="Th" for v in valid_life)/len(valid_life) if valid_life else None,
        "broad_therophyte_share_including_mixed_Th_lifeforms":
            sum("Th" in v for v in valid_life)/len(valid_life) if valid_life else None,
        "strict_annual_share_of_known_lifespans":
            sum(v=="annual" for v in valid_year)/len(valid_year) if valid_year else None,
        "broad_annual_share_including_mixed_lifespans":
            sum("annual" in v for v in valid_year)/len(valid_year) if valid_year else None,
        "explicitly_labelled_alien_share_of_all_recorded_taxa":
            sum(v is not None and "alien" in v.lower() for v in alien)/len(alien),
        "fraction_unclassified_Nativeness_status_of_all_recorded_taxa":
            sum(v is None for v in alien)/len(alien),
        "no_trait_missing_taxa_set_to_zero":True,
    }


def pair_source(geo,zones,sp,traits):
    paired=[]
    for (lid,zid),one in zones.items():
        if one["area"]!="BS":continue
        if zid not in ("LBS","RBS"):
            raise ValueError("Original railway shoulder side zone id not LBS/RBS")
        side=zid[0]
        other=zones.get((lid,side+"BT"))
        if other is None:continue
        if other["area"]!="BT" or other.get("area_side")!=one.get("area_side"):
            raise ValueError("Not the same railway location and line side for inner/outer ballast")
        item={
            "location_id":lid,
            "track_side":side,
            "section_number":geo[lid]["track_section_number"],
            "track_type":geo[lid]["track_type"],
            "survey_year":(datetime(1899,12,30)+timedelta(
                days=float(geo[lid]["date"]))).year,
            "BS_cover_source_numeric_not_guaranteed_percentage":one["cover"],
            "BT_cover_source_numeric_not_guaranteed_percentage":other["cover"],
            "source_cover_value_above_100_on_either_side":bool(
                (one["cover"] is not None and one["cover"]>100)
                or (other["cover"] is not None and other["cover"]>100)),
            "BS":trait_shares(sp[(lid,zid)],traits),
            "BT":trait_shares(sp[(lid,side+"BT")],traits),
        }
        if item["survey_year"] not in (2022,2023):
            raise ValueError("Actual railway vegetation source year outside original 2022/23")
        paired.append(item)
    if len(paired)!=173 or len({p["location_id"] for p in paired})!=89:
        raise ValueError("Actual same-side BS/BT source matched pairs and locations drifted")
    if len({p["section_number"] for p in paired})!=15:
        raise ValueError("Original railway segments not represented in paired BS/BT source")
    return paired


MEASURES=(
    ("cover","BS_cover_source_numeric_not_guaranteed_percentage",
     "BT_cover_source_numeric_not_guaranteed_percentage"),
    ("richness","species_richness_per_original_zone_not_area_standardized",None),
    ("therophyte_strict","strict_therophyte_share_of_known_lifeforms",None),
    ("annual_strict","strict_annual_share_of_known_lifespans",None),
    ("therophyte_broad","broad_therophyte_share_including_mixed_Th_lifeforms",None),
    ("annual_broad","broad_annual_share_including_mixed_lifespans",None),
    ("explicit_alien_label","explicitly_labelled_alien_share_of_all_recorded_taxa",None),
)


def pair_measure(pair,measure):
    f=next(x for x in MEASURES if x[0]==measure)
    if f[2] is not None:
        a=pair[f[1]]
        b=pair[f[2]]
    else:
        a=pair["BS"].get(f[1])
        b=pair["BT"].get(f[1])
    if a is None or b is None:return None
    if type(a) not in (int,float) or type(b) not in (int,float):
        raise ValueError("Original railway microhabitat source not quantitative")
    return float(a-b)


def locations_independent(pairs):
    byloc=defaultdict(list)
    for p in pairs:
        byloc[p["location_id"]].append(p)
    frame=[]
    for lid,side_pairs in sorted(byloc.items()):
        a=side_pairs[0]
        if any(z["section_number"]!=a["section_number"] or
               z["survey_year"]!=a["survey_year"] for z in side_pairs):
            raise ValueError("Railway source location sides have contradictory sections/years")
        row={
            "rail_location_id":lid,
            "section_number":a["section_number"],
            "track_type":a["track_type"],
            "survey_year":a["survey_year"],
            "n_matched_same_side_BS_BT_pairs":len(side_pairs),
            "has_any_original_cover_above100":any(
                z["source_cover_value_above_100_on_either_side"] for z in side_pairs),
        }
        for name,_,_ in MEASURES:
            sample=[pair_measure(p,name) for p in side_pairs]
            nonnull=[v for v in sample if v is not None]
            row["BS_minus_BT_"+name]=mean(nonnull) if nonnull else None
        frame.append(row)
    if len(frame)!=89:
        raise ValueError("Roadside railway source location cohort changed")
    return frame


def location_summary(locations,measure):
    col="BS_minus_BT_"+measure
    src=[x[col] for x in locations if x[col] is not None]
    return {
        "n_original_field_locations_with_measure":len(src),
        "mean_BS_minus_BT":round(mean(src),6) if src else None,
        "median_BS_minus_BT":round(median(src),6) if src else None,
        "n_locations_BS_greater":sum(v>0 for v in src),
        "n_locations_BS_equal":sum(v==0 for v in src),
        "n_locations_BT_greater":sum(v<0 for v in src),
        "no_independent_region_p_values":True
    }


def section_summary(locations,measure,min_locations=MIN_SOURCE_LOCATIONS_PER_SEGMENT):
    by=defaultdict(list)
    for site in locations:
        y=site["BS_minus_BT_"+measure]
        if y is not None:by[site["section_number"]].append(y)
    rows=[{"section_number":str(k),"n_field_locations":len(vals),
           "mean_BS_minus_BT":round(mean(vals),6)}
          for k,vals in sorted(by.items(),key=lambda q:str(q[0]))
          if len(vals)>=min_locations]
    return {
        "min_source_locations_per_rail_segment":min_locations,
        "n_eligible_original_rail_segments_NOT_countries":len(rows),
        "n_segments_BS_greater":sum(r["mean_BS_minus_BT"]>0 for r in rows),
        "n_segments_BS_equal":sum(r["mean_BS_minus_BT"]==0 for r in rows),
        "n_segments_BT_greater":sum(r["mean_BS_minus_BT"]<0 for r in rows),
        "median_rail_segment_equal_mean_BS_minus_BT":round(
            median(r["mean_BS_minus_BT"] for r in rows),6) if rows else None,
        "rail_segments":rows
    }


def run(original):
    loc,zone,sp,trait,codebook,integrity=validate_original_keys(original)
    p=pair_source(loc,zone,sp,trait)
    locations=locations_independent(p)
    result={
        "complete":True,
        "UTC":datetime.now(timezone.utc).isoformat(),
        "independent_original_source_DOI":SOURCE_DOI,
        "underlying_original_published_RAIL_field_surveys_not_MIREN_roadside":True,
        "year_source":"2022–2023 after 2021 Slovenia herbicide ban, not repeated pre/post intervention",
        "original_source_integrity":integrity,
        "same_location_same_side_matched_ballast_shoulder_vs_toe":{},
        "n_same_side_BS_BT_zone_pairs_NOT_independent_field_locations":len(p),
        "n_distinct_actual_railway_field_locations":len(locations),
        "n_actual_original_rail_sections":15,
        "section_level_support":{},
        "track_type_source_sensitivity":{},
        "survey_year_source_sensitivity":{},
        "cover_original_gt100_exclusion_sensitivity":{},
        "DO_NOT_INFER_actual_train_wind_mowing_viable_seed_or_adaptation":True,
    }
    for name,_,_ in MEASURES:
        result["same_location_same_side_matched_ballast_shoulder_vs_toe"][name]=location_summary(
            locations,name)
        result["section_level_support"][name]=section_summary(locations,name)
        for typ in ("SINGLE","DOUBLE"):
            result["track_type_source_sensitivity"].setdefault(typ,{})[name]=location_summary(
                [x for x in locations if x["track_type"]==typ],name)
        for year in (2022,2023):
            result["survey_year_source_sensitivity"].setdefault(str(year),{})[name]=location_summary(
                [x for x in locations if x["survey_year"]==year],name)
        result["cover_original_gt100_exclusion_sensitivity"][name]=location_summary(
            [x for x in locations if not x["has_any_original_cover_above100"]],name)
    result["limitations"]=[
        "This is genuine RAILWAY floristic field data, NOT direct measured seed carriage from freight cars or trackside seed traps; 2026 paper original authors own survey observations.",
        "BS/BT are similar named ~1–2m ballast zones over 100m transects but real measured widths/areas are not equalized, so raw per-zone richness is not species/100m2 comparable to MIREN ROAD plots.",
        "The source 'cover' field contains one value 200; we treat it as reported COVER UNITS not unconditionally a 0–100 percentage, and provide exclusion sensitivity.",
        "No plant-source confirmed native status: blank Nativeness.status is unknown, not native. Only explicitly tagged alien taxa are classified; all other taxa not called natives.",
        "Only 89 surveyed RAIL locations across 15 nonrandomly selected sections, from one country; two sides of track per location are not independent 173 sites.",
        "2022–2023 is after 2021 herbicide policy change but no true before-treatment survey, mowing/chemical application or slipstream velocity per plot; cannot identify policy/stress effect.",
        "Trait 'annual' and 'Th' are plant species-level categories, not measured plasticity/fitness or evolved site-specific genotype.",
        "Do not naively compare Slovenian RAIL percentages and MIREN global ROAD 100m2 richness as if they were matched identical experimental plots.",
    ]
    return result


def main():
    ap=argparse.ArgumentParser()
    for name in ("source-xlsx","verified-v102","output"):
        ap.add_argument("--"+name,required=True)
    a=ap.parse_args()
    original=extract_original_xlsx(a.source_xlsx,a.verified_v102)
    report=run(original)
    path=Path(a.output);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps({
        "complete":True,
        "source":report["original_source_integrity"],
        "n_same_track_side_pairs":report["n_same_side_BS_BT_zone_pairs_NOT_independent_field_locations"],
        "n_field_locations":report["n_distinct_actual_railway_field_locations"],
        "location_means":report["same_location_same_side_matched_ballast_shoulder_vs_toe"],
        "rail_segment_signs":{k:{key:v[key] for key in (
            "n_eligible_original_rail_segments_NOT_countries",
            "n_segments_BS_greater","n_segments_BT_greater",
            "median_rail_segment_equal_mean_BS_minus_BT")}
            for k,v in report["section_level_support"].items()},
        "track_type":{k:{t:v[t]["mean_BS_minus_BT"] for t in (
            "annual_strict","therophyte_strict","cover")}
            for k,v in report["track_type_source_sensitivity"].items()},
        "no_rail_freight_seed_or_mowing_gust_causal_claim":True
    },ensure_ascii=False),flush=True)


if __name__=="__main__":
    main()
