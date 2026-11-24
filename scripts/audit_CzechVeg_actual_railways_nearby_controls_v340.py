#!/usr/bin/env python3
"""v340: source-certified CzechVeg railway substratum vs nearby-only labels.

Real source v1.1 original full CzechVeg headers SHA pinned; no source plot
species values or original coordinates rehosted. This is a matching
FEASIBILITY gate, not a biological effect estimate or valid causal control.

Only H(A)BITAT text can certify likely rail substrate; LOCALITY-only rail
mentions are landmarks. Categories must be audited by original papers before
independent biological observations can be asserted.
"""
import argparse,csv,hashlib,json,math,re,unicodedata
from collections import Counter,defaultdict
from pathlib import Path
from statistics import median
from audit_CzechVeg_railway_headers_source_v339 import (
    EXPECTED_MD5,EXPECTED_PLOTS,ascii_cz,classify,good_number)

def rail_certainty(text):
    h=ascii_cz(text)
    status,_=classify(text,"")
    if status!="rail_habitat_text_candidate":
        return "no_explicit_rail_habitat"
    # Strong counterexamples where railway is a physical landmark but the
    # plot is plainly outside the railway substrate/management unit.
    nearby=re.search(
        r"\b\d+\s*m\s+od\s+zelez|\bpod\s+zelez|\bpod\s+trati\b|"
        r"\bnad\s+zelez|\bnad\s+kolej|\bnad\s+trati|"
        r"\bokol[iy]\s+stavby|\bokol[iy]\s+zelezn|"
        r"\bhromad\w*\s+hlin|\bvilkova\s+ctvrt|"
        r"\bchodniku|\brakosina\s+pod\b|\boperna\s+cihlova\s+zed",
        h)
    if nearby:
        return "railway_nearby_or_below_not_confirmed"
    if re.search(r"\bnebyla\s+nikdy\b|\bnevyuzivan\w*",h):
        return "unused_railway_structure"
    if re.search(r"\bkolejist\w*|\bmezi\s+kolej\w*|\bzelez\w*\s+most\w*|\bzelez\w*\s+viadukt\w*",h):
        return "railway_track_yard_or_structure"
    if re.search(r"\bzelez\w*\s+nasyp\w*|\bzelez\w*\s+prikop\w*|"
                 r"\bzelez\w*\s+zarez\w*|\bzelez\w*\s+vleck\w*|"
                 r"\bsvah\s+zelez",h):
        return "railway_embankment_ditch_or_cut"
    if re.search(r"\bnadraz\w*",h):
        return "station_property_not_necessarily_track"
    return "railway_related_habitat_ambiguous"

CERTIFIED={"railway_track_yard_or_structure","railway_embankment_ditch_or_cut"}
# Two intentionally different observation design gates, fixed before source read.
MATCH_GATES={
 "G1_same_source_high_precision":dict(radius_km=.5,year_delta=3,area_ratio=2,
                                     max_precision_m=100,same_reference=True),
 "G2_same_source_extended":dict(radius_km=2,year_delta=7,area_ratio=4,
                                max_precision_m=1000,same_reference=True),
 "G3_cross_source_high_precision":dict(radius_km=.5,year_delta=3,area_ratio=2,
                                      max_precision_m=100,same_reference=False),
 "G4_cross_source_extended":dict(radius_km=2,year_delta=7,area_ratio=4,
                                 max_precision_m=1000,same_reference=False),
}

def src_metadata(row):
    lat=good_number(row.get("DEG_LAT"))
    lon=good_number(row.get("DEG_LON"))
    year=good_number(row.get("YEAR"))
    area=good_number(row.get("SURF_AREA"))
    precision=good_number(row.get("PRECISION"))
    alt=good_number(row.get("ALTITUDE"))
    if lat is not None and lon is not None:
        if not (48.5<=lat<=51.5 and 12<=lon<=19):
            lat=lon=None
    return {"lat":lat,"lon":lon,"year":int(year) if year is not None and
         1880<=year<=2026 and float(year).is_integer() else None,
         "area":area if area is not None and area>0 else None,
         "precision":precision if precision is not None and precision>=0 else None,
         "altitude":alt if alt is not None and -500<alt<3000 else None,
         "reference":str(row.get("REFERENCE") or ""),
         "cov_trees":good_number(row.get("COV_TREES")),
         "cov_shrubs":good_number(row.get("COV_SHRUBS")),
         "plot_id":str(row.get("RELEVE_NR") or "")}

def distance_km(a,b):
    rad=math.radians
    lat1,lat2=rad(a["lat"]),rad(b["lat"])
    x=math.sin((lat2-lat1)/2)**2+math.cos(lat1)*math.cos(lat2)*math.sin(
      (rad(b["lon"])-rad(a["lon"]))/2)**2
    return 12742.0176*math.asin(min(1.,math.sqrt(x)))

def in_bounds(a,b,gate):
    if any(a[k] is None or b[k] is None
           for k in ("lat","lon","year","area","precision")):
        return False,None
    if max(a["precision"],b["precision"])>gate["max_precision_m"]:
        return False,None
    if abs(a["year"]-b["year"])>gate["year_delta"]:
        return False,None
    if max(a["area"],b["area"])/min(a["area"],b["area"])>gate["area_ratio"]:
        return False,None
    if gate["same_reference"] and (not a["reference"] or
                                   a["reference"]!=b["reference"]):
        return False,None
    d=distance_km(a,b)
    return (d<=gate["radius_km"]),d

def short_area_summary(rows):
    areas=[r["area"] for r in rows if r["area"] is not None]
    years=[r["year"] for r in rows if r["year"] is not None]
    prec=[r["precision"] for r in rows if r["precision"] is not None]
    return {"n":len(rows),"area_known":len(areas),
            "area_median":median(areas) if areas else None,
            "year_known":len(years),"year_min":min(years) if years else None,
            "year_max":max(years) if years else None,
            "precision_known":len(prec),
            "precision_le100m":sum(p<=100 for p in prec),
            "precision_le1000m":sum(p<=1000 for p in prec),
            "n_references":len({r["reference"] for r in rows if r["reference"]}),
            "plausible_coordinate_pairs":sum(r["lat"] is not None for r in rows)}

def audit(path):
    raw=hashlib.md5()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda:f.read(1<<21),b""):raw.update(chunk)
    if raw.hexdigest()!=EXPECTED_MD5:
        raise ValueError("Publisher v1.1 117739 Czech plots checksum mismatch")
    type_stats=Counter()
    habitats=defaultdict(Counter)
    selected=[]
    all_rows=[]
    n=0
    with Path(path).open("r",encoding="utf-8-sig",newline="") as fh:
        reader=csv.DictReader(fh)
        required={"RELEVE_NR","HABITAT","LOCALITY","SURF_AREA","YEAR","DEG_LAT",
                  "DEG_LON","PRECISION","REFERENCE"}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError("Original Czech source metadata missing precise matching fields")
        for row in reader:
            n+=1
            status,ind=classify(row["HABITAT"],row["LOCALITY"])
            typ=(rail_certainty(row["HABITAT"])
                 if status=="rail_habitat_text_candidate"
                 else "other_or_locality_only")
            type_stats[typ]+=1
            if status=="rail_habitat_text_candidate":
                habitats[typ][ascii_cz(row["HABITAT"])[:150]]+=1
            meta=src_metadata(row)
            if typ in CERTIFIED:
                selected.append(meta)
            if status=="no_rail_text_evidence":
                all_rows.append(meta)
    if n!=EXPECTED_PLOTS:raise ValueError("Original 117739 source row count changed")
    # Coarse latitude/longitude index for all NON-rail reference plots; bin ~2km.
    grid=defaultdict(list)
    for v in all_rows:
        if v["lat"] is not None:
            key=(int(v["lat"]/.02),int(v["lon"]/.02))
            grid[key].append(v)
    matches={}
    for gate_name,gate in MATCH_GATES.items():
        n_valid_rail=0
        n_with_control=0
        counts=[]
        closest=[]
        same_ref=Counter()
        for r in selected:
            if any(r[k] is None for k in ("lat","lon","year","area","precision")):
                continue
            if r["precision"]>gate["max_precision_m"]:continue
            n_valid_rail+=1
            ix=int(r["lat"]/.02);iy=int(r["lon"]/.02)
            neigh=[]
            for i in range(ix-2,ix+3):
                for j in range(iy-3,iy+4):
                    neigh.extend(grid.get((i,j),[]))
            valid=[]
            for x in neigh:
                ok,d=in_bounds(r,x,gate)
                if ok:valid.append(d)
            if valid:
                n_with_control+=1
                counts.append(len(valid));closest.append(min(valid))
                if gate["same_reference"]:same_ref[r["reference"]]+=1
        matches[gate_name]={"gate_definition":gate,
            "n_rail_plots_with_all_meta_and_precision_gate":n_valid_rail,
            "n_rail_plots_with_at_least_one_plausible_control":n_with_control,
            "n_source_control_candidates_for_matched_rail_plot":{
                "mean":sum(counts)/len(counts) if counts else None,
                "min":min(counts) if counts else None,
                "max":max(counts) if counts else None},
            "median_nearest_potential_control_geodesic_km":
                  median(closest) if closest else None,
            "n_distinct_reference_sources_with_matched_rail_plots":
                  len(same_ref) if gate["same_reference"] else None}
    return {
        "version":"v340","doi":"10.5281/zenodo.19231152",
        "source_sha_md5":EXPECTED_MD5,
        "source_total":n,
        "source_railway_certainty_counts":dict(type_stats),
        "strict_rail_source_metadata":short_area_summary(selected),
        "certification_rules":["Only explicit on-track/rail-yard or actual rail embankment/ditch/cutting HABITAT descriptions considered certified candidates",
            "HABITAT merely near/above/below railway, or station property, not certified as true track habitat",
            "Original HABITAT cannot prove active line, substrate management, or actual sampling orientation",
            "Railway-locality landmarks rejected as habitat signals"],
        "candidate_habitat_labels_by_certainty":{
          cat:[{"normalized_habitat_label":label,"n":value}
               for label,value in freqs.most_common(30)]
          for cat,freqs in habitats.items()},
        "potential_control_feasibility":matches,
        "data_not_yet_species_linked_or_temperature_tested":True,
        "data_release_safety":["No source GPS/plot identity, no row-level botanical plots or original site descriptions exported",
              "Spatially close plots may come from different vegetation studies and habitats; cannot assert matched causal counterfactual",
              "Even same-reference controls need author-study, area, physiognomy and date matching",
              "No original EIVE or source botanical species table used in this source classification audit"]}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source-headers",required=True)
    p.add_argument("--output",required=True)
    a=p.parse_args()
    result=audit(a.source_headers)
    out=Path(a.output);out.parent.mkdir(exist_ok=True,parents=True)
    out.write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding="utf-8")
    print(json.dumps({
        "version":result["version"],"counts":result["source_railway_certainty_counts"],
        "strict_source_summary":result["strict_rail_source_metadata"],
        "potential_control_feasibility":result["potential_control_feasibility"],
        "candidate_habitat_labels_by_certainty":result["candidate_habitat_labels_by_certainty"]
        },ensure_ascii=False),flush=True)
if __name__=="__main__":main()
