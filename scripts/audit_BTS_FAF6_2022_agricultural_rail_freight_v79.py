#!/usr/bin/env python3
"""v7.9: real public BTS FAF6.0 STATE 2022 freight-origin/destination audit.

BTS FAF6 provides estimated REGIONAL commodity flow weights by domestic
transport mode and SCTG2 category. PREDECLARED railway-vs-truck exposure
and potential crop-linked propagule supply are not actual shipments of
viable non-native seeds, path-specific train usage, seed establishment,
railway wind, mowing, plasticity or local evolution.

Only aggregate state-origin/state-destination groupings and national
source measures leave this script; do not claim that regional OD
tonnage passed any particular 50km mapped railway or FRA crossing.
"""
from __future__ import annotations
import argparse,csv,io,json,math,zipfile
from collections import Counter,defaultdict
from datetime import datetime,timezone
from hashlib import sha256
from pathlib import Path
import requests

URL="https://faf.ornl.gov/faf6/data/Download_Files/ZIP/FAF6.0_State.zip"
OFFICIAL_PAGE="https://www.bts.gov/faf"
VERSION="FAF6.0 / 2022 BASE YEAR / STATE OD CSV"
MAX_DOWNLOAD_BYTES=80000000
MAX_UNCOMPRESSED_BYTES=380000000
MAX_CSV_ROWS=4500000
MODES={"1":"TRUCK","2":"RAIL"}
TRADE={"1":"domestic","2":"imports","3":"exports"}
CROPLINKED={"02","03"}
RESOURCE_CONTROLS={"10","11","12"}
REQUIRED={"dms_orig","dms_dest","sctg2","dms_mode","trade_type","tons_2022"}


def normalize_code(value):
    s=str(value or "").strip()
    if s=="":
        return None
    try: v=int(s)
    except ValueError:return None
    return f"{v:02d}" if 0<=v<=999 else None


def source_float(value):
    try:val=float(str(value).strip())
    except (TypeError,ValueError):
        raise ValueError("Malformed official FAF6.0 tonnage source field")
    if not math.isfinite(val) or val<0:
        raise ValueError("Invalid negative/nonfinite FAF6 flow tonnage")
    return val


def records_from_faf6_zip(data):
    if len(data)<250 or len(data)>MAX_DOWNLOAD_BYTES or not zipfile.is_zipfile(io.BytesIO(data)):
        raise ValueError("Official BTS FAF6 STATE source ZIP missing/truncated/too large")
    with zipfile.ZipFile(io.BytesIO(data)) as arc:
        candidates=[f for f in arc.infolist() if f.filename.lower().endswith(".csv")
                    and "faf6" in f.filename.lower()
                    and not f.is_dir()]
        if len(candidates)!=1:
            raise ValueError("Expected one genuine FAF6 STATE CSV inside official source ZIP")
        src=candidates[0]
        if src.file_size<=100 or src.file_size>MAX_UNCOMPRESSED_BYTES:
            raise ValueError("Official FAF6 CSV source size outside safe bounded range")
        with arc.open(src) as payload:
            with io.TextIOWrapper(payload,encoding="utf-8-sig",newline="") as stream:
                reader=csv.DictReader(stream)
                if reader.fieldnames is None or not REQUIRED.issubset(set(reader.fieldnames)):
                    raise ValueError("FAF6.0 STATE CSV missing official mode/commodity/year fields. Actual original columns: "+repr(reader.fieldnames))
                for i,row in enumerate(reader):
                    if i>=MAX_CSV_ROWS:
                        raise ValueError("FAF6 source exceeds predeclared max rows")
                    if not isinstance(row,dict) or None in row:
                        raise ValueError("FAF6 CSV source malformed field count")
                    yield row


def summarize_records(records):
    all_rows=0
    coded=Counter()
    sums=defaultdict(float)
    rowcounts=Counter()
    origin=defaultdict(lambda:defaultdict(float))
    destination=defaultdict(lambda:defaultdict(float))
    per_sctg=defaultdict(float)
    total_counted=0
    unique_state_codes=set()
    for row in records:
        all_rows+=1
        mode=normalize_code(row.get("dms_mode"))
        tradetype=normalize_code(row.get("trade_type"))
        sctg=normalize_code(row.get("sctg2"))
        o=normalize_code(row.get("dms_orig"))
        d=normalize_code(row.get("dms_dest"))
        tons=source_float(row.get("tons_2022"))
        if o is None or d is None or sctg is None:
            coded["incomplete_origin_dest_or_commodity"]+=1
            continue
        unique_state_codes.update((o,d))
        if tradetype not in {"01","02","03"}:
            coded["unsupported_trade_category"]+=1
            continue
        tr=TRADE[str(int(tradetype))]
        if mode not in ("01","02"):
            coded["other_mode_not_in_rail_truck_comparator"]+=1
            continue
        m=MODES[str(int(mode))]
        group=("crop_codes_02_03" if sctg in CROPLINKED else
               "resource_codes_10_11_12" if sctg in RESOURCE_CONTROLS else
               "all_other_categories")
        source_key=(tr,group,m)
        sums[source_key]+=tons
        rowcounts[source_key]+=1
        total_counted+=1
        if tr=="domestic":
            origin[(o,group)][m]+=tons
            destination[(d,group)][m]+=tons
            per_sctg[(sctg,m)]+=tons
    if all_rows<500 or total_counted<100 or not unique_state_codes:
        raise ValueError("Official FAF6 source population too small/incomplete")
    if not sums.get(("domestic","crop_codes_02_03","RAIL"),0)>0 or not sums.get(
            ("domestic","crop_codes_02_03","TRUCK"),0)>0:
        raise ValueError("No actual 2022 domestic crop-linked freight under rail AND truck")
    def ratio(v):
        den=v.get("RAIL",0)+v.get("TRUCK",0)
        return round(v.get("RAIL",0)/den,6) if den>0 else None
    def values(v):
        return {"rail_thousand_short_tons":round(v.get("RAIL",0),6),
                "truck_thousand_short_tons":round(v.get("TRUCK",0),6),
                "rail_fraction_among_only_rail_and_truck":ratio(v),
                "both_modes_positive":v.get("RAIL",0)>0 and v.get("TRUCK",0)>0}
    national={}
    for trade in TRADE.values():
        national[trade]={}
        for group in ("crop_codes_02_03","resource_codes_10_11_12","all_other_categories"):
            subset={mode:sums.get((trade,group,mode),0) for mode in ("RAIL","TRUCK")}
            national[trade][group]={**values(subset),
                "n_actual_source_rows":{
                    mode:rowcounts.get((trade,group,mode),0) for mode in ("RAIL","TRUCK")}}
    def asdict(grouped):
        result={}
        for (state,group),value in sorted(grouped.items()):
            result.setdefault(state,{})[group]=values(value)
        return result
    states_o=asdict(origin);states_d=asdict(destination)
    crop_states=[k for k,v in states_o.items()
                 if v.get("crop_codes_02_03",{}).get("both_modes_positive")]
    return {
        "n_original_official_FAF6_state_OD_commodity_mode_rows":all_rows,
        "n_rows_rail_or_truck_comparator":total_counted,
        "n_distinct_state_code_categories_in_source":len(unique_state_codes),
        "source_not_classified_reason_counts":dict(coded),
        "national_trade_type_x_commodity_group_2022":national,
        "n_state_origin_codes_with_domestic_crop_rail_and_truck_both":len(crop_states),
        "by_source_state_origin_code_DOMESTIC_ONLY":states_o,
        "by_source_state_destination_code_DOMESTIC_ONLY":states_d,
        "by_domestic_crop_commodity_code_02_03":{
            code:{m:round(per_sctg.get((code,m),0),6) for m in ("RAIL","TRUCK")}
            for code in sorted(CROPLINKED)},
        "source_mode_codes":{"1":"Truck","2":"Rail"},
        "source_crop_potential_codes":{"02":"Cereal grains including seed","03":"Other agricultural products except feed"},
        "source_resource_comparator_codes":sorted(RESOURCE_CONTROLS),
        "NOT_actual_verified_seed_propagules_or_track_section_freight_volume":True,
    }


def run(fetch=requests.get):
    r=fetch(URL,timeout=(30,180),headers={"User-Agent":"tetsu-US-FAF6-real-state-crop-railroad-source-v79"})
    r.raise_for_status()
    raw=r.content
    agg=summarize_records(records_from_faf6_zip(raw))
    return {
        "complete":True,
        "read_at_utc":datetime.now(timezone.utc).isoformat(),
        "official_source_url":URL,"official_publisher_page":OFFICIAL_PAGE,
        "FAF_release":VERSION,
        "official_published_CSV_ZIP_sha256":sha256(raw).hexdigest(),
        "source_compressed_bytes":len(raw),
        "measurement_unit":"thousands of short tons of 2022 estimated commodity flows",
        "data":agg,
        "not_tonnes_assigned_to_50km_rail_or_FRA_crossing_or_viable_seeds":True,
        "limitations":[
            "FAF6 is ESTIMATED commodity origin-destination flow, not train-level seed load measurements.",
            "Rail mode is FAF DOMESTIC mode for the OD shipment, not actual railroad track or route taken; do not overlay OD tons onto a shortest route as verified traffic.",
            "Commodity SCTG2=02 includes cereals, some seed; SCTG2=03 includes other agriculture, but neither measures viable non-native propagule contamination.",
            "This is 2022 one-year state-scale mode/commodity data, whereas FRA reported crossing counts are 2018-2026 and GBIF botanical records are 2018-2025 in preselected 50km cells.",
            "No direct causal seed delivery, railway verge establishment, mowing, wind or genetic adaptation effects can be inferred.",
            "Other transport modes are excluded from the RAIL/(RAIL+TRUCK) ratio denominator; this is NOT fraction of all freight modes.",
            "State code categories are source-declared economic zones, not geographic intersection with previously studied EPSG:6933 50km railway cells.",
        ]
    }


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--output",required=True)
    args=ap.parse_args()
    out=run()
    dst=Path(args.output);dst.parent.mkdir(parents=True,exist_ok=True)
    dst.write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps({
       "complete":True,"source_rows":out["data"]["n_original_official_FAF6_state_OD_commodity_mode_rows"],
       "states_with_both_crop_modes":out["data"]["n_state_origin_codes_with_domestic_crop_rail_and_truck_both"],
       "domestic_crop_modes":out["data"]["national_trade_type_x_commodity_group_2022"]["domestic"]["crop_codes_02_03"],
       "NOT_a_train_seed_dispersal_or_stress_adaptation_effect":True
    },ensure_ascii=False),flush=True)


if __name__=="__main__":
    main()
