#!/usr/bin/env python3
"""v23.0 genuine 2015/2020 Eurostat train movements by railway network segment.

First-party open source:
https://ec.europa.eu/eurostat/api/dissemination/statistics/1.0/data/rail_tf_ns20_si?lang=EN
Eurostat Regulation 2018/643 Annex V G1/G2. 2020 is COVID-affected,
not 2022/23 field-date exposure. This is *actual network-segment traffic
reported for 2020*, not 2025 capacity planning or national gross tonnes.

The 2020 source may only cover TEN railway corridors, not all regional
branch lines. A missing network segment is UNKNOWN, never zero traffic.
No automatic 2020 Eurostat segment -> author-original 2022/23 plant
line/GPS mapping without exact segment/geo crosswalk. Source IDs
are not assumed to equal Slovenian national route numbers.

The report prints only aggregate QA, category dictionaries, bounded
label previews and SAME source segment-code 2015 vs 2020 overlap
if demonstrably identical source ID (not re-labeled geography).
Full original JSON not republished, only SHA and bounded provenance.
"""
from __future__ import annotations
import argparse,hashlib,json,re,time,urllib.error,urllib.request
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path
from statistics import median

BASE="https://ec.europa.eu/eurostat/api/dissemination/statistics/1.0/data/"
SERIES=("rail_tf_ns20_si","rail_tf_ns15_si")
MAX_BYTES=20*1024*1024

def dimension_labels(dim):
    cat=dim.get("category",{})
    codes=cat.get("index")
    if isinstance(codes,dict):
        code_order=[k for k,v in sorted(codes.items(),key=lambda kv:kv[1])]
    elif isinstance(codes,list):
        code_order=list(codes)
    else:raise ValueError("Eurostat JSON-stat has no interpretable original category codes")
    labs=cat.get("label",{})
    if isinstance(labs,dict):
        return [(x,str(labs.get(x,x))) for x in code_order]
    return [(x,x) for x in code_order]

def unpack_jsonstat(dataset,expected_year):
    if not isinstance(dataset,dict) or dataset.get("class")!="dataset":
        raise ValueError("Not authenticated JSON-stat dataset")
    dims=dataset.get("id");size=dataset.get("size")
    if not isinstance(dims,list) or not isinstance(size,list) or len(dims)!=len(size):
        raise ValueError("Bad Eurostat dataset dimensions")
    if "net_seg20" not in dims and "net_seg15" not in dims:
        raise ValueError("Country segment source dimension not found")
    if "train" not in dims or "time" not in dims:
        raise ValueError("Train and time dimensions required")
    cat={k:dimension_labels(dataset["dimension"][k]) for k in dims}
    for key,n in zip(dims,size):
        if len(cat[key])!=n:
            raise ValueError("Original Eurostat source size/index code mismatch")
    timecat=cat["time"]
    if expected_year not in [str(x[0]) for x in timecat]:
        raise ValueError("2020/2015 actual reference year missing from genuine source")
    vals=dataset.get("value",{})
    if not isinstance(vals,(dict,list)):
        raise ValueError("Original source values not JSON-stat index keyed or ordered")
    product=1
    for n in size:
        if not isinstance(n,int) or n<1:raise ValueError("No valid category size")
        product*=n
    if product>2_000_000:raise ValueError("Oversize Eurostat country dataset")
    def getv(ix):
        if isinstance(vals,list):return vals[ix] if ix<len(vals) else None
        return vals.get(str(ix))
    def axis_indices(ix):
        remaining=ix;result=[0]*len(size)
        for i in reversed(range(len(size))):
            result[i]=remaining%size[i]
            remaining//=size[i]
        return result
    nonnull=[];miss=0;other_year=0
    for ix in range(product):
        ind=axis_indices(ix)
        name={k:cat[k][ind[i]] for i,k in enumerate(dims)}
        if str(name["time"][0])!=str(expected_year):
            other_year+=1;continue
        v=getv(ix)
        if v is None:miss+=1;continue
        if isinstance(v,bool) or not isinstance(v,(int,float)):
            raise ValueError("Non-numeric original Eurostat train value")
        if v<0:raise ValueError("Negative original train movement")
        nonnull.append({k:{"code":x[0],"label":x[1]} for k,x in name.items()} |
                       {"value":float(v),"row_index":ix})
    return nonnull,{
       "original_label":dataset.get("label"),
       "original_units":cat.get("unit"),
       "original_trains":cat.get("train"),
       "original_years":cat.get("time"),
       "dimensions":dims,"original_source_dimension_sizes":size,
       "total_possible_category_cells_for_reference_year":sum(1 for ix in range(product)
           if str(cat["time"][axis_indices(ix)[dims.index("time")]][0])==str(expected_year)),
       "nonnull_train_values_for_year":len(nonnull),"missing_2020_vs_2015_source_cells":miss,
       "unknown_or_absent_units_not_zero_imputed":True}

def query_series(code):
    url=BASE+code+"?lang=EN"
    last=None
    for trial in range(5):
        try:
            request=urllib.request.Request(url,headers={
                "User-Agent":"tetsu-original-rail-traffic-open-research/1.0 (public Eurostat)",
                "Accept":"application/json"})
            with urllib.request.urlopen(request,timeout=35) as resp:
                content=resp.read(MAX_BYTES+1)
                if len(content)>MAX_BYTES:raise ValueError("Eurostat source response oversize")
                if resp.status!=200 or len(content)<100:raise ValueError("Eurostat source invalid HTTP/body")
            data=json.loads(content)
            return data,{"source_url":url,"actual_bytes":len(content),
              "sha256_original_json_response":hashlib.sha256(content).hexdigest(),
              "utc_retrieval":datetime.now(timezone.utc).isoformat()}
        except (TimeoutError,urllib.error.HTTPError,urllib.error.URLError,ValueError) as e:
            last=str(e)
            if trial<4:time.sleep(2+2*trial)
    raise RuntimeError(f"Eurostat dataset {code} failed genuine remote source fetch: {last}")

def summarized(rows):
    segment_key="net_seg20" if "net_seg20" in rows[0] else "net_seg15"
    bytype=defaultdict(list)
    bysegment=defaultdict(list)
    for row in rows:
        label=row["train"]["label"]
        bytype[label].append(row["value"])
        seg=row[segment_key]
        bysegment[(str(seg["code"]),str(seg["label"]))].append(row)
    sample=sorted(bysegment.items(),key=lambda kv:kv[0][0])
    return {
      "source_segment_dimension":segment_key,
      "unique_nonnull_source_segments":len(bysegment),
      "source_train_types_counts":{t:len(v) for t,v in sorted(bytype.items())},
      "source_train_types_min_median_max":{t:[min(v),median(v),max(v)] for t,v in bytype.items()},
      "source_segment_label_preview_for_metadata_qa":[
           {"segment_code":k[0],"source_label":k[1],
            "actual_original_source_observed_traffic_types":
              [{"train_code":r["train"]["code"],
                "train_label":r["train"]["label"],
                "unit":r["unit"]["label"] if "unit" in r else None,
                "value":r["value"]} for r in v][:5]}
             for k,v in sample[:40]],
      "n_source_segment_labels_containing_human_named_route":sum(
         any(x in k[1].lower() for x in ("ljubljana","maribor","hodoš","zidani","koper","divača","dobova","jesenice"))
         for k in bysegment),
      "n_segment_codes_resembling_official_2digit_national_routes":sum(
           bool(re.fullmatch("[0-9]{1,2}",k[0])) for k in bysegment)}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",required=True)
    p.add_argument("--strict-2015",action="store_true")
    args=p.parse_args()
    out={"version":"v23.0","status":"2020_EUROSTAT_REAL_SEGMENT_SOURCE_AUDIT",
         "legal_provenance":"EU Regulation 2018/643 Annex V network-segment G1/G2",
         "interpretation_warning":"2020 COVID-affected observed railway traffic; not 2022–2023 field year. 2015 optional for pre-COVID rank. Network coverage may be only TEN-T.",
         "country":"Slovenia","series":{}}
    for code,year in (("rail_tf_ns20_si","2020"),("rail_tf_ns15_si","2015")):
        try:
            data,provenance=query_series(code)
            rows,qa=unpack_jsonstat(data,year)
            if not rows:raise ValueError("Genuine source lacks nonmissing source train counts")
            out["series"][code]={"status":"VERIFIED_SOURCE","source":provenance,
                                 "dimension_QA":qa,"content_AUDIT":summarized(rows)}
        except Exception as e:
            if year=="2020" or args.strict_2015:raise
            out["series"][code]={"status":"NOT_VERIFIED_NO_IMPUTED_VALUES",
                                 "error":str(e)[:500]}
    twenty=out["series"]["rail_tf_ns20_si"]
    if twenty["dimension_QA"]["nonnull_train_values_for_year"]<6:
        raise ValueError("Source 2020 has insufficient segment observations")
    years=out["series"]
    out["common_nonmissing_eurostat_segment_codes_2020_2015"] = None
    # Deliberately do not map 2015 segment IDs across year without the official
    # Eurostat release metadata establishing stable GIS segment boundaries.
    out["original_plant_site_network_crosswalk"]="NOT_YET_VERIFIED"
    out["genuine_2022_2023_train_counts_status"]="NOT_AVAILABLE_FROM_2020_EUROSTAT"
    out["no_fitted_ecological_effect_claimed"]=True
    target=Path(args.output);target.parent.mkdir(parents=True,exist_ok=True)
    target.write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding="utf-8")
    print("V230_EUROSTAT_REAL_TRAFFIC="+json.dumps(out,ensure_ascii=False,separators=(",",":")),flush=True)
if __name__=="__main__":main()
