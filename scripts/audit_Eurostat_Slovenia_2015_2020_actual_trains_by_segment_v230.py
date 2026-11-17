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
import argparse,hashlib,json,math,re,time,urllib.error,urllib.request
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
      "source_all_observed_segment_label_and_traffic_metadata":[
           {"segment_code":k[0],"source_label":k[1],
            "actual_original_source_observed_traffic_types":
              [{"train_code":r["train"]["code"],
                "train_label":r["train"]["label"],
                "unit":r["unit"]["label"] if "unit" in r else None,
                "value":r["value"]} for r in v][:5]}
             for k,v in sample],
      "n_source_segment_labels_containing_human_named_route":sum(
         any(x in k[1].lower() for x in ("ljubljana","maribor","hodoš","zidani","koper","divača","dobova","jesenice"))
         for k in bysegment),
      "n_segment_codes_resembling_official_2digit_national_routes":sum(
           bool(re.fullmatch("[0-9]{1,2}",k[0])) for k in bysegment)}

def source_segment_count_rows(original_rows):
    """Source identifies endpoint-labelled physical segments, NOT national line codes."""
    arr=defaultdict(dict)
    for x in original_rows:
        key=(x['net_seg20']['code'] if 'net_seg20' in x else x['net_seg15']['code'],
             x['net_seg20']['label'] if 'net_seg20' in x else x['net_seg15']['label'])
        typ=x['train']['code']
        if typ in arr[key]:raise ValueError('Duplicated original Eurostat segment×train')
        arr[key][typ]=x['value']
    return arr

def spearman_ties_pairs(pairs):
    """Average rank on original segment observations, ties averaged."""
    if len(pairs)<8:return None
    def rank(v):
        perm=sorted(range(len(v)),key=lambda i:v[i])
        result=[0.]*len(v)
        lo=0
        while lo<len(v):
            hi=lo+1
            while hi<len(v) and v[perm[hi]]==v[perm[lo]]:hi+=1
            avg=0.5*(lo+hi-1)+1.
            for i in perm[lo:hi]:result[i]=avg
            lo=hi
        return result
    x=rank([a for a,b in pairs]);y=rank([b for a,b in pairs])
    mx=sum(x)/len(x);my=sum(y)/len(y)
    cov=sum((a-mx)*(b-my) for a,b in zip(x,y))
    den=math.sqrt(sum((a-mx)**2 for a in x)*sum((b-my)**2 for b in y))
    return cov/den if den else None

def exact_same_segment_year_comparison(year2020,year2015):
    a=source_segment_count_rows(year2020)
    b=source_segment_count_rows(year2015)
    common=set(a)&set(b)
    diff_label={k for k in set(k[0] for k in a)&set(k[0] for k in b)
                if k not in {x[0] for x in common}}
    if len(common)<15:
        raise ValueError('Fewer than 15 identical codes AND labels between Eurostat original years')
    out={'source_2015_nonmissing_segments':len(b),
         'source_2020_nonmissing_segments':len(a),
         'source_identical_endpoint_and_segment_code_in_both_years':len(common),
         'source_code_present_both_years_but_endpoint_label_different_n':len(diff_label),
         'not_assume_same_named_segment_has_unchanged_railway_geometry':True,
         'null_and_implicit_missing_not_zero_filled':True,
         'per_train_type':{}}
    for typ in ('TRN_GD','TRN_PAS'):
        pairs=[(b[k][typ],a[k][typ]) for k in common if typ in a[k] and typ in b[k]]
        med2015=median(x for x,y in pairs) if pairs else None
        med2020=median(y for x,y in pairs) if pairs else None
        ratios=[y/x for x,y in pairs if x>0]
        out['per_train_type'][typ]={
           'n_same_source_segment_code_label_with_both_numbers':len(pairs),
           'spearman_observed_2015_to_2020_frequency_rank':spearman_ties_pairs(pairs),
           'median_2015_annual_trains_on_matched_segments':med2015,
           'median_2020_annual_trains_on_matched_segments':med2020,
           'median_of_individual_2020_2015_ratio_excluding_2015zero':median(ratios) if ratios else None,
           'n_original_2015_true_zero':sum(x==0 for x,y in pairs),
           'n_original_2020_true_zero':sum(y==0 for x,y in pairs)}
    return out

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",required=True)
    p.add_argument("--strict-2015",action="store_true")
    args=p.parse_args()
    out={"version":"v23.0","status":"2020_EUROSTAT_REAL_SEGMENT_SOURCE_AUDIT",
         "legal_provenance":"EU Regulation 2018/643 Annex V network-segment G1/G2",
         "interpretation_warning":"2020 COVID-affected observed railway traffic; not 2022–2023 field year. 2015 optional for pre-COVID rank. Network coverage may be only TEN-T.",
         "country":"Slovenia","series":{}}
    raw_rows={}
    for code,year in (("rail_tf_ns20_si","2020"),("rail_tf_ns15_si","2015")):
        try:
            data,provenance=query_series(code)
            rows,qa=unpack_jsonstat(data,year)
            raw_rows[year]=rows
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
    if "2015" in raw_rows:
        out["same_original_segment_code_and_identical_name_2015_2020_comparison"]=(
            exact_same_segment_year_comparison(raw_rows["2020"],raw_rows["2015"]))
    else:
        out["same_original_segment_code_and_identical_name_2015_2020_comparison"]=None
    # Exact code and endpoint labels do not guarantee physically unmodified
    # railway segment geometry through time.
    out["original_plant_site_network_crosswalk"]="NOT_YET_VERIFIED"
    out["genuine_2022_2023_train_counts_status"]="NOT_AVAILABLE_FROM_2020_EUROSTAT"
    out["no_fitted_ecological_effect_claimed"]=True
    target=Path(args.output);target.parent.mkdir(parents=True,exist_ok=True)
    target.write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding="utf-8")
    print("V230_EUROSTAT_REAL_TRAFFIC="+json.dumps(out,ensure_ascii=False,separators=(",",":")),flush=True)
if __name__=="__main__":main()
