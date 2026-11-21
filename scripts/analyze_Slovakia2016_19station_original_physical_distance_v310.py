#!/usr/bin/env python3
"""v310 geographic sensitivity with ACTUAL printed author DMS station positions.

This is post-hoc sensitivity after ordinal-station source result. Original
genuine 318 taxa×19 stations, five richness checks and nine universal source
taxa remain hard prerequisites. Source published station positions, not
running trains/actual train OD/travel distance.
"""
from __future__ import annotations
import argparse,json,math,re,random
from collections import Counter
from pathlib import Path
import numpy as np
import fitz
from analyze_Slovakia2016_19station_original_taxon_incidence_v310 import (
    source_list_clauses,parse_clauses,validate_original_source,station_sets,beta_slope,
    swap_states,DRAWS_PER_CHAIN,BURN,THIN,SEEDS,jaccard,
)
from audit_Slovakia2016_original_19station_318taxa_pdf_preflight_v310 import preflight

GPS_PATTERN=re.compile(
 r"(?<!\d)(?P<id>[1-9]|1[0-9])\.\s+"
 r"(?P<name>[A-ZÀ-ž][^,\n]{2,80})\s*,\s*"
 r"(?P<latd>\d{2})°\s*(?P<latm>\d{1,2})[´′']\s*(?P<lats>\d{1,2})[˝″\"]"
 r"\s*,\s*"
 r"(?P<lond>\d{2})°\s*(?P<lonm>\d{1,2})[´′']\s*(?P<lons>\d{1,2})[˝″\"]"
 r"\s*,\s*(?P<alt>\d{2,4})\s*m\s*n\.\s*m\.",
 re.UNICODE)

def parse_source_dms(text):
    coords={}
    for m in GPS_PATTERN.finditer(text):
        d=m.groupdict();idx=int(d["id"])
        if idx in coords:raise ValueError("Original Slovak source station DMS repeated")
        lat=int(d["latd"])+int(d["latm"])/60+int(d["lats"])/3600
        lon=int(d["lond"])+int(d["lonm"])/60+int(d["lons"])/3600
        altitude=int(d["alt"])
        if not (48<=lat<51 and 18<=lon<23 and 300<=altitude<=1200):
            raise ValueError("Original station source DMS outside Slovakia valley bounds")
        coords[idx]=(lat,lon,altitude,d["name"].strip())
    if len(coords)!=19 or sorted(coords)!=list(range(1,20)):
        raise ValueError("Original station coordinate source incomplete: "+str(sorted(coords)))
    return coords

def geodist(a,b):
    lat1,lon1=a[:2];lat2,lon2=b[:2]
    phi1,phi2=map(math.radians,(lat1,lat2))
    dphi=phi2-phi1;dlam=math.radians(lon2-lon1)
    h=math.sin(dphi/2)**2+math.cos(phi1)*math.cos(phi2)*math.sin(dlam/2)**2
    return 6371.0088*2*math.asin(min(1,math.sqrt(h)))

def calc_geographic_original(pdf:Path,draws=DRAWS_PER_CHAIN):
    audit=preflight(pdf)
    d=fitz.open(str(pdf))
    # Station gazetteer is printed across pages 38,39 (PDF pages2,3)
    source="\n".join(d[i].get_text("text",sort=True) for i in (1,2))
    d.close()
    coords=parse_source_dms(source)
    entries,_=source_list_clauses(pdf)
    obs,_=parse_clauses(entries)
    val=validate_original_source(obs)
    source_sites=station_sets(obs)
    x=np.array([geodist(coords[i],coords[j]) for i in range(1,20) for j in range(i+1,20)])
    if not (0<float(x.min())<float(x.max())<100):
        raise ValueError("Original station physical separation invalid")
    b,mean=beta_slope(source_sites,lambda i,j:geodist(coords[i],coords[j]))
    bo,_=beta_slope(source_sites,lambda i,j:abs(i-j))
    if b>=0:sign="positive"
    else:sign="negative"
    # Whole physical station holdout; each leaves all dependent pairs involving
    # that physical station, NOT one of the 171 dependent dyads.
    loo=[]
    for excluded in range(1,20):
        pairs=[(i,j) for i in range(1,20) for j in range(i+1,20) if excluded not in (i,j)]
        y=np.array([jaccard(source_sites[i],source_sites[j]) for i,j in pairs])
        X=np.array([[1.,geodist(coords[i],coords[j])] for i,j in pairs])
        loo.append(float(np.linalg.lstsq(X,y,rcond=None)[0][1]))
    source_counts={i:len(source_sites[i]) for i in source_sites}
    baseline_freq=Counter(t for r in source_sites.values() for t in r)
    null=[];chain=[]
    for ix,seed in enumerate(SEEDS):
        rr=random.Random(seed+123)
        work={k:set(v) for k,v in source_sites.items()}
        accepted=swap_states(work,rr,BURN)
        z=[]
        for _ in range(draws):
            accepted+=swap_states(work,rr,THIN)
            if any(len(work[i])!=n for i,n in source_counts.items()) or Counter(t for r in work.values() for t in r)!=baseline_freq:
                raise ValueError("Source site richness/occupancy invalid under physical-distance swap null")
            v,_=beta_slope(work,lambda i,j:geodist(coords[i],coords[j]))
            null.append(v);z.append(v)
        chain.append({"chain":ix+1,"acceptance":accepted/(BURN+draws*THIN),
                      "mean_null_geodesic_km_slope":float(np.mean(z))})
    ns=np.array(null)
    # Maintain source nonindependence and don't call technical chain tail a field p value.
    return {
       "version":"v310_physical_distance_posthoc_sensitivity",
       "original_source_SHA256":audit["pdf_sha256"],
       "original_stations_DMS_true_in_author_source":len(coords),
       "original_botanical_taxa_strict_original_incidence_qa":val["n_pub_source_total_taxa"],
       "original_taxa_x_true_station_presence":val["n_original_station_taxon_presence"],
       "source_extremes_elevation_m":[min(v[2] for v in coords.values()),max(v[2] for v in coords.values())],
       "source_stations_reported_transect_length_km":60,
       "pairwise_geographic_distance_km_minmax":[float(x.min()),float(x.max())],
       "n_source_dyads_depend_on_only_19_stations":len(x),
       "original_observed_Jaccard_mean":mean,
       "source_original_geodesic_distance_km_slope":b,
       "source_original_station_order_index_slope_reproduced":bo,
       "leave_one_entire_original_station_out_geodesic_slope_minmax":[min(loo),max(loo)],
       "n_leave_source_station_out_negative_geodesic_slope":sum(i<0 for i in loo),
       "original_source_taxa_frequency_and_site_richness_fixed_null_mean":float(ns.mean()),
       "observed_minus_source_fixed_null_slope":float(b-ns.mean()),
       "technical_null_2p5_97p5_NOT_field_CI":list(map(float,np.quantile(ns,[.025,.975]))),
       "technical_lower_rank_NOT_independent_field_p":float((1+sum(ns<=b))/(len(ns)+1)),
       "technical_draws_NOT_independent_plots":len(ns),
       "strict_source_original_19_sites_318_taxa_margins_maintained":True,
       "original_publisher_station_DMS_site_coordinates_NOT_republished":True,
       "original_freight_and_passenger_train_count_data_not_recorded":True,
       "source_railway_station_botanical_survey_period_2012_14_not_official_live_train_year":True,
       "source_reference_caveat":"same river valley along one rail line; direct distance correlates with elevation/climate and local regional species pools, not a seed transport experiment",
       "technical_shuffle_chain_diagnostics":chain
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--pdf",required=True,type=Path)
    p.add_argument("--out",required=True,type=Path)
    args=p.parse_args()
    out=calc_geographic_original(args.pdf)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(out,ensure_ascii=False,indent=2)+"\n",encoding="utf8")
    print(json.dumps(out,ensure_ascii=False,indent=2))
if __name__=="__main__":main()
