#!/usr/bin/env python3
"""v22.3 original 89 Slovenian vegetation sites × 2024 OSM main/branch line tags.

Real author 2022/23 vegetation XLSX source must be SHA-verified through v10.2.
Real OSM Slovenia archive 2024-01-01 PBF must match v15.0 cache SHA.
Source 2023 official operator distinguishes main/regional but OSM usage
is crowdsourced proxy, not guaranteed to reproduce official classification.
This is a covariate support check and descriptive observational pilot;
scheduled departure counts, actual through train counts, freight OD,
mowing, wind speed, and rail-mediated seed flow are not observed.
"""
from __future__ import annotations
import argparse,hashlib,json,math
from collections import Counter,defaultdict
from pathlib import Path
from statistics import mean,median

from analyze_Slovenia_rail_microhabitats_v103 import (
    extract_original_xlsx,validate_original_keys,pair_source,locations_independent
)
from analyze_Slovenia_source_rail_functional_taxonomic_identity_v141 import (
    source_sites_and_guilds,functional_overlap,jaccard
)
from analyze_Slovenia_rail_spatial_matched_taxon_identity_v142 import (
    original_gps_coordinates,haversine_km
)

ORIGINAL_PBF_SHA="b572258fe000d6dfc76dc63afe5456121c5f463995e132a5dcd79957d819ae39"
N_SITES=89
N_SECTIONS=15
TAGGED_MAX_KM=.100
MAX_CLASS_AMBIGUITY_KM=.050
VALID_CLASSES=("main","branch")
DIST_CUTS=(0,5,10,20,40,80,160,500)

def load_classified_way_nodes(pbf):
    try:import osmium
    except ImportError as e:raise RuntimeError("Install osmium for genuine historic PBF") from e
    tags=Counter();nodes={}
    class ClassRail(osmium.SimpleHandler):
        def way(self,w):
            if w.tags.get("railway") not in ("rail","narrow_gauge"):return
            if w.tags.get("construction") or w.tags.get("disused"):return
            tags["all_source_active_rail_ways"]+=1
            use=w.tags.get("usage") or "unclassified"
            service=w.tags.get("service") or ""
            if service in ("yard","siding","spur","crossover"):return
            if use not in VALID_CLASSES:
                tags["railway_without_explicit_main_or_branch_tag"]+=1
                return
            tags["explicit_"+use+"_ways"]+=1
            for n in w.nodes:
                if not n.location.valid():continue
                ref=int(n.ref)
                v=nodes.get(ref)
                coords=(float(n.lat),float(n.lon))
                if v is None:
                    nodes[ref]=(coords,{use})
                else:
                    if any(abs(a-b)>1e-7 for a,b in zip(v[0],coords)):
                        raise ValueError("One OSM node has inconsistent geographic location")
                    v[1].add(use)
    ClassRail().apply_file(str(pbf),locations=True,idx="sparse_mem_array")
    if len(nodes)<1000 or sum(tags.values())<1000:
        raise ValueError("OSM railway way tags absent or archive no longer valid")
    return [(lat,lon,tuple(sorted(k))) for (lat,lon),k in nodes.values()],dict(tags)

def classify_geo_nodes(pts,sites,maxdist=TAGGED_MAX_KM,ambiguity=MAX_CLASS_AMBIGUITY_KM):
    import numpy as np
    from scipy.spatial import cKDTree
    km_lon=111.320*math.cos(math.radians(46));km_lat=111.132
    arr=np.asarray([[a*km_lat,b*km_lon] for a,b,_ in pts],dtype=float)
    tree=cKDTree(arr)
    out=[]
    for site in sites:
        query=(site["lat"]*km_lat,site["lon"]*km_lon)
        idx=tree.query_ball_point(query,r=maxdist*1.4)
        if not idx:
            out.append({"class":None,"reason":"no_explicit_main_branch_within_radius",
                        "nearest_explicit_tag_m":None});continue
        cand=sorted((haversine_km(site["lat"],site["lon"],pts[i][0],pts[i][1]),i) for i in idx)
        cand=[x for x in cand if x[0]<=maxdist]
        if not cand:
            out.append({"class":None,"reason":"no_explicit_main_branch_within_exact_km",
                        "nearest_explicit_tag_m":None});continue
        min_by_class={c:min((d for d,i in cand if c in pts[i][2]),default=float("inf")) for c in VALID_CLASSES}
        main,branch=min_by_class["main"],min_by_class["branch"]
        best=min(main,branch)
        if best>maxdist:
            out.append({"class":None,"reason":"too_far","nearest_explicit_tag_m":best*1000});continue
        if abs(main-branch)<ambiguity:
            out.append({"class":None,"reason":"two_main_branch_tags_too_close",
                        "nearest_explicit_tag_m":best*1000});continue
        cl="main" if main<branch else "branch"
        out.append({"class":cl,"reason":"one_unambiguous_OSM_usage_tag",
                    "nearest_explicit_tag_m":round(best*1000,1)})
    return out

def class_sections(rows):
    d=defaultdict(list)
    for r in rows:
        d[r["section"]].append(r["class"])
    ans={}
    for name,vals in d.items():
        count=Counter(v for v in vals if v in VALID_CLASSES)
        n=sum(count.values())
        strongest=count.most_common(1)[0] if count else (None,0)
        # Do not recode an internally mixed route using majority station count.
        cat=(strongest[0] if n>=3 and strongest[1]/n>=.8 else None)
        ans[name]={"class":cat,"n_source_sites":len(vals),
                   "classified_site_n":n,
                   "main":count["main"],"branch":count["branch"],
                   "unclassified_source_site_n":len(vals)-n}
    return ans

def pilot(original,pbf):
    raw=Path(pbf).read_bytes()
    digest=hashlib.sha256(raw).hexdigest()
    if digest!=ORIGINAL_PBF_SHA:
        raise ValueError("Authentic Geofabrik 2024 OSM original PBF SHA drift")
    del raw
    locs,zones,sp,traits,cb,src=validate_original_keys(original)
    paired=pair_source(locs,zones,sp,traits)
    loc=locations_independent(paired)
    fields,group,meta=source_sites_and_guilds(original)
    gps,gpsqa=original_gps_coordinates(original,fields)
    lookup={str(x["rail_location_id"]):x for x in loc}
    nodes,tags=load_classified_way_nodes(pbf)
    classes=classify_geo_nodes(nodes,gps)
    frame=[]
    for site,c in zip(gps,classes):
        p=lookup.get(str(site["site"]))
        if p is None or str(p["section_number"])!=site["section"]:
            raise ValueError("Source original BS/BT section and GPS site identifiers differ")
        frame.append({"section":site["section"],"year":site["year"],"class":c["class"],
             "classification_reason":c["reason"],
             "nearest_tagged_node_m":c["nearest_explicit_tag_m"],
             "track_type_source":str(site["track_type"]),
             "annual_delta":p["BS_minus_BT_annual_strict"],
             "therophyte_delta":p["BS_minus_BT_therophyte_strict"],
             "cover_delta":p["BS_minus_BT_cover"],
             "alien_delta":p["BS_minus_BT_explicit_alien_label"],
             "site_species":site["known"],
             "field":site})
    assert len(frame)==N_SITES
    bysection=class_sections(frame)
    mainsec={k for k,v in bysection.items() if v["class"]=="main"}
    branchsec={k for k,v in bysection.items() if v["class"]=="branch"}
    summary={}
    for kind in VALID_CLASSES:
        sub=[r for r in frame if r["class"]==kind]
        def meas(k):
            xs=[float(r[k]) for r in sub if r[k] is not None and math.isfinite(float(r[k]))]
            return {"n_sites":len(xs),"BS_minus_BT_mean":round(mean(xs),7) if xs else None,
                    "BS_minus_BT_median":round(median(xs),7) if xs else None}
        summary[kind]={"n_89_original_sites":len(sub),
          "n_original_sections":len({r["section"] for r in sub}),
          "n_sections_with_3plus_classified_sites":sum(
             v["class"]==kind for v in bysection.values()),
          "by_original_track_type":dict(Counter(r["track_type_source"] for r in sub)),
          "annual":meas("annual_delta"),"therophyte":meas("therophyte_delta"),
          "cover_source_units":meas("cover_delta"),"alien_source_status":meas("alien_delta")}
    # Deliberately descriptive pair sample, cannot be independent site pair evidence.
    pool=defaultdict(list)
    for i in range(len(frame)):
        a=frame[i]
        if a["class"] is None:continue
        for j in range(i+1,len(frame)):
            b=frame[j]
            if b["class"] is None:continue
            di=haversine_km(a["field"]["lat"],a["field"]["lon"],b["field"]["lat"],b["field"]["lon"])
            bucket=next((f"{lo}-{hi}" for lo,hi in zip(DIST_CUTS[:-1],DIST_CUTS[1:]) if lo<=di<hi),None)
            if bucket is None:continue
            k="both_main" if a["class"]==b["class"]=="main" else (
                "both_branch" if a["class"]==b["class"]=="branch" else "mixed_class")
            pool[k,bucket].append((jaccard(a["site_species"],b["site_species"]),di))
    pairs={}
    for (kind,d),v in sorted(pool.items()):
        pairs.setdefault(kind,{})[d]={"dependent_pairs":len(v),
              "jaccard_mean":round(mean(x for x,_ in v),6),
              "geographic_distance_mean_km":round(mean(x for _,x in v),4)}
    return {
      "source_scope":"authentic Slovenia 2022/23 publisher original railway vegetation + 2024-01 OSM published GIS",
      "source_original_rail_survey_locations":N_SITES,
      "n_original_sections":N_SECTIONS,
      "original_vegetation_source_89_locations_checksum":meta["original_source_supplement_SHA256"],
      "OSM_archive_sha256_checked":digest,
      "OSM_rail_way_tags":tags,
      "n_unique_tagged_main_branch_rail_nodes":len(nodes),
      "n_original_sites_by_tagged_class":dict(Counter(x["class"] or "unknown" for x in frame)),
      "classification_failure_causes":dict(Counter(x["classification_reason"] for x in frame if x["class"] is None)),
      "median_nearest_explicit_tag_m":median(x["nearest_tagged_node_m"] for x in frame if x["nearest_tagged_node_m"] is not None),
      "by_section_counts_without_coordinates":bysection,
      "group_descriptive_BS_minus_BT_annual_and_cover":summary,
      "dependent_pairwise_taxon_jaccard_by_distance":pairs,
      "regression_can_assess_main_vs_branch_only_if_both_support_ge3_sections":
        len(mainsec)>=3 and len(branchsec)>=3,
      "no_evidence_of_train_departure_count_actual_freight_OD_or_local_management":True,
      "do_not_interpret_PBF_usage_as_verified_2023_official_main_regional_or_train_count":True,
      "train_site_sections_not_random_replicates":True,
      "sampling_areas_BS_BT_not_confirmed_equal":True,
      "source_2023_official_main_regional_map":"https://infrastruktura.sz.si/wp-content/uploads/2023/02/NETWORK-STATEMENT-2023.pdf",
      "OSM_tag_guide":"https://wiki.openstreetmap.org/wiki/Key:usage"}
def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source-xlsx",required=True)
    p.add_argument("--source-audit",required=True)
    p.add_argument("--osm-pbf",required=True)
    p.add_argument("--output",required=True)
    a=p.parse_args()
    original=extract_original_xlsx(a.source_xlsx,a.source_audit)
    summary=pilot(original,a.osm_pbf)
    dst=Path(a.output);dst.parent.mkdir(parents=True,exist_ok=True)
    dst.write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding="utf-8")
    print("V223_SUMMARY="+json.dumps(summary,ensure_ascii=False,separators=(",",":")))
if __name__=="__main__":main()
