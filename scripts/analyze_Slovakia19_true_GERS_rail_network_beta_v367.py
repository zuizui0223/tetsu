#!/usr/bin/env python3
"""v367 Slovak 2016 original Poprad–Plaveč ONE railway line:
19 original author DMS locations with actual 318-taxon station floras,
true Overture 2026 source GERS shortest-path km vs geographic km.

Study-wide species-frequency+site-richness fixed null, 19 whole-station
heldout prediction, NO same-line-vs-other-line contrast inside one line.
Year of source GIS not matched to original 2012-14 field operations.
"""
from __future__ import annotations
import argparse,json,math,random,statistics
from collections import Counter
from itertools import combinations
from pathlib import Path

from audit_original_CZ_SK_station_GERS_topology_v366 import (
    original_SK_sites,original_overture_rail_source,coordinate_source_window,
    source_rail_gps_coverage,SK_SHA
)
from recompute_Slovenia_GERS_true_geodesic_network_v364 import (
    corrected_physical_graph,corrected_source_attachments
)
from analyze_Slovenia_geod_and_rail_distance_matched_v359 import physical_network_distances

SK_SOURCE_SITE_N=19
SK_SOURCE_TAXA_N=318
SK_SOURCE_SITE_TAXON_INCIDENCES=1763
EXPECTED_GERS_CONNECTED_PAIRS=171
EXPECTED_LEGACY_GEO_SLOPE=-0.00200218
LAMBDA=.05
DRAWS=199
BURN=18000
THIN=650
SEEDS=(310901,310902)

def source_linear_slope(x,y):
    import numpy as np
    x=np.asarray(x,dtype=float);y=np.asarray(y,dtype=float)
    if len(x)!=EXPECTED_GERS_CONNECTED_PAIRS or len(y)!=len(x):
        raise ValueError("Source original Slovak 171 station pairs incomplete")
    if not np.all(np.isfinite(x)) or not np.all(np.isfinite(y)):
        raise ValueError("Original source railway distance or botanical Jaccard invalid")
    if x.max()-x.min()<5:
        raise ValueError("Source railway distance range too narrow for original study")
    return float(np.linalg.lstsq(np.stack([np.ones(len(x)),x],axis=1),y,rcond=None)[0][1])

def source_physical_distance_triangle_check(sites,attached,network):
    """Check source train-network metres cannot shortcut physical point geometry.

    Nearest source station original DMS point may be ~tens of metres off
    source mapped track. Triangle inequality allows the two source GPS→rail
    offsets, but no arbitrary km-large path shorter than straight distance.
    This QC uses NO original taxon identity or Jaccard score.
    """
    from pyproj import Geod
    geod=Geod(ellps="WGS84")
    geodelta=[];limits=[]
    rows=Counter()
    for i,j in combinations(range(len(sites)),2):
        net=network.get((i,j))
        if net is None:
            rows["source_GERS_physical_network_path_unavailable"]+=1
            continue
        src=abs(geod.inv(sites[i]["lon"],sites[i]["lat"],
                         sites[j]["lon"],sites[j]["lat"])[2])/1000
        short=src-net
        u=(attached[i].get("original_source_track_offset_m") or 0.) + (
           attached[j].get("original_source_track_offset_m") or 0.)
        rows["source_rail_pairs_with_path"]+=1
        rows["source_path_shorter_than_GPS_chord_above50m"]+=short>.05
        rows["source_path_shorter_than_GPS_chord_above100m"]+=short>.10
        rows["source_path_below_GPS_chord_minus_OWN_TRACK_offsets_plus5m"]+=short*1000>u+5.
        geodelta.append(short*1000)
        limits.append(u)
    return {
       "original_source_CZECH_SLOVAK_network_triangulation_QA_not_botany":True,
       "original_source_station_pair_QA_counts":dict(rows),
       "source_max_original_geodesic_distance_minus_GERS_railpath_m":max(geodelta) if geodelta else None,
       "source_min_original_geodesic_distance_minus_GERS_railpath_m":min(geodelta) if geodelta else None,
       "source_site_to_rail_GPS_displacement_bound_max_sum_m":max(limits) if limits else None,
       "source_same_real_171_original_Slovak_pair_support":rows["source_rail_pairs_with_path"]==171,
       "can_claim_original_network_geodesic_physical_distance_consistency":
           rows["source_path_below_GPS_chord_minus_OWN_TRACK_offsets_plus5m"]==0
    }


def station_pairs(sites,coords,network):
    from pyproj import Geod
    from analyze_Slovakia2016_19station_original_taxon_incidence_v310 import jaccard
    geo=Geod(ellps="WGS84")
    pairs=[]
    for i,j in combinations(range(19),2):
        n=network.get(i,j)
        if n is None or not 0<n<300:
            raise ValueError("Original author Slovakia ALL 171 stations must have real GERS path")
        meters=abs(float(geo.inv(sites[i]["lon"],sites[i]["lat"],
                                 sites[j]["lon"],sites[j]["lat"])[2]))
        km=meters/1000.
        # Geometrically valid source network route may be shorter than the
        # original station DMS-point chord only by source point-to-rail offsets.
        # Original per-station offsets were validated before this function.
        common=jaccard(coords[i+1],coords[j+1])
        pairs.append({"i":i,"j":j,"geo":km,"network":n,"jaccard":common})
    if len(pairs)!=EXPECTED_GERS_CONNECTED_PAIRS:
        raise ValueError("Original 19 Slovak station pairs source incomplete")
    return pairs

def fit_ridge(train,test,use_network):
    """Ridge with same fixed hyperparameter, training fold only standardization."""
    import numpy as np
    def getfeatures(x):
        g=math.log1p(x["geo"])
        vals=[g,g*g]
        if use_network:vals.append(math.log1p(x["network"]))
        return vals
    x=np.asarray([getfeatures(z) for z in train])
    y=np.asarray([z["jaccard"] for z in train])
    xt=np.asarray([getfeatures(z) for z in test])
    mu=x.mean(axis=0);sd=x.std(axis=0);sd=np.where(sd<1e-10,1.,sd)
    x=np.column_stack((np.ones(len(x)),(x-mu)/sd))
    xt=np.column_stack((np.ones(len(xt)),(xt-mu)/sd))
    w=np.ones(len(train))/len(train)
    pen=np.eye(x.shape[1]);pen[0,0]=0.
    beta=np.linalg.solve(x.T@(w[:,None]*x)+LAMBDA*pen,
                         x.T@(w*y))
    return np.clip(xt@beta,0,1)

def source_station_holdout(pairs):
    import numpy as np
    folds=[]
    for s in range(19):
        train=[r for r in pairs if s not in (r["i"],r["j"])]
        test=[r for r in pairs if s in (r["i"],r["j"])]
        if len(train)!=153 or len(test)!=18:
            raise ValueError("Original Slovakia leave whole station out must have 153 training & 18 test pairs")
        y=np.asarray([r["jaccard"] for r in test])
        a=fit_ridge(train,test,False)
        b=fit_ridge(train,test,True)
        folds.append({"source_test_M0_geography_MSE":float(np.mean((y-a)**2)),
                      "source_test_M1_geo_plus_network_MSE":float(np.mean((y-b)**2)),
                      "source_M0_minus_M1_test_MSE":float(np.mean((y-a)**2)-np.mean((y-b)**2))})
    gain=statistics.mean(x["source_M0_minus_M1_test_MSE"] for x in folds)
    return {"source_original_entire_station_holdout_folds_n":19,
            "source_holdout_pairs_train_per_fold":153,
            "source_holdout_pairs_test_per_fold":18,
            "source_geo_only_MSE":statistics.mean(x["source_test_M0_geography_MSE"] for x in folds),
            "source_geo_plus_network_MSE":statistics.mean(x["source_test_M1_geo_plus_network_MSE"] for x in folds),
            "source_independent_true_GERS_network_gain_M0_minus_M1_MSE":gain,
            "source_heldout_stations_where_network_improves_MSE":sum(x["source_M0_minus_M1_test_MSE"]>0 for x in folds),
            "source_ridge_lambda_prespecified":LAMBDA,
            "source_pairwise_test_not_19_independent_network_railways":True}

def original_species_frequency_null(taxa,pairs):
    from analyze_Slovakia2016_19station_original_taxon_incidence_v310 import (
        swap_states,station_sets,jaccard
    )
    allsites=station_sets(taxa)
    if len(allsites)!=19 or sum(len(x) for x in allsites.values())!=SK_SOURCE_SITE_TAXON_INCIDENCES:
        raise ValueError("Original PDF original source 19 station flora incidence changed")
    degrees=Counter(x for vals in allsites.values() for x in vals)
    richness={i:len(v) for i,v in allsites.items()}
    x_geo=[r["geo"] for r in pairs]
    x_net=[r["network"] for r in pairs]
    original=[r["jaccard"] for r in pairs]
    observed={
        "network":source_linear_slope(x_net,original),
        "geographic":source_linear_slope(x_geo,original)
    }
    results={"network":[],"geographic":[]};chains=[]
    for seed in SEEDS:
        draw={i:set(allsites[i]) for i in allsites}
        rng=random.Random(seed+3671124)
        acc=swap_states(draw,rng,BURN)
        local={"network":[],"geographic":[]}
        for _ in range(DRAWS):
            acc+=swap_states(draw,rng,THIN)
            if (Counter(x for z in draw.values() for x in z)!=degrees or
                any(len(draw[i])!=n for i,n in richness.items())):
                raise ValueError("Original 19 station species source presence exact site and species margins not preserved")
            y=[jaccard(draw[r["i"]+1],draw[r["j"]+1]) for r in pairs]
            for metric,v in (("network",x_net),("geographic",x_geo)):
                slope=source_linear_slope(v,y)
                results[metric].append(slope);local[metric].append(slope)
        chains.append({"source_seed":seed,"source_site_richness_and_species_frequency_preserved":True,
            "source_checkerboard_acceptance":acc/(BURN+DRAWS*THIN),
            "source_null_network_slope_mean":statistics.mean(local["network"]),
            "source_null_geographic_slope_mean":statistics.mean(local["geographic"])})
    out={}
    for metric,vals in results.items():
        z=sorted(vals)
        def q(p):
            v=p*(len(z)-1);i=int(v)
            return z[i] if i==len(z)-1 else z[i]*(1-v+i)+z[i+1]*(v-i)
        out[metric]={"source_observed_original_Jaccard_per_distance_km_slope":observed[metric],
           "original_site_and_taxon_frequency_fixed_null_mean":statistics.mean(vals),
           "source_observed_minus_null_slope":observed[metric]-statistics.mean(vals),
           "technical_null_2p5_97p5_not_field_CI":[q(.025),q(.975)],
           "technical_lower_rank_NOT_independent_field_p":
                (1+sum(z<=observed[metric]+1e-12 for z in vals))/(len(vals)+1),
           "technical_correlated_null_draws_NOT_independent_field":len(vals)}
    return {"source_original_survey_station_n":19,
            "original_source_margin_preservation_exact":True,
            "null_slopes":out,"source_null_chains":chains}

def original_station_LOO_slopes(pairs):
    import numpy as np
    out=[]
    for i in range(19):
        keep=[p for p in pairs if i not in (p["i"],p["j"])]
        if len(keep)!=153:raise ValueError("Original one-station leave-out does not retain exactly 153 source pairs")
        out.append(source_linear_slope_on_slice(keep))
    return {"source_full_station_Loo_true_network_slope_range":[min(out),max(out)],
            "original_stations_with_negative_GERS_network_distance_decay":sum(x<0 for x in out),
            "original_station_LOO_slope_n":len(out)}

def source_linear_slope_on_slice(pairs):
    import numpy as np
    a=np.asarray([r["network"] for r in pairs]);y=np.asarray([r["jaccard"] for r in pairs])
    return float(np.linalg.lstsq(np.stack([np.ones(len(a)),a],axis=1),y,rcond=None)[0][1])

def run(pdf,rail_factory=None):
    import fitz
    from pyproj import Transformer
    from analyze_Slovakia2016_19station_original_taxon_incidence_v310 import (
        source_list_clauses,parse_clauses,validate_original_source
    )
    from analyze_Slovakia2016_19station_original_physical_distance_v310 import (
        parse_source_dms,geodist
    )
    from audit_original_CZ_SK_station_GERS_topology_v366 import original_overture_rail_source
    if rail_factory is None:rail_factory=original_overture_rail_source
    source,meta=original_SK_sites(pdf)
    book=fitz.open(str(pdf))
    original_DMS_text="\n".join(book[i].get_text("text",sort=True) for i in (1,2))
    book.close()
    DMS=parse_source_dms(original_DMS_text)
    parts,_=source_list_clauses(Path(pdf))
    taxa,_=parse_clauses(parts)
    qa=validate_original_source(taxa)
    if qa["n_original_station_taxon_presence"]!=SK_SOURCE_SITE_TAXON_INCIDENCES:
        raise ValueError("Source 2016 Slovak 19 original stations and 318 taxa incompatible")
    bbox=coordinate_source_window(source)
    source_rail,railQA=rail_factory(bbox)
    T=Transformer.from_crs("EPSG:4326","EPSG:32634",always_xy=True)
    _,G,segments,topqa,metric=corrected_physical_graph(source_rail,T)
    if len(segments)!=745 or G.number_of_nodes()!=1768 or G.number_of_edges()!=1936:
        raise ValueError("v366 original Slovak 745 segment source and railway connector graph drift")
    enhanced,source_attachments,aqa=corrected_source_attachments(G,segments,source,T,threshold=180.)
    network=physical_network_distances(enhanced,source_attachments)
    cover=source_rail_gps_coverage(source,source_attachments,network)
    if (cover["source_original_true_physical_rail_matched_station_n"]!=19 or
        cover["source_real_network_connected_original_station_pairs_within_300km"]!=171):
        raise ValueError("Original v366 Slovak 19/19 and 171/171 physical mapped rail source network support drift")
    from analyze_Slovakia2016_19station_original_taxon_incidence_v310 import station_sets
    flora=station_sets(taxa)
    triangle=source_physical_distance_triangle_check(source,source_attachments,network)
    print(json.dumps({"original_Slovak_GERS_source_path_geometric_QA_BEFORE_BOTANY":triangle},
                     ensure_ascii=False),flush=True)
    if not triangle["can_claim_original_network_geodesic_physical_distance_consistency"]:
        raise ValueError("Original GERS path triangle exceeds observed GPS-to-track original source offset; no source biological network slope fitted")
    pairs=station_pairs(source,flora,network)
    vals=[p["jaccard"] for p in pairs]
    slope_network=source_linear_slope([p["network"] for p in pairs],vals)
    slope_geo=source_linear_slope([p["geo"] for p in pairs],vals)
    geo_previous=[geodist(DMS[i],DMS[j]) for i,j in combinations(range(1,20),2)]
    old_slope=source_linear_slope(geo_previous,vals)
    if abs(old_slope-EXPECTED_LEGACY_GEO_SLOPE)>1e-5:
        raise ValueError("Original publication v310 geographical distance slope changed")
    detours=[p["network"]/p["geo"] for p in pairs]
    result={"version":"v367-SK-genuine-original-19-station-one-railway-route-GERS-network-beta",
        "source_original_Slovak_author_PDF_SHA256":SK_SHA,
        "source_original_19_station_DMS_and_318_taxa_1763_incidences":True,
        "source_GERS_rail_physical_distance_consistency_QA":triangle,
        "source_actual_phyiscal_original_author_railway_routes_n":1,
        "source_overtyre_release":"2026-11-18.0",
        "source_true_19_station_GERS_connected_pairs_n":len(pairs),
        "original_field_Jaccard_mean":statistics.mean(vals),
        "original_source_true_network_km_minmax":[min(p["network"] for p in pairs),max(p["network"] for p in pairs)],
        "original_source_geodesic_km_minmax":[min(p["geo"] for p in pairs),max(p["geo"] for p in pairs)],
        "source_true_GERS_network_to_geodesic_distance_ratio_median":statistics.median(detours),
        "original_true_network_Jaccard_per_network_km_slope":slope_network,
        "original_WGS84_geodetic_Jaccard_per_geographic_km_slope":slope_geo,
        "source_prior_v310_original_approx_geographic_km_slope_reproduced":old_slope,
        "source_publisher_19_entire_stations_not_171_independent_source_botanical_transects":True,
        "source_entire_original_station_leave_one_out_network_distance_slope":original_station_LOO_slopes(pairs),
        "source_19_station_heldout_geo_and_network_prediction":source_station_holdout(pairs),
        "source_true_GERS_network_and_geographic_distance_species_frequency_null":
             original_species_frequency_null(taxa,pairs),
        "source_original_botanical_station_distances_are_NOT_train_2012_2014_OD":True,
        "original_physically_connected_graph_is_not_an_independent_between_line_comparison":True,
        "no_source_taxon_names_original_station_DMS_or_rail_graph_geometries_republished":True}
    return result

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--pdf",required=True)
    p.add_argument("--output",required=True)
    a=p.parse_args()
    data=run(a.pdf)
    dest=Path(a.output);dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(data,ensure_ascii=False,indent=2),encoding="utf8")
    print(json.dumps(data,ensure_ascii=False),flush=True)
if __name__=="__main__":main()
