#!/usr/bin/env python3
"""v359b: exact same original 27 Slovenian geo+rail network matched source
site triplets tested against technical fixed-site-richness and species-
frequency-preserving source null, including same year×track layout null.
Real original site species sets (89) source verified. Not a field p-value.
"""
from __future__ import annotations
import argparse,json,math
from pathlib import Path

from analyze_Slovenia_rail_microhabitats_v103 import extract_original_xlsx
from analyze_Slovenia_rail_route_nearest_plant_microhabitat_v277 import to_original_sites
from analyze_Slovenia_Overture_connector_network_botany_v358 import (
    real_overture,build_rail_graph,site_attach,SOURCE_UTM,SOURCE_ORIGINAL_SHA,RELEASE
)
from analyze_Slovenia_geod_and_rail_distance_matched_v359 import (
    physical_network_distances,joint_geo_and_network_triples,
    botanical_comparison,SOURCE_V358_MATCHED,SOURCE_V358_AMBIGUOUS
)
from audit_Slovenia_same_line_geodistance_matched_frequency_null_v281 import (
    source_groups,original_matched_gap_for_rows,experiment
)

EXPECTED_ORIGINAL_SITES=89
EXPECTED_ORIGINAL_SECTIONS=15
EXPECTED_NETWORK_MATCHED=(86,3)
EXPECTED_MATCHED_NETWORK_TAXONOMIC_GAP=0.027033311999075067
EXPECTED_GEO_NETWORK_ANCHORS=27
EXPECTED_GEO_NETWORK_AUTHOR_SECTIONS=7
PANEL_NULLS=("all_original_sites","within_original_survey_year_and_track_layout")
DRAWS=99
BURN=3000
THIN=350
SEED=20261124+359_000

def shape_original_for_null(sites,triads):
    pairs=[]
    for x in triads:
        pairs.append({"anchor_index":x["anchor"],
                      "within_index":x["same"],
                      "other_index":x["different"]})
    rows=[set(s["BS"]|s["BT"]) for s in sites]
    if len(rows)!=EXPECTED_ORIGINAL_SITES:
        raise ValueError("Unexpected original source flora site × taxa incidence population")
    baseline,sections=original_matched_gap_for_rows(rows,sites,pairs)
    if len(pairs)!=EXPECTED_GEO_NETWORK_ANCHORS or len(sections)!=EXPECTED_GEO_NETWORK_AUTHOR_SECTIONS:
        raise ValueError("Original v359 true physical network source matching support drift")
    if abs(baseline-EXPECTED_MATCHED_NETWORK_TAXONOMIC_GAP)>1e-10:
        raise ValueError("Source original v359 Jaccard baseline drift under network+geo matched identical 27 sites")
    return pairs,baseline

def run(original,rail_factory=real_overture):
    from pyproj import Transformer
    bysection,metadata=to_original_sites(original)
    sites=[s for s in sorted(bysection) for s in bysection[s]]
    if len(sites)!=EXPECTED_ORIGINAL_SITES or len(bysection)!=EXPECTED_ORIGINAL_SECTIONS:
        raise ValueError("Original botanical source n=89/15 no longer authentic")
    b=(min(x["lon"] for x in sites)-.12,min(x["lat"] for x in sites)-.12,
       max(x["lon"] for x in sites)+.12,max(x["lat"] for x in sites)+.12)
    original_rails,railqa=rail_factory(b)
    transformer=Transformer.from_crs("EPSG:4326",SOURCE_UTM,always_xy=True)
    graph,segments,geoQA=build_rail_graph(original_rails,transformer)
    network,attached=site_attach(graph,segments,sites,transformer)
    nmatch=sum(x["status"]=="matched_to_real_rail_connector_network" for x in attached)
    namb=sum(x["status"]=="ambiguous_multiple_disconnected_source_tracks" for x in attached)
    if (nmatch,namb)!=EXPECTED_NETWORK_MATCHED:
        raise ValueError("Source original v359 geographic track source intersection drift")
    network_paths=physical_network_distances(network,attached)
    triplets,selectqa=joint_geo_and_network_triples(sites,network_paths,"geo_and_network")
    bqa=botanical_comparison(sites,triplets)
    if not bqa["source_qualifies_predefined_minimum_12_anchors_4_sections"]:
        raise ValueError("Actual network+geodesic matched original botanical panel not identifiable")
    pairs,observed=shape_original_for_null(sites,triplets)
    if abs(bqa["source_habitat_comparisons"]["BS_and_BT"]["source_author_section_equal_same_minus_other_Jaccard"]-observed)>1e-10:
        raise AssertionError("Original source 27/7 Jaccard in v359/technical null disagree")
    checks={}
    for i,family in enumerate(PANEL_NULLS):
        groups=source_groups(sites,family)
        out=experiment(
            sites,pairs,groups,seed=SEED+99991*i,
            draws=DRAWS,burnin=BURN,thin=THIN)
        if abs(out["source_matched_observed_line_gap"]-observed)>1e-10:
            raise AssertionError("Same original physical source pairs cannot change observed Jaccard")
        checks[family]=out
    return {"version":"v359b-source-authentic-27-anchors-7-sections-network-geo-fixed-plant-frequency-null",
        "original_Slovenia_publisher_DOI":"10.1007/s44473-026-00155-z",
        "original_publisher_source_sha256":SOURCE_ORIGINAL_SHA,
        "original_Overture_release":RELEASE,
        "original_Overture_source_status":railqa,
        "original_publisher_plant_site_n":len(sites),
        "original_publisher_sections_n":len(bysection),
        "original_global_physical_connector_source_rail_lines_in_region":len(segments),
        "original_true_railway_matched_site_n":nmatch,
        "original_ambiguous_source_site_n":namb,
        "original_geo_and_rail_network_matched_source_anchors":len(pairs),
        "original_author_sections_with_strict_matched_source_data":EXPECTED_GEO_NETWORK_AUTHOR_SECTIONS,
        "original_v359_matched_source_Jaccard_same_minus_other":observed,
        "source_frequency_nulls_same_original_27_anchor_triples":checks,
        "fixed_source_sampled_botanical_incidence_89_sites_not_27_only":True,
        "no_source_gps_raw_taxa_railway_shapes_named_source_line_ID_or_user_records_published":True,
        "technical_frequency_rank_is_not_randomized_field_p_value":True,
        "conclusions_cannot_establish_train_seed_transport_or_global_railway_network_homogenization":True}

def main():
    p=argparse.ArgumentParser()
    for name in ("original-xlsx","source-audit","output"):
        p.add_argument("--"+name,required=True)
    a=p.parse_args()
    original=extract_original_xlsx(a.original_xlsx,a.source_audit)
    obj=run(original)
    out=Path(a.output);out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(obj,ensure_ascii=False,indent=2),encoding="utf8")
    print(json.dumps(obj,ensure_ascii=False),flush=True)

if __name__=="__main__":main()
