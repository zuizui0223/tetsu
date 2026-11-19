#!/usr/bin/env python3
"""v29.8 source-agnostic strict joint conditioning on true tree clade AND
real original study source species occupancy-fifth.

NOT phylogenetic variance explained, no row/species/trait imputations, no
invented biological effect when source clade×occupancy strata cannot shuffle.

The occupation quintiles are assigned across the unchanged ORIGINAL source
taxa and its actual site observations. For MIREN ROAD versus INLAND the SAME
frequency bin and SAME randomized paired 2D height/seed species vector is
used in BOTH habitats. Phylogenetic clades come from the genuine first-party
V.PhyloMaker2 75%-depth tree, NOT taxonomic genus heuristics.
"""
from __future__ import annotations
from collections import Counter
import math
from typing import Any
import numpy as np

MIN_SOURCE_EXCHANGEABLE_JOINT_GROUPS=5
MIN_SOURCE_EXCHANGEABLE_ORIGINAL_INCIDENCE=.25
MIN_SOURCE_GEO_OR_ROAD_STRATA=2

def original_global_taxon_occupancy_quintiles(P:np.ndarray)->list[np.ndarray]:
    """Freeze source occupancy; use same tie-stable v29.6/v29.7 formula."""
    P=np.asarray(P)
    if P.ndim!=2 or P.shape[1]<10 or P.shape[0]<2:
        raise ValueError("Original source site×taxon matrix invalid for five frequency bins")
    if not np.all((P==0)|(P==1)):
        raise ValueError("Source occupancy must be original 0/1 incidence, not CWM matrix")
    order=np.argsort(P.sum(axis=0),kind="stable")
    return [x.astype(int) for x in np.array_split(order,5)]

def joint_true_clade_and_occupancy_groups(
    phylo_clades:list[np.ndarray],occupancy_quintiles:list[np.ndarray],
    n_original_taxa:int)->list[np.ndarray]:
    """Mathematical intersections partition original species *exactly*.

    Intersections with one species remain fixed in any trait label null.
    Do not merge sparse phylogenetic groups after seeing biological outcomes.
    """
    if not phylo_clades or len(occupancy_quintiles)!=5:
        raise ValueError("Source actual tree clades and all five fixed original frequency bins needed")
    if n_original_taxa<10:raise ValueError("Original taxon pool insufficient")
    def validate(partition):
        flat=[int(y) for group in partition for y in group]
        if len(flat)!=n_original_taxa or sorted(flat)!=list(range(n_original_taxa)):
            raise ValueError("Original source species are not exactly uniquely partitioned")
        if any(len(g)<1 for g in partition):raise ValueError("Empty source clade or occupancy bin")
    validate(phylo_clades)
    validate(occupancy_quintiles)
    strata=[]
    for clade in phylo_clades:
        s=set(map(int,clade))
        for occbin in occupancy_quintiles:
            matched=sorted(s.intersection(map(int,occbin)))
            if matched:strata.append(np.asarray(matched,dtype=int))
    validate(strata)
    return strata

def original_source_joint_exchangeability(
    P:np.ndarray,groups:list[np.ndarray],
    original_source_strata:np.ndarray,
    habitat_Ps:dict[str,np.ndarray]|None=None,
    original_road_strata:dict[str,list[np.ndarray]]|None=None)->dict:
    """Data support gate only; never view observed CWM effect signs."""
    P=np.asarray(P)
    if P.ndim!=2 or P.shape[1]<10:
        raise ValueError("Expected original site×source species incidence")
    n=P.shape[1]
    flat=sorted([int(i) for a in groups for i in a])
    if flat!=list(range(n)):raise ValueError("Original source phylogenetic frequency strata do not partition taxa")
    if original_source_strata.shape!=(P.shape[0],):
        raise ValueError("Original study source geography label per real site missing")
    mobile=np.asarray(sorted(set(int(i) for g in groups if len(g)>=2 for i in g)),dtype=int)
    frac=float(P[:,mobile].sum()/P.sum()) if len(mobile) and P.sum() else 0.
    source_group_fractions={}
    for field in sorted(set(original_source_strata)):
        sub=P[original_source_strata==field]
        source_group_fractions[str(field)]=(float(sub[:,mobile].sum()/sub.sum())
                    if len(mobile) and sub.sum() else 0.)
    habitat_fraction={}
    if habitat_Ps:
        for key,mat in habitat_Ps.items():
            q=np.asarray(mat)
            if q.ndim!=2 or q.shape[1]!=n or not np.all((q==0)|(q==1)):
                raise ValueError("Original paired road & inland incidence matrix invalid")
            habitat_fraction[key]=float(q[:,mobile].sum()/q.sum()) if len(mobile) and q.sum() else 0.
    road_pair_ok=0
    road_mobility={}
    if original_road_strata:
        for road_label, pair in original_road_strata.items():
            if len(pair)!=2:raise ValueError("Expected paired genuine field road-side and inland source mats")
            ra,rb=(np.asarray(x) for x in pair)
            vals=[]
            for q in (ra,rb):
                if q.ndim!=2 or q.shape[1]!=n or not np.all((q==0)|(q==1)):
                    raise ValueError("Original source physical roadside/inland sites not frozen")
                vals.append(float(q[:,mobile].sum()/q.sum()) if len(mobile) and q.sum() else 0.)
            road_mobility[str(road_label)]=vals
            road_pair_ok+=min(vals)>=MIN_SOURCE_EXCHANGEABLE_ORIGINAL_INCIDENCE
    multi=sum(len(g)>=2 for g in groups)
    n_geos=sum(x>=MIN_SOURCE_EXCHANGEABLE_ORIGINAL_INCIDENCE
               for x in source_group_fractions.values())
    allowed=(multi>=MIN_SOURCE_EXCHANGEABLE_JOINT_GROUPS and
             frac>=MIN_SOURCE_EXCHANGEABLE_ORIGINAL_INCIDENCE and
             n_geos>=MIN_SOURCE_GEO_OR_ROAD_STRATA and
             all(x>=MIN_SOURCE_EXCHANGEABLE_ORIGINAL_INCIDENCE
                 for x in habitat_fraction.values()) and
             (not original_road_strata or
              road_pair_ok>=MIN_SOURCE_GEO_OR_ROAD_STRATA))
    return {
        "n_original_phylo_clade_times_occupancy_intersection_strata":len(groups),
        "n_source_phylo_clade_x_occupancy_strata_size_ge2":multi,
        "n_original_trait_known_source_taxa_in_exchangeable_joint_groups":len(mobile),
        "fraction_original_source_site_species_incidence_exchangeable":frac,
        "n_original_source_geographical_field_groups_with_ge25pct_mobility":n_geos,
        "source_geographical_group_incidence_exchangeability":source_group_fractions,
        "source_habitat_exchangeable_incidence_fraction":habitat_fraction,
        "n_original_roadside_and_inland_field_roads_with_ge25pct_mobility":road_pair_ok,
        "source_original_road_pair_mobile_fraction":road_mobility,
        "predeclared_min_multispecies_joint_strata":MIN_SOURCE_EXCHANGEABLE_JOINT_GROUPS,
        "predeclared_min_original_site_taxon_incidence_fraction":MIN_SOURCE_EXCHANGEABLE_ORIGINAL_INCIDENCE,
        "predeclared_min_original_geographies_or_paired_roads":MIN_SOURCE_GEO_OR_ROAD_STRATA,
        "joint_phylogenetic_clade_AND_source_occupancy_null_identifiable":bool(allowed),
        "support_status":(
          "IDENTIFIABLE_PHYLO75_AND_OCCUPANCY_JOINT_LABEL_NULL" if allowed
          else "NOT_IDENTIFIABLE_PHYLO_OCCUPANCY_JOINT_LABEL_NULL"),
        "did_not_change_original_species_site_membership_or_traits":True,
        "technical_label_swap_cannot_be_called_independent_field_replicate":True,
    }

def occupation_group_source_labels(assignment:dict[str,int],groups:list[np.ndarray],names:list[str])->dict:
    """Optional zero-values-free audit only; never export names in source QA."""
    if len(set(names))!=len(names):raise ValueError("Original botanical species name duplicates")
    return {"source_total_original_taxa":len(names),
            "total_plant_species_partitioned":
                sum(len(g) for g in groups),
            "n_original_effectively_exchangeable_strata":
                sum(len(g)>=2 for g in groups)}
