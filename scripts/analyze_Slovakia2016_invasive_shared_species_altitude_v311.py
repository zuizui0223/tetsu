#!/usr/bin/env python3
"""v311 original Slovak 2016 author-listed invasive taxa spatial sharedness.

Exploratory after v309/v310. Author source has 318 taxa/19 stations and
names 14 of 15 reported invasive taxa in prose; we NEVER invent the 15th.
For dyad i,j, invasive contribution = shared original named invader taxa /
ALL observed taxon union. Source exact split reproduces all-taxa Jaccard.
No train count exists in original author dataset.
"""
from __future__ import annotations
import argparse, hashlib, json, math, random, re
from collections import Counter,defaultdict
from pathlib import Path
import numpy as np

from audit_Slovakia2016_original_19station_318taxa_pdf_preflight_v310 import preflight
from analyze_Slovakia2016_19station_original_taxon_incidence_v310 import (
    source_list_clauses,parse_clauses,validate_original_source,station_sets,
    beta_slope,swap_states,jaccard,BURN,THIN,SEEDS
)
from analyze_Slovakia2016_19station_original_physical_distance_v310 import (
    parse_source_dms,geodist
)
import fitz

PUBLISHED_NAMED_INVASIVES=(
    "Amaranthus retroflexus",
    "Conyza canadensis",
    "Galinsoga parviflora",
    "Galinsoga urticifolia",
    "Matricaria discoidea",
    "Bidens frondosa",
    "Epilobium ciliatum",
    "Impatiens parviflora",
    "Negundo aceroides",
    "Rumex patientia",
    "Solidago canadensis",
    "Solidago gigantea",
    "Apera spica-venti",
    "Echinochloa crus-galli",
)
N_DECLARED=14
EXPECTED_ORIGINAL_GEO_BETA=-0.0020021839482198068
SEEDS_311=(311681,311682)
DRAWS=199
BURNIN=18000
SPACING=650
GROUPS=("all","author_named_invasive_14","others_not_claimed_native")

def match_original_author_names(observed:dict)->dict[str,str]:
    """Only exact author binomial or source scientific entry with distinct suffix.

    An unmatched species is NOT completed by GBIF taxonomic guessing.
    """
    if len(PUBLISHED_NAMED_INVASIVES)!=N_DECLARED:
        raise ValueError("Original author's exact list did not retain 14 predeclared names")
    mapping={}
    for desired in PUBLISHED_NAMED_INVASIVES:
        hits=[x for x in observed if x==desired or x.startswith(desired+" ")]
        if len(hits)!=1:
            raise ValueError("Strict Slovak original source taxon name cannot be uniquely matched: "+
                            desired+" candidates="+str(hits))
        mapping[desired]=hits[0]
    if len(set(mapping.values()))!=N_DECLARED:
        raise ValueError("Two source invasive names map to same author-original taxon")
    return mapping

def source_covariates(coords:dict):
    k=sorted(coords)
    if k!=list(range(1,20)):raise ValueError("Original station IDs changed")
    pairs=[(i,j) for i in k for j in k if i<j]
    km=np.array([geodist(coords[i],coords[j]) for i,j in pairs],float)
    alt=np.array([abs(coords[i][2]-coords[j][2])/100 for i,j in pairs],float)
    meanalt=np.array([(coords[i][2]+coords[j][2])/200 for i,j in pairs],float)
    if len(km)!=171 or min(km)<0.1 or max(km)>100:
        raise ValueError("True original 19 field stations have invalid geographic support")
    return pairs,km,alt,meanalt

def dyad_components(site_sets:dict,source_invasives:set,pairs:list)->np.ndarray:
    rows=np.zeros((len(pairs),3),float)
    for n,(i,j) in enumerate(pairs):
        a,b=site_sets[i],site_sets[j]
        u=a|b
        if not u:raise ValueError("Original real floristic station pair union empty")
        shared=a&b
        inv=len(shared&source_invasives)/len(u)
        others=(len(shared)-len(shared&source_invasives))/len(u)
        rows[n,:]=[len(shared)/len(u),inv,others]
    if not np.allclose(rows[:,0],rows[:,1]+rows[:,2],atol=1e-12,rtol=0):
        raise ValueError("Source invaded/not-explicit-invaded decomposition invalid")
    return rows

def design(km,alt=None,meanalt=None):
    predictors=[km]
    if alt is not None:predictors.append(alt)
    if meanalt is not None:predictors.append(meanalt)
    return np.column_stack([np.ones(len(km))]+predictors)

def fit(X,C):
    if len(X)!=len(C) or np.linalg.matrix_rank(X)!=X.shape[1]:
        raise ValueError("Original station dyadic predictor design rank deficient")
    co=np.linalg.lstsq(X,C,rcond=None)[0]
    if max(abs(co[:,0]-co[:,1]-co[:,2]))>1e-9:
        raise ValueError("Ecological group coefficients do not add up")
    return co

def standard_condition(km,alt,meanalt=None):
    X=np.column_stack([km,alt]+([meanalt] if meanalt is not None else []))
    X=(X-X.mean(axis=0))/X.std(axis=0)
    mat=np.column_stack([np.ones(len(km)),X])
    return float(np.linalg.cond(mat)),float(np.corrcoef(km,alt)[0,1])

def source_elevation_tertiles(coords):
    by=sorted(range(1,20),key=lambda i:(-coords[i][2],i))
    groups=[by[:6],by[6:12],by[12:]]
    if sorted(x for g in groups for x in g)!=list(range(1,20)) or list(map(len,groups))!=[6,6,7]:
        raise ValueError("Original elevation tertiles need original 19 observed sites")
    return groups

def original_traffic_free_analysis(pdf:Path,draws=DRAWS):
    audit=preflight(pdf)
    if audit["pdf_sha256"]!="99465bdd7ed56985e65b1163f16eca4f702fc821619d84c89bf986d4d3736ae1":
        raise ValueError("2016 original author paper source SHA drift")
    clauses,_=source_list_clauses(pdf)
    obs,_=parse_clauses(clauses)
    qa=validate_original_source(obs)
    source_g=match_original_author_names(obs)
    invasive=set(source_g.values())
    original=station_sets(obs)
    d=fitz.open(str(pdf))
    coords=parse_source_dms("\n".join(d[i].get_text("text",sort=True) for i in (1,2)))
    d.close()
    pairs,km,alt,meanalt=source_covariates(coords)
    C=dyad_components(original,invasive,pairs)
    X=design(km)
    joint=design(km,alt)
    expanded=design(km,alt,meanalt)
    raw=fit(X,C)
    if abs(raw[1,0]-EXPECTED_ORIGINAL_GEO_BETA)>1e-10:
        raise ValueError("Source v310 genuinely observed Jaccard geographical beta not reproduced")
    fitjoint=fit(joint,C)
    fitfull=fit(expanded,C)
    cond,corr=standard_condition(km,alt)
    condfull,_=standard_condition(km,alt,meanalt)
    loo=[]
    for excluded in range(1,20):
        keep=[ind for ind,(i,j) in enumerate(pairs) if excluded not in (i,j)]
        if len(keep)!=153:raise ValueError("Whole source station deletion did not remove 18 pairs")
        b=fit(joint[keep],C[keep])[1,:]
        loo.append({"excluded":excluded,"distance_slope_all":float(b[0]),
                    "distance_slope_invasive":float(b[1])})
    orig_freq=Counter(t for xx in original.values() for t in xx)
    orig_rich={k:len(x) for k,x in original.items()}
    mode_results={}
    for mode in ("all_19_stations","source_elevation_tertiles"):
        strata=[list(range(1,20))] if mode=="all_19_stations" else source_elevation_tertiles(coords)
        bystratum=[Counter(t for i in g for t in original[i]) for g in strata]
        nulldistance=[]
        nulljoint=[]
        chain=[]
        for chain_id,seed in enumerate(SEEDS_311):
            rr=random.Random(seed+(0 if mode=="all_19_stations" else 10000))
            work={k:set(v) for k,v in original.items()}
            accepted=0
            def advance(n):
                nonlocal accepted
                for _ in range(n):
                    g=strata[rr.randrange(len(strata))]
                    i,j=rr.sample(g,2)
                    a,b=work[i]-work[j],work[j]-work[i]
                    if not a or not b:continue
                    a1=rr.choice(sorted(a));b1=rr.choice(sorted(b))
                    work[i].remove(a1);work[i].add(b1)
                    work[j].remove(b1);work[j].add(a1)
                    accepted+=1
            advance(BURNIN)
            partial=[]
            for step in range(draws):
                advance(SPACING)
                if Counter(t for st in work.values() for t in st)!=orig_freq:
                    raise ValueError("Taxon original total 19-site prevalence not fixed")
                if any(len(work[k])!=n for k,n in orig_rich.items()):
                    raise ValueError("Original author field station taxa richness not fixed")
                if mode!="all_19_stations":
                    for g,prior in zip(strata,bystratum):
                        if Counter(t for i in g for t in work[i])!=prior:
                            raise ValueError("Source species elevation-tertile within-stratum prevalences changed")
                z=dyad_components(work,invasive,pairs)
                bd=fit(X,z)[1,:]
                bj=fit(joint,z)[1,:]
                nulldistance.append(bd)
                nulljoint.append(bj)
                partial.append(bd)
            chain.append({"chain":chain_id+1,"accepted_fraction":accepted/(BURNIN+draws*SPACING),
                          "mean_jaccard_distance_null":[float(x) for x in np.mean(partial,axis=0)]})
        null1=np.array(nulldistance)
        null2=np.array(nulljoint)
        if null1.shape!=(2*draws,3) or null2.shape!=(2*draws,3):
            raise ValueError("398 null draws not present")
        slopes={}
        for n,g in enumerate(GROUPS):
            obsval=float(raw[1,n])
            v=null1[:,n]
            if not np.isfinite(v).all():
                raise ValueError("Original constrained null contains nonfinite slope")
            slopes[g]={
                "observed_geodesic_per_km_slope":obsval,
                "conditional_source_null_mean":float(v.mean()),
                "observed_minus_conditional_null_mean":obsval-float(v.mean()),
                "source_null_2p5_97p5_technical_NOT_ECOLOGICAL_CI":list(map(float,np.quantile(v,[.025,.975]))),
                "technical_lower_rank_NOT_FIELD_p":float((1+sum(v<=obsval))/(len(v)+1)),
                "technical_upper_rank_NOT_FIELD_p":float((1+sum(v>=obsval))/(len(v)+1)),
                "km_slope_AFTER_adjusting_original_altitude_gap":float(fitjoint[1,n]),
                "null_altitude_adjusted_km_slope_mean":float(null2[:,n].mean()),
                "km_slope_adjusted_altitude_observed_minus_null":float(fitjoint[1,n]-null2[:,n].mean())
            }
        if abs(slopes["all"]["observed_minus_conditional_null_mean"]-
               slopes["author_named_invasive_14"]["observed_minus_conditional_null_mean"]-
               slopes["others_not_claimed_native"]["observed_minus_conditional_null_mean"])>1e-10:
            raise ValueError("True group source taxa frequency null exact additive property not preserved")
        mode_results[mode]={
            "strata_station_counts":[len(g) for g in strata],
            "preserves_original_field_site_richness":True,
            "preserves_original_taxon_incidence_frequency_within_elevation_tertiles":mode=="source_elevation_tertiles",
            "n_correlated_technical_null_draws_NOT_independent_field_studies":len(null1),
            "original_source_taxon_group_spatial_slopes":slopes,
            "chain_source_margins_and_acceptance":chain
        }
    return {
        "version":"v311", "source_original_author_2016_PDF_SHA":audit["pdf_sha256"],
        "source_source_historical_railway_plant_field_units":qa,
        "n_source_author_invasive_names_explicitly_listed":len(source_g),
        "original_author_reports_15_total_invaders_but_14_named":True,
        "named_invasives_not_ALL_98_original_alien_species":True,
        "source_author_invasive_names_exact_original_matched":list(source_g),
        "source_country":"Slovakia", "n_source_stations":19,
        "n_source_station_pairs_which_are_dependent":len(pairs),
        "source_pair_geodesic_km_minmax":[float(km.min()),float(km.max())],
        "source_pair_absolute_elevation_gap_m_minmax":[float(alt.min()*100),float(alt.max()*100)],
        "source_geographic_distance_and_absolute_altitude_difference_correlation":corr,
        "standardized_predictor_condition_number_distance_plus_altitude":cond,
        "standardized_predictor_condition_number_distance_plus_altitude_plus_mean_altitude":condfull,
        "source_altitude_adjusted_distance_coefficients_by_group":dict(zip(GROUPS,map(float,fitjoint[1,:]))),
        "source_elevation_gap_coefficients_by_group":dict(zip(GROUPS,map(float,fitjoint[2,:]))),
        "source_mean_altitude_sensitivity_distance_coefficients_by_group":dict(zip(GROUPS,map(float,fitfull[1,:]))),
        "source_original_jaccard_mean":float(np.mean(C[:,0])),
        "observed_source_distance_only_coefficients_by_group":dict(zip(GROUPS,map(float,raw[1,:]))),
        "n_original_stations_with_altitude_adjusted_all_Jaccard_distance_slope_negative":sum(x["distance_slope_all"]<0 for x in loo),
        "n_original_stations_with_altitude_adjusted_invasive_sharedness_distance_slope_positive":sum(x["distance_slope_invasive"]>0 for x in loo),
        "source_19_whole_station_leave_one_out_altitude_adjusted_distance_all_minmax":[min(x["distance_slope_all"] for x in loo),max(x["distance_slope_all"] for x in loo)],
        "source_19_whole_station_leave_one_out_altitude_adjusted_distance_invasive_minmax":[min(x["distance_slope_invasive"] for x in loo),max(x["distance_slope_invasive"] for x in loo)],
        "null_modes":mode_results,
        "v309_2015_2020_Slovenia_train_volume_NOT_tested_here":True,
        "actual_2012_2014_passenger_freight_traffic_source_records":0,
        "field_plot_nonrail_native_reference":False,
        "no_railway_transport_or_causal_management_claims":True,
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--pdf",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    a=p.parse_args()
    out=original_traffic_free_analysis(a.pdf)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(out,ensure_ascii=False,indent=2,allow_nan=False)+"\n",encoding="utf8")
    print(json.dumps({
        "source_original_taxa":out["source_source_historical_railway_plant_field_units"]["n_pub_source_total_taxa"],
        "source_invader_taxa_matched":out["n_source_author_invasive_names_explicitly_listed"],
        "true_station_pairs_not_independent":out["n_source_station_pairs_which_are_dependent"],
        "geo_vs_alt_pair_correlation":out["source_geographic_distance_and_absolute_altitude_difference_correlation"],
        "altitude_adjusted_geo_slope_total":out["source_altitude_adjusted_distance_coefficients_by_group"]["all"],
        "altitude_adjusted_geo_slope_invaders":out["source_altitude_adjusted_distance_coefficients_by_group"]["author_named_invasive_14"],
        "invader_source_frequency_excess":{k:v["original_source_taxon_group_spatial_slopes"]["author_named_invasive_14"]["observed_minus_conditional_null_mean"] for k,v in out["null_modes"].items()},
    },ensure_ascii=False,indent=2))
if __name__=="__main__":main()
