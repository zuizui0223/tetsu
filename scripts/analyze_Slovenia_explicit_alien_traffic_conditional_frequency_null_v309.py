#!/usr/bin/env python3
"""v309: Is explicit-ALIEN source shared-species freight beta above strict year×layout null?

Exploratory after v305/v306. This computes the exact Jaccard attribution
to publisher-native-status explicit alien taxa and several preselected guilds.
The outcome is a source-identified same-taxon portion of *all-taxon Jaccard*,
not an alien-only Jaccard (which changes the original model denominator).
Technical constrained-checkerboard draws are not independent ecological tests.
"""
from __future__ import annotations
import argparse, hashlib, json, math, random
from collections import Counter, defaultdict
from pathlib import Path
import numpy as np
from analyze_Slovenia_rail_microhabitats_v103 import extract_original_xlsx
from analyze_Slovenia_historical_train_frequency_local_species_coupling_v275 import mapped_segment_data
from analyze_Slovenia_traffic_shared_species_v306 import (
    site_cohort,site_pair_dyads,traits_and_names,zoneset,
    metric_and_contributions,exposure_design,fit_slope,
    swap_communities,finite_original_within_line_traffic_labels,
    SEEDS,DRAWS,BURN,SPACING,ZONES,
)
from analyze_Slovenia_rail_route_nearest_plant_microhabitat_v277 import jaccard_similarity

BASELINE_SHA="5d9cb0ba983d5718c2b25ad2498484e337087ae00f58a1d9e10fe99860725a1e"
FROZEN_ORIGINAL_SHA="26be7ca9513dd54136d4ec499d8413c1cb288045f87fe3cf8814731f2edac99f"
GROUPS=("all","explicit_alien","not_explicit_alien","strict_annual",
        "wind_Epilobium_Phragmites","Bidens_archetype","missing_nativeness")
MODES=("source_year_track_layout","global")
TOL=1e-10

def make_guild_sets(original_sets:list[set],traits:dict)->dict[str,set]:
    source=set().union(*original_sets)
    if not source:raise ValueError("Empty genuinely observed original botanical taxa")
    alien={t for t in source if isinstance(traits.get(t,{}).get("Nativeness.status"),str)
           and "alien" in traits[t]["Nativeness.status"].lower()}
    annual={t for t in source if traits.get(t,{}).get("Lifespan.categories")=="annual"}
    wind={t for t in source if traits.get(t,{}).get("Diaspore.strategy") in ("Epilobium","Phragmites")}
    bidens={t for t in source if traits.get(t,{}).get("Diaspore.strategy")=="Bidens"}
    missing={t for t in source if traits.get(t,{}).get("Nativeness.status") is None}
    return {
        "all":source,
        "explicit_alien":alien,
        "not_explicit_alien":source-alien,
        "strict_annual":annual,
        "wind_Epilobium_Phragmites":wind,
        "Bidens_archetype":bidens,
        "missing_nativeness":missing
    }

def dyad_contribution_vectors(site_sets:list[set],dyads:list[dict],guildsets:dict)->np.ndarray:
    """All taxa weighted by |site-pair species union|, never group-only Jaccard.

    Vector columns are in predeclared GROUPS order and share identical
    geography/traffic support. all == explicit_alien + not_explicit_alien.
    """
    data=np.zeros((len(dyads),len(GROUPS)),float)
    masksets=[guildsets[g] for g in GROUPS]
    for i,row in enumerate(dyads):
        for a,b in row["members"]:
            common=site_sets[a]&site_sets[b]
            all_species=site_sets[a]|site_sets[b]
            if not all_species:raise ValueError("Real source site pair lacks plant taxa")
            fac=1/(len(all_species)*len(row["members"]))
            data[i,:]+=np.array([len(common&mask) for mask in masksets],float)*fac
    if not np.allclose(data[:,0],data[:,1]+data[:,2],atol=1e-12,rtol=0):
        raise ValueError("Original total Jaccard not partitioned by alien source status")
    return data

def observed_and_technical_null(
    sites:list[dict],dyads:list[dict],zone:str,traits:dict,baseline_zone:dict,
    original_groups:dict,keys:list[str],draws:int=DRAWS
)->dict:
    if zone not in ZONES or draws!=DRAWS:raise ValueError("Deviation from frozen source null")
    original=[zoneset(s,zone) for s in sites]
    group_sets=make_guild_sets(original,traits)
    observed_vec=dyad_contribution_vectors(original,dyads,group_sets)
    X=exposure_design(dyads)
    if np.linalg.matrix_rank(X)!=3:raise ValueError("Original dyad traffic rank invalid")
    w=np.linalg.pinv(X)[1,:] # exact linear derivative d(beta_freight)/d(site-pair Jaccard)
    obs=w@observed_vec
    if abs(obs[0]-baseline_zone["no_urban_observed_freight_slope"])>TOL:
        raise ValueError("Original no-urban v306 real dyad freight regression differs")
    lookup={
        "explicit_alien":"explicit_alien",
        "not_explicit_alien":"not_explicitly_alien_INCLUDES_native_unknown_crop",
        "strict_annual":"strict_annual",
        "wind_Epilobium_Phragmites":"author_wind_Epilobium_Phragmites",
        "Bidens_archetype":"author_Bidens_not_measured_rail_attachment",
        "missing_nativeness":"source_missing_nativeness",
    }
    for g,k in lookup.items():
        j=GROUPS.index(g)
        v=baseline_zone["source_defined_guilds"][k]["signed_contribution_sum"]
        if abs(float(obs[j])-v)>TOL:
            raise ValueError(f"Genuine publisher v306 source group {zone}/{g} drift: {obs[j]} vs {v}")
    all_taxa=sorted(group_sets["all"])
    source_y, source_taxon_C=metric_and_contributions(original,dyads,all_taxa)
    if not np.allclose(source_y,observed_vec[:,0],rtol=0,atol=1e-12):
        raise ValueError("Source taxon-level actual Jaccard does not equal group attribution")
    null={}
    for mode in MODES:
        if mode=="global":
            strata=[list(range(len(sites)))]
        else:
            d=defaultdict(list)
            for i,s in enumerate(sites):d[(s["year"],s["layout"])].append(i)
            strata=[x for x in d.values() if len(x)>=2]
        original_size=[len(x) for x in original]
        original_freq=Counter(x for row in original for x in row)
        orig_freq_by_stratum={str(k):Counter(x for i in inds for x in original[i]) for k,inds in enumerate(strata)}
        results=[]
        chains=[]
        for chain,seed in enumerate(SEEDS):
            m=[set(r) for r in original]
            rand=random.Random(seed+101*ZONES.index(zone)+(1 if mode!="global" else 0))
            accepted=swap_communities(m,rand,BURN,strata)
            partial=[]
            for _ in range(draws):
                accepted+=swap_communities(m,rand,SPACING,strata)
                if [len(x) for x in m]!=original_size:
                    raise ValueError("Original individual source site richness changed under null")
                if Counter(x for row in m for x in row)!=original_freq:
                    raise ValueError("Original source taxon incidence frequency changed under null")
                if mode!="global":
                    for k,inds in enumerate(strata):
                        if Counter(x for i in inds for x in m[i])!=orig_freq_by_stratum[str(k)]:
                            raise ValueError("Conditional source-year x track-layout taxon occupancy changed")
                q=w@dyad_contribution_vectors(m,dyads,group_sets)
                if abs(q[0]-q[1]-q[2])>TOL:
                    raise ValueError("Draw group contributions do not sum")
                results.append(q);partial.append(q)
            chains.append({
                "chain":chain+1,
                "accepted_source_incidence_checkerboard_fraction":accepted/(BURN+draws*SPACING),
                "mean_source_freight_slope_by_group":dict(zip(GROUPS,map(float,np.mean(partial,axis=0))))
            })
        vals=np.asarray(results,float)
        if vals.shape!=(2*draws,len(GROUPS)):
            raise ValueError("Original technical null draw total failed")
        parts={}
        for i,g in enumerate(GROUPS):
            z=vals[:,i]
            score=float(obs[i])
            parts[g]={
                "actual_source_original_freight_coefficient":score,
                "original_89_site_x_source_taxon_occurrence_NULL_mean":float(np.mean(z)),
                "actual_minus_source_NULL":float(score-np.mean(z)),
                "source_technical_null_2p5_97p5_NOT_population_CI":list(map(float,np.quantile(z,[.025,.975]))),
                "source_technical_upper_rank_NOT_independent_field_p":float((1+sum(z>=score))/(len(z)+1)),
                "source_technical_lower_rank_NOT_independent_field_p":float((1+sum(z<=score))/(len(z)+1)),
            }
        if abs(parts["all"]["actual_minus_source_NULL"]-
               (parts["explicit_alien"]["actual_minus_source_NULL"]+
                parts["not_explicit_alien"]["actual_minus_source_NULL"]))>TOL:
            raise ValueError("Source group NULL excess does not decompose original all-species NULL excess")
        return_null={
            "source_null_mode":mode,
            "n_correlated_technical_source_draws_NOT_field_studies":len(vals),
            "fixed_all_89_site_source_richness":True,
            "fixed_original_taxon_89_site_occupancy":True,
            "fixed_taxon_source_year_x_track_layout_frequency":mode!="global",
            "source_guild_coefficient_results":parts,
            "chain_diagnostics":chains
        }
        null[mode]=return_null
    historical=finite_original_within_line_traffic_labels(
        keys,original_groups,dyads,observed_vec[:,1])
    historical_all=finite_original_within_line_traffic_labels(
        keys,original_groups,dyads,observed_vec[:,0])
    return {
        "zone":zone,
        "n_all_source_taxa":len(group_sets["all"]),
        "n_source_taxa_by_group":{g:len(v) for g,v in group_sets.items()},
        "source_decomposition_uses_full_flora_union_in_denominator":True,
        "real_source_shared_species_group_freight_slopes":dict(zip(GROUPS,map(float,obs))),
        "original_v306_source_and_technical_full_taxa_exact_replay":True,
        "source_2020_official_48_within_line_freight_reassignments_EXPLICIT_ALIEN":historical,
        "source_2020_official_48_within_line_freight_reassignments_ALL_SPECIES":historical_all,
        "source_technical_taxon_incidence_margin_nulls":null,
    }

def main():
    parser=argparse.ArgumentParser()
    for name in ("source-xlsx","source-audit","station-pdf","v306-source","output"):
        parser.add_argument("--"+name,required=True,type=Path)
    a=parser.parse_args()
    if hashlib.sha256(a.source_xlsx.read_bytes()).hexdigest()!=FROZEN_ORIGINAL_SHA:
        raise ValueError("True CC-BY source botanist original SHA changed")
    if hashlib.sha256(a.v306_source.read_bytes()).hexdigest()!=BASELINE_SHA:
        raise ValueError("Replayed frozen v306 source result artifact SHA changed")
    original=extract_original_xlsx(a.source_xlsx,a.source_audit)
    traffic,tqa=mapped_segment_data(original,a.station_pdf)
    sites,keys,original_groups,sqa=site_cohort(original,traffic)
    traits,codebook,traitqa=traits_and_names(original)
    dyads=site_pair_dyads(sites,keys,original_groups)
    v306=json.loads(a.v306_source.read_text(encoding='utf8'))
    if (v306["n_true_original_botanical_sites"],v306["n_mapped_traffic_sites"],
        v306["n_true_traffic_segments"],v306["n_DEpendent_segment_pairs_not_independent"])!=(89,27,9,36):
        raise ValueError("Genuine original v306 source biological sample changed")
    out={
        "version":"v309",
        "design":"EXPLORATORY after v306; exact source-trait fraction of same-species Jaccard traffic effect; no urban covariate in this model",
        "source_original_Slovenia_author_botanical_sha256":FROZEN_ORIGINAL_SHA,
        "previous_v306_original_source_result_SHA256":BASELINE_SHA,
        "n_original_botanical_sites":len(sites),
        "n_sites_with_true_official_traffic":sum(map(len,original_groups.values())),
        "n_true_national_2020_traffic_segments":len(keys),
        "n_source_original_national_lines_with_true_traffic":len({x[0][1]["section"] for x in original_groups.values()}),
        "n_dependent_dyads_NOT_independent":len(dyads),
        "true_original_site_x_species_source_data_NOT_republished":True,
        "v306_original_source_traffic_QA":tqa,
        "zones":{}
    }
    for zone in ZONES:
        out["zones"][zone]=observed_and_technical_null(
            sites,dyads,zone,traits,v306["zones"][zone],original_groups,keys)
    target=a.output
    target.parent.mkdir(parents=True,exist_ok=True)
    target.write_text(json.dumps(out,indent=2,ensure_ascii=False,allow_nan=False)+"\n",encoding="utf8")
    print(json.dumps({
        "status":"v309_REAL_AUTHENTIC_SOURCE_COMPLETE",
        "true_stations":len(sites),
        "source_n_traffic_segments":len(keys),
        "results":{z:{
            "real_alien_freight_slope":x["real_source_shared_species_group_freight_slopes"]["explicit_alien"],
            "strict_original_year_x_track_freq_null_alien_excess":x["source_technical_taxon_incidence_margin_nulls"]["source_year_track_layout"]["source_guild_coefficient_results"]["explicit_alien"]["actual_minus_source_NULL"],
            "strict_null_alien_technical_upper_rank":x["source_technical_taxon_incidence_margin_nulls"]["source_year_track_layout"]["source_guild_coefficient_results"]["explicit_alien"]["source_technical_upper_rank_NOT_independent_field_p"],
            "48_original_within_line_traffic_label_positive_rank":x["source_2020_official_48_within_line_freight_reassignments_EXPLICIT_ALIEN"]["technical_one_sided_positive_rank_NOT_field_p"]
        } for z,x in out["zones"].items()}
    },ensure_ascii=False,indent=2))

if __name__=="__main__":main()
