#!/usr/bin/env python3
"""v28.4 — original CZECH 35 railway stations × Pladias dispersal/lifeform
on fixed geodistance-matched SAME author stretch vs DIFFERENT stretch pairs.

Source Kutlvašr et al (2024) Preslia original S1/S2, 2020–22, actual
39 railway stations / 5 source study stretches / 761 recoverable current
taxa. Exactly 35 whole-station vegetation assemblages match a THIRD-PARTY
2026 registry centroid by full station-name equality. Centroids are NOT
2020–22 original plant-plot GPS. Pladias dispersal strategy & Raunkiaer
life form matched to 544 of 761 actual taxa (5081/6941 site observations
across original 39; exact 35-site coverage computed in this audit).

Source Sádlo 'Epilobium' and 'Phragmites' strategy archetypes are used
as a WIND-ASSOCIATED DISPERSAL STRATEGY PROXY, NOT verified measured
individual dispersal distance/vehicle attachment/aerodynamic speed.
All other strategy archetypes comprise 'other'. Original therophyte is
a life-form PROXY, not observed field lifespan or resistance to mowing.
Missing traits omitted from all proxy groups, NEVER non-wind-imputed.

Same original source physical station anchor i; another station j from
same author 5 original stretches and k from DIFFERENT source stretch,
chosen only using original name-matched WGS84 registry station-centroid
geography, predeclared <=160 km and source |dij-dik| <=10 km.
All comparisons reuse the EXACT 33 source anchor triples from v28.2;
species or trait identities are not used to choose candidate sites.
Main primary contrast:
  [same−different plant Jaccard for wind-associated taxa] −
  [same−different Jaccard for all other known-strategy taxa]
on identical anchors where both taxa groups have defined Jaccard.
Four supplementary source guilds:
  therophyte×wind, therophyte×other, nontherophyte×wind,
  nontherophyte×other.
Avoid pseudo-independent same/other pair counts: anchor mean ->
true source STRETCH mean -> study total equal five stretch weights.
Gate each guild >=4 source stretches and >=12 anchors else NOT
IDENTIFIABLE. Score 544-trait-known taxon-only group as reference,
and 761 all-source taxa separately as pre-registered comparator.

For each source guild independently, 2×2 checkerboard source
randomizations preserve EXACT each source station×guild species richness
and each source taxon's station occupancy frequency. Two nulls:
  - all original 35 stations global source species prevalence;
  - original 5 source stretches' species occupancy also fixed.
Latter may OVERCONTROL true line-selective dispersal, so report as
a sensitivity not a causal negative. Two independent seeded chains,
79 draws/chain, correlated Markov samples are not independent
field experimental p values. Each group defined from true published
source, not post hoc positive species labels. No source raw data,
site names/centroids, source taxon-by-station rows or guild taxonomy
rows republished; output only aggregate figures and source QA counts.

Do not claim natural seed dispersal, actual train seed carriage,
freight network flow, mowing/wind adaptive selection, or heritability.
"""
from __future__ import annotations

import argparse,json,math,random
from collections import Counter,defaultdict
from pathlib import Path
from statistics import mean

from analyze_Czech_2024_real_35station_centroid_geography_identity_v154 import (
    original_35_station_sites
)
from analyze_Czech_2024_39station_functional_identity_v143 import (
    GROUPS
)
from audit_Czech_2024_true_35_stations_stretch_geodistance_matched_v282 import (
    candidates,aggregate
)
from audit_RAIL_ROAD_source_frequency_fixed_reference_null_v270 import (
    swaps,source_freq
)

STRATEGIES=("wind_strategy_proxy","other_strategy_proxy")
GUILDS=tuple(GROUPS)
LABELS=("all_544_both_traits",)+STRATEGIES+GUILDS
MAIN="wind_minus_other_source_section_similarity_gap"
NULLS=("all_35_station_source_frequency",
       "within_original_5_author_stretches_source_frequency")
NCHAINS=2
DRAWS=79
BURNIN=2000
THIN=280
SEED=284118
MIN_STRETCHES=4
MIN_ANCHORS=12
EXPECTED_MATCHES=33
EXPECTED_ALL_MATCHED_TRAITS=544


def original_station_traits(s1:Path,s2:Path,disp:Path,life:Path,sr70:Path)->tuple[list[dict],dict,dict,list[dict],dict]:
    sites,guild,qa=original_35_station_sites(s1,s2,disp,life,sr70)
    if len(sites)!=35 or len(guild)!=EXPECTED_ALL_MATCHED_TRAITS:
        raise ValueError("Original Czech 35 located 544 exact published Pladias matched taxa drift")
    geom=[{
        "source_stretch":s["section"],"lat":s["lat"],"lon":s["lon"],
        "area_ha":s["area_ha"],"taxa":s["all_original_taxa"],
    } for s in sites]
    triads,match_qa=candidates(geom,160.,10.,None)
    support=aggregate(triads,match_qa)
    if not support["identified"] or len(triads)!=EXPECTED_MATCHES or support["source_original_stretches_with_real_geographic_match"]!=5:
        raise ValueError("Prior independent Czech v28.2 original 33 station source geographic triads changed")
    # Numeric source trait coverage should be recovered without imputing
    # any missing traits, even if unknown 4 stations were omitted.
    cover={
        "source_unique_observed_2020_22_35site_taxa":len(set.union(*(s["all_original_taxa"] for s in sites))),
        "source_35site_trait_known_incidence":sum(len(s["known"]) for s in sites),
        "source_35site_all_actual_incidence":sum(len(s["all_original_taxa"]) for s in sites),
        "source_35site_original_plant_lifeform_dispersal_trait_coverage_fraction":
            sum(len(s["known"]) for s in sites)/sum(len(s["all_original_taxa"]) for s in sites),
        "original_35_sites_known_taxa_min":min(len(s["known"]) for s in sites),
        "Pladias_Sadlo_wind_archetype_is_not_empirical_seed_dispersal_rate":True,
        "Raunkiaer_therophyte_not_strict_annual_lifespan":True,
        "original_source_study_stretches_not_active_train_OD":True
    }
    return sites,guild,qa,triads,cover


def guild_source_matrix(sites:list[dict],guild:dict)->dict:
    result={label:[set() for _ in sites] for label in GUILDS}
    for i,s in enumerate(sites):
        if not set(s["known"])<=set(guild):
            raise ValueError("Original published Pladias unknown source taxa assigned to guild")
        for taxon in s["known"]:
            result[guild[taxon]][i].add(taxon)
    if any(set.union(*(result[g][i] for g in GUILDS))!=s["known"]
           for i,s in enumerate(sites)):
        raise ValueError("No source 35 station Pladias trait taxa may be duplicated or dropped")
    return result


def group_set_at(rows:dict,label:str,i:int)->set:
    if label in GUILDS:return rows[label][i]
    if label=="all_544_both_traits":return set.union(*(rows[g][i] for g in GUILDS))
    if label=="wind_strategy_proxy":return rows["th_wind"][i]|rows["nonth_wind"][i]
    if label=="other_strategy_proxy":return rows["th_other"][i]|rows["nonth_other"][i]
    raise ValueError("Undeclared Pladias strategy group")


def jaccard(a:set,b:set)->float|None:
    union=a|b
    return len(a&b)/len(union) if union else None


def matched_source_metrics(rows:dict,sites:list[dict],triads:list[dict])->dict:
    by_label={lab:defaultdict(list) for lab in LABELS}
    source_common=defaultdict(list)
    for p in triads:
        i,j,k=(p[q] for q in ("source_station_i","source_station_same","source_station_other"))
        source_stretch=sites[i]["section"]
        matched={}
        for label in LABELS:
            a,b,c=(group_set_at(rows,label,v) for v in (i,j,k))
            near=jaccard(a,b); far=jaccard(a,c)
            if near is None or far is None:
                continue # source plot has zero taxa for both inputs; NOT 0 similarity
            value=near-far
            by_label[label][source_stretch].append(value)
            matched[label]=value
        # This difference-in-differences uses identical ORIGINAL anchor
        # triads, not the two guilds' potentially distinct study supports.
        if all(label in matched for label in STRATEGIES):
            source_common[source_stretch].append(
                matched["wind_strategy_proxy"]-
                matched["other_strategy_proxy"])
    values={}
    for label in LABELS:
        lines={line:mean(s) for line,s in by_label[label].items()}
        n=sum(len(x) for x in by_label[label].values())
        values[label]={
            "identified":len(lines)>=MIN_STRETCHES and n>=MIN_ANCHORS,
            "source_usable_station_anchors":n,
            "source_distinct_author_stretches":len(lines),
            "source_stretch_equal_same_minus_other_Jaccard":mean(lines.values()) if lines else None,
            "_source_line_means":lines
        }
    lines={s:mean(z) for s,z in source_common.items()}
    n=sum(len(z) for z in source_common.values())
    values[MAIN]={
        "identified":len(lines)>=MIN_STRETCHES and n>=MIN_ANCHORS,
        "source_usable_station_anchors":n,
        "source_distinct_author_stretches":len(lines),
        "source_stretch_equal_wind_minus_other_matched_gap":mean(lines.values()) if lines else None,
        "_source_line_means":lines
    }
    return values


def exact_source_margins(rows:dict,groups:list[list[int]])->dict:
    return {guild:[{
        "original_each_station_guild_richness":[len(rows[guild][i]) for i in indices],
        "original_each_taxon_station_occupancy":source_freq([rows[guild][i] for i in indices])
    } for indices in groups] for guild in GUILDS}


def source_strata(sites:list[dict],kind:str)->list[list[int]]:
    if kind==NULLS[0]:return [list(range(len(sites)))]
    if kind==NULLS[1]:
        grouped=defaultdict(list)
        for i,s in enumerate(sites):grouped[s["section"]].append(i)
        if len(grouped)!=5:raise ValueError("Source original Czech 5 stretches lost")
        return [indices for _,indices in sorted(grouped.items())]
    raise ValueError("No source selected null after looking at observed group effects")


def draw_source_matrix(rows:dict,strata:list[list[int]],rng:random.Random,attempts:int)->int:
    total=0
    for guild in GUILDS:
        for indices in strata:
            if len(indices)<2:continue
            # Pointers refer to the same source sets. No row recreation.
            sets=[rows[guild][i] for i in indices]
            total+=swaps(sets,rng,attempts)
    return total


def technical_null(sites:list[dict],guild:dict,triads:list[dict],
                   kind:str,seed:int,draws:int=DRAWS,
                   burnin:int=BURNIN,thin:int=THIN)->dict:
    if draws<30:raise ValueError("Too few technical proposals to characterize source null")
    initial=guild_source_matrix(sites,guild)
    strata=source_strata(sites,kind)
    margins=exact_source_margins(initial,strata)
    original=matched_source_metrics(initial,sites,triads)
    if not original[MAIN]["identified"]:
        raise ValueError("Original Pladias wind-vs-other source support too limited for main interaction")
    preserved={lab:(original[lab]["identified"],original[lab]["source_usable_station_anchors"],
                    original[lab]["source_distinct_author_stretches"]) for lab in LABELS+(MAIN,)}
    nulls={lab:[] for lab in LABELS+(MAIN,) if original[lab]["identified"]}
    null_stretch=defaultdict(list)
    chains=[]
    for ch in range(NCHAINS):
        matrices={g:[set(t) for t in initial[g]] for g in GUILDS}
        rng=random.Random(seed+ch*7121)
        accept=draw_source_matrix(matrices,strata,rng,burnin)
        run_means=defaultdict(list)
        for _ in range(draws):
            accept+=draw_source_matrix(matrices,strata,rng,thin)
            q=matched_source_metrics(matrices,sites,triads)
            if {lab:(q[lab]["identified"],q[lab]["source_usable_station_anchors"],
                     q[lab]["source_distinct_author_stretches"]) for lab in LABELS+(MAIN,)}!=preserved:
                raise ValueError("Source original matched station×guild support drift under exact richness-preserving checkerboard")
            for label in nulls:
                observed_key=("source_stretch_equal_wind_minus_other_matched_gap" if label==MAIN
                              else "source_stretch_equal_same_minus_other_Jaccard")
                nulls[label].append(q[label][observed_key])
                run_means[label].append(q[label][observed_key])
                for sec,val in q[label]["_source_line_means"].items():
                    null_stretch[label,sec].append(val)
        if exact_source_margins(matrices,strata)!=margins:
            raise ValueError("True source station guild species richness or taxon frequency not conserved")
        chains.append({
            "technical_chain":ch+1,
            "n_accepted_checkerboard_swaps":accept,
            "n_proposals_per_source_guild_per_source_stratum":burnin+draws*thin,
            "source_main_wind_minus_other_mean":
                mean(run_means[MAIN])
        })
    def quantile(vals:list[float],q:float)->float:
        a=sorted(vals);n=len(a);pos=(n-1)*q;low=int(pos)
        return a[low] if low==n-1 else a[low]*(1-(pos-low))+a[low+1]*(pos-low)
    computed={}
    for label in LABELS+(MAIN,):
        q=original[label]
        if not q["identified"]:
            computed[label]={
                "identified":False,
                "source_actual_station_anchors_with_two_defined_species_assemblages":q["source_usable_station_anchors"],
                "source_original_stretches":q["source_distinct_author_stretches"],
                "NOT_zero_biological_effect":True}
            continue
        metric=("source_stretch_equal_wind_minus_other_matched_gap" if label==MAIN
                else "source_stretch_equal_same_minus_other_Jaccard")
        observed=q[metric]
        null=nulls[label]
        line=q["_source_line_means"]
        computed[label]={
            "identified":True,
            "original_station_anchor_sites":q["source_usable_station_anchors"],
            "original_author_stretches":q["source_distinct_author_stretches"],
            "observed_equal_source_stretch_Jaccard_gap":observed,
            "original_guild_frequency_fixed_null_mean":mean(null),
            "observed_minus_frequency_fixed_null":observed-mean(null),
            "technical_correlated_null_q025":quantile(null,.025),
            "technical_correlated_null_q975":quantile(null,.975),
            "source_stretches_observed_exceed_technical_null":sum(
                val>mean(null_stretch[label,sec]) for sec,val in line.items()
            ),
            "technical_upper_rank_NOT_field_experiment_p":
                (1+sum(v>=observed-1e-12 for v in null))/(len(null)+1),
            "source_original_stretch_LOO_range":[
                min(mean([v for i,v in enumerate(line.values()) if i!=k]) for k in range(len(line))),
                max(mean([v for i,v in enumerate(line.values()) if i!=k]) for k in range(len(line)))
            ],
        }
    return {
        "source_original_35_station_5_stretch_fixed_trait_richness_taxon_frequency_null_kind":kind,
        "source_technical_correlated_null_draws":draws*NCHAINS,
        "source_two_independent_seeds_acceptance":chains,
        "predeclared_full_trait_four_guild_and_wind_other_results":computed,
        "frozen_original_true_matched_anchor_site_indices_not_exported":True,
        "site_guild_richness_and_taxon_occurrence_fixed_exact":True,
        "not_a_randomized_train_seed_field_experiment":True
    }


def run(s1:Path,s2:Path,disp:Path,life:Path,sr70:Path,
        draws:int=DRAWS,burnin:int=BURNIN,thin:int=THIN)->dict:
    sites,guild,source,matched,coverage=original_station_traits(s1,s2,disp,life,sr70)
    originals=guild_source_matrix(sites,guild)
    q=matched_source_metrics(originals,sites,matched)
    outputs={kind:technical_null(sites,guild,matched,kind,
        seed=SEED+17017*k,draws=draws,burnin=burnin,thin=thin)
        for k,kind in enumerate(NULLS)}
    original_results={label:{
        k:v for k,v in q[label].items() if not k.startswith("_")
        } for label in LABELS+(MAIN,)}
    return {
        "version":"v28.4",
        "source_original_Czech_Preslia_DOI":"10.23855/preslia.2024.247",
        "genuine_Preslia_original_35_station_5_stretch_plus_2026_SR70_and_Pladias_source_QA":source,
        "genuine_matched_Pladias_source_taxonomic_coverage_QA":coverage,
        "original_33_geographical_source_matched_anchor_plants_by_trait_group":original_results,
        "all_guild_conditioned_genuine_original_source_nulls":outputs,
        "Pladias_Sadlo_Epilobium_Phragmites_archetypes_NOT_measured_wind_transport_rate":True,
        "same_source_station_anchor_triples_from_original_761_species_independent_geo_matching":True,
        "source_matched_stretches_original_real_countries_not_independent_replicates":True,
        "Czech_2020_22_station_entire_flora_not_source_exact_BS_BT_100m_Slovenia":True,
        "no_rail_freight_OD_or_mowing_wind_adaptation_or_seed_arrival_measured":True,
        "raw_original_plant_taxon_station_geometry_and_source_trait_records_NOT_republished":True
    }


def main():
    p=argparse.ArgumentParser()
    for name in ("s1","s2","dispersal","lifeform","sr70","output"):
        p.add_argument("--"+name,required=True,type=Path)
    p.add_argument("--draws",type=int,default=DRAWS)
    p.add_argument("--burnin",type=int,default=BURNIN)
    p.add_argument("--thin",type=int,default=THIN)
    args=p.parse_args()
    report=run(args.s1,args.s2,args.dispersal,args.lifeform,args.sr70,
               args.draws,args.burnin,args.thin)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    print(json.dumps({"version":"v28.4",
        "original_source_35_station_coverage":report["genuine_matched_Pladias_source_taxonomic_coverage_QA"],
        "source_actual_guild_group_gap":report["original_33_geographical_source_matched_anchor_plants_by_trait_group"]
    },ensure_ascii=False))
    for kind,panel in report["all_guild_conditioned_genuine_original_source_nulls"].items():
        print(json.dumps({
            "version":"v28.4","source_null":kind,
            "real_original_guild_results":{
                lab:{
                    k:v for k,v in val.items()
                    if k in ("identified","original_station_anchor_sites","original_author_stretches",
                    "observed_equal_source_stretch_Jaccard_gap","original_guild_frequency_fixed_null_mean",
                    "observed_minus_frequency_fixed_null","source_stretches_observed_exceed_technical_null",
                    "technical_upper_rank_NOT_field_experiment_p")
                } for lab,val in panel["predeclared_full_trait_four_guild_and_wind_other_results"].items()
            }},ensure_ascii=False))


if __name__=="__main__":
    main()
