#!/usr/bin/env python3
"""v14.1: genuine 89-location Slovenian RAIL vegetation field study
functional-profile-constrained taxonomic identity sharing.

Author-original Šabić, Grašič & Rovanšek (2026), DOI
10.1007/s44473-026-00155-z; published 401279-byte 2026 CC BY
supplement SHA256 26be7ca9513dd54136d4ec499d8413c1cb288045f87fe3cf8814731f2edac99f.

Uses same EXACT original BS/BT source zones on matched track sides.
Side pairs are unioned WITHIN 89 independent railway field locations.
Locs nested within 15 ORIGINAL railway route SECTIONS.
Only taxa with actual BOTH 'Diaspore.strategy' and 'Lifespan.categories'
known enter primary functional+taxonomic shared-species analysis
(261 taxa in full original 89-location source).
Missing trait codes never set to 'other'; 312 original full taxa
site-incidence taxonomic similarity is also reported as sensitivity.

This FIRST real-plot same study pilot checks whether two plots
belonging to same source section share more taxonomic species
identities than a fixed site × functional-guild and
species × regional-frequency null would predict.

4 source-defined guilds: strict-annual vs other known lifespan,
and published Epilobium+Phragmites wind strategy archetypes
vs other original known dispersal strategies. Guild names are
STRATEGY ARCHETYPES, not species genera or actual seed transport.

Null: 2×2 checkerboard swaps within each guild, conserving EACH
original field site's exact guild taxon richness and EACH source
taxon's study-wide observed original site occupancy.
Sensitivity restricts swaps to same original YEAR × TRACK-TYPE
stratum, additionally conserving source species prevalence within
those strata. Markov chains, their swap counts and mixing noted.
No false independent 3916-pair sample sizes: this is descriptive
within-one-country/15-section permutation and leave-section-out.
No false physical rail connectivity predictor or rail freight
direction; SAME SOURCE SECTION is a section membership indicator
and both nearby climate and local regional pools are unadjusted.
No non-rail habitat controls, so even perfect same-section
similarity could NOT establish 'rail-specific functional convergence'.
No train-carriage seed flux, mowing dates, gusts or adaptation observed.

No original XLSX redistributed; source CI re-downloads previously
hash-audited original publisher XLSX.
"""
from __future__ import annotations
import argparse, copy, json, math, random
from collections import Counter, defaultdict
from datetime import datetime, timedelta
from itertools import combinations
from pathlib import Path
from statistics import mean, median

from analyze_Slovenia_rail_microhabitats_v103 import (
    SOURCE_DOI, SOURCE_SHA256, extract_original_xlsx, validate_original_keys
)
from analyze_Slovenia_rail_species_occupancy_permutation_v109 import (
    original_location_matched_sets, source_guild
)

SEED=20261115
DRAWS_PER_CHAIN=180
N_CHAINS=2
BURNIN_PROPOSALS=32000
PROPOSALS_BETWEEN_DRAWS=6000
FUNCTIONAL_CATEGORIES=(
    "strict_annual_wind",
    "strict_annual_other",
    "other_known_wind",
    "other_known_other",
)
SWAP_MODES=("all_years_both_track_types","same_survey_year_and_track_type")
CATEGORICAL_SOURCE_NOT_SPATIALLY_CALIBRATED=True


def source_functional_label(source_trait):
    d=source_guild(source_trait)
    if d is None:return None
    if d["stratum"]=="strict_annual":
        return "strict_annual_wind" if d["wind"] else "strict_annual_other"
    return "other_known_wind" if d["wind"] else "other_known_other"


def source_sites_and_guilds(original):
    matched, traits, join_audit=original_location_matched_sets(original)
    original_loc, z, sp, tr, cb, other_audit=validate_original_keys(original)
    labels={code:source_functional_label(t) for code,t in traits.items()}
    known={k:v for k,v in labels.items() if v is not None}
    sites=[]
    for key,x in sorted(matched.items(),key=lambda it:str(it[0])):
        all_taxa=set(x["BS"])|set(x["BT"])
        keep=all_taxa&set(known)
        field=original_loc[key]
        year=(datetime(1899,12,30)+timedelta(days=float(field["date"]))).year
        if year not in (2022,2023):
            raise ValueError("Original source field survey year must be 2022 or 2023")
        section=str(x["section_number"])
        sites.append({
            "site":str(key),
            "section":section,
            "year":year,
            "track_type":str(field["track_type"]),
            "known":set(keep),
            "all_original_taxa":set(all_taxa),
            "n_original_matched_track_sides":x["n_matched_sides"]
        })
    if len(sites)!=89 or len({x["section"] for x in sites})!=15:
        raise ValueError("Original Slovenian real 89-source-site, 15 rail-section cohort changed")
    if len(known)!=261:
        # Exact source-known population must be logged, not silently
        # weakened to a more favorable posthoc functional sample.
        raise ValueError("Original Slovenian dual-known source functional taxa !=261")
    if not all(len(x["known"])>=1 for x in sites):
        raise ValueError("An original site has no jointly strategy/lifespan-known observed taxa")
    return sites, known, {
        "original_authors_field_study_DOI":SOURCE_DOI,
        "original_source_supplement_SHA256":SOURCE_SHA256,
        "original_sampled_locations":len(sites),
        "original_route_sections":15,
        "original_matched_side_BS_BT_pairs":sum(s["n_original_matched_track_sides"] for s in sites),
        "original_strategy_lifespan_known_taxa":len(known),
        "all_published_source_trait_taxa":len(traits),
        "original_qualifying_site_taxa_observations":sum(len(x["known"]) for x in sites),
        "original_all_site_taxa_observations_deduplicated_per_location":sum(
            len(x["all_original_taxa"]) for x in sites),
        "known_fraction_site_median":median(len(x["known"])/len(x["all_original_taxa"])
                                            for x in sites),
        "known_fraction_site_min":min(len(x["known"])/len(x["all_original_taxa"])
                                         for x in sites),
        "source_site_count_by_year":dict(Counter(x["year"] for x in sites)),
        "source_site_count_by_track_type":dict(Counter(x["track_type"] for x in sites)),
        "source_sampled_rail_sections_have_no_verified_actual_train_OD_or_freight_flow":True,
        "original_geographical_coordinates_source_columns":sorted(original_loc[next(iter(original_loc))]),
    }


def jaccard(x,y):
    union=len(x|y)
    return len(x&y)/union if union else None


def functional_overlap(x,y,traits):
    a=Counter(traits[k] for k in x)
    b=Counter(traits[k] for k in y)
    if not x or not y:return None
    return 1-0.5*sum(abs(a[k]/len(x)-b[k]/len(y))
                     for k in FUNCTIONAL_CATEGORIES)


def pair_summary(sites, groups):
    n=len(sites)
    if n<4:raise ValueError("Require real independent source sites, not single pair")
    within,between=[],[]
    fullwithin,fullbetween=[],[]
    fwithin,fbetween=[],[]
    for i,j in combinations(range(n),2):
        a=sites[i];b=sites[j]
        same=(a["section"]==b["section"])
        y=jaccard(a["known"],b["known"])
        full=jaccard(a["all_original_taxa"],b["all_original_taxa"])
        fs=functional_overlap(a["known"],b["known"],groups)
        (within if same else between).append(y)
        (fullwithin if same else fullbetween).append(full)
        (fwithin if same else fbetween).append(fs)
    if not within or not between:
        raise ValueError("Original 15 sections must include both within- and between-section site pairs")
    def means(a,b):
        return {"mean_same_rail_section":mean(a),
                "mean_different_rail_section":mean(b),
                "same_minus_different":mean(a)-mean(b)}
    return {
        "n_source_sites":n,
        "n_dependent_within_section_site_pairs_NOT_independent_studies":len(within),
        "n_dependent_between_section_site_pairs_NOT_independent_studies":len(between),
        "taxonomic_Jaccard_similarity_dual_trait_known_same_taxa":means(within,between),
        "taxonomic_Jaccard_similarity_ALL_312_original_source_taxa_sensitivity":
            means(fullwithin,fullbetween),
        "functional_4guild_proportional_overlap_on_dual_known_taxa":
            means(fwithin,fbetween),
        "original_site_pair_dependency_not_pseudoreplicated":True,
    }


def make_matrices(sites,guild):
    allsites=[x["known"] for x in sites]
    return {group:[x & {code for code,g in guild.items() if g==group}
                   for x in allsites] for group in FUNCTIONAL_CATEGORIES}


def row_col_margins(matrix):
    bygroup={}
    for g,rows in matrix.items():
        cols=Counter(code for row in rows for code in row)
        bygroup[g]={
            "site_guild_counts":tuple(len(r) for r in rows),
            "taxon_source_occupancies":dict(cols),
        }
    return bygroup


def swap_proposals(matrices,site_keys,rng,nproposals,stratify):
    if stratify:
        strata=defaultdict(list)
        for i,site in enumerate(site_keys):
            strata[(site["year"],site["track_type"])].append(i)
        allowable=[x for x in strata.values() if len(x)>=2]
    else:allowable=[list(range(len(site_keys)))]
    if not allowable:raise ValueError("Original source strata cannot supply two sites to any null swaps")
    accepted=0
    groups=tuple(FUNCTIONAL_CATEGORIES)
    for _ in range(nproposals):
        indices=rng.choice(allowable)
        i,j=rng.sample(indices,2)
        rows=matrices[rng.choice(groups)]
        left=rows[i]-rows[j]
        right=rows[j]-rows[i]
        if not left or not right:continue
        a=rng.choice(tuple(left));b=rng.choice(tuple(right))
        rows[i].remove(a);rows[i].add(b)
        rows[j].remove(b);rows[j].add(a)
        accepted+=1
    return accepted


def matrix_taxon_lists(matrices,n):
    result=[set() for _ in range(n)]
    for _,rows in matrices.items():
        for i,group_set in enumerate(rows):
            result[i].update(group_set)
    return result


def taxonomic_gap(sets,sections):
    within_total=between_total=0.
    within_n=between_n=0
    for i,j in combinations(range(len(sets)),2):
        val=jaccard(sets[i],sets[j])
        if val is None:continue
        if sections[i]==sections[j]:
            within_total+=val;within_n+=1
        else:
            between_total+=val;between_n+=1
    if not within_n or not between_n:
        raise ValueError("Cannot estimate same-vs-other source-section taxonomic gap")
    return within_total/within_n-between_total/between_n


def null_diagnostics(sites,guild,mode,draws_per_chain=DRAWS_PER_CHAIN,
                     nchains=N_CHAINS,burnin=BURNIN_PROPOSALS,
                     spacing=PROPOSALS_BETWEEN_DRAWS,seed=SEED):
    if mode not in SWAP_MODES:raise ValueError("Undeclared source null family")
    stratified=(mode=="same_survey_year_and_track_type")
    originals=make_matrices(sites,guild)
    margins=row_col_margins(originals)
    observed=[s["known"] for s in sites]
    groups=[s["section"] for s in sites]
    obs=taxonomic_gap(observed,groups)
    null=[]
    perchain=[]
    for chain in range(nchains):
        m=copy.deepcopy(originals)
        r=random.Random(seed+1009*chain+(1 if stratified else 0))
        accepted=swap_proposals(m,sites,r,burnin,stratified)
        for _ in range(draws_per_chain):
            accepted+=swap_proposals(m,sites,r,spacing,stratified)
            if row_col_margins(m)!=margins:
                raise ValueError("Source guild-margin-preserving swap invalid")
            draw=matrix_taxon_lists(m,len(sites))
            if any(len(x)!=len(y) for x,y in zip(draw,observed)):
                raise ValueError("Original site richness per source functional stratum not retained")
            null.append(taxonomic_gap(draw,groups))
        nprop=burnin+draws_per_chain*spacing
        perchain.append({
           "chain":chain+1,"proposed_source_checkerboard_swaps":nprop,
           "accepted_source_checkerboard_swaps":accepted,
           "accepted_fraction":accepted/nprop})
    if len(null)<20:raise ValueError("Source null requires enough stochastic observations")
    sorted_null=sorted(null)
    def q(p):
        j=p*(len(null)-1);l=math.floor(j);h=math.ceil(j)
        return sorted_null[l]+(sorted_null[h]-sorted_null[l])*(j-l)
    p_plus=(1+sum(z>=obs for z in null))/(len(null)+1)
    return {
        "randomization_family":mode,
        "observed_source_section_taxonomic_similarity_gap":
            obs,
        "expected_source_gap_under_strict_guild_and_species_occupancy_margin_null":mean(null),
        "source_excess_similarity_gap_above_guild_frequency_null":obs-mean(null),
        "source_chain_null_gap_percentile_2_5":q(0.025),
        "source_chain_null_gap_percentile_97_5":q(0.975),
        "source_guild_and_species_occupancy_fixed_margin_upper_tail_fraction_plus_one":
            p_plus,
        "MonteCarlo_draws_NOT_independent_field_studies":len(null),
        "separate_chains_and_within_chain_autocorrelation_not_exact_permutation":perchain,
        "note_about_P_value":"Exploratory source fixed-margin swap reference diagnostic, NOT valid train/dispersal causal p-value; MCMC samples correlated",
        "original_year_track_type_stratum_species_counts_fixed":stratified,
        "year_stratified_swap_still_does_not_control_geographic_distance_or_local_species_pools":True,
    }


def section_leaveout(sites,guild):
    original_sections=sorted({s["section"] for s in sites})
    rows=[]
    for section in original_sections:
        sel=[s for s in sites if s["section"]!=section]
        if len({x["section"] for x in sel})<2:continue
        stat=taxonomic_gap([s["known"] for s in sel],[s["section"] for s in sel])
        rows.append({"excluded_source_section":section,
                     "n_remaining_actual_original_sites":len(sel),
                     "same_minus_different_taxonomic_similarity":stat})
    return {
        "original_source_rail_sections_left_out":len(rows),
        "n_leave_one_source_section_out_positive_section_clustering":
            sum(x["same_minus_different_taxonomic_similarity"]>0 for x in rows),
        "min_same_minus_different":min(x["same_minus_different_taxonomic_similarity"] for x in rows),
        "max_same_minus_different":max(x["same_minus_different_taxonomic_similarity"] for x in rows),
        "source_section_leaveout_list":rows,
    }


def run(original):
    sites,guild,source=source_sites_and_guilds(original)
    observed=pair_summary(sites,guild)
    test={
       mode:null_diagnostics(sites,guild,mode)
       for mode in SWAP_MODES
    }
    return {
        "complete":True,"version":"v14.1",
        "original_source_provenance":source,
        "original_primary_89_site_CW_BS_BT_source_beta_summary":observed,
        "original_guild_frequency_fixed_margins_stochastic_reference":test,
        "section_cluster_robustness_NOT_independent_country_replication":section_leaveout(sites,guild),
        "biological_limits":[
          "Original author selected plant-visible Slovenian locations; 15 source railway segments, not randomly representative of railway freight services.",
          "Section membership is NOT a graph shortest path, railway service schedule, train-specific seed carriage, section freight volume or a track-seed dispersal event.",
          "Four broad source-lifespan and source-diaspore guilds are coarse traits; functional overlap is a compositional proxy NOT Gower multivariate functional beta hypervolume.",
          "Original shared species patterns may arise from geographic proximity, source geology/climate, habitat, regional source pools, and nonrandom line selection.",
          "Degree-preserving null fully fixes source site guild richness and original taxon frequency; it cannot rule out unmeasured environment or local common supply hubs.",
          "MCMC checkerboard swaps are autocorrelated, not independent null experimental field sites or exact p-value randomization.",
          "Site pairs share field locations, so 89 choose 2 pairs are NOT independent samples, and the 15 original line sections are nested in only ONE country.",
          "Only actual dual strategy+lifespan known species are in the primary taxonomic/functional paired comparisons; full taxon source set is sensitivity, not imputed guilds.",
          "Year and track type stratification controls only those two original metadata categories; it does NOT adjust physical rail/road geographic paths or temporal operations.",
          "No local nonrail comparator in source; same-section functional proximity cannot establish rail-specific environmental convergence, and taxonomic surplus is not direct train-seed dispersal evidence.",
          "No same-species within-rail height, common garden, trackside cutting actual dates, phenology response, germination, fitness, plasticity or genetically evolved adaptation."
        ],
        "status":"REAL author-original plant community composition pilot, NOT a transport-network causal test or railway meta-analysis"
    }


def main():
    p=argparse.ArgumentParser()
    for name in ("source-xlsx","source-audit102","output"):
        p.add_argument("--"+name,required=True)
    a=p.parse_args()
    original=extract_original_xlsx(a.source_xlsx,a.source_audit102)
    result=run(original)
    dest=Path(a.output);dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding="utf-8")
    summary={
        "complete":True,"version":"v14.1",
        "sample_source":result["original_source_provenance"],
        "observed":result["original_primary_89_site_CW_BS_BT_source_beta_summary"],
        "swap_null":result["original_guild_frequency_fixed_margins_stochastic_reference"],
        "section_leaveout":{k:v for k,v in result[
            "section_cluster_robustness_NOT_independent_country_replication"].items()
            if k!="source_section_leaveout_list"},
        "no_train_seed_causality_or_cross_country_replication":True
    }
    print(json.dumps(summary,ensure_ascii=False),flush=True)


if __name__=="__main__":main()
