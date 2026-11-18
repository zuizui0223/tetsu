#!/usr/bin/env python3
"""v28.5: four ORIGINAL railway flora studies, exact taxonomic plant
assemblages FROZEN, permute paired continuous functional traits across
their observed species. A serious independent null distinct from v24.5/24.6
source taxon-incidence margin checkerboards.

Four independently published SOURCE studies, five European countries:
  CZ Kutlvašr 2024, 39 actual entire-station floras (5 stretches)
  IT Adorni/Ghillani 2020, 28 actual entire-station floras (plain/other)
  HR+BA Jasprica 2017, 15 actual entire-station floras (three climates),
        HR and BA are ONE paper not two independently sampled studies
  SI Šabić et al 2026, 89 true 100-m ballast microhabitat physical
        site unions BS+BT, 15 railway national-line sections. SI is
        explicitly not an 89-STATION flora.

Taxon identities / their exact site distributions / every site species
richness / every study original geographical block / original observed
taxonomic turnover are all fully FROZEN. Only swap 2D paired LEDA
z[log canopy height(m), log seed mass(mg)] species values between
source species occurring in that study, preserving correlation between
the two traits and each study's actual observed trait distribution.
The Czech/Italy/HRBA common z score axes were frozen in v24.5 before
Slovenia and preserved for Slovenia as reference exactly as v24.6.
Never re-normalize each country's original 2D axes.

Compare:
 A. each original study's paired 2D observed trait vector shuffled
    uniformly among its ACTUAL trait-known source flora;
 B. original observed taxa stratified into FIVE mutually exclusive
    species source-site frequency quantiles; shuffle within quantiles
    (approximate prevalence–trait relationship retained, but not exact).
Quantile bins are defined by original species occupancy frequency,
NOT ecological outcomes. Within species pair vectors moved together.
Original species matrix itself is NEVER changed.

Metrics evaluated over same original field sites:
- mean pairwise site CWM2D Euclidean distance, all site-pairs and
  within original field-study geography group;
- 1D height and seed-mass CWM differences independently;
- within original site trait spread trace, a distinct alpha statistic.

A positive observed-minus-null for between-SITE CWM distance is
consistent with original species identities being located nonrandomly
relative to their source trait values. It does NOT alone prove
rail-specific environmental selection or adaptation: source local flora,
spatial environment, phylogeny, historical occupancy, trait sampling
missingness and sample area can all generate it.
A negative contrast implies site-level trait convergence conditional on
observed taxonomic source assemblages, NOT railway-vs-nonrail.
These four field source studies are not randomly selected globally and
never create 171 equivalent railway stations or independent 4-country
train experiments.

Every source downloaded and SHA256 pinned in GitHub Actions: Czech
publisher S1 and corrected S2, Italian original 54p PDF, HR+BA original
33p PDF, Slovenia author original XLSX, original University Oldenburg
LEDA canopy/seed files, and v24.5/v24.6 validated aggregate
reference metric artifacts. No original plant records, exact locality
or source trait tables exported.
"""
from __future__ import annotations
import argparse,hashlib,json,math
from pathlib import Path
from statistics import mean
import numpy as np

from validate_Europe_3studies_82stations_original_LEDA_2D_trait_beta_v245 import (
    publisher_czech,publisher_italy,publisher_hr_ba,
    traits_join,source_cohort,get_metrics,LEDASHA
)
from validate_SI_89_original_track_sites_Europe82_trait_axis_v246 import (
    genuine_slovenian_89_microhabitat_unions
)
from audit_SI_CZ_original_continuous_trait_coverage_v220 import read_leda
from analyze_Slovenia_rail_microhabitats_v103 import extract_original_xlsx

STUDIES=("CZ_Czechia","IT_Italy_Parma","HR_BA_2017_one_study","SI_Slovenia_100m_track")
N={"CZ_Czechia":39,"IT_Italy_Parma":28,
   "HR_BA_2017_one_study":15,"SI_Slovenia_100m_track":89}
MODES=("permute_whole_observed_study_trait_pool",
       "permute_within_original_species_occupancy_quintiles")
N_CHAINS=2
DRAWS=299
SEED=285118
SOURCE_REQUIRED_STUDIES=4
# Exact successful v24.5/v24.6 first-party source-audited checks. Earlier
# aggregate artifact retrieval can fail with GitHub cross-run credentials;
# instead regenerate the actual 3-study axis from ALL original SHA-gated
# publisher materials and verify every scalar against these frozen values.
FROZEN_3_STUDY_LOG_MEANS=(-0.5802365754696401,0.2936113186962433)
FROZEN_3_STUDY_LOG_SDS=(1.3801000545693773,2.1657660559995824)
FROZEN_SOURCE_CWM_MEANS={
    "CZ_Czechia":0.16501788905535814,
    "IT_Italy_Parma":0.2014200675821101,
    "HR_BA_2017_one_study":0.2232461322998969,
    "SI_Slovenia_100m_track":0.35003008736045554,
}
METRICS=("cwm2D_mean_all_source_pairs",
         "cwm2D_mean_within_original_geography",
         "cwm_height_absolute_mean_within_original_geography",
         "cwm_seed_mass_absolute_mean_within_original_geography",
         "source_site_within_2D_trait_variance_trace")


def frozen_reference(originals:dict,raw_z:dict)->tuple[dict,dict,dict]:
    studies,original_qa=originals
    traits,traitqa=traits_join(studies,raw_z["height"],raw_z["seed"])
    sourcecohort,cov=source_cohort(studies,traits)
    if len(sourcecohort)!=3 or {k:len(x) for k,x in sourcecohort.items()}!={k:N[k] for k in N if k!="SI_Slovenia_100m_track"}:
        raise ValueError("Original 82 complete Preslia, Parma, Croatia/Bosnia railway stations disappeared")
    if traitqa["source_original_LEDA_both_canopy_height_seed_mass_species_n"]!=622:
        raise ValueError("Original v24.5 European fixed 622 species trait axis changed")
    mu=traitqa["species_equal_global_log_mean"]
    sd=traitqa["species_equal_global_log_SD"]
    if (any(abs(x-y)>1e-12 for x,y in zip(mu,FROZEN_3_STUDY_LOG_MEANS))
        or any(abs(x-y)>1e-12 for x,y in zip(sd,FROZEN_3_STUDY_LOG_SDS))):
        raise ValueError("Original 2026 Nov 17 source 3-publisher European 2D reference z axis changed")
    # Exactly the held-out Slovenian source transformation frozen in v24.6,
    # calculated solely from the 3 other independently published source
    # studies, NOT ever scaled anew using Slovenian plant observations.
    si_all={taxon:((math.log(raw_z["height"][taxon])-mu[0])/sd[0],
                   (math.log(raw_z["seed"][taxon])-mu[1])/sd[1])
            for taxon in set(raw_z["height"])&set(raw_z["seed"])}
    si_sites,sl_qa=genuine_slovenian_89_microhabitat_unions(
        original_qa["slovenia"],si_all)
    if len(si_sites)!=89:
        raise ValueError("Slovenian original source 89 100m railway plots disappeared")
    sourcecohort["SI_Slovenia_100m_track"]=si_sites
    source_trait={}
    observed={}
    for k in STUDIES:
        source_trait[k]=si_all if k=="SI_Slovenia_100m_track" else traits
        if k not in sourcecohort or len(sourcecohort[k])!=N[k]:
            raise ValueError(f"One original study source site population changed: {k}")
        baseline=get_metrics(sourcecohort[k],source_trait[k])[
            "functional_CWM_centroid_distance_all_stationpairs"]
        if abs(baseline-FROZEN_SOURCE_CWM_MEANS[k])>1e-10:
            raise ValueError(f"Reanalyzed genuine original {k} source CWM does not match frozen v24.5/v24.6 author-source audit")
        observed[k]=baseline
    return sourcecohort,source_trait,{
        "original_3_study_v245_author_source":{k:v for k,v in original_qa.items() if k!="slovenia"},
        "original_3_study_trait_axis":traitqa,
        "original_3_study_trait_site_coverage":cov,
        "original_Slovenian_89_100m_source":sl_qa,
        "v245_AND_v246_prior_true_original_observed_CWM_baselines_exactly_reproduced":True,
        "v245_v246_prior_genuine_source_frozen_observed_means":observed,
        "v245_refit_from_original_3_study_source_data_not_private_action_artifact":True,
        "source_four_independently_published_railway_field_studies":True
    }


def build_matrix(sites:list[dict],traits:dict)->tuple[np.ndarray,np.ndarray,list[str],dict]:
    species=sorted(set().union(*(site["taxa"] for site in sites)))
    if any(t not in traits for t in species):
        raise ValueError("Original site taxon with no verified University Oldenburg traits found")
    matrix=np.zeros((len(sites),len(species)),dtype=np.float64)
    idx={s:i for i,s in enumerate(species)}
    for i,site in enumerate(sites):
        for sp in site["taxa"]:matrix[i,idx[sp]]=1.
    counts=matrix.sum(axis=1)
    if np.any(counts<5) or len(species)<25:
        raise ValueError("Source authentic original taxon coverage below survey unit threshold")
    trait=np.array([traits[t] for t in species],dtype=np.float64)
    if trait.shape!=(len(species),2) or not np.isfinite(trait).all():
        raise ValueError("Original matched source height/mass values nonfinite")
    W=matrix/counts[:,None]
    geo=[x["group"] for x in sites]
    ii,jj=np.triu_indices(len(sites),1)
    within=np.array([geo[i]==geo[j] for i,j in zip(ii,jj)],dtype=bool)
    if not within.any() or within.all():
        raise ValueError("Independent source study original geography strata do not have within and between sites")
    qa={
        "source_original_sites":len(sites),
        "source_actual_plant_species_with_both_traits":len(species),
        "original_source_site_known_trait_species_min_max":[int(counts.min()),int(counts.max())],
        "genuine_original_geography_group_count":len(set(geo)),
        "source_independent_sitepairs_not_inference_replicates":int(len(ii)),
        "source_within_original_geography_pair_count_not_independent":int(within.sum()),
        "original_each_site_richness_and_species_occupancy_NEVER_changed":True,
        "original_2D_trait_vector_pairs_permuted_together":True
    }
    return W,trait,geo,{"i":ii,"j":jj,"within":within,"species_site_occurrence":matrix.sum(axis=0),**qa}


def metric_values(W:np.ndarray,T:np.ndarray,ii:np.ndarray,jj:np.ndarray,
                  within:np.ndarray)->dict:
    if len(W.shape)!=2 or T.shape!=(W.shape[1],2):
        raise ValueError("Source incidence and trait shapes not conformable")
    C=W@T
    d=C[ii]-C[jj]
    norm=np.sqrt(np.sum(d*d,axis=1))
    trace=np.sum(W@(T*T)-C*C,axis=1)
    if np.any(trace<-1e-10):
        raise ValueError("Original site trait variance trace negative after bivariate paired shuffling")
    result={
        "cwm2D_mean_all_source_pairs":float(norm.mean()),
        "cwm2D_mean_within_original_geography":float(norm[within].mean()),
        "cwm_height_absolute_mean_within_original_geography":float(np.abs(d[within,0]).mean()),
        "cwm_seed_mass_absolute_mean_within_original_geography":float(np.abs(d[within,1]).mean()),
        "source_site_within_2D_trait_variance_trace":float(trace.mean()),
    }
    if not all(math.isfinite(x) for x in result.values()):
        raise ValueError("Source genuine two-trait functional outcome nonfinite")
    return result


def occupancy_quintile_strata(occupancy:np.ndarray)->list[np.ndarray]:
    # Global deterministic species order by source occurrence frequency,
    # stable index for taxa tied in original source occupancy.
    order=np.argsort(occupancy,kind="mergesort")
    return [q for q in np.array_split(order,5) if len(q)>=2]


def permutation(W:np.ndarray,T:np.ndarray,ii:np.ndarray,jj:np.ndarray,
                within:np.ndarray,strata:list[np.ndarray],seed:int,draws:int=DRAWS)->dict:
    if draws<50:raise ValueError("Insufficient technical independently permuted source trait assignments")
    initial=metric_values(W,T,ii,jj,within)
    source_sorted=np.sort(T,axis=0)
    rng=np.random.default_rng(seed)
    vec={k:[] for k in METRICS}
    shuffled_counts=[]
    for d in range(draws):
        labels=np.arange(len(T))
        for group in strata:
            labels[group]=rng.permutation(group)
        candidate=T[labels] # 2D height/seed vectors stay paired
        if d==0:
            if not np.array_equal(np.sort(candidate,axis=0),source_sorted):
                raise AssertionError("Study source 2D functional coordinate pool changed")
            if not np.allclose(np.corrcoef(candidate.T),np.corrcoef(T.T),rtol=1e-12,atol=1e-12):
                raise AssertionError("Source height/seed pair covariance altered by trait shuffle")
        metrics=metric_values(W,candidate,ii,jj,within)
        for k in METRICS:vec[k].append(metrics[k])
    return {
        "observed":initial,"draws":vec,"source_permutation_draws":draws,
        "no_site_species_matrix_altered":True,
        "bivariate_height_mass_vectors_moved_together":True,
        "source_trait_values_multiset_and_covariance_exactly_preserved":True
    }


def quantile(xs:list[float],p:float)->float:
    a=sorted(xs);q=(len(a)-1)*p;lo=int(q)
    return a[lo] if lo==len(a)-1 else a[lo]*(1-(q-lo))+a[lo+1]*(q-lo)


def source_study_null(sites:list[dict],traits:dict,study:str,draws:int=DRAWS)->dict:
    W,T,geo,qa=build_matrix(sites,traits)
    ii,jj,within=qa.pop("i"),qa.pop("j"),qa.pop("within")
    n=int(W.shape[1])
    occupancy=qa.pop("species_site_occurrence")
    groups={
        MODES[0]:[np.arange(n)],
        MODES[1]:occupancy_quintile_strata(occupancy)
    }
    output={}
    for k,mode in enumerate(MODES):
        pool=groups[mode]
        if sum(map(len,pool))!=n:
            raise ValueError("Original species trait stratification dropped taxa")
        allvals={name:[] for name in METRICS}
        chainmeans=[]
        obs=None
        for ch in range(N_CHAINS):
            draw=permutation(W,T,ii,jj,within,pool,
                seed=SEED+17911*list(STUDIES).index(study)+401*ch+19031*k,
                draws=draws)
            if obs is None:obs=draw["observed"]
            elif any(abs(obs[m]-draw["observed"][m])>1e-12 for m in METRICS):
                raise ValueError("Identical authentic source study plant data changed across trait permutations")
            for metric in METRICS:allvals[metric].extend(draw["draws"][metric])
            chainmeans.append({
                "chain":ch+1,
                "n_technical_independent_label_permutations":draws,
                "source_CWM_within_geography_null_mean":
                    mean(draw["draws"]["cwm2D_mean_within_original_geography"])
            })
        values={}
        for metric in METRICS:
            v=allvals[metric];actual=obs[metric]
            values[metric]={
                "observed":actual,
                "bivariate_trait_label_permutation_source_null_mean":mean(v),
                "original_observed_minus_trait_randomized_null":actual-mean(v),
                "technical_q025_NOT_field_CI":quantile(v,.025),
                "technical_q975_NOT_field_CI":quantile(v,.975),
                "technical_source_upper_tail_rank_NOT_field_p":
                    (1+sum(x>=actual for x in v))/(1+len(v)),
                "technical_source_lower_tail_rank_NOT_field_p":
                    (1+sum(x<=actual for x in v))/(1+len(v))
            }
        output[mode]={
            "source_actual_species":n,
            "n_source_species_prevalence_bands":len(pool),
            "source_site_species_community_matrix_is_exact_original":True,
            "source_taxon_observed_prevalence_not_permutes":True,
            "source_2D_trait_vector_correlation_preserved":True,
            "n_technical_non_field_trait_shuffles":len(allvals[METRICS[0]]),
            "chains":chainmeans,
            "metrics":values
        }
    return {"source_study":study,"source_QA":qa,"source_null_results":output}


def run(cz_s1:Path,cz_s2:Path,it_pdf:Path,hr_pdf:Path,
        si_xlsx:Path,si_audit:Path,
        height:Path,seed:Path,draws:int=DRAWS)->dict:
    for file,dig in ((height,LEDASHA["height"]),(seed,LEDASHA["seed"])):
        if hashlib.sha256(file.read_bytes()).hexdigest()!=dig:
            raise ValueError("Author-original Oldenburg LEDA reference traits source SHA changed")
    h,_=read_leda(height,"height_m")
    s,_=read_leda(seed,"seed_mass_mg")
    published={
        "CZ_Czechia":publisher_czech(cz_s1,cz_s2),
        "IT_Italy_Parma":publisher_italy(it_pdf),
        "HR_BA_2017_one_study":publisher_hr_ba(hr_pdf)
    }
    studies={k:v[0] for k,v in published.items()}
    qc={k:v[1] for k,v in published.items()}
    original=extract_original_xlsx(si_xlsx,si_audit)
    cohorts,traits,gate=frozen_reference((studies,{**qc,"slovenia":original}),
        {"height":h,"seed":s})
    results={k:source_study_null(cohorts[k],traits[k],k,draws)
             for k in STUDIES}
    positives={mode:{
        "independent_original_studies_with_observed_functional_CWM_differentiation_above_traitshuffle_null":
            sum(r["source_null_results"][mode]["metrics"]["cwm2D_mean_within_original_geography"][
                "original_observed_minus_trait_randomized_null"]>0 for r in results.values()),
        "independent_original_studies_with_observed_alpha_trait_trace_below_traitshuffle_null":
            sum(r["source_null_results"][mode]["metrics"]["source_site_within_2D_trait_variance_trace"][
                "original_observed_minus_trait_randomized_null"]<0 for r in results.values()),
        "original_independent_study_count_not_field_sites":len(results)
    } for mode in MODES}
    if len(results)!=SOURCE_REQUIRED_STUDIES or sum(r["source_QA"]["source_original_sites"] for r in results.values())!=171:
        raise ValueError("Four-source independent survey sample units drift")
    return {
        "version":"v28.5",
        "genuine_original_4_publisher_studies_5_countries_and_different_plot_grains_QA":gate,
        "four_independent_field_publication_group_trait_label_null":results,
        "two_distinct_nulls_descriptive_study_coherence":positives,
        "three_station_studies_plus_slovenian_89_true_100m_rail_line_study_not_171_rail_stations":True,
        "does_not_estimate_external_nonrail_regional_pool_railway_specific_convergence":True,
        "source_2D_trait_pairs_shuffled_together_taxa_and_original_geography_sites_FROZEN":True,
        "source_original_LEDA_missingness_nonrandom_and_study_grain_not_identical":True,
        "n_technical_trait_label_shuffles_not_randomized_field_experiments":N_CHAINS*draws,
        "no_direct_freight_seed_transport_mowing_wind_adaptation_genetic_heritability":True,
        "no_original_source_site_species_GPS_or_trait_value_table_exported":True
    }


def main():
    p=argparse.ArgumentParser()
    names=("cz_s1","cz_s2","it_pdf","hr_pdf","si_xlsx","si_audit",
           "height","seed","output")
    for name in names:
        p.add_argument("--"+name.replace("_","-"),required=True,type=Path)
    p.add_argument("--draws",type=int,default=DRAWS)
    a=p.parse_args()
    report=run(a.cz_s1,a.cz_s2,a.it_pdf,a.hr_pdf,a.si_xlsx,a.si_audit,
               a.height,a.seed,a.draws)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(report,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    for study,data in report["four_independent_field_publication_group_trait_label_null"].items():
        for mode,results in data["source_null_results"].items():
            print(json.dumps({
                "version":"v28.5","source_original_publication":study,
                "n_actual_sites":data["source_QA"]["source_original_sites"],
                "source_matched_taxa":data["source_QA"]["source_actual_plant_species_with_both_traits"],
                "trait_permutation_mode":mode,
                "CWM_2D_within_geography":results["metrics"]["cwm2D_mean_within_original_geography"],
                "height_component":results["metrics"]["cwm_height_absolute_mean_within_original_geography"],
                "seed_mass_component":results["metrics"]["cwm_seed_mass_absolute_mean_within_original_geography"],
                "site_alpha_trait_variance_trace":results["metrics"]["source_site_within_2D_trait_variance_trace"],
            },ensure_ascii=False))


if __name__=="__main__":
    main()
