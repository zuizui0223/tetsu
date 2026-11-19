#!/usr/bin/env python3
"""v29.7 Genuine MIREN mountain ROADS counterfactual for European RAIL v29.6.

Source-author MIREN study: 2007-2023 original independent mountain ROAD and
same-site same-area inland plant data, 167804 original observation rows,
original Zenodo MD5. NOT any RAILWAY observation or fifth rail publication.

Keep original Region×Road×Site first eligible year and BOTH roadside Plot1
and interior Plot3 incidence FROZEN. Original EU calibrated paired LEDA
canopy-height/seed-mass coordinates match the frozen four railway source
studies. Reassign bivariate trait vectors only among original species in
the SAME original V.PhyloMaker2 tree clade (root-to-tip depth .25/.50/.75).
Also run frequency-quintile NULL on the same phylo-tip-only species pool.

Original genus source study, original plot cohort, source ecology and site
identity never regenerated; unmatched phylo species are omitted ONLY after
prespecified source-support gates, never grafted. No biological effects are
created when true source plot richness/phylo exchangeability gate fails.
No original plant site×taxon or private coordinates exported.
"""
from __future__ import annotations
import argparse, hashlib, json, math
from collections import Counter, defaultdict
from pathlib import Path
from statistics import mean
import numpy as np

from audit_SI_CZ_original_continuous_trait_coverage_v220 import read_leda
from analyze_global_MIREN_road_interior_LEDA_2D_functional_beta_v249 import (
    read_source_miren, traits_from_source, regional_comparison, MIREN_MD5, MIREN_DOI
)
from validate_MIREN_global_ROAD_2D_LEDA_coverage_sensitivity_v250 import valid_region_panels
from validate_MIREN_ROAD_original_site_frozen_trait_label_counterfactual_v291 import (
    PANEL_DEFINITIONS, PINNED_PRIOR_ROAD_RESULT
)
from audit_Europe_4independent_railstudies_fixed_community_trait_label_null_v285 import (
    FROZEN_3_STUDY_LOG_MEANS, FROZEN_3_STUDY_LOG_SDS
)
from validate_Europe_3studies_82stations_original_LEDA_2D_trait_beta_v245 import LEDASHA
from audit_VPhyloMaker2_RData_python_no_R_environment_v296 import inspect_tree_rdata,TREES,COMMIT
from analyze_four_rail_actual_3phylo_clade_null_python_v296 import (
    original_root_depths, original_true_clade_groups, permute_paired_traits, source_json_native
)

MODES=("same_phylogeny_tip_original_region_species_occupancy_quintile",
       "original_megatree_monophyletic_rootdepth_cut_0.25",
       "original_megatree_monophyletic_rootdepth_cut_0.50",
       "original_megatree_monophyletic_rootdepth_cut_0.75")
FRACTIONS=(None,.25,.50,.75)
SEEDS=(297119,297120)
DRAWS_PER_CHAIN=199
MIN_DIRECT_TIP_UNIQUE_TAXON_FRACTION=.60
MIN_DIRECT_TIP_ORIGINAL_PLOT_INCIDENCE_FRACTION=.60
MIN_ORIGINAL_SITE_RETENTION_AFTER_PHYLO=.80
MIN_ROADS=2
MIN_ORIGINAL_SITES_EACH_ROAD=4
MIN_EXCHANGEABLE_CLADES=5
MIN_EXCHANGEABLE_INCIDENCE_EACH_ZONE=.25
MIN_ORIGINAL_ROADS_WITH_EXCHANGEABLE_SITE_TAXA=2
MIN_REGION_CODES_FOR_DESCRIPTIVE_AVERAGE=3
SCIENTIFIC_METRICS=('joint','height','seed')
REAL_ZLOG_REFERENCE_THREE_RAIL_ORIGINAL_STUDIES=True

def trait_pool_original_matched_regions(records,raw_LEDA_h,raw_LEDA_s):
    """Verify Europe z-axis and MIREN original source taxa original values."""
    x,qa=traits_from_source(records,raw_LEDA_h,raw_LEDA_s)
    if len(x)!=1096:raise ValueError('Original author MIREN exact dual LEDA taxa count drift')
    trait={name:(
        (math.log(raw_LEDA_h[name])-FROZEN_3_STUDY_LOG_MEANS[0])/FROZEN_3_STUDY_LOG_SDS[0],
        (math.log(raw_LEDA_s[name])-FROZEN_3_STUDY_LOG_MEANS[1])/FROZEN_3_STUDY_LOG_SDS[1],
    ) for name in x}
    if set(x)!=set(trait):raise ValueError('Original source matched taxa set changed with Europe axes')
    return x,trait,qa

def construct_original_region_tree_panel(original_roads, original_zlog_traits, tip_lookup,site_species_threshold):
    """Filter one source REGION; never split paired roadside/interior same-source site.

    Keep original Road×Site first eligible source year and exact source original
    100m² plot exposure. If phylogeny matching loses sites/roads, report loss.
    Do not use original site taxon or private coordinates in output QA.
    """
    original_source_total=sum(len(sites) for road,sites in original_roads)
    pool={s for _,sites in original_roads for row in sites
          for z in ('roadside','interior') for s in row[z]}
    exact_pool=pool & set(tip_lookup)
    total_incidence=sum(len(row[z]) for _,sites in original_roads for row in sites
                        for z in ('roadside','interior'))
    source_tip_incidence=sum(len(row[z]&exact_pool) for _,sites in original_roads
                             for row in sites for z in ('roadside','interior'))
    retained=[]
    n_qualified_2zone_sites=0
    original_years=set()
    for road,sites in original_roads:
        chosen=[]
        for row in sites:
            original_years.add(int(row['year']))
            a=row['roadside']&exact_pool
            b=row['interior']&exact_pool
            if min(len(a),len(b))<site_species_threshold:continue
            n_qualified_2zone_sites+=1
            chosen.append({'roadside':a,'interior':b,'original_year':row['year']})
        if len(chosen)>=MIN_ORIGINAL_SITES_EACH_ROAD:
            retained.append((road,chosen))
    kept_sites=sum(len(sites) for _,sites in retained)
    ratio_taxa=len(exact_pool)/len(pool) if pool else 0.
    ratio_occ=source_tip_incidence/total_incidence if total_incidence else 0.
    ratio_sites=kept_sites/original_source_total if original_source_total else 0.
    status=(len(retained)>=MIN_ROADS and
            ratio_taxa>=MIN_DIRECT_TIP_UNIQUE_TAXON_FRACTION and
            ratio_occ>=MIN_DIRECT_TIP_ORIGINAL_PLOT_INCIDENCE_FRACTION and
            ratio_sites>=MIN_ORIGINAL_SITE_RETENTION_AFTER_PHYLO)
    qa={
        'original_source_roads':len(original_roads),
        'original_source_true_spatial_sites':original_source_total,
        'original_original_trait_known_species_pool':len(pool),
        'original_direct_tree_tip_matched_species_pool':len(exact_pool),
        'source_taxa_exact_tree_tip_fraction':ratio_taxa,
        'source_roadside_and_interior_original_observation_tip_fraction':ratio_occ,
        'original_first_selected_field_survey_years_count':len(original_years),
        'original_site_pairs_with_phylogenetic_traits_above_min_after_tree_cut':n_qualified_2zone_sites,
        'retained_original_roads_with_4_or_more_sites':len(retained),
        'retained_original_source_site_pairs':kept_sites,
        'retained_original_source_paired_site_fraction':ratio_sites,
        'source_min_direct_tip_species_each_plot':site_species_threshold,
        'predeclared_min_unique_species_tip_fraction':MIN_DIRECT_TIP_UNIQUE_TAXON_FRACTION,
        'predeclared_min_observed_plant_plot_tip_fraction':MIN_DIRECT_TIP_ORIGINAL_PLOT_INCIDENCE_FRACTION,
        'predeclared_min_original_site_retention':MIN_ORIGINAL_SITE_RETENTION_AFTER_PHYLO,
        'direct_tip_road_inland_original_field_cohort_support_pass':status,
        'original_source_first_year_roadside_interior_pairs_NEVER_split':True,
        'source_original_railway_plots_n':0,
        'source_original_miren_roads_NOT_RAILWAYS':True,
    }
    if not status:
        return None,qa
    # Only retained true sites/roads should define eligible species in null.
    names=sorted({s for _,sites in retained for row in sites
                   for z in ('roadside','interior') for s in row[z]})
    trait=np.asarray([original_zlog_traits[s] for s in names],dtype=np.float64)
    position={taxon:i for i,taxon in enumerate(names)}
    road_design=[]
    occ=np.zeros(len(names),dtype=np.int64)
    for road,sites in retained:
        matrices={}
        for zone in ('roadside','interior'):
            P=np.zeros((len(sites),len(names)),dtype=np.uint8)
            for i,row in enumerate(sites):
                for taxon in row[zone]:P[i,position[taxon]]=1
            if np.any(P.sum(axis=1)<site_species_threshold):
                raise ValueError('Original source matched site species richness not preserved after tree')
            occ+=P.sum(axis=0)
            matrices[zone]=P/P.sum(axis=1)[:,None]
        road_design.append(matrices)
    if np.any(occ<=0) or not np.all(np.isfinite(trait)):
        raise ValueError('Source phylogenetic trait pool includes unobserved or invalid plant trait')
    return {'roads':road_design,'species':names,
            'source_trait_vectors':trait,'original_occupancy':occ},qa

def road_beta_three_metrics(W:np.ndarray,T:np.ndarray)->dict:
    C=W@T
    ii,jj=np.triu_indices(len(W),1)
    if len(ii)==0:raise ValueError("Original road has fewer than 2 independent field site locations")
    delta=C[ii]-C[jj]
    return {'height':float(np.abs(delta[:,0]).mean()),
            'seed':float(np.abs(delta[:,1]).mean()),
            'joint':float(np.linalg.norm(delta,axis=1).mean())}

def region_beta_for_two_zones(panel,T)->dict:
    outputs={z:{k:[] for k in SCIENTIFIC_METRICS}
             for z in ('roadside','interior')}
    for road in panel['roads']:
        for zone in ('roadside','interior'):
            d=road_beta_three_metrics(road[zone],T)
            for met,val in d.items():outputs[zone][met].append(val)
    result={zone:{metric:mean(vals) for metric,vals in data.items()}
            for zone,data in outputs.items()}
    result['roadside_minus_interior']={metric:
        result['roadside'][metric]-result['interior'][metric]
        for metric in SCIENTIFIC_METRICS}
    return result

def occupancy_quintile_groups(incidence:np.ndarray)->list[np.ndarray]:
    order=np.argsort(incidence,kind='stable')
    return [v for v in np.array_split(order,5) if len(v)>0]

def region_tree_clade_exchangeability(panel,groups)->dict:
    mobile=sorted(set(i for grp in groups if len(grp)>=2 for i in grp))
    n_groups=sum(len(grp)>=2 for grp in groups)
    total={z:0. for z in ('roadside','interior')}
    mob={z:0. for z in ('roadside','interior')}
    n_roads_mobile=0
    for road in panel['roads']:
        valid=True
        for zone in ('roadside','interior'):
            # Original W is 1 / n source taxa per site, not raw incidence. We
            # need true source species×site incidence fractions, recovered
            # exactly by counting nonzero cells rather than using centroids.
            W=road[zone]
            all_occ=float(np.count_nonzero(W))
            n_mob=float(np.count_nonzero(W[:,mobile])) if mobile else 0.
            total[zone]+=all_occ
            mob[zone]+=n_mob
            valid &= n_mob>=MIN_EXCHANGEABLE_INCIDENCE_EACH_ZONE*all_occ
        if valid:n_roads_mobile+=1
    fracs={z:mob[z]/total[z] if total[z] else 0. for z in total}
    eligible=(n_groups>=MIN_EXCHANGEABLE_CLADES and
              min(fracs.values())>=MIN_EXCHANGEABLE_INCIDENCE_EACH_ZONE and
              n_roads_mobile>=MIN_ORIGINAL_ROADS_WITH_EXCHANGEABLE_SITE_TAXA)
    return {
        'source_tree_clades_with_two_or_more_joint_LEDA_species':n_groups,
        'source_trait_known_roadside_incidence_exchangeable_fraction':fracs['roadside'],
        'source_trait_known_interior_incidence_exchangeable_fraction':fracs['interior'],
        'n_original_roads_roadside_and_interior_have_25percent_exchangeable':n_roads_mobile,
        'predeclared_phylogenetic_source_null_support':eligible,
        'source_true_site_plant_incidence_and_raw_author_road_ids_unchanged':True,
    }

def run_region_null(panel,groups,mode,tree_index,region_index,draws):
    support=region_tree_clade_exchangeability(panel,groups)
    if mode==MODES[0]:
        support['predeclared_phylogenetic_source_null_support']=True
        status='SOURCE_ORIGINAL_PHYLOTIP_OCCUPANCY_QUINTILE_REFERENCE'
    else:
        status=('SOURCE_PHYLO_CLADE_NULL_IDENTIFIABLE' if support[
            'predeclared_phylogenetic_source_null_support'] else
            'NOT_IDENTIFIABLE_ROAD_INTERIOR_PHYLOGENETIC_CLADE_NULL')
    if status.startswith('NOT_IDENTIFIABLE'):
        return {'status':status,'source_support':support,'effects':None}
    if draws<5:raise ValueError('Insufficient technical trait-label null permutations')
    values=panel['source_trait_vectors']
    obs=region_beta_for_two_zones(panel,values)
    examples={z:{key:[] for key in SCIENTIFIC_METRICS}
              for z in ('roadside','interior','roadside_minus_interior')}
    # A single trait assignment is shared by both physically paired habitats,
    # not independently randomized within the two vegetation zones.
    for seed in SEEDS:
        rng=np.random.default_rng(seed+17011*tree_index+5009*region_index+
                                  sum(ord(x) for x in mode))
        for _ in range(draws):
            p=permute_paired_traits(values,groups,rng)
            trial=region_beta_for_two_zones(panel,p)
            for z in examples:
                for k in SCIENTIFIC_METRICS:
                    examples[z][k].append(trial[z][k])
    results={}
    for z in examples:
        results[z]={}
        for metric,v in examples[z].items():
            u=mean(v);actual=obs[z][metric]
            results[z][metric]={
              'source_original_observed':actual,
              'original_tree_clade_trait_label_null_mean':u,
              'observed_minus_trait_label_source_null':actual-u,
              'technical_upper_rank_NOT_field_p':
                (1+sum(x>=actual for x in v))/(1+len(v))
            }
    return {'status':status,'source_support':support,
            'source_original_region_two_zone_effects':results,
            'n_technical_shuffles_not_field_studies':2*draws,
            'source_paired_2D_height_seed_original_trait_values_preserved':True,
            'original_site_species_plant_roadside_and_interior_FIXED':True}

def run(miren:Path,height:Path,seed:Path,tree_dir:Path,output:Path,
        draws:int=DRAWS_PER_CHAIN)->dict:
    for p,expected in ((height,LEDASHA['height']),(seed,LEDASHA['seed'])):
        if hashlib.sha256(p.read_bytes()).hexdigest()!=expected:
            raise ValueError('Oldenburg LEDA original source SHA mismatch')
    records,source_qa=read_source_miren(miren)
    source_local_traits,global_z_traits,trait_qa=trait_pool_original_matched_regions(
        records,*[read_leda(p,kind)[0] for p,kind in (
           (height,'height_m'),(seed,'seed_mass_mg'))])
    if source_qa['original_source_record_QA']['original_source_species_observation_rows']!=167804:
        raise ValueError('MIREN source 167804 original floristic observations changed')
    true_panels={}
    for panel_key,(source_min,share,minyear,nreg,nroad,nsite) in PANEL_DEFINITIONS.items():
        regions,support=valid_region_panels(records,source_local_traits,source_min,share,minyear)
        if (len(regions),sum(len(v) for v in regions.values()),
            sum(len(s) for v in regions.values() for _,s in v))!=(nreg,nroad,nsite):
            raise ValueError('Original first-year within-road MIREN field site cohort source changed')
        if panel_key=='primary_5species_25percent':
            old=mean(regional_comparison(roads,source_local_traits,reg)[
                'source_2D_CWM_beta_road_minus_interior'] for reg,roads in regions.items())
            if abs(old-PINNED_PRIOR_ROAD_RESULT)>1e-9:
                raise ValueError('Previously verified v24.9 MIREN source habitat 2D beta drift')
        true_panels[panel_key]={'regions':regions,'original_field_support':support,
                                'source_plot_min_direct_tip':source_min}
    all_tree={}
    for tree_index,tree in enumerate(TREES):
        original,tree_qa=inspect_tree_rdata(tree_dir/('GBOTB.extended.'+tree+'.rda'),tree)
        tips=[str(v).replace('_',' ') for v in np.asarray(original['tip.label']).reshape(-1)]
        if len(set(tips))!=len(tips):raise ValueError('Original published phylo normalized tip labels duplicate')
        lookup={t:i+1 for i,t in enumerate(tips)}
        parent,depth,depth_qa=original_root_depths(
            np.asarray(original['edge']),np.asarray(original['edge.length']),
            len(tips),int(np.asarray(original['Nnode']).reshape(-1)[0]))
        all_tree[tree]={
            'original_published_phylogenetic_tree_QA':tree_qa,
            'original_published_tree_branch_length_QA':depth_qa,
            'panels':{}
        }
        for panel_key,panel in true_panels.items():
            results={}
            for region_index,(region,source_roads) in enumerate(sorted(panel['regions'].items())):
                matched,cohort_qa=construct_original_region_tree_panel(
                    source_roads,global_z_traits,lookup,panel['source_plot_min_direct_tip'])
                if matched is None:
                    results[region]={'source_tree_support':cohort_qa,
                                     'status':'NOT_IDENTIFIABLE_ORIGINAL_DIRECT_TIP_FIELD_COHORT',
                                     'conditional_trait_label_nulls':None}
                    continue
                source_taxa=matched['species']
                original_positions=[lookup[t] for t in source_taxa]
                occupation=occupancy_quintile_groups(matched['original_occupancy'])
                modes={}
                for i,mode in enumerate(MODES):
                    fraction=FRACTIONS[i]
                    groups=(occupation if i==0 else original_true_clade_groups(
                        original_positions,parent,depth,
                        fraction*depth_qa['original_total_root_to_tip_distance_median']))
                    modes[mode]=run_region_null(matched,groups,mode,tree_index,region_index,draws)
                results[region]={'source_tree_support':cohort_qa,
                                 'status':'ORIGINAL_DIRECT_TIP_MIREN_FIELD_COHORT_SUPPORTED',
                                 'conditional_trait_label_nulls':modes}
            all_tree[tree]['panels'][panel_key]={'source_original_field_support':panel['original_field_support'],
                                                'original_source_region_results':results}
    summaries={}
    for tree,x in all_tree.items():
        summaries[tree]={}
        for panel_key,dat in x['panels'].items():
            mode_reports={}
            for mode in MODES:
                signed={z:{k:{} for k in SCIENTIFIC_METRICS}
                        for z in ('roadside','interior','roadside_minus_interior')}
                unsupported=[]
                for region,source in dat['original_source_region_results'].items():
                    modes=source['conditional_trait_label_nulls']
                    if modes is None or modes[mode]['effects'] is None:
                        unsupported.append(region);continue
                    for zone in signed:
                        for metric in SCIENTIFIC_METRICS:
                            signed[zone][metric][region]=modes[mode][
                                'source_original_region_two_zone_effects'][zone][metric][
                                'observed_minus_trait_label_source_null']
                n=len(signed['roadside']['height'])
                desc=None
                if n>=MIN_REGION_CODES_FOR_DESCRIPTIVE_AVERAGE:
                    desc={z:{m:{
                        'original_region_code_equal_mean_source_excess':mean(signed[z][m].values()),
                        'original_source_regions_signed_positive':sum(v>0 for v in signed[z][m].values()),
                        'source_original_region_codes_with_evaluable_phylogeny':n,
                        'by_original_source_region_code':signed[z][m],
                        'region_leave_one_out_minmax_not_global_CI':[
                            min(mean(v for key,v in signed[z][m].items() if key!=drop)
                                for drop in signed[z][m]),
                            max(mean(v for key,v in signed[z][m].items() if key!=drop)
                                for drop in signed[z][m])
                        ],
                    } for m in SCIENTIFIC_METRICS} for z in signed}
                mode_reports[mode]={
                    'source_original_region_codes_with_results':n,
                    'source_original_region_codes_without_phylogeny_support':unsupported,
                    'region_equal_descriptive_effects_if_3plus':desc,
                    'not_independent_field_publications_n':1,
                }
            summaries[tree][panel_key]=mode_reports
    out={
        'version':'v29.7',
        'science':'MIREN original independent NON-RAIL mountain road field program, roadside vs same-site interior source trait-geography benchmark to 4 EU railway study phylogenetic null',
        'source_original_MIREN':source_qa,
        'source_original_LEDA_coverage_QA':trait_qa,
        'source_phylogenetic_backbones_original_author_git_commit':COMMIT,
        'original_tree_nomenclatures_not_three_independent_field_studies':list(TREES),
        'source_Europe_three_publications_fixed_log_z_axis_used':True,
        'source_original_first_year_same_site_same_100m2_roadside_interior_field_sampling_frozen':True,
        'original_MIREN_independent_field_programs_n':1,
        'source_panels_not_independent_field_programs_n':len(PANEL_DEFINITIONS),
        'source_predeclared_direct_tip_min_taxa_per_original_plot_and_retention_applied':True,
        'source_phylo_clade_source_eligibility_independent_of_effect_sign':True,
        'original_tree_real_bivariate_height_seed_source_unchanged':True,
        'original_source_plant_taxon_raw_values_individual_site_IDs_GPS_NOT_REEXPORTED':True,
        'not_original_railway_surveys':True,
        'no_local_evolution_or_traffic_seed_dispersal_selection_mechanism_identified':True,
        'source_tree_original_per_region_support_and_effects':all_tree,
        'source_region_equal_descriptive_summaries_only_if_supported':summaries,
    }
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(out,indent=2,ensure_ascii=False,default=source_json_native),encoding='utf8')
    print(json.dumps({
        'source_original_provenance':'SUCCESS_ROADS_NOT_RAIL',
        'source_original_field_8regions_22roads_272site_vs_strict_5regions_13roads_192site':True,
        'by_tree_source_phylo_support_and_stature':{
            tree:{
              panel:{mode:{
                'n_region_codes':data['source_original_region_codes_with_results'],
                'roadside_height':data['region_equal_descriptive_effects_if_3plus']['roadside']['height']['original_region_code_equal_mean_source_excess']
                   if data['region_equal_descriptive_effects_if_3plus'] else None,
                'interior_height':data['region_equal_descriptive_effects_if_3plus']['interior']['height']['original_region_code_equal_mean_source_excess']
                   if data['region_equal_descriptive_effects_if_3plus'] else None,
                'road_minus_interior_height_excess':data['region_equal_descriptive_effects_if_3plus']['roadside_minus_interior']['height']['original_region_code_equal_mean_source_excess']
                   if data['region_equal_descriptive_effects_if_3plus'] else None}
                   for mode,data in modes.items()}
              for panel,modes in panels.items()}
              for tree,panels in summaries.items()}
        },indent=2,ensure_ascii=False,default=source_json_native))
    return out

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--miren',type=Path,required=True)
    p.add_argument('--height',type=Path,required=True)
    p.add_argument('--seed',type=Path,required=True)
    p.add_argument('--tree-dir',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--draws-per-chain',type=int,default=DRAWS_PER_CHAIN)
    args=p.parse_args()
    run(args.miren,args.height,args.seed,args.tree_dir,args.output,args.draws_per_chain)
if __name__=='__main__':
    main()
