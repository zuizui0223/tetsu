#!/usr/bin/env python3
"""v25.0 robust source-year/site-specific MIREN mountain ROAD functional
beta sensitivity to dramatic global LEDA height/seed missingness.

Source: authentic 2025 MIREN Zenodo 15749536 global 25+ original
region codes, 2007-2023 2×50m actual roadside plot1 and interior
plot3, true 167804 plant records, true source MD5. Unique first survey
per source Region×Road×Site. NOT railway field observations.

v24.9 primary: 1096 / 5739 MIREN source Latin taxa have joint
original LEDA data (19.1%), so 8 source regions / 22 roads / 272 sites
at min>=5 joint-LEDA species each and trait coverage >=25%; 7/8 regions
show negative road-minus-interior functional CWM beta. This is sample
selection risk, not a 25-region global effect.

This sensitivity reuses same original actual source CSV and LEDA,
compares all nested metadata-only support panels:
  A: >=5 species and >=25% in each original plot (v24.9)
  B: >=5 and >=50%
  C: >=10 and >=25%
  D: >=10 and >=50%
  E: >=5 >=25% but genuine survey year >=2015 (no extrapolation
     from 2007 baseline to 2025 traffic intensities).
Require each regional source cohort >=2 distinct roads with
>=4 separately surveyed spatial sites per road.
Don't swap dates, don't include roads/rail in the same study class.
For each strictly original region, compare same road actual site-pairs
roadside CWM distance vs interior, same 2D species median log trait
axes and all original taxon presence; rank site pairs not independent
field replicates. No positive-significance fishing, no formal p-value.

Report at least original study region support and n signed regions,
and each region leave-one-out descriptive sign (not a field CI).
The preexisting v24.9 result uses local 59 repeated rarefaction
draws, NOT independent site bootstrap. Here thresholding is its
separate sensitivity, not a new experimental dataset.
"""
from __future__ import annotations
import argparse,hashlib,json,math
from collections import Counter,defaultdict
from pathlib import Path
from statistics import mean,median

from audit_SI_CZ_original_continuous_trait_coverage_v220 import read_leda
from validate_Europe_3studies_82stations_original_LEDA_2D_trait_beta_v245 import LEDASHA
from analyze_global_MIREN_road_interior_LEDA_2D_functional_beta_v249 import (
   read_source_miren,traits_from_source,regional_comparison,MIREN_DOI,MIREN_STUDY_DOI
)

PANELS={
    "primary_5taxa_25pct":(5,.25,2007),
    "strict_5taxa_50pct":(5,.50,2007),
    "strict_10taxa_25pct":(10,.25,2007),
    "strict_10taxa_50pct":(10,.50,2007),
    "contemporary_5taxa_25pct_2015plus":(5,.25,2015),
}
MIN_ROADS_PER_REGION=2
MIN_SITES_PER_ROAD=4

def valid_region_panels(records,trait_values,min_taxa,min_share,earliest_year):
    by=defaultdict(dict)
    for (region,road,site,year,zone),taxa in records.items():
        if year<earliest_year or zone not in (1,3):continue
        by[region,road,site,year][zone]=taxa
    earliest={}
    source_year_pairs=Counter()
    qualified=Counter()
    validtraits=set(trait_values)
    for (reg,road,site,year),z in sorted(by.items()):
        if 1 not in z or 3 not in z:continue
        source_year_pairs[reg]+=1
        one=z[1]&validtraits;three=z[3]&validtraits
        if len(one)<min_taxa or len(three)<min_taxa:continue
        if len(one)/len(z[1])<min_share or len(three)/len(z[3])<min_share:continue
        key=(reg,road,site)
        if key not in earliest or year<earliest[key]['year']:
            earliest[key]={'region':reg,'road':road,'site':site,'year':year,
                 'roadside':one,'interior':three}
        qualified[reg]+=1
    byroad=defaultdict(list)
    for row in earliest.values():
        byroad[(row['region'],row['road'])].append(row)
    eligible=defaultdict(list)
    for (region,road),items in byroad.items():
        if len(items)>=MIN_SITES_PER_ROAD:
            eligible[region].append((road,items))
    selected={region:roads for region,roads in eligible.items()
              if len(roads)>=MIN_ROADS_PER_REGION}
    return selected,{
      'threshold_joint_LEDA_species_each_plot':min_taxa,
      'threshold_joint_LEDA_fraction_each_plot':min_share,
      'source_first_allowed_survey_year':earliest_year,
      'original_source_regions_with_any_true_1vs3_site_year':len(source_year_pairs),
      'source_region_true_plot_years_with_both_1_3':dict(source_year_pairs),
      'source_region_siteyears_passed_trait_coverage':dict(qualified),
      'source_original_true_unique_site_earliest_eligible_count':len(earliest),
      'source_regional_2road_4site_support_n':len(selected),
      'source_eligible_distinct_roads':sum(len(v) for v in selected.values()),
      'source_eligible_spatial_sites':sum(len(x) for gr in selected.values() for r,x in gr),
      'original_source_region_codes_eligible':sorted(selected)}

def describe_panel(regions,values):
    if not regions:
        return {'identified':False,'status':'NO_INDEPENDENT_SOURCE_REGION_SUPPORT','not_zero_functional_road_effect':True}
    observed={reg:regional_comparison(roads,values,reg) for reg,roads in sorted(regions.items())}
    d={reg:r['source_2D_CWM_beta_road_minus_interior'] for reg,r in observed.items()}
    shared={reg:r['source_species_Jaccard_road_minus_interior'] for reg,r in observed.items()}
    if len(d)>=4:
        leave={reg:mean(x for k,x in d.items() if k!=reg) for reg in d}
    else:leave={}
    return {
       'identified':True,
       'n_source_original_region_codes':len(d),
       'source_original_n_true_roads':sum(v['n_original_roads'] for v in observed.values()),
       'source_original_unique_spatial_sites':sum(v['n_original_spatial_sites'] for v in observed.values()),
       'region_equal_weight_original_CWM_roadside_minus_interior_mean':mean(d.values()),
       'median_region_roadside_minus_interior_CWM':median(d.values()),
       'n_region_codes_roadside_MORE_functionally_converged_negative':sum(x<0 for x in d.values()),
       'n_region_codes_roadside_MORE_functionally_differentiated_positive':sum(x>0 for x in d.values()),
       'source_original_region_functional_beta_effects':d,
       'source_tax_Jaccard_roadside_minus_interior_region_equal_mean':mean(shared.values()),
       'original_region_leave_one_region_out_mean_range':[
         min(leave.values()),max(leave.values())] if leave else None,
       'all_original_region_LOO_means_negative':
         bool(leave) and all(x<0 for x in leave.values()),
       'small_source_study_region_codes_not_global_population_random_sample':True}

def run(miren,height,seed):
    plant,source=read_source_miren(miren)
    values,traits=traits_from_source(plant,height,seed)
    panels={}
    for panel,th in PANELS.items():
        groups,support=valid_region_panels(plant,values,*th)
        effects=describe_panel(groups,values)
        panels[panel]={'source_original_study_support_QA':support,'descriptive_road_functional_beta':effects}
    primary=panels['primary_5taxa_25pct']
    if primary['source_original_study_support_QA']['source_eligible_spatial_sites']!=272 or primary[
       'source_original_study_support_QA']['source_regional_2road_4site_support_n']!=8:
        raise ValueError('Genuine v24.9 original same-source baseline 272 sites/8 regions drifted')
    actual=primary['descriptive_road_functional_beta']['region_equal_weight_original_CWM_roadside_minus_interior_mean']
    if abs(actual+0.10381093795211047)>1e-9:
        raise ValueError('True original v24.9 global MIREN road source CWM result drifted')
    return {'version':'v25.0',
       'status':'SOURCE_VERIFIED_GLOBAL_MIREN_MOUNTAIN_ROADS_2D_TRAIT_MISSINGNESS_SENSITIVITY_NOT_RAILWAYS',
       'source_original_MIREN_Zenodo_v3_MD5':source['original_csv_MD5'],
       'source_original_MIREN_DOI':MIREN_DOI,
       'source_original_2026_MIREN_publication_DOI':MIREN_STUDY_DOI,
       'n_original_source_regions_in_csv':source['n_source_original_regions'],
       'fraction_original_MIREN_source_species_with_both_LEDA_traits':traits['MIREN_overall_species_coverage_fraction'],
       'n_source_exact_LEDA_joint_species':traits['original_source_LEDA_2D_joint_species'],
       'source_region_specific_2D_functional_beta_support_and_sign_by_predeclared_threshold':panels,
       'not_reclassify_MIREN_roadside_plots_as_railway_plots':True,
       'explanatory_limitations':[
         'Original 25-region MIREN publication and 26 region codes in CSV must not be silently treated as 26 independent ecosystems',
         'Source 1096/5739 species =19 percent 2D LEDA coverage, due geographic Europe-biased traits and taxonomic availability',
         'Roadside vs interior are 100m2 source field plots, not track ballast vs nonrail, even when same T transect',
         'Road and interior source presence data selected for both metadata/surviving LEDA taxa; later years do not equal independent source observations',
         'Strongest sensitivity may lose most tropical/cold study regions; no global inference from selected subset',
         'Original regional road/source site nestedness, each road pair dependence and MCMC absent do not allow independent pairwise p-values',
         'No direct railway train movement, seed rain, growing height plasticity or adaptation tested']}
def main():
    p=argparse.ArgumentParser()
    for key in ('miren','height','seed','output'):
        p.add_argument('--'+key,required=True)
    a=p.parse_args()
    for f,s in ((a.height,LEDASHA['height']),(a.seed,LEDASHA['seed'])):
        if hashlib.sha256(Path(f).read_bytes()).hexdigest()!=s:
            raise ValueError('Original Oldenburg LEDA source SHA changed')
    height,hqa=read_leda(a.height,'height_m');seed,sqa=read_leda(a.seed,'seed_mass_mg')
    report=run(a.miren,height,seed)
    out=Path(a.output);out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(report,indent=2,ensure_ascii=False),encoding='utf8')
    compact={'status':report['status'],
      'source_n_region_codes':report['n_original_source_regions_in_csv'],
      'joint_LEDA_fraction':report['fraction_original_MIREN_source_species_with_both_LEDA_traits'],
      'panels':{name:{
         'n_regions':p['source_original_study_support_QA']['source_regional_2road_4site_support_n'],
         'n_roads':p['source_original_study_support_QA']['source_eligible_distinct_roads'],
         'n_sites':p['source_original_study_support_QA']['source_eligible_spatial_sites'],
         'results':p['descriptive_road_functional_beta']}
         for name,p in report['source_region_specific_2D_functional_beta_support_and_sign_by_predeclared_threshold'].items()}}
    print('V250_MIREN_ROAD_TRUE_TRAIT_MISSINGNESS_SENSITIVITY='+json.dumps(compact,ensure_ascii=False,separators=(',',':')),flush=True)
if __name__=='__main__':main()
