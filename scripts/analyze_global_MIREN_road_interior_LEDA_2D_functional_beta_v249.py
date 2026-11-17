#!/usr/bin/env python3
"""v24.9: global SAME-PROTOCOL 25-region ROAD vs interior benchmark for
transport corridor vegetation; intentionally NOT a fifth railway study.

Actual source MIREN survey vascular plant species in mountain ROADS:
Seipel, Haider & MIREN consortium, Zenodo v3 15749536, original source
CSV MIREN_road_survey_(2007-2023)_zen_v3_002507.csv MD5
b640195cc44f252ccb6ddb28eb6fff13 (17.9 MB).
2026 GEB paper DOI 10.1111/geb.70284, 25 regions globally,
~170,000 species records, 3364 plots, 6854 taxa, 2007-2023.

Each true MIREN site contains exactly 2 × 50 m roadside plot 1,
intermediate plot 2 and interior plot 3; individual source plant
taxa in each plot and survey year. The method is related to (but NOT
identical to) Rashid et al 2021 India Himalayan 31 Railway T transects.
Do NOT imply MIREN roads are railways or railway-specific effects.

Source quality: original Region/Road/Site/Plot/Year, Accepted.Name.MIREN,
Plot_id and Status; preserve original site ID & year. Select FIRST
source year per stable (Region,Road,Site) in which both 1 and 3
true plant lists have sufficient JOINT LEDA 2D (canopy height,
seed mass), avoiding repeated site-year pseudoreplication.
Original roadside (zone1) vs interior (zone3) same 100m2 area, but
different vegetation; first usable year may vary across real sites.

Within each region, SAME ROAD and SAME PAIRS OF SOURCE SITES
estimate mean CWM pairwise Euclidean distance of roadside plot1
minus interior plot3. Average within original road, then region;
regional replicates are units, not thousands of site pairs.
Positive=roadside more functionally differentiated, negative=roadside
relative convergence. Also trace within source plot1 minus plot3,
and source species Jaccard between sites, separately.
No inference to plant evolution or train transport.

Sensitivity: original source raw richness matched by technical
species rarefaction of joint-LEDA species in each SAME site to
min(roadside known, interior known), 59 technical draws, then
repeat same regional paired-beta summaries. These are sampling
conditional draws, NOT independent ecological confidence intervals.

Strict >=5 species each (source cover fraction >=0.25), regional
support >=2 genuine roads each with >=4 unique original sites;
report missingness explicitly, never fill source plots with zero
or impute unknown LEDA values. No source coordinates, species×plot
rows or raw MIREN table redistributed in GitHub.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,math,random,re
from collections import Counter,defaultdict
from itertools import combinations
from pathlib import Path
from statistics import mean,median

from audit_SI_CZ_original_continuous_trait_coverage_v220 import read_leda,binomial
from analyze_Slovenia_original_LEDA_height_seed_breadth_v224 import components,quart
from validate_Europe_3studies_82stations_original_LEDA_2D_trait_beta_v245 import LEDASHA

MIREN_MD5="b640195cc44f252ccb6ddb28eb6fff13"
MIREN_DOI="10.5281/zenodo.15749536"
MIREN_STUDY_DOI="10.1111/geb.70284"
SOURCE_FIELDS=("Plot_id","Region","Road","Site","Plot","Year","Accepted.Name.MIREN")
MIN_LEDA=5
MIN_TRAIT_SHARE=0.25
MIN_ROADS_PER_REGION=2
MIN_SITES_PER_ROAD=4
RARE_DRAWS=59
SEED=20261117+249

def read_source_miren(source):
    raw=Path(source).read_bytes()
    if hashlib.md5(raw).hexdigest()!=MIREN_MD5:
        raise ValueError("MIREN true original 2025 Zenodo source CSV MD5 changed")
    records=defaultdict(set);sourceids=defaultdict(set);qa=Counter()
    country_groups=defaultdict(Counter)
    regions=set();rawspecies=set();siteyears=set();yearvalues=set()
    with open(source,encoding="utf-8-sig",newline="") as f:
        rdr=csv.DictReader(f)
        absent=set(SOURCE_FIELDS)-set(rdr.fieldnames or [])
        if absent:raise ValueError("MIREN published original CSV lacks columns "+str(sorted(absent)))
        for line in rdr:
            qa["original_source_species_observation_rows"]+=1
            region=line["Region"].strip();road=line["Road"].strip()
            site=line["Site"].strip();year=line["Year"].strip()
            zone=line["Plot"].strip();pid=line["Plot_id"].strip()
            if not all((region,road,site,year,zone,pid)):
                qa["missing_original_geographical_or_plot_year_source_codes"]+=1;continue
            if not year.isdigit() or not 2007<=int(year)<=2023:
                qa["invalid_true_original_survey_year"]+=1;continue
            # Verified MIREN protocol has 1 roadside, 2 intermediate,
            # 3 interior, and suffix A indicates replaced source plot.
            m=re.fullmatch(r"([123])([A-Za-z]?)",zone)
            if not m:
                qa["invalid_source_plot_type_or_unexpected_suffix"]+=1;continue
            typ=int(m.group(1))
            if not pid.endswith("."+zone):
                qa["source_plot_id_not_consistent_with_author_plot_column"]+=1
                continue
            name=line["Accepted.Name.MIREN"].strip()
            species=binomial(name)
            if species is None:
                qa["source_nonbinomial_species_name_not_assigned_trait"]+=1
                continue
            key=(region,road,site,int(year),typ)
            records[key].add(species)
            sourceids[key].add(pid)
            rawspecies.add(species);regions.add(region);yearvalues.add(int(year))
            siteyears.add((region,road,site,int(year)))
            country_groups[region][typ]+=1
    if not 80000<=qa["original_source_species_observation_rows"]<=300000:
        raise ValueError("True original MIREN global 170k rows not present")
    if len(regions)<20 or len(rawspecies)<3000:
        raise ValueError("Source MIREN genuine global 25 regions/plant taxa coverage unexpectedly sparse")
    # A duplicate true source row of same plot-year/taxon is collapsed
    # to presence; two source plot IDs of same 1/2/3 zone in a year
    # are ambiguous and discarded rather than falsely combined.
    multi={k for k,v in sourceids.items() if len(v)!=1}
    for key in multi:records.pop(key,None)
    return records,{
       'verified_source_DOI':MIREN_DOI,'original_csv_MD5':MIREN_MD5,
       'source_2026_PAPER_DOI':MIREN_STUDY_DOI,
       'source_original_csv_bytes':len(raw),
       'original_source_record_QA':dict(qa),
       'n_source_original_regions':len(regions),
       'n_source_unique_original_plant_binomials':len(rawspecies),
       'n_source_true_site_years':len(siteyears),
       'n_source_possible_plot_year_species_sets':len(sourceids),
       'n_source_ambiguous_duplicate_same_plot_type_year':len(multi),
       'source_original_year_range':[min(yearvalues),max(yearvalues)],
       'n_region_original_siteyear_observations':dict(Counter(k[0] for k in siteyears)),
       'source_plot_1_2_3_are_ROADS_NOT_RAILWAYS':True}

def traits_from_source(records,h,s):
    names=set().union(*records.values())
    known={species:(math.log(h[species]),math.log(s[species]))
           for species in names if species in h and species in s}
    if len(known)<250:raise ValueError('Too few unambiguous LEDA 2D road study taxa')
    mu=[mean(v[k] for v in known.values()) for k in (0,1)]
    sd=[math.sqrt(mean((v[k]-mu[k])**2 for v in known.values())) for k in (0,1)]
    if min(sd)<.02:raise ValueError('Non-usable actual 2D global source trait axis')
    return {t:((v[0]-mu[0])/sd[0],(v[1]-mu[1])/sd[1]) for t,v in known.items()},{
       'original_MIREN_unique_plant_binomials':len(names),
       'original_source_LEDA_2D_joint_species':len(known),
       'MIREN_overall_species_coverage_fraction':len(known)/len(names),
       'region_independent_same_2D_axes':True,
       'source_trait_zlog_mean':mu,'source_trait_zlog_std':sd}

def complete_paired_first_surveys(records,values):
    by=defaultdict(dict)
    for key,taxa in records.items():
        region,road,site,year,plot=key
        if plot in (1,3):by[region,road,site,year][plot]=taxa
    possible=Counter();qual=Counter();selected={}
    for (region,road,site,year),zones in sorted(by.items()):
        if not all(z in zones for z in (1,3)):continue
        possible[region]+=1
        a=zones[1]&set(values);b=zones[3]&set(values)
        if len(a)<MIN_LEDA or len(b)<MIN_LEDA:continue
        if len(a)/len(zones[1])<MIN_TRAIT_SHARE or len(b)/len(zones[3])<MIN_TRAIT_SHARE:continue
        base=(region,road,site)
        if base not in selected or year<selected[base]['year']:
            selected[base]={'region':region,'road':road,'site':site,
                'year':year,'roadside':set(a),'interior':set(b),
                'observed_source_species_roadside':len(zones[1]),
                'observed_source_species_interior':len(zones[3])}
        qual[region]+=1
    supported=defaultdict(list)
    for row in selected.values():
        supported[row['region'],row['road']].append(row)
    geo=defaultdict(list)
    for (region,road),items in supported.items():
        if len(items)>=MIN_SITES_PER_ROAD:geo[region].append((road,items))
    primary={reg:groups for reg,groups in geo.items()
             if len(groups)>=MIN_ROADS_PER_REGION}
    qa={'original_region_year_matched_zone1_zone3_opportunities':dict(possible),
        'source_region_year_trait_joint_qualified':dict(qual),
        'qualified_unique_original_stable_road_site_first_year_only':len(selected),
        'n_distinct_source_regions_with_at_least_two_4site_roads':len(primary),
        'source_region_qualifying_road_counts':{
          reg:len(gr) for reg,gr in sorted(primary.items())},
        'n_main_roadside_vs_interior_paired_site_records':sum(
          len(sites) for groups in primary.values() for road,sites in groups),
        'original_first_eligible_year_used_once_per_region_road_site':True}
    return primary,qa

def source_profile(taxa,traits):
    obj=components(taxa,traits)
    return {'centroid':(obj['height_mean'],obj['seed_mean']),
            'trace':obj['joint_trace'],'taxa':taxa}

def regional_comparison(sites,values,region_name=None,rng=None,rarefy=False):
    # sites = list (road, matching original site rows) for one REGION,
    # NEVER compare sites from different roads/years without controls.
    roadout=[];zone_rich=[]
    for road,original in sites:
        src=[]
        for row in original:
            a=row['roadside'];b=row['interior']
            if rarefy:
                n=min(len(a),len(b))
                a=set(rng.sample(sorted(a),n))
                b=set(rng.sample(sorted(b),n))
            src.append((source_profile(a,values),source_profile(b,values)))
            zone_rich.append((len(a),len(b)))
        d={'road_beta_distance':[],'interior_beta_distance':[],
           'road_tax_Jaccard':[],'interior_tax_Jaccard':[]}
        for (ai,bi),(aj,bj) in combinations(src,2):
            d['road_beta_distance'].append(math.dist(ai['centroid'],aj['centroid']))
            d['interior_beta_distance'].append(math.dist(bi['centroid'],bj['centroid']))
            d['road_tax_Jaccard'].append(len(ai['taxa']&aj['taxa'])/len(ai['taxa']|aj['taxa']))
            d['interior_tax_Jaccard'].append(len(bi['taxa']&bj['taxa'])/len(bi['taxa']|bj['taxa']))
        roadout.append({
          'sitewise_trace_road_minus_interior':mean(a['trace']-b['trace'] for a,b in src),
          'roadside_between_sites_beta':mean(d['road_beta_distance']),
          'interior_between_sites_beta':mean(d['interior_beta_distance']),
          'beta_road_minus_interior':mean(d['road_beta_distance'])-mean(d['interior_beta_distance']),
          'jaccard_road_minus_interior':mean(d['road_tax_Jaccard'])-mean(d['interior_tax_Jaccard']),
          'n_true_source_site_pairs_within_road':len(d['road_beta_distance'])})
    return {
       'n_original_roads':len(roadout),
       'n_original_spatial_sites':sum(len(s) for road,s in sites),
       'source_2D_CWM_beta_road_minus_interior':mean(x['beta_road_minus_interior'] for x in roadout),
       'source_alpha_2D_trace_road_minus_interior':mean(x['sitewise_trace_road_minus_interior'] for x in roadout),
       'source_species_Jaccard_road_minus_interior':mean(x['jaccard_road_minus_interior'] for x in roadout),
       'roadsite_cwm_distance':mean(x['roadside_between_sites_beta'] for x in roadout),
       'interior_cwm_distance':mean(x['interior_between_sites_beta'] for x in roadout),
       'source_site_LEDA_species_median_road':median(p[0] for p in zone_rich),
       'source_site_LEDA_species_median_interior':median(p[1] for p in zone_rich),
       'source_dependent_site_pairs_dont_count_as_independent_road_experiments':sum(
           q['n_true_source_site_pairs_within_road'] for q in roadout)}

def run(miren,h,seed,draws=RARE_DRAWS):
    rec,qsource=read_source_miren(miren)
    traits,traitQA=traits_from_source(rec,h,seed)
    regions,qa=complete_paired_first_surveys(rec,traits)
    base={reg:regional_comparison(group,traits,reg)
          for reg,group in sorted(regions.items())}
    if not base:
        return {'version':'v24.9','status':'GLOBAL_ROAD_SOURCE_TRAIT_COVERAGE_NOT_ENOUGH',
                'original_source_QA':qsource,'source_trait_QA':traitQA,
                'source_roadside_interior_paired_QA':qa,
                'not_evidence_that_railways_lack_effect':True}
    signed=[z['source_2D_CWM_beta_road_minus_interior'] for z in base.values()]
    rng=random.Random(SEED)
    drawsamples=[]
    for draw in range(draws):
        rr={reg:regional_comparison(group,traits,reg,rng,rarefy=True)
            for reg,group in regions.items()}
        drawsamples.append({
            'road_beta_mean_equal_region':mean(v['source_2D_CWM_beta_road_minus_interior'] for v in rr.values()),
            'alpha_trace_mean_equal_region':mean(v['source_alpha_2D_trace_road_minus_interior'] for v in rr.values()),
            'jaccard_mean_equal_region':mean(v['source_species_Jaccard_road_minus_interior'] for v in rr.values()),
            'n_regions_beta_negative':sum(x['source_2D_CWM_beta_road_minus_interior']<0 for x in rr.values())})
    return {'version':'v24.9','status':'SOURCE_VERIFIED_GLOBAL_25REGION_ROAD_MIREN_FUNCTIONAL_BETA_NOT_RAILWAY',
       'original_MIREN_road_study_DOI':MIREN_STUDY_DOI,
       'original_Zenodo_road_data_DOI':MIREN_DOI,
       'raw_source_QA':qsource,'MIREN_true_height_seed_traits':traitQA,
       'original_exact_spatial_pairs_and_trait_coverage_QA':qa,
       'n_regions_in_primary_functional_beta_source_support':len(base),
       'original_road_groups_as_independent_within_region_units':sum(v['n_original_roads'] for v in base.values()),
       'source_independent_spatial_site_n':sum(v['n_original_spatial_sites'] for v in base.values()),
       'descriptive_original_road_vs_interior_functional_beta_by_region':base,
       'unrarefied_original_equal_region_summary':{
         'n_roadside_more_differentiated_regions':sum(x>0 for x in signed),
         'n_roadside_more_converged_regions':sum(x<0 for x in signed),
         'mean_road_minus_interior_CWM_beta':mean(signed),
         'median_road_minus_interior_CWM_beta':median(signed),
         'mean_road_minus_interior_alpha_trace':mean(v['source_alpha_2D_trace_road_minus_interior'] for v in base.values()),
         'mean_road_minus_interior_taxonomic_Jaccard':mean(v['source_species_Jaccard_road_minus_interior'] for v in base.values())},
       'technical_same_site_species_richness_rarefaction_by_equal_region':{
         k:{'mean':mean(z[k] for z in drawsamples),
            'q025_technical_NOT_region_CI':quart([z[k] for z in drawsamples],.025),
            'q975_technical_NOT_region_CI':quart([z[k] for z in drawsamples],.975)}
         for k in drawsamples[0]},
       'source_roads_are_not_authentic_railroads':True,
       'source_land_interior_plot_is_not_independent_region_wide_natural_pool':True,
       'limitations':[
         'Global MIREN mountain ROADS not railway substrate; cannot add 25 regions as independent railway field studies',
         'Only first survey per road/site included; selection by LEDA joint trait availability not a random regional sample',
         'Same region original road sites paired only to their own road interior 2x50m subplot; station pairs and roads within region dependent',
         'Outside Europe and tropics may have poor European Oldenburg species median trait coverage; no genus-level imputation',
         'Nonnative plants alone may dominate transport effects; future need nativestatus separate if supported',
         'Source original elevations/climate/years not regressed out beyond same-site same-year paired road vs interior design',
         'Technical rarefaction draws hold original observed regional plants fixed and are not new independent biological confidence intervals',
         'No train frequency, freight seed dispersal, railway mowing or within-species railway adaptation']}
def main():
    p=argparse.ArgumentParser()
    for key in ('miren','height','seed','output'):
        p.add_argument('--'+key,required=True)
    p.add_argument('--draws',type=int,default=RARE_DRAWS)
    args=p.parse_args()
    for file,dig in ((args.height,LEDASHA['height']),(args.seed,LEDASHA['seed'])):
        if hashlib.sha256(Path(file).read_bytes()).hexdigest()!=dig:
            raise ValueError('Original university LEDA source checksum changed')
    h,qa=read_leda(args.height,'height_m')
    s,qb=read_leda(args.seed,'seed_mass_mg')
    data=run(args.miren,h,s,args.draws)
    out=Path(args.output);out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(data,ensure_ascii=False,indent=2),encoding="utf8")
    compact={
      'status':data['status'],
      'raw_source_QA':data.get('raw_source_QA'),
      'traitQA':data.get('MIREN_true_height_seed_traits') or data.get('source_trait_QA'),
      'paired_QA':data.get('original_exact_spatial_pairs_and_trait_coverage_QA') or data.get('source_roadside_interior_paired_QA'),
      'n_regions':data.get('n_regions_in_primary_functional_beta_source_support'),
      'n_roads':data.get('original_road_groups_as_independent_within_region_units'),
      'n_source_sites':data.get('source_independent_spatial_site_n'),
      'actual_region_results':data.get('descriptive_original_road_vs_interior_functional_beta_by_region'),
      'unrarefied':data.get('unrarefied_original_equal_region_summary'),
      'rarefied':data.get('technical_same_site_species_richness_rarefaction_by_equal_region')}
    print('V249_TRUE_GLOBAL_MIREN_ROAD_NOT_RAILWAY='+json.dumps(compact,ensure_ascii=False,separators=(',',':')),flush=True)
if __name__=='__main__':main()
