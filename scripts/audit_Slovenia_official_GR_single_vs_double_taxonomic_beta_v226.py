#!/usr/bin/env python3
"""v22.6: G/R railway vs SINGLE/DOUBLE confounding and line-pair Jaccard.

Existing authenticated original SI source and LEDA SHA as in v22.4/v22.5.
Prevents naive attribution of national main/regional difference to train
departures when the regional source sites are all single track.
"""
from __future__ import annotations
import argparse,hashlib,json,math
from collections import defaultdict
from itertools import combinations
from pathlib import Path
from statistics import mean
from audit_SI_CZ_original_continuous_trait_coverage_v220 import read_leda,binomial
from analyze_Slovenia_rail_microhabitats_v103 import (
    extract_original_xlsx,validate_original_keys,pair_source,locations_independent
)
from select_two_original_rail_20_species_for_inat_v164 import original_slovenian_EPPO_codebook_binomial_mapping
from analyze_Slovenia_source_rail_functional_taxonomic_identity_v141 import (
    source_sites_and_guilds,jaccard
)
from analyze_Slovenia_rail_spatial_matched_taxon_identity_v142 import (
    original_gps_coordinates,haversine_km
)
from analyze_Slovenia_original_LEDA_height_seed_breadth_v224 import (
    HEIGHT_SHA,SEED_SHA,all_original_matched_taxa,ztransform,create_sites,
    site_effects,section_estimates,DIMENSIONS
)
from analyze_Slovenia_gov_official_main_regional_traits_v225 import (
    official_class,categorical_contrast,KNOWN,CUTS
)

def equals_stratum(lines,units,m):
    relevant=[x for x in units if str(x['section']) in lines and x.get(m) is not None]
    by=defaultdict(list)
    for x in relevant:by[x['section']].append(float(x[m]))
    per={s:mean(q) for s,q in by.items()}
    return {'n_locations':len(relevant),'n_national_source_lines':len(per),
            'n_G_main_lines':sum(official_class(s)=='main_G' for s in per),
            'n_R_regional_lines':sum(official_class(s)=='regional_R' for s in per),
            'overall_line_equal_mean':mean(per.values()) if per else None,
            'G_vs_R':categorical_contrast(per)}

def line_pair_geo_taxonomy(original):
    sites,guild,meta=source_sites_and_guilds(original)
    attached,gpsqa=original_gps_coordinates(original,sites)
    unit=defaultdict(list)
    for a,b in combinations(attached,2):
        if a['section']==b['section']:continue
        kinds=(official_class(a['section']),official_class(b['section']))
        group=('G_G' if kinds==('main_G','main_G') else
               'R_R' if kinds==('regional_R','regional_R') else 'G_R')
        dist=haversine_km(a['lat'],a['lon'],b['lat'],b['lon'])
        bucket=next((f'{lo}-{hi}' for lo,hi in zip(CUTS[:-1],CUTS[1:]) if lo<=dist<hi),None)
        if bucket is None:continue
        sectionpair=tuple(sorted((a['section'],b['section']),key=lambda s:int(s)))
        unit[(group,bucket,sectionpair)].append(jaccard(a['known'],b['known']))
    by_group=defaultdict(list)
    for (group,bucket,linepair),values in unit.items():
        by_group[(group,bucket)].append(mean(values))
    report={}
    for (group,bucket),scores in sorted(by_group.items()):
        report.setdefault(bucket,{})[group]={
          'n_unique_original_line_pairs':len(scores),
          'mean_of_line_pair_mean_Jaccard':mean(scores),
          'min_line_pair_mean':min(scores),
          'max_line_pair_mean':max(scores)}
    supported={}
    for band,grp in report.items():
        # Require 3 distinct ORIGINAL line pairs in each of same G and R
        if all(q in grp and grp[q]['n_unique_original_line_pairs']>=3 for q in ('G_G','R_R')):
            supported[band]={'G_G_minus_R_R_source_linepair_equal_Jaccard':
                grp['G_G']['mean_of_line_pair_mean_Jaccard']-
                grp['R_R']['mean_of_line_pair_mean_Jaccard'],
                'n_G_G_source_line_pairs':grp['G_G']['n_unique_original_line_pairs'],
                'n_R_R_source_line_pairs':grp['R_R']['n_unique_original_line_pairs']}
    return report,supported

def run(original,height,seed):
    loc,zone,sp,tr,cb,qc=validate_original_keys(original)
    whole=locations_independent(pair_source(loc,zone,sp,tr))
    if len(whole)!=89 or {str(x['section_number']) for x in whole}!=set(KNOWN):
        raise ValueError('True SI source 89-site, official 15 national line ids mismatch')
    real_full=[{'section':str(x['section_number']),
                'track':str(x['track_type']),
                'annual':x['BS_minus_BT_annual_strict'],
                'therophyte':x['BS_minus_BT_therophyte_strict'],
                'cover':x['BS_minus_BT_cover']} for x in whole]
    full_single=[x for x in real_full if x['track']=='SINGLE']
    full_double=[x for x in real_full if x['track']=='DOUBLE']
    by_full={}
    for name,rows in (('all',real_full),('single_only',full_single),('double_only',full_double)):
        by_full[name]={m:equals_stratum(set(KNOWN),rows,m) for m in ('annual','therophyte','cover')}
    names,cbq=original_slovenian_EPPO_codebook_binomial_mapping(original)
    src=all_original_matched_taxa(zone,sp)
    vals,qa=ztransform({t:binomial(names.get(t)) for t in src},height,seed)
    if len(vals)!=224:raise ValueError('v22.4 expected joint species 224 source drift')
    pairs,allsides=create_sites(original,vals,min_taxa=5)
    sites=site_effects(pairs,vals)
    full_by_site={str(x['rail_location_id']):x for x in whole}
    data=[]
    for site,r in sites.items():
        true=full_by_site[site]
        for m in DIMENSIONS:r[m]=r['delta'][m]
        data.append({'section':r['section'],'track':str(true['track_type']),
                     **{m:r['delta'][m] for m in DIMENSIONS}})
    if len(data)!=82:raise ValueError('v22.4 82 source sites changed')
    main_single={s for s in KNOWN if official_class(s)=='main_G'}
    out={}
    for label,rows in (('all_joint_5species',data),
                       ('source_single_only', [x for x in data if x['track']=='SINGLE']),
                       ('source_double_only',[x for x in data if x['track']=='DOUBLE'])):
        out[label]={m:equals_stratum(set(KNOWN),rows,m) for m in DIMENSIONS}
    # Fixed local support, >=10 known taxa on BOTH matched zones:
    # strict sensitivity differs from main 5+ and must not be silently pooled.
    filt={k:[p for p in prs if p['nBS']>=10 and p['nBT']>=10] for k,prs in pairs.items()}
    filt={k:v for k,v in filt.items() if v}
    sc=site_effects(filt,vals)
    strict=[]
    for site,record in sc.items():
        strict.append({'section':record['section'],
                       'track':str(full_by_site[site]['track_type']),
                       **record['delta']})
    if len(strict)!=40 or sum(map(len,filt.values()))!=56:
        raise ValueError('v22.5 original 40 strict sites / 56 pairs changed')
    out['strict_10species_each_zone']={m:equals_stratum(set(KNOWN),strict,m) for m in DIMENSIONS}
    out['strict_single_only']={m:equals_stratum(set(KNOWN),[x for x in strict if x['track']=='SINGLE'],m)
                               for m in DIMENSIONS}
    geography,matched=line_pair_geo_taxonomy(original)
    return {'status':'v22.6_SOURCE_SI_OFFICIAL_G_R_SINGLE_TRACK_CONFOUNDER_AND_PAIRWISE_GEOGRAPHY',
        'source_original_survey_years':'2022/2023',
        'original_89_annual_and_cover_by_track':by_full,
        'original_82_LEDA_functional_breadth_by_track_and_quality':out,
        'original_full_source_train_route_track_cross_count':{
           'G_main_SINGLE_original_sites':sum(x['track']=='SINGLE' and official_class(x['section'])=='main_G' for x in real_full),
           'G_main_DOUBLE_original_sites':sum(x['track']=='DOUBLE' and official_class(x['section'])=='main_G' for x in real_full),
           'R_regional_SINGLE_original_sites':sum(x['track']=='SINGLE' and official_class(x['section'])=='regional_R' for x in real_full),
           'R_regional_DOUBLE_original_sites':sum(x['track']=='DOUBLE' and official_class(x['section'])=='regional_R' for x in real_full)},
        'original_site_pair_taxon_Jaccard_by_distance_and_national_G_R':geography,
        'geographic_distance_bins_with_at_least_3_distinct_original_line_pairs_for_GG_and_RR':matched,
        'strict_cautions':['comparison within single track has very few G lines and is spatially unbalanced',
           'all pairwise source Jaccard estimates are dependent and only one nation/study',
           'national G/R line category is not observed passenger/freight daily train counts',
           'geographic binning only cannot eliminate regional species pool or environment confounding',
           'smaller strict trait-richness source cohort selects different sites and sides post hoc',
           'data do not prove any selection, transport or allelic convergence']}
def main():
    p=argparse.ArgumentParser()
    for n in ('source_xlsx','source_audit','height','seed','output'):
        p.add_argument('--'+n.replace('_','-'),required=True)
    a=p.parse_args()
    source=extract_original_xlsx(a.source_xlsx,a.source_audit)
    if hashlib.sha256(Path(a.height).read_bytes()).hexdigest()!=HEIGHT_SHA:raise ValueError('LEDA original height source SHA drift')
    if hashlib.sha256(Path(a.seed).read_bytes()).hexdigest()!=SEED_SHA:raise ValueError('LEDA original seed source SHA drift')
    height,hqa=read_leda(a.height,'height_m')
    seed,sqa=read_leda(a.seed,'seed_mass_mg')
    report=run(source,height,seed)
    path=Path(a.output);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf8')
    print('V226_SI_SINGLE_TRACK_RESULTS='+json.dumps(report,ensure_ascii=False,separators=(',',':')),flush=True)
if __name__=='__main__':main()
