#!/usr/bin/env python3
"""v24.6: genuinely independent fourth railway-flora original study (Slovenia 2026).

NOT pooled with station-wide Czech 39, Italy 28, HR+BA 15. Slovenian
2022–23 field data are 89 100-m trackside original author locations,
173 same-side BS-vs-BT ballast microhabitat pairs nested in 15 national
route sections, not 89 station-wide repeated surveys.

Use the SAME source-authenticated Europe CZ+IT+HR/BA 622-species
z(log LEDA canopy height[m]), z(log seed mass[mg]) axes from v24.5
original audit (GitHub Actions artifact, SHA and source count verified).
Source Slovenia independently from Šabić et al 2026 genuine publisher
401279 byte XLSX / source SHA, original 89 sites in 15 route sections;
each site is union of original BS+BT species across matched sides, not
wild-sampled station boundary or estimated "non-rail" species.

Within 15 source lines, preserve EACH field site's exact dual-LEDA-
species richness and each species' frequency within its SOURCE
route section via local checkerboard swaps. Sections with just one
site are kept unchanged rather than sampled with replacement. The
15 geographical route sections are NOT randomized natural regions.

Compute CWM Euclidean distances within- and between-section, site
2D population variance trace and source-matched Latin-species Jaccard.
2 chains * 100 technical draws; do not treat 3916 site pairs as
independent field replicates or source internal null p-values.

This is a fourth *observational study* across a different spatial
grain, not a fourth independent 39-station-type dataset. A positive
deviation of CWM beta from the conditional railway pool is NOT
evidence rail-specific homogenization or railway plant adaptation.
"""
from __future__ import annotations
import argparse,hashlib,json,math,random
from collections import Counter,defaultdict
from pathlib import Path
from statistics import mean,median
from analyze_Slovenia_rail_microhabitats_v103 import (
    extract_original_xlsx,validate_original_keys
)
from analyze_Slovenia_source_rail_functional_taxonomic_identity_v141 import source_sites_and_guilds
from select_two_original_rail_20_species_for_inat_v164 import original_slovenian_EPPO_codebook_binomial_mapping
from audit_SI_CZ_original_continuous_trait_coverage_v220 import read_leda,binomial
from analyze_Slovenia_original_LEDA_height_seed_breadth_v224 import (
    HEIGHT_SHA,SEED_SHA,quart
)
from validate_Europe_3studies_82stations_original_LEDA_2D_trait_beta_v245 import (
    get_metrics,regional_margin,METRICS
)

N_SITES=89
N_SECTIONS=15
KNOWN_SI_BOTH_TRAITS=224
SEED=20261117+246
CHAINS=2
DRAWS=100
BURNIN=18000
GAP=1700
SOURCE_V245_STUDY_COUNT=3
SOURCE_V245_STATION_N=82
SOURCE_V245_GLOBAL_TRAIT_N=622

def source_european_trait_axis(qa,h,seed):
    if qa.get('status')!='SOURCE_VERIFIED_CZ_IT_HR_BA_82_INDEPENDENT_RAIL_STATIONS_2D_LEDA_INTERNAL_TRAIT_NULL':
        raise ValueError('Source European original genuine v24.5 validation status drift')
    s=qa.get('source_all_three_study_taxa_shared_trait_axis_QA') or {}
    signs=qa.get('original_three_study_descriptive_sign_coherence') or {}
    if s.get('source_original_LEDA_both_canopy_height_seed_mass_species_n')!=SOURCE_V245_GLOBAL_TRAIT_N:
        raise ValueError('Real Europe fixed 622-source-species global continuous trait pool drift')
    if signs.get('number_original_independent_publisher_studies')!=SOURCE_V245_STUDY_COUNT or (
        signs.get('original_station_number_across_three_studies')!=SOURCE_V245_STATION_N):
        raise ValueError('European 82 source original stations no longer independently verified')
    mu=s['species_equal_global_log_mean'];sd=s['species_equal_global_log_SD']
    if len(mu)!=2 or len(sd)!=2 or min(sd)<=0:
        raise ValueError('Original v24.5 global scaling malformed')
    # This 3-country original source pooled scaling is FROZEN BEFORE
    # computing Slovenian response values: nonleaky fourth-study test.
    overlap=set(h)&set(seed)
    return {taxon:((math.log(h[taxon])-mu[0])/sd[0],
                   (math.log(seed[taxon])-mu[1])/sd[1])
            for taxon in overlap},{
       'source_original_reference_study_count':3,
       'original_reference_rail_station_count':82,
       'original_reference_trait_unique_species_count':622,
       'no_slovenian_species_used_to_refit_reference_trait_axis':True,
       'frozen_v245_reference_mean':mu,'frozen_v245_reference_sd':sd}

def genuine_slovenian_89_microhabitat_unions(original,all_trait_values):
    sites,guild,source=source_sites_and_guilds(original)
    names,cbq=original_slovenian_EPPO_codebook_binomial_mapping(original)
    raw_local={k:len(s['all_original_taxa']) for k,s in enumerate(sites)}
    used=[];all_joint=set()
    for site in sites:
        codes={binomial(names.get(code)) for code in site['all_original_taxa']}
        codes.discard(None)
        valid=codes&set(all_trait_values)
        # We need sufficient real original trait coverage per original
        # independent track transect before calculating variance.
        if len(valid)<5:
            continue
        all_joint|=valid
        used.append({'station':site['site'],
            'group':site['section'],
            'year':site['year'],
            'track_type':site['track_type'],
            'taxa':valid,
            'source_original_taxon_n':len(codes),
            'source_original_matched_sides':site['n_original_matched_track_sides']})
    if len(used)<80 or len({x['group'] for x in used})<13:
        raise ValueError('Genuine original Slovenia matched 89 locations 15 route sections lacks source 2D coverage')
    original_all={binomial(names.get(code)) for x in sites
                  for code in x['all_original_taxa']}
    original_all.discard(None)
    if len(original_all)<260:
        raise ValueError('Source original Slovenian full 89-site flora has too few author taxa')
    return used,{
       'source_original_100m_matched_BS_BT_locations':len(sites),
       'source_original_route_sections':len({site['section'] for site in sites}),
       'real_source_100m_locations_ge5_joint_LEDA_species':len(used),
       'source_real_route_sections_ge5_joint_LEDA_species':len({x['group'] for x in used}),
       'matched_source_unique_taxa_in_genuine_original_plant_sites':len(original_all),
       'source_unique_taxa_with_both_external_LEDA_median_traits':len(all_joint),
       'source_captured_original_author_survey_2022_2023':True,
       'one_source_location_is_not_a_whole_railway_station':True,
       'original_mean_matched_dual_species_per_field_location':
           mean(len(x['taxa']) for x in used),
       'original_min_median_max_dual_species_per_field_location':[
           min(len(x['taxa']) for x in used),
           median(len(x['taxa']) for x in used),
           max(len(x['taxa']) for x in used)],
       'source_publisher_route_community_sample_size':dict(Counter(x['group'] for x in used)),
       'source_provenance_QA':source,
       'original_taxonomic_EPPO_to_binomial_QA':cbq}

def source_slovenian_stratified_null(sites,traits,draws=DRAWS,
                                   burnin=BURNIN,gap=GAP,seed=SEED):
    rg=defaultdict(list)
    for i,s in enumerate(sites):rg[s['group']].append(i)
    active=[z for z in rg.values() if len(z)>=2]
    if len(active)<7:raise ValueError('Insufficient original source railway route sections with >=2 independent 100m trackside sites')
    nzero=len(rg)-len(active)
    origmargin=regional_margin(sites)
    richness=[len(x['taxa']) for x in sites]
    rnd=random.Random(seed)
    state=[set(x['taxa']) for x in sites]
    accepted=0
    def step(proposals):
        nonlocal accepted
        for _ in range(proposals):
            pool=rnd.choice(active)
            i,j=rnd.sample(pool,2)
            aa=state[i]-state[j];bb=state[j]-state[i]
            if not aa or not bb:continue
            a=rnd.choice(sorted(aa));b=rnd.choice(sorted(bb))
            state[i].remove(a);state[i].add(b)
            state[j].remove(b);state[j].add(a)
            accepted+=1
    step(burnin)
    out=[]
    for trial in range(draws):
        step(gap)
        changed=[{**s,'taxa':state[i]} for i,s in enumerate(sites)]
        if [len(x) for x in state]!=richness or regional_margin(changed)!=origmargin:
            raise ValueError('Source regional and location plant richness marginal constraints changed')
        out.append(get_metrics(changed,traits))
    if accepted<1000:raise ValueError('Original Slovenian source route stratified checkerboard swaps failed to mix')
    return out,{
      'source_n_route_sections_eligible_for_within_line_swaps':len(active),
      'source_n_route_sections_with_single_site_unchanged':nzero,
      'total_accepted_source_swaps':accepted,
      'all_source_original_site_richness_and_route_frequency_held_fixed':True,
      'n_technical_draws_not_independent_field':len(out),
      'n_technical_swap_proposals':burnin+draws*gap}

def run(original,v245,height,seed,draws=DRAWS):
    traits,qaxis=source_european_trait_axis(v245,height,seed)
    sites,qa=genuine_slovenian_89_microhabitat_unions(original,traits)
    obs=get_metrics(sites,traits)
    chains=[];nulls=[]
    for n in range(CHAINS):
        q,mix=source_slovenian_stratified_null(
            sites,traits,seed=SEED+n*17051,draws=draws)
        chains.append(mix);nulls.extend(q)
    comparison={}
    for key in METRICS:
        val=obs[key]
        vec=[x[key] for x in nulls]
        avg=mean(vec)
        comparison[key]={
          'observed_genuine_100m_author_rail_sites':val,
          'within_original_route_section_site_richness_and_taxon_frequency_fixed_null_mean':avg,
          'observed_minus_source_internal_null':val-avg,
          'technical_null_q025_NOT_field_CI':quart(vec,.025),
          'technical_null_q975_NOT_field_CI':quart(vec,.975),
          'technical_lower_tail_rank_plusone_NOT_field_p':
              (1+sum(x<=val for x in vec))/(1+len(vec)),
          'no_comparison_to_natural_nonrail_regional_pool':True}
    return {'version':'v24.6',
       'status':'AUTHENTIC_4TH_ORIGINAL_STUDY_SLOVENIA_89_100M_FIELD_SITES_SAME_EUROPE_TRAIT_AXIS',
       '2026_slovenian_field_source_QA':qa,
       'same_2026_Europe_3studies_82stations_source_trait_scaling_QA':qaxis,
       'true_source_number_source_100m_sites':len(sites),
       'n_original_route_sections':len({x['group'] for x in sites}),
       'source_conditional_null_QA':chains,
       'Slovenian_4th_study_same_two_tratis_under_different_field_sampling_grain':comparison,
       'not_an_equal_station_survey_study':True,
       'not_a_natural_reference_control':True,
       'interpretation_stops':[
         'Slovenian 2022–23 true original 89 field stations are author 100m ballast shoulder+toe trackside location unions, NOT a station-wide 2020 original railway station flora',
         'Same global zlog trait scale from Europe 3 original station studies was fixed before Slovenia original response looked at',
         'Source original section nestedness and 15 sections not 15 independent studies, source site pairs not independent n=3916',
         'Source local seed/height species-median traits not in situ phenotype, genotype or adaptation',
         'Conditional null fixes the whole original railway site source taxon occurrences, not independent natural regional flora',
         'No railway traffic counts matched to 2022/23 direct seed movement and freight origin destination',
         'Differences in BS vs BT source sampled area, occurrence detection and spatial survey protocol not controlled',
         'Never pool 89 track locations with true 82 station-wide floras to claim 171 equal-area station replicates']}
def main():
    p=argparse.ArgumentParser()
    for f in ('source_xlsx','source_audit','v245_audit','height','seed','output'):
        p.add_argument('--'+f.replace('_','-'),required=True)
    p.add_argument('--draws',type=int,default=DRAWS)
    a=p.parse_args()
    for file,dig in ((a.height,HEIGHT_SHA),(a.seed,SEED_SHA)):
        if hashlib.sha256(Path(file).read_bytes()).hexdigest()!=dig:
            raise ValueError('Original Oldenburg LEDA public SHA differs')
    original=extract_original_xlsx(a.source_xlsx,a.source_audit)
    v245=json.loads(Path(a.v245_audit).read_text(encoding='utf8'))
    h,hq=read_leda(a.height,'height_m');s,sq=read_leda(a.seed,'seed_mass_mg')
    result=run(original,v245,h,s,a.draws)
    path=Path(a.output);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf8')
    summary={k:result[k] for k in ('status','true_source_number_source_100m_sites',
       'n_original_route_sections','2026_slovenian_field_source_QA','source_conditional_null_QA',
       'Slovenian_4th_study_same_two_tratis_under_different_field_sampling_grain')}
    print('V246_SLOVENIA_FOURTH_TRUE_SOURCE_GRAIN='+json.dumps(summary,ensure_ascii=False,separators=(',',':')),flush=True)
if __name__=='__main__':main()
