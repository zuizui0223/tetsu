#!/usr/bin/env python3
"""v22.0: Original-source SI/CZ continuous-trait coverage gate (not an effect test).

SI: source-verified 2026 92-site XLSX, paired BS/BT ballast microhabitats.
CZ: source-verified 2024 S1 + S2 erratum, whole-station contemporary flora.
LEDA: canopy height [m] and seed mass [mg] directly from University Oldenburg.
Only aggregate coverage is released; no publisher/LEDA raw data redistributed.
"""
from __future__ import annotations
import argparse,csv,hashlib,io,json,math,re
from collections import Counter,defaultdict
from pathlib import Path
from statistics import median
from analyze_Slovenia_rail_microhabitats_v103 import extract_original_xlsx,validate_original_keys
from select_two_original_rail_20_species_for_inat_v164 import original_slovenian_EPPO_codebook_binomial_mapping
from analyze_Czech_2024_39station_functional_identity_v143 import read_original_czech_publisher
from analyze_Slovenia_original_guild_vs_species_microhabitat_15section_predict_v172 import author_source_rows

LEDA={
 'height_m':('canopy_height.txt','[m]'),
 'seed_mass_mg':('seed_mass.txt','[mg]'),
}
BASE='https://uol.de/f/5/inst/biologie/ag/landeco/download/LEDA/Data_files/'

def binomial(raw):
    if not isinstance(raw,str):return None
    v=raw.strip()
    m=re.match(r'^([A-Z][A-Za-z-]+)\s+([a-z][a-z-]+)(?=\s|$)',v)
    if not m or m.group(2) in ('sp','spp','agg','cf','aff','indet'):return None
    return m.group(1)+' '+m.group(2)

def read_leda(path,which):
    raw=Path(path).read_bytes()
    if len(raw)<50000:raise ValueError('LEDA source suspiciously small')
    text=None;enc=None
    for e in ('utf-8-sig','utf-16','cp1252','latin-1'):
        try:
            candidate=raw.decode(e)
            if '\x00' not in candidate:
                text=candidate;enc=e;break
        except UnicodeDecodeError:continue
    if text is None:raise ValueError('Could not decode original LEDA data')
    lines=text.splitlines()
    header=next((i for i,s in enumerate(lines[:40]) if 'SBS.name' in s and ';' in s),None)
    if header is None or header<2:raise ValueError('LEDA SBS.name header not detected; genuine first lines='+repr(lines[:10])[:2500])
    reader=csv.DictReader(io.StringIO('\n'.join(lines[header:])),delimiter=';')
    cols=reader.fieldnames or []
    namecol=next((s for s in cols if 'sbs.name' in s.lower()),None)
    valuecol=next((s for s in cols if 'single value' in s.lower() and LEDA[which][1] in s.lower()),None)
    if namecol is None or valuecol is None:
        raise ValueError('LEDA source trait column missing: '+str(cols[:25]))
    numeric=defaultdict(list);qc=Counter()
    for r in reader:
        qc['rows']+=1
        name=binomial((r.get(namecol) or '').strip(' "'))
        if name is None:qc['invalid_binomial']+=1;continue
        rawvalue=(r.get(valuecol) or '').strip(' "')
        try:v=float(rawvalue)
        except ValueError:qc['missing_or_non_numeric']+=1;continue
        if not math.isfinite(v) or v<=0:
            qc['nonpositive_or_nonfinite']+=1;continue
        numeric[name].append(v)
    if len(numeric)<1000:raise ValueError('Original LEDA table yields unexpectedly few species')
    return {n:median(a) for n,a in numeric.items()},{
      'original_public_source':BASE+LEDA[which][0],
      'sha256_actual_download':hashlib.sha256(raw).hexdigest(),
      'source_file_bytes':len(raw),'file_encoding':enc,'header_at_zero_based_line':header,
      'original_name_column':namecol,'original_measure_column':valuecol,
      'unique_positive_trait_species':len(numeric),
      'n_original_positive_trait_measurements':sum(map(len,numeric.values())),
      'median_across_repeated_measurements':True,'QC':dict(qc)}

def summarize(taxa,units,trait_db):
    names={t:binomial(name) for t,name in taxa.items()}
    ident={t:n for t,n in names.items() if n is not None}
    n_events=sum(len(u['codes']) for u in units)
    if not taxa or n_events<=0:raise ValueError('No real original plant presences')
    out={'taxon_codes_n':len(taxa),
         'source_taxa_with_unambiguous_name_n':len(ident),
         'unique_species_binomials_n':len(set(ident.values())),
         'n_real_source_sites_or_zones':len(units),
         'n_real_site_taxon_incidence':n_events,'traits':{}}
    matched={}
    for tr,values in trait_db.items():
        codes={t for t,n in ident.items() if n in values}
        matched[tr]=codes
        per=[len(u['codes']&codes) for u in units]
        out['traits'][tr]={
            'n_source_taxa_joined':len(codes),
            'n_source_taxa_total':len(taxa),
            'fraction_source_taxa_joined':len(codes)/len(taxa),
            'incidence_joined':sum(per),
            'incidence_total':n_events,
            'fraction_incidence_joined':sum(per)/n_events,
            'n_sites_ge5_for_one_trait_variance':sum(x>=5 for x in per),
            'n_sites_ge10':sum(x>=10 for x in per),
            'median_site_trait_species_count':median(per)}
    joint=set.intersection(*matched.values())
    per=[len(u['codes']&joint) for u in units]
    out['joint_height_seed_mass']={
        'n_source_taxa_joined':len(joint),
        'fraction_source_taxa_joined':len(joint)/len(taxa),
        'incidence_joined':sum(per),
        'fraction_incidence_joined':sum(per)/n_events,
        'n_sites_ge5_joint':sum(x>=5 for x in per),
        'n_sites_ge10_joint':sum(x>=10 for x in per)}
    return out,{**matched,'both':joint}

def audit_si(path,audit,traits):
    original=extract_original_xlsx(path,audit)
    # Replicate v17.2 original validated 89 sites, 15 sections, 261 guild taxa, 2729 events.
    known,auth=author_source_rows(original)
    if len(known)!=2729 or len({x['author_original_section'] for x in known})!=15:
        raise ValueError('Slovenia genuine 89-site cohort changed')
    loc,zones,sp,alltr,book,qc=validate_original_keys(original)
    names,codebook=original_slovenian_EPPO_codebook_binomial_mapping(original)
    units=[];pairs=[]
    for (lid,zid),record in zones.items():
        if record.get('area')!='BS' or zid not in ('LBS','RBS'):continue
        other=(lid,zid[0]+'BT')
        if other not in zones or zones[other].get('area')!='BT':continue
        A=set(sp.get((lid,zid),set()))
        B=set(sp.get(other,set()))
        if not A or not B:continue
        units.extend([{'codes':A},{'codes':B}])
        pairs.append({'site':lid,'section':str(loc[lid]['track_section_number']),
                      'BS':A,'BT':B})
    sites={p['site'] for p in pairs};sections={p['section'] for p in pairs}
    if len(sites)!=89 or len(sections)!=15:
        raise ValueError('Slovenia same-side BS/BT genuine matched-source support changed')
    source_taxa=set().union(*(u['codes'] for u in units))
    report,matches=summarize({t:names.get(t) for t in source_taxa},units,traits)
    report['original_census']={'author_locations_92_source':len(loc),
       'QC_real_matched_locations':len(sites),'line_sections':len(sections),
       'matched_same_side_BS_BT_pairs':len(pairs),
       'original_two_guild_known_taxa':auth['n_actual_taxa_with_both_author_lifespan_and_diaspore_source_codes'],
       'original_source_site_taxon_events_union_sides':len(known),
       'original_source_EPPO_species_codebook_QA':codebook}
    report['paired_diversity_support']={}
    for tr,codes in matches.items():
        valid=[p for p in pairs if len(p['BS']&codes)>=5 and len(p['BT']&codes)>=5]
        robust=[p for p in pairs if len(p['BS']&codes)>=10 and len(p['BT']&codes)>=10]
        report['paired_diversity_support'][tr]={
         'pairs_ge5_each_zone':len(valid),
         'distinct_sites_ge5_each_zone':len({p['site'] for p in valid}),
         'distinct_sections_ge5_each_zone':len({p['section'] for p in valid}),
         'pairs_ge10_each_zone':len(robust),
         'median_BS_trait_species':median(len(p['BS']&codes) for p in pairs),
         'median_BT_trait_species':median(len(p['BT']&codes) for p in pairs),
         'predeclared_support_gate_ge60_pairs_and_ge10_sections':
              len(valid)>=60 and len({p['section'] for p in valid})>=10}
    report['physical_zone_contrast_real']='BS vs BT at same original location and track side'
    return report

def audit_cz(s1,s2,traits):
    floras,meta,auth=read_original_czech_publisher(s1,s2)
    units=[{'codes':s} for s in floras.values()]
    taxa=set().union(*floras.values())
    if len(floras)!=39 or len(taxa)!=761 or sum(map(len,floras.values()))!=6941:
        raise ValueError('Original Czech station source incidence drifted')
    report,_=summarize({x:x for x in taxa},units,traits)
    report['original_census']={'stations':39,
      'unique_taxa':761,'original_station_taxon_presences':6941,
      'original_railway_stretches':len({x['stretch'] for x in meta.values()}),
      'publisher_original_source_QA':auth}
    report['physically_matched_BS_BT_comparison_possible']=False
    return report

def main():
    p=argparse.ArgumentParser()
    for arg in ('si_original','si_audit','cz_s1','cz_s2','height','seed','output'):
        p.add_argument('--'+arg.replace('_','-'),required=True)
    z=p.parse_args()
    traits={};provenance={}
    for label,path in (('height_m',z.height),('seed_mass_mg',z.seed)):
        traits[label],provenance[label]=read_leda(path,label)
    si=audit_si(z.si_original,z.si_audit,traits)
    cz=audit_cz(z.cz_s1,z.cz_s2,traits)
    result={'version':'v22.0','status':'SOURCE_VERIFIED_TRAIT_COVERAGE_ONLY',
      'data_provider':provenance,'Slovenia':si,'Czechia':cz,
      'source_inference_limits':[
        'LEDA species median is not an observed railway individual phenotype',
        'SI 89 independent locations but side pairs and taxa repeat within location',
        'CZ station-wide incidence has three habitat types pooled; no BS/BT response',
        'No independent local natural species pools; no regional-filtering effect yet',
        'No train counts, herbicide pressure, mowing exposure or genotype measured',
        'No fuzzy matches, genus-to-species substitution or imputed missing traits']}
    dst=Path(z.output);dst.parent.mkdir(parents=True,exist_ok=True)
    dst.write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf8')
    print(json.dumps({'provenance':provenance,
                      'SI':{'coverage':si['traits'],'joint':si['joint_height_seed_mass'],
                            'paired':si['paired_diversity_support'],'n_sites':si['original_census']['QC_real_matched_locations']},
                      'CZ':{'coverage':cz['traits'],'joint':cz['joint_height_seed_mass']}},
                     ensure_ascii=False,indent=2))
if __name__=='__main__':main()
