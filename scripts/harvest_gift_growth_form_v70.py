#!/usr/bin/env python3
"""Fetch genuine public GIFT 1.2.1 growth form for existing six-country GBIF GRIIS taxa.

Frozen plant taxon opportunity first. Exact accepted species names only. This is
NOT a measured mowing/wind exposure, a railway adaptation, or new dispersal.
"""
from __future__ import annotations
import argparse
from collections import Counter
from datetime import datetime, timezone
import json, math, re, time
from pathlib import Path
import requests

BLOCKS_SIX=('us_midwest','fr_burgundy','ca_ontario','gb_yorkshire','de_saxony','au_southeast')
PERIODS=('2010_2017','2018_2025')
GIFT='https://gift.uni-goettingen.de/api/extended/index.php'
GBIF='https://api.gbif.org/v1/species'
TRAIT_ID='1.2.1'  # GIFT growth form, official R vignette
CHUNK=10000
MAX_TRAIT_PAGES=45
MAX_TAXA=12000
MAX_GIFT_SPECIES_PAGES=6
MIN_AGREEMENT=.66


def checked_rows(data, columns):
    """Validate actual public API JSON without inventing missing values."""
    if not isinstance(data,list):
        raise ValueError('GIFT API expected JSON array, got missing/changed response contract')
    for row in data:
        if not isinstance(row,dict) or not set(columns).issubset(row):
            raise ValueError('GIFT API row schema drift or missing mandatory fields')
    return data


def species_binomial(value):
    """Reject fuzzy and infraspecific string matches; exact genus+epithet only."""
    if not isinstance(value,str): return None
    parts=value.strip().split()
    if len(parts)!=2 or not re.fullmatch(r'[A-Z][A-Za-z\-]+',parts[0]) or not re.fullmatch(r'[a-z][a-z\-]+',parts[1]):
        return None
    return ' '.join(parts).casefold()


def growth_group(value):
    """Conservative, explicit morphology group; climbers and mixed forms unknown."""
    if not isinstance(value,str):return None
    x=re.sub(r'[_\-]+',' ',value.strip().lower())
    if x in {'herb','herbs','herbaceous','graminoid','forb','graminoids','forbs'}:
        return 'herbaceous'
    if x in {'tree','trees','shrub','shrubs','subshrub','subshrubs','woody shrub','woody tree'}:
        return 'woody'
    return None


def make_source_keyset(flora):
    if flora.get('complete') is not True or flora.get('source_errors') or not set(BLOCKS_SIX).issubset(flora.get('regions',{})):
        raise ValueError('Original six-country GRIIS source not complete')
    per={}
    keys=set()
    for block in BLOCKS_SIX:
        region=flora['regions'][block]
        source=set()
        for cid,cell in region['cells'].items():
            if not cell.get('source_pre_outcome_both_period_recording_eligible'):
                continue
            for period in PERIODS:
                p=cell.get('periods',{}).get(period,{})
                if not p.get('full_facet_source_pagination_completed'):
                    continue
                counts=p.get('observed_alien_specieskey_counts',{})
                if not isinstance(counts,dict):raise ValueError('Invalid GBIF species facet count')
                for k,v in counts.items():
                    if not str(k).isdigit() or not isinstance(v,int) or isinstance(v,bool) or v<=0:
                        raise ValueError('GBIF taxon key or count source invalid')
                    source.add(str(k))
        per[block]=sorted(source,key=int)
        keys.update(source)
    if not keys or len(keys)>MAX_TAXA:
        raise ValueError('Unexpected six-country GBIF taxon count')
    return sorted(keys,key=int),per


def source_get(client,url,params=None):
    r=client.get(url,params=params,timeout=(15,120))
    r.raise_for_status()
    try: return r.json()
    except ValueError as exc:raise ValueError('Non-JSON external source response') from exc


def GIFT_sources(client):
    meta=checked_rows(source_get(client,GIFT,{'query':'traits_meta'}),('Lvl3','Trait2','count'))
    rows=[r for r in meta if str(r['Lvl3'])==TRAIT_ID]
    if len(rows)!=1 or 'growth' not in str(rows[0]['Trait2']).lower() or 'form' not in str(rows[0]['Trait2']).lower():
        raise ValueError('Official GIFT 1.2.1 growth-form identity changed; stop instead of mixing traits')
    count=int(rows[0]['count'])
    if count<1 or count>MAX_TRAIT_PAGES*CHUNK:
        raise ValueError('GIFT meta trait count out of bounded source range')
    ntrait_pages=math.ceil(count/CHUNK)
    # The GIFT package itself retrieves species in six 100k pages; same API/offset.
    species={}
    for i in range(MAX_GIFT_SPECIES_PAGES):
        records=checked_rows(source_get(client,GIFT,{'query':'species','startat':i*100000}),('work_ID','work_species'))
        for row in records:
            wid=str(row['work_ID'])
            nm=species_binomial(row['work_species'])
            if nm:
                if wid in species and species[wid]!=nm:
                    raise ValueError('Conflicting GIFT species ID')
                species[wid]=nm
        if not records:break
    if not species:raise ValueError('Empty GIFT actual species table')
    work={}
    classification=Counter()
    for i in range(ntrait_pages):
        recs=checked_rows(source_get(client,GIFT,{'query':'traits','traitid':TRAIT_ID,
                       'biasref':0,'biasderiv':0,'startat':i*CHUNK}),('work_ID','trait_value','agreement'))
        for r in recs:
            wid=str(r['work_ID'])
            a=r.get('agreement')
            try: agreement=float(a)
            except (ValueError,TypeError):
                classification['missing_agreement']+=1;continue
            if not math.isfinite(agreement) or agreement<MIN_AGREEMENT or agreement>1:
                classification['below_threshold_or_invalid']+=1;continue
            group=growth_group(r['trait_value'])
            if group is None:
                classification['unclassified_value']+=1;continue
            if wid in work and work[wid]!=group:
                work[wid]='ambiguous'
            elif wid not in work:
                work[wid]=group
        time.sleep(.02)
    byname={}
    for wid,cat in work.items():
        nm=species.get(wid)
        if nm is None or cat=='ambiguous':continue
        if nm in byname and byname[nm]!=cat:
            byname[nm]='ambiguous'
        elif nm not in byname:
            byname[nm]=cat
    return byname,{'GIFT_meta_ID':TRAIT_ID,'GIFT_meta_label':rows[0]['Trait2'],
        'GIFT_source_meta_count':count,'GIFT_source_pages':ntrait_pages,
        'GIFT_species_workIDs_names_parsed':len(species),
        'GIFT_unclassified_qa':dict(classification),'GIFT_biasref':False,
        'GIFT_biasderiv':False,'GIFT_categorical_min_agreement':MIN_AGREEMENT,
        'GIFT_public_API':GIFT}


def GBIF_exact_names(client,keys):
    names={}; excluded=Counter(); seen={}
    for i,k in enumerate(keys):
        data=source_get(client,f'{GBIF}/{k}')
        if not isinstance(data,dict):raise ValueError('GBIF species API changed')
        accepted=str(data.get('acceptedKey',k))
        if (str(data.get('rank','')).upper()!='SPECIES' or
            str(data.get('kingdom','')).lower()!='plantae' or
            str(data.get('taxonomicStatus',data.get('status',''))).upper() not in {'ACCEPTED','PROVISIONALLY_ACCEPTED','PROVISIONALLY ACCEPTED'} or
            accepted!=k):
            excluded['not_strict_current_accepted_species']+=1;continue
        nm=species_binomial(data.get('canonicalName'))
        if nm is None:
            excluded['invalid_or_nonbinomial_canonical_name']+=1;continue
        names[k]=nm
        seen.setdefault(nm,[]).append(k)
        if i%80==0:time.sleep(.08)
    # Species concepts mapping many source keys to one name must not silently
    # become 2 independently labeled source taxa.
    ambiguous={k for ks in seen.values() if len(ks)>1 for k in ks}
    for k in ambiguous:names.pop(k)
    excluded['duplicate_gbif_accepted_names_removed']=len(ambiguous)
    return names,dict(excluded)


def analyze(flora,client):
    keys,countrykeys=make_source_keyset(flora)
    byname,source=GIFT_sources(client)
    names,excluded=GBIF_exact_names(client,keys)
    mapped={k:byname[nm] for k,nm in names.items() if byname.get(nm) in {'herbaceous','woody'}}
    status={}
    for country,cc in countrykeys.items():
        ccset=set(cc)
        match=ccset & set(mapped)
        status[country]={
            'n_distinct_source_GRIIS_specieskeys':len(ccset),
            'n_exact_GBIF_current_accepted_name':len(ccset&set(names)),
            'n_GIFT_growth_form_mapped':len(match),
            'n_herbaceous':sum(mapped[k]=='herbaceous' for k in match),
            'n_woody':sum(mapped[k]=='woody' for k in match),
            'source_GIFT_coverage_species_fraction':round(len(match)/len(ccset),6) if ccset else None,
        }
    return {'complete':True,'utc':datetime.now(timezone.utc).isoformat(),
        'source':source,'n_union_GBIF_original_source_specieskeys':len(keys),
        'n_exact_GBIF_current_accepted_name':len(names),
        'n_growth_form_accepted_taxon_matches':len(mapped),
        'gbif_taxonomy_and_unmatched_qa':excluded,'by_block':status,
        'growth_form_by_source_gbif_species_key':mapped,
        'limitations':['GIFT species growth form is not measured mowing or slipstream stress tolerance.',
          'GBIF index species keys are source records, not a species census or observed absence.',
          'Exact accepted binomial match only; taxonomy drift, synonyms, ambiguous growth forms and missing traits excluded.',
          'No actual train service, wind, cutting dates, genetic adaptation or seed transfer measured.'],
        'NO_causal_rail_mowing_wind_adaptation_claim':True}


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--flora',required=True)
    p.add_argument('--output',required=True)
    a=p.parse_args()
    flora=json.loads(Path(a.flora).read_text(encoding='utf-8'))
    client=requests.Session();client.headers.update({'User-Agent':'tetsu-GIFT-growth-form-GRIIS-macro/7.0 (research)'})
    result=analyze(flora,client)
    dst=Path(a.output);dst.parent.mkdir(parents=True,exist_ok=True)
    dst.write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps({'complete':True,'source_taxa':result['n_union_GBIF_original_source_specieskeys'],
        'trait_mapped':result['n_growth_form_accepted_taxon_matches'],'country_coverage':result['by_block'],
        'not_evolution_or_seed_transport':True}),flush=True)

if __name__=='__main__':main()