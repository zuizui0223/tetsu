#!/usr/bin/env python3
"""v369q independent original source GZIP 100k record QUICK quality grade audit.

Never download full 13.7GB or publish observation IDs, exact latlon,
individual taxa or observer IDs. Verify actual original 2026-10-27 iNat
photo-linked metadata schema and quality grade enum, early record year
mixture and GPS accuracies before a full-world original stream run.
"""
from __future__ import annotations
import argparse,csv,gzip,hashlib,io,json
from collections import Counter
from pathlib import Path

from scan_fullworld_2025_iNat_plant_macro_grid_v369 import (
    S3,REGION,OBS_ETAG,OBS_BYTES,source_client
)
FIRST_N=100000
SOURCE_CSV_KEY='observations.csv.gz'
def run():
    client=source_client()
    meta=client.head_object(Bucket=S3,Key=SOURCE_CSV_KEY)
    if meta['ETag'].strip('"')!=OBS_ETAG or int(meta['ContentLength'])!=OBS_BYTES:
        raise ValueError('Original public iNat monthly source changed or not pinned')
    result=Counter()
    response=client.get_object(Bucket=S3,Key=SOURCE_CSV_KEY,
         IfMatch='"'+OBS_ETAG+'"')
    with io.TextIOWrapper(gzip.GzipFile(fileobj=response['Body']),
          encoding='utf-8',newline='',errors='strict') as f:
        rows=csv.DictReader(f,delimiter='\t',quoting=csv.QUOTE_NONE)
        names=rows.fieldnames
        if names is None or not {'observed_on','quality_grade','taxon_id','positional_accuracy','latitude','longitude'}<=set(names):
            raise ValueError('Original public iNat observed-source header changed')
        for ix,r in enumerate(rows):
            if ix>=FIRST_N:break
            label=str(r['quality_grade']).strip().lower().replace(' ','_').replace('-','_')
            result['source_quality_'+label]+=1
            result['original_first_year_'+str(r['observed_on'])[:4]]+=1
            result['source_first_records_accuracy_present']+=bool(r['positional_accuracy'])
            result['source_first_records_georeference_present']+=bool(r['latitude'] and r['longitude'])
            result['source_first_records_plants_2025_quality_research_grade']+=str(r['observed_on']).startswith('2025-') and label in ('research','research_grade')
    return {'version':'v369q-pinned-first100k-2026-10-27-observation-real-schema-grade-source-audit',
     'source_first_n_original_observation_rows':FIRST_N,
     'source_original_gzip_ETag':OBS_ETAG,
     'source_original_observation_records_quality_distribution':{
       k:v for k,v in result.items() if k.startswith('source_quality_')},
     'source_sampled_first100k_original_date_year_counts':{
       k:v for k,v in result.items() if k.startswith('original_first_year_')},
     'source_first100k_nonunique_quality_2025_and_accuracy_statistics':{
       k:v for k,v in result.items() if not (k.startswith('source_quality_') or k.startswith('original_first_year_'))},
     'original_source_grade_labels_case_normalized':True,
     'original_source_precise_obs_GPS_taxon_UUID_or_observer_NOT_republished':True,
     'not_a_global_source_census':True}
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True);a=p.parse_args()
 obj=run();o=Path(a.output);o.parent.mkdir(parents=True,exist_ok=True)
 o.write_text(json.dumps(obj,indent=2,ensure_ascii=False),encoding='utf8')
 print(json.dumps(obj,ensure_ascii=False),flush=True)
