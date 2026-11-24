#!/usr/bin/env python3
"""v369: scan the COMPLETE worldwide iNaturalist monthly photos-linked
public observations stream (13.70GB gzip) exactly ONCE, build 2025
research-grade wild PLANT SPECIES at global 100km and 200km grids.

NO original individual observations, observers, sensitive species GPS or
cell-by-species identities are uploaded. Source does not contain geoprivacy:
100km/200km outcomes are macroecological observation-opportunity maps
ONLY, never evidence for actual railway ballast plants within 80m.
"""
from __future__ import annotations
import argparse,csv,gzip,hashlib,io,json,math,time
from collections import Counter,defaultdict
from pathlib import Path

OBS_ETAG='4539360677c8190d9aa12dbf6982ee72-27'
TAXA_ETAG='0a6a162c57bb2f24ea4f00c76c20080d'
OBS_BYTES=13697505636
TAXA_BYTES=40184871
S3='inaturalist-open-data'
REGION='us-east-1'
ROOT_PLANT_TAXON=47126
YEAR='2025'
OBS_QUALITY='research_grade'
X_START=-17_400_000
Y_START=-7_400_000
X_END=17_400_000
Y_END=7_400_000
CELL_SIZE_M=100_000
SOURCE_REGION_M=(100_000,200_000)
MAX_UNCERT_M=1000
PROGRESS_EVERY_SOURCE_ROWS=10_000_000

def source_client():
    import boto3
    from botocore import UNSIGNED
    from botocore.config import Config
    return boto3.client('s3',region_name=REGION,config=Config(
        signature_version=UNSIGNED,region_name=REGION,
        connect_timeout=30,read_timeout=120,retries={'max_attempts':5}))

class HashReader:
    """gzip reads source compressed original S3 bytes in one pass, hashes in transit."""
    def __init__(self,body):
        self.body=body
        self.hash=hashlib.sha256()
        self.bytes=0
    def read(self,n=-1):
        b=self.body.read(n)
        self.hash.update(b)
        self.bytes+=len(b)
        return b
    def close(self):self.body.close()

def open_original_source(client,key,etag,expected_bytes):
    head=client.head_object(Bucket=S3,Key=key)
    if head['ETag'].strip('"')!=etag or int(head['ContentLength'])!=expected_bytes:
        raise ValueError('Original 2026-10-27 monthly iNaturalist S3 source ETag or length changed: '+key)
    resp=client.get_object(Bucket=S3,Key=key,IfMatch='"'+etag+'"')
    if int(resp['ContentLength'])!=expected_bytes:
        raise ValueError('Downloaded source S3 object compressed length changed')
    wrapped=HashReader(resp['Body'])
    text=io.TextIOWrapper(gzip.GzipFile(fileobj=wrapped),encoding='utf-8',
                        newline='',errors='strict')
    rows=csv.DictReader(text,delimiter='\t',quoting=csv.QUOTE_NONE)
    return rows,wrapped,text

def source_ancestry_ids(value):
    if value is None:return []
    return [k for k in str(value).strip().replace(chr(92),"/").split("/") if k]


def taxonomy_source(client):
    rows,raw,stream=open_original_source(client,'taxa.csv.gz',TAXA_ETAG,TAXA_BYTES)
    try:
        wanted={'taxon_id','ancestry','rank','active','name'}
        if not wanted<=set(rows.fieldnames or []):
            raise ValueError('Original public iNaturalist taxon hierarchy columns missing')
        taxa={};qa=Counter()
        for row in rows:
            qa['original_global_iNat_all_taxonomy_records']+=1
            tid=int(row['taxon_id'])
            if tid in taxa:raise ValueError('Original global taxon ID duplicated')
            rank=row['rank'].strip().lower()
            is_active=str(row['active']).strip().lower() in ('t','true','1')
            # Original iNaturalist open metadata contains slash-separated
            # ancestry in the current 2026 source despite README describing
            # a backslash. Support both; never classify 1.3m species as
            # nonplants because an original taxonomy delimiter changed.
            ancestors=source_ancestry_ids(row['ancestry'])
            if not ancestors:ancestors=[]
            # Compressed original taxonomy ancestor IDs, not original names
            taxa[tid]=(rank,is_active,ancestors)
        qa['original_source_taxon_bytes_decompressed'] = None
    finally:stream.close()
    if raw.bytes!=TAXA_BYTES:raise ValueError('Incomplete original full taxon compressed S3 source')
    plant_species=set()
    for tid,(rank,active,anc) in taxa.items():
        if rank=='species' and active and str(ROOT_PLANT_TAXON) in anc:
            plant_species.add(tid)
    # Public photo-linked open data taxa need not include 100,000
    # independently observed active plant species. Diagnose the ACTUAL
    # source rank/ancestry coverage before interpreting a smaller count.
    raw_rank_count=Counter(rank for rank,active,anc in taxa.values()
                           if rank=='species')
    raw_ancestry_plant_substring_count=sum(
        '47126' in anc for rank,active,anc in taxa.values() if rank=='species')
    print(json.dumps({'SOURCE_TAXONOMY_PREFLIGHT_ONLY':True,
        'n_original_taxonomy_rows':len(taxa),
        'n_original_all_species_rank':raw_rank_count['species'],
        'n_original_species_with_plant_ancestor_47126':raw_ancestry_plant_substring_count,
        'n_original_ACTIVE_plant_species':len(plant_species),
        'source_first_five_plant_taxon_ancestry_only_NO_species_names':[
            taxa[k][2][:6] for k in list(plant_species)[:5]]
        }),flush=True)
    if len(plant_species)<10_000:raise ValueError(
        'Original global active plant species taxon tree unexpectedly small; correct taxonomy source parser before continuing')
    canonical={k:k for k in plant_species}
    infraspecies=('subspecies','variety','form','subvariety')
    for tid,(rank,active,anc) in taxa.items():
        if not active or rank not in infraspecies or str(ROOT_PLANT_TAXON) not in anc:continue
        hit=next((int(x) for x in reversed(anc) if x and int(x) in plant_species),None)
        if hit is not None:canonical[tid]=hit
    # Taxonomic tree parent Genus / Family from the ORIGINAL dynamic source,
    # not a substitute for LIFE-FORM or plant growth/dispersal functional traits.
    genus={};family={}
    for tid in plant_species:
        _,_,ancestry=taxa[tid]
        for ancestor in reversed(ancestry):
            if not ancestor:continue
            parent=taxa.get(int(ancestor))
            if parent is None:continue
            if tid not in genus and parent[0]=='genus':genus[tid]=int(ancestor)
            if tid not in family and parent[0]=='family':family[tid]=int(ancestor)
            if tid in genus and tid in family:break
    qa.update({
        'original_plant_species_taxa_active':len(plant_species),
        'original_plant_species_or_infrafamily_taxa_to_species_ID_map_n':len(canonical),
        'source_species_with_original_genus_ID':len(genus),
        'source_species_with_original_family_ID':len(family)})
    return canonical,genus,family,dict(qa),raw.hash.hexdigest()

def source_coarse_cell(lon,lat,transform,grid_m):
    if not all(math.isfinite(v) for v in (lon,lat)) or not(-180<=lon<=180 and -90<=lat<=90):
        return None
    x,y=transform.transform(lon,lat)
    if not (X_START<=x<X_END and Y_START<=y<Y_END):return None
    return int((x-X_START)//grid_m),int((y-Y_START)//grid_m)

def cell_contents():
    return {'original_2025_public_photo_observations':0,
            'source_plant_species':set(),'source_observers':set()}

MIN_PAIR_CELL_OBS=50
MIN_PAIR_CELL_OBSERVERS=5
MIN_PAIR_CELL_PLANT_SPECIES=20

def source_local_taxonomic_beta_pairs(cell_data,cell_m):
    """Global coarsely georeferenced adjacent source flora comparison.

    Both cells must have independently adequate source public plant
    sampling effort. This Jaccard describes the ENTIRE 100km/200km
    landscape, NOT railway ballast vegetation. The raw 2025 taxon IDs,
    rare-species sensitive lat/long and observer IDs do not leave CI.
    """
    pairs=[]
    valid={}
    for key,entry in cell_data.items():
        if (entry['original_2025_public_photo_observations']>=MIN_PAIR_CELL_OBS
            and len(entry['source_observers'])>=MIN_PAIR_CELL_OBSERVERS
            and len(entry['source_plant_species'])>=MIN_PAIR_CELL_PLANT_SPECIES):
            valid[key]=entry
    for (i,j),left in valid.items():
        for dx,dy in ((1,0),(0,1)):
            other=valid.get((i+dx,j+dy))
            if other is None:continue
            a=left['source_plant_species'];b=other['source_plant_species']
            union=len(a|b)
            if union<=0:raise ValueError('Original world macro beta zero union despite source richness gate')
            shared=len(a&b)
            pairs.append({
              'source_grid_width_m':cell_m,
              'source_cell_A_original_grid_ix_iy':[i,j],
              'source_cell_B_original_grid_ix_iy':[i+dx,j+dy],
              'source_biological_Jaccard_similarity_between_ENTIRE_landscape_cells':round(shared/union,8),
              'source_shared_plants_species_n':shared,
              'source_union_plants_species_n':union,
              'source_grid_A_2025_photo_observations_n':left['original_2025_public_photo_observations'],
              'source_grid_B_2025_photo_observations_n':other['original_2025_public_photo_observations'],
              'source_grid_A_distinct_original_observers_n':len(left['source_observers']),
              'source_grid_B_distinct_original_observers_n':len(other['source_observers']),
              'source_grid_A_original_plant_species_n':len(a),
              'source_grid_B_original_plant_species_n':len(b),
              'source_original_full_landscape_beta_NOT_railroad_field_plot_beta':True})
    return pairs,{"n_source_cells_meeting_effort_and_taxon_source_requirements":len(valid),
        "n_original_adjacent_source_landscape_beta_pairs":len(pairs),
        "sampling_gate_min_original_2025_photo_observations_per_cell":MIN_PAIR_CELL_OBS,
        "sampling_gate_min_original_distinct_observers_per_cell":MIN_PAIR_CELL_OBSERVERS,
        "sampling_gate_min_original_source_species_per_cell":MIN_PAIR_CELL_PLANT_SPECIES,
        "source_pair_records_dependent_and_NOT_independent_biological_studies":True}


def count_original_2025_plants(client,canonical,genus,family):
    from pyproj import Transformer
    transform=Transformer.from_crs('EPSG:4326','EPSG:6933',always_xy=True)
    rows,raw,stream=open_original_source(client,'observations.csv.gz',
                                        OBS_ETAG,OBS_BYTES)
    columns=set(rows.fieldnames or [])
    required={'observer_id','latitude','longitude','positional_accuracy',
              'taxon_id','quality_grade','observed_on'}
    if not required<=columns:
        raise ValueError('Global original iNaturalist observations source field drift')
    if {'geoprivacy','coordinates_obscured','taxon_geoprivacy'}&columns:
        raise ValueError('Source geoprivacy schema changed, revalidate preregistration')
    qa=Counter()
    cells={m:defaultdict(cell_contents) for m in SOURCE_REGION_M}
    nbad=0;t0=time.monotonic()
    try:
        for row in rows:
            qa['original_global_iNat_photo_linked_rows_read']+=1
            if (qa['original_global_iNat_photo_linked_rows_read']%PROGRESS_EVERY_SOURCE_ROWS)==0:
                print(json.dumps({'source_original_global_photo_rows_read':qa[
                    'original_global_iNat_photo_linked_rows_read'],
                    'source_2025_valid_research_plant_rows':qa[
                    'source_year_2025_public_research_grade_georeferenced_plant_obs'],
                    'source_seconds_elapsed':round(time.monotonic()-t0,1)}),flush=True)
            date=str(row['observed_on'])
            if not date.startswith(YEAR+'-'):
                continue
            qa['original_2025_photo_linked_all_taxa_source_records']+=1
            quality=str(row['quality_grade']).strip().lower().replace(' ','_').replace('-','_')
            if quality!=OBS_QUALITY:
                qa['original_2025_plant_source_failed_research_grade_or_other_taxa']+=1
                continue
            sid=row['taxon_id']
            try:species=canonical.get(int(sid))
            except (ValueError,TypeError):species=None
            if species is None:
                qa['original_2025_research_grade_not_confirmed_plant_species_taxon']+=1
                continue
            try:
                uncertainty=float(row['positional_accuracy'])
                lon=float(row['longitude']);lat=float(row['latitude'])
            except (TypeError,ValueError):
                qa['source_plant_observation_no_finite_uncertainty_or_WGS84']+=1
                continue
            if (not math.isfinite(uncertainty) or uncertainty<0
                or uncertainty>MAX_UNCERT_M):
                qa['source_plant_observation_excluded_original_accuracy_gt1000m']+=1
                continue
            large=source_coarse_cell(lon,lat,transform,200_000)
            small=source_coarse_cell(lon,lat,transform,100_000)
            if large is None or small is None:
                qa['source_plant_outside_full_global_fixed_equalarea_grid']+=1
                continue
            qa['source_year_2025_public_research_grade_georeferenced_plant_obs']+=1
            observer=row['observer_id']
            if not observer:
                observer=None
            for k,cell in ((100_000,small),(200_000,large)):
                val=cells[k][cell]
                val['original_2025_public_photo_observations']+=1
                val['source_plant_species'].add(species)
                if observer is not None:
                    val['source_observers'].add(int(observer))
    finally:stream.close()
    if raw.bytes!=OBS_BYTES:
        raise ValueError('FULL ORIGINAL 13.7GB S3 global plant observation snapshot NOT completely streamed')
    if qa['original_global_iNat_photo_linked_rows_read']<50_000_000:
        raise ValueError('Source original global observation file unexpectedly truncated')
    result={}
    for km,cmap in cells.items():
        items=[]
        for (ix,iy),val in sorted(cmap.items()):
            taxa=val['source_plant_species']
            items.append({'source_equal_area_grid_ix':ix,'source_equal_area_grid_iy':iy,
                'source_grid_width_m':km,
                'source_year_2025_photo_linked_plants_observation_count':
                    val['original_2025_public_photo_observations'],
                'source_2025_plant_species_richness':len(taxa),
                'source_2025_plant_genus_richness':len({genus[k] for k in taxa if k in genus}),
                'source_2025_plant_family_richness':len({family[k] for k in taxa if k in family}),
                'source_2025_distinct_observers_with_photo_linked_plant_records':
                    len(val['source_observers'])})
        original_pairs,pairQA=source_local_taxonomic_beta_pairs(cmap,km)
        result[str(km)]={'n_global_occupied_observation_grid_cells':len(items),
                'source_full_world_100km_or_200km_gridded_plant_effort_and_taxonomic_counts':items,
                'source_fullworld_adjacent_100km_200km_ORIGINAL_LANDSCAPE_botanical_Jaccard_pairs':original_pairs,
                'source_original_macro_landscape_beta_opportunity_QA':pairQA}
    return result,dict(qa),raw.hash.hexdigest()

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--output-summary',required=True)
    p.add_argument('--output-grid',required=True)
    args=p.parse_args()
    client=source_client()
    species,genus,family,tqa,tsha=taxonomy_source(client)
    cells,obsqa,osha=count_original_2025_plants(client,species,genus,family)
    header={'version':'v369-original-2026-10-27-global-photo-linked-plant-2025-100km-and-200km-source',
      'original_source_iNat_obs_S3_ETag':OBS_ETAG,
      'source_iNat_original_2026_10_27_observations_gzip_SHA256_full':osha,
      'original_source_iNat_taxa_S3_ETag':TAXA_ETAG,
      'source_iNat_original_2026_10_27_taxa_gzip_SHA256_full':tsha,
      'original_total_source_observation_compressed_bytes':OBS_BYTES,
      'source_year_2025':YEAR,
      'source_source_plant_rank':'active species plus canonical original infraspecies mapped to parent species',
      'source_unobscured_original_observations_NOT_confirmed':True,
      'original_coarse_grid_100km_and_200km_are_OBSERVATION_OPPORTUNITY_not_causal_rail_species_presence':True,
      'source_no_full_site_plot_species_absences':True,
      'source_research_grade_photo_linked_NOT_all_global_plant_observations':True,
      'source_global_taxa_QA':tqa,
      'source_global_photo_observations_QA':obsqa,
      'source_grid_n_by_width':{x:len(z['source_full_world_100km_or_200km_gridded_plant_effort_and_taxonomic_counts'])
              for x,z in cells.items()},
      'original_taxonomic_beta_is_entire_adjacent_landscape_flora_NOT_rail_specific_plot':True,
      'source_original_macro_landscape_beta_pair_QA_by_grid_width':{
          x:z['source_original_macro_landscape_beta_opportunity_QA'] for x,z in cells.items()},
      'original_specific_plant_observer_IDs_species_IDs_and_GPS_NOT_PUBLISHED':True}
    dest=Path(args.output_grid);dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps({'source_summary':header,'source_grids':cells},
          ensure_ascii=False,indent=2),encoding='utf8')
    short=Path(args.output_summary);short.parent.mkdir(parents=True,exist_ok=True)
    short.write_text(json.dumps(header,ensure_ascii=False,indent=2),encoding='utf8')
    print(json.dumps(header,ensure_ascii=False),flush=True)

if __name__=='__main__':main()
