"""Rebuild Australian analysis inputs from hash-verified public source files."""
import argparse
import hashlib
import json
import time
import zipfile
from pathlib import Path
import pandas as pd
import requests
from pyproj import Transformer
from extract_australian_rail_reserve_v303 import extract


def prepare(folder):
    folder=Path(folder);folder.mkdir(parents=True,exist_ok=True)
    root=Path(__file__).resolve().parents[1]
    manifest=json.loads((root/'docs/tetsu_v303_source_manifest.json').read_text())
    for item in manifest:
        if item['file'] not in ['wa_rail_reserve.pdf','austraits-7.0.0.zip']:continue
        target=folder/item['file']
        if not target.exists():
            with requests.get(item['url'],stream=True,timeout=90) as response:
                response.raise_for_status()
                with target.with_suffix(target.suffix+'.part').open('wb') as f:
                    for chunk in response.iter_content(1024*1024):f.write(chunk)
            target.with_suffix(target.suffix+'.part').replace(target)
        if hashlib.sha256(target.read_bytes()).hexdigest()!=item['sha256']:
            raise ValueError('Source checksum mismatch: '+str(target))
    data=extract(folder/'wa_rail_reserve.pdf')
    (folder/'wa_extracted.json').write_text(json.dumps(data,ensure_ascii=False,indent=2),encoding='utf-8')
    # Latin-1 is a lossless byte decoding for this release's mixed text encoding.
    # Selected taxon names, trait keys and numeric columns are ASCII; no errors are ignored.
    with zipfile.ZipFile(folder/'austraits-7.0.0.zip') as z:
        parts=[]
        for chunk in pd.read_csv(z.open('austraits-7.0.0/traits.csv'),dtype=str,encoding='latin1',chunksize=100000):
            parts.append(chunk[chunk.trait_name.str.contains('plant_height',na=False)])
        pd.concat(parts).to_csv(folder/'austraits_height.csv',index=False)
    transformer=Transformer.from_crs(28350,4326,always_xy=True)
    points=[]
    for m in data['metadata']:
        assert m['zone']==50
        lon,lat=transformer.transform(m['easting'],m['northing'])
        assert -25<lat<-20 and 118<lon<122
        points.append(dict(site=m['site'],lat=lat,lon=lon))
    for start in range(0,len(points),70):
        subset=points[start:start+70]
        response=requests.get('https://api.opentopodata.org/v1/srtm30m',params={
            'locations':'|'.join(f"{r['lat']},{r['lon']}" for r in subset)},timeout=60)
        response.raise_for_status();body=response.json()
        assert body['status']=='OK' and len(body['results'])==len(subset)
        for row,item in zip(subset,body['results']):
            assert abs(row['lat']-item['location']['lat'])<1e-6
            assert abs(row['lon']-item['location']['lng'])<1e-6
            assert item['elevation'] is not None
            row['elevation_m']=item['elevation']
        (folder/f'wa_dem_response_{start}.json').write_text(json.dumps(body,indent=2))
        if start+70<len(points):time.sleep(1.1)
    (folder/'wa_dem.json').write_text(json.dumps(dict(crs_source='EPSG:28350 GDA94 MGA zone 50',
        crs_target='EPSG:4326',source='OpenTopoData SRTM30m v3',points=points),indent=2))
    print('Verified public sources; rebuilt matrix, trait subset, and 138 DEM lookups.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--data',required=True);a=p.parse_args();prepare(a.data)
