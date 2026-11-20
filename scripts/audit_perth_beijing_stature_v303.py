"""Descriptive three-site panels; not independent tests of elevation effects."""
import argparse
import hashlib
import io
import json
import re
import zipfile
from pathlib import Path
import numpy as np
import openpyxl
import requests
from audit_SI_CZ_original_continuous_trait_coverage_v220 import read_leda, binomial


def run(supplement,leda):
    raw=Path(supplement).read_bytes()
    assert hashlib.sha256(raw).hexdigest()=='baea2386655fd1ddcd02b9d2dfd3fe79fc2402e2ee4d488f2ae44fbd1d31ea14'
    assert hashlib.sha256(Path(leda).read_bytes()).hexdigest()=='cbf5e761390d78a843512dfc1139acfc1154d46d6c1995a827354840a9f78a92'
    height,_=read_leda(leda,'height_m')
    z=zipfile.ZipFile(io.BytesIO(raw));files=[n for n in z.namelist() if n.endswith('.xlsx')]
    assert len(files)==1
    book=openpyxl.load_workbook(io.BytesIO(z.read(files[0])),data_only=True);rows=[]
    for col,city,year,cc,start in [(0,'Perth',2025,0,4),(2,'Perth',2024,0,9),(7,'Beijing',2025,5,4),(9,'Beijing',2024,5,9)]:
        current=None
        for row in book.worksheets[1].values:
            value=row[col]
            if not isinstance(value,str):continue
            if value.startswith('Site '):
                site=int(re.search(r'Site\s+(\d)',value)[1])
                lat,lon=map(float,book.worksheets[0].cell(start+site-1,cc+1).value.split(','))
                current=dict(city=city,year=year,site=site,typology=value.split('—')[-1].strip(),lat=lat,lon=lon,taxa=[])
                rows.append(current)
            elif current is not None and re.match(r'^[A-Z][a-z]+\s+[a-z]',value):current['taxa'].append(value.strip())
    assert len(rows)==12
    response=requests.get('https://api.opentopodata.org/v1/srtm30m',params={
        'locations':'|'.join(f"{r['lat']},{r['lon']}" for r in rows)},timeout=60)
    response.raise_for_status();body=response.json();assert body['status']=='OK' and len(body['results'])==12
    for row,item in zip(rows,body['results']):
        assert abs(item['location']['lat']-row['lat'])<1e-6
        assert abs(item['location']['lng']-row['lon'])<1e-6
        taxa=set(row['taxa']);known={binomial(t) for t in taxa if binomial(t) in height}
        row.update(elevation_m=item['elevation'],n_taxa=len(taxa),known=len(known),coverage=len(known)/len(taxa),
            mean_log_leda=float(np.mean([np.log(height[t]) for t in known])))
    panels=[]
    for city in ['Perth','Beijing']:
        for year in [2024,2025]:
            selected=[r for r in rows if r['city']==city and r['year']==year]
            x=np.array([r['elevation_m'] for r in selected]);y=np.array([r['mean_log_leda'] for r in selected])
            panels.append(dict(city=city,year=year,n=3,elevation_range_m=[float(x.min()),float(x.max())],
                descriptive_slope_per_100m=float(np.polyfit(x,y,1)[0]*100),
                n_5species_50pct_supported=sum(r['known']>=5 and r['coverage']>=.5 for r in selected)))
    moves=[]
    for city in ['Perth','Beijing']:
        for site in [1,2,3]:
            pair=[r for r in rows if r['city']==city and r['site']==site]
            if (pair[0]['lat'],pair[0]['lon'])!=(pair[1]['lat'],pair[1]['lon']):moves.append(dict(city=city,site=site))
    return dict(source_studies=1,panels=panels,changed_coordinates=moves,
        warning='Descriptive only: three sites, typology confounding, limited trait coverage; years not all same locations.',
        site_summaries=[{k:v for k,v in r.items() if k not in ['lat','lon','taxa']} for r in rows])


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--supplement',required=True);p.add_argument('--leda',required=True);p.add_argument('--output',required=True)
    a=p.parse_args();r=run(a.supplement,a.leda);Path(a.output).write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r['panels'],indent=2))
