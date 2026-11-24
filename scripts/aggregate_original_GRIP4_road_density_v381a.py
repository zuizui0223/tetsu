#!/usr/bin/env python3
"""v381a original GRIP4 road density, exact overlap with frozen global equal-area cells."""
import argparse, hashlib,json,math,tempfile,zipfile
from pathlib import Path
import numpy as np
from scipy.sparse import csr_matrix
from pyproj import Transformer
import rasterio
from scan_fullworld_2025_iNat_plant_macro_grid_v369 import X_START,Y_START
MD5="24d3ed5c079f7a898fbc4a5d68be337a"
W=100000
NX=348
NY=148
def weights(edges,targets):
    edges=np.asarray(edges,dtype=float);targets=np.asarray(targets,dtype=float)
    if np.any(np.diff(targets)<=0):raise ValueError("Target edges not increasing")
    lo=np.minimum(edges[:-1],edges[1:]);hi=np.maximum(edges[:-1],edges[1:])
    if not np.isfinite(edges).all() or np.any(lo>=hi):raise ValueError("Invalid source edges")
    rows=[];cols=[];values=[]
    for j,(a,b) in enumerate(zip(targets[:-1],targets[1:])):
        k=np.flatnonzero((lo<b)&(hi>a))
        v=np.minimum(hi[k],b)-np.maximum(lo[k],a)
        rows.extend([j]*len(k));cols.extend(k.tolist());values.extend(v.tolist())
    return csr_matrix((values,(rows,cols)),shape=(len(targets)-1,len(lo)))
def area_aggregate(values,valid,xw,yw):
    if values.shape!=valid.shape or values.shape!=(yw.shape[1],xw.shape[1]):
        raise ValueError("Source area matrix shape mismatch")
    v=np.where(valid,values,0.)
    area=(xw@(yw@valid.astype(float)).T).T
    numerator=(xw@(yw@v).T).T
    if np.any(area<0) or np.any(numerator<0):raise ValueError("Invalid source area")
    return numerator,area
def grid_rows(a,n,w):
    rows=[]
    for y,x in zip(*np.where(a>0)):
        valid=float(a[y,x]);num=float(n[y,x])
        if valid>w*w*(1+1e-9):raise ValueError("Raster overlap area exceeds grid area")
        rows.append({"source_equal_area_grid_ix":int(x),"source_equal_area_grid_iy":int(y),
            "source_grid_width_m":w,"source_GRIP4_valid_area_m2":valid,
            "source_GRIP4_valid_fraction":valid/(w*w),
            "source_GRIP4_observed_road_density_m_per_km2":num/valid,
            "source_GRIP4_implied_road_km_in_valid_area":num/1e9})
    return rows
def roll(a):return a.reshape(NY//2,2,NX//2,2).sum(axis=(1,3))
def run(archive):
    with archive.open("rb") as f:
        if hashlib.file_digest(f,"md5").hexdigest()!=MD5:raise ValueError("Publisher ZIP MD5 differs")
    with archive.open("rb") as f:sha=hashlib.file_digest(f,"sha256").hexdigest()
    with zipfile.ZipFile(archive) as z,tempfile.TemporaryDirectory() as td:
        members=[s for s in z.namelist() if s.lower().endswith(".asc")]
        print(json.dumps({"source_zip_members":z.namelist()}),flush=True)
        if len(members)!=1:raise ValueError("Expected one original density ASCII grid")
        dest=Path(td)/"original_total.asc"
        dest.write_bytes(z.read(members[0]))
        with rasterio.open(dest) as r:
            t=r.transform
            meta={"source_rows":r.height,"source_columns":r.width,"source_bounds":list(r.bounds),
                  "source_transform":list(t),"source_nodata":r.nodata,"source_crs":str(r.crs)}
            print(json.dumps(meta),flush=True)
            if r.count!=1 or r.width<4000 or r.height<1500 or r.width>4320 or r.height>2160:
                raise ValueError("Unexpected global 5 arcmin raster size")
            if r.crs is not None and r.crs.to_epsg()!=4326:
                raise ValueError("Unexpected source geographic CRS")
            if not np.allclose([t.a,t.e,t.b,t.d],[1/12,-1/12,0,0],rtol=0,atol=1e-7):
                raise ValueError("Not original 5 arcmin raster")
            lon=t.c+np.arange(r.width+1)*t.a;lat=t.f+np.arange(r.height+1)*t.e
            if lon.min()<-180.0001 or lon.max()>180.0001 or lat.min()<-90.0001 or lat.max()>90.0001:
                raise ValueError("Source extent out of world")
            tf=Transformer.from_crs(4326,6933,always_xy=True)
            x,_=tf.transform(np.clip(lon,-180,180),np.zeros(lon.size))
            _,y=tf.transform(np.zeros(lat.size),np.clip(lat,-90,90))
            xw=weights(x,X_START+np.arange(NX+1)*W)
            yw=weights(y,Y_START+np.arange(NY+1)*W)
            raw=r.read(1,masked=True)
            valid=~np.ma.getmaskarray(raw)
            v=np.asarray(raw.data,dtype=np.float64)
            if np.any(~np.isfinite(v[valid])) or np.any(v[valid]<0):
                raise ValueError("Negative/nonfinite original valid road density")
            if valid.sum()<100000:raise ValueError("Source raster support unexpectedly sparse")
            n,a=area_aggregate(v,valid,xw,yw)
            qa={"source_valid_pixels":int(valid.sum()),"source_zero_valid_pixels":int((valid&(v==0)).sum()),
                "source_missing_pixels":int((~valid).sum()),"source_valid_range":[float(v[valid].min()),float(v[valid].max())]}
    a2,n2=roll(a),roll(n)
    if not math.isclose(float(a.sum()),float(a2.sum()),rel_tol=1e-12) or not math.isclose(float(n.sum()),float(n2.sum()),rel_tol=1e-12):
        raise ValueError("Original road area/length not conserved")
    return {"version":"v381a-publisher-original-GRIP4-2018-global-road-density",
        "source_doi":"10.5281/zenodo.6420961","source_zip_MD5":MD5,"source_zip_SHA256":sha,
        "source_density_units":"m of road per km2 original pixel area","source_meta":meta,"source_QA":qa,
        "source_valid_projected_area_m2":float(a.sum()),"source_implied_road_km_in_valid_projected_area":float(n.sum()/1e9),
        "source_zero_density_NOT_imputed_for_nodata":True,
        "source_CRS_if_missing_assigned_from_publisher_5arcmin_WGS84_documentation":True,
        "source_historical_road_inventory_NOT_2026_full_road_network":True,
        "source_density_assumed_uniform_within_original_5arcmin_pixel":True,
        "source_global_GRIP4_road_density_by_grid":{"100000":grid_rows(a,n,W),"200000":grid_rows(a2,n2,2*W)}}
def main():
    p=argparse.ArgumentParser();p.add_argument("--archive",type=Path,required=True);p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    o=run(a.archive);a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(o,allow_nan=False))
    print(json.dumps({k:v for k,v in o.items() if k!="source_global_GRIP4_road_density_by_grid"}),flush=True)
    print(json.dumps({k:len(v) for k,v in o["source_global_GRIP4_road_density_by_grid"].items()}),flush=True)
if __name__=="__main__":main()
