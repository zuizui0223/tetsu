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
def aggregate_land_density(density,valid,land_km2,land_valid,x,y,xw,yw):
    # Density denominator is publisher land area, NOT geographical pixel area.
    # Allocate land (and implied road length) proportional to exact
    # intersection with each target grid, assuming within-pixel uniformity.
    pixel=np.abs(np.diff(y))[:,None]*np.abs(np.diff(x))[None,:]
    land=np.where(land_valid,land_km2*1e6,0.)
    if np.any(~np.isfinite(land)) or np.any(land<0) or np.any(land>1.05*pixel+1):
        raise ValueError("Publisher land area inconsistent with equal-area source pixel geometry")
    proportion=land/pixel
    included=land_valid&valid&(land>0)
    all_land,_=area_aggregate(proportion,land_valid,xw,yw)
    known_land,_=area_aggregate(proportion,included,xw,yw)
    weighted_road,_=area_aggregate(np.where(included,density,0.)*proportion,included,xw,yw)
    return all_land,known_land,weighted_road
def grid_rows(land,a,n,w):
    rows=[]
    for y,x in zip(*np.where(land>0)):
        total=float(land[y,x]);valid=float(a[y,x]);num=float(n[y,x])
        if total>w*w*1.05 or valid>total+1e-3 or valid<0 or num<0:
            raise ValueError("Raster overlap land area invalid")
        rows.append({"source_equal_area_grid_ix":int(x),"source_equal_area_grid_iy":int(y),
            "source_grid_width_m":w,"source_GRIP4_total_land_area_m2":total,
            "source_GRIP4_land_fraction_of_cell":total/(w*w),
            "source_GRIP4_valid_road_land_area_m2":valid,
            "source_GRIP4_road_coverage_of_source_land":valid/total,
            "source_GRIP4_observed_road_density_m_per_km2_land":num/valid if valid else None,
            "source_GRIP4_implied_road_km_in_valid_land":num/1e9 if valid else None})
    return rows
def roll(a):return a.reshape(NY//2,2,NX//2,2).sum(axis=(1,3))
def run(archive):
    with archive.open("rb") as f:
        if hashlib.file_digest(f,"md5").hexdigest()!=MD5:raise ValueError("Publisher ZIP MD5 differs")
    with archive.open("rb") as f:sha=hashlib.file_digest(f,"sha256").hexdigest()
    with zipfile.ZipFile(archive) as z,tempfile.TemporaryDirectory() as td:
        members={Path(s).name.lower():s for s in z.namelist() if not s.endswith("/")}
        print(json.dumps({"source_zip_members":z.namelist()}),flush=True)
        required=("grip4_total_dens_m_km2.asc","grip4_area_land_km2.asc","readme.txt")
        if not set(required)<=set(members):raise ValueError("Publisher original density/land/readme bundle incomplete")
        print(z.read(members["readme.txt"]).decode("utf8",errors="replace"),flush=True)
        dest=Path(td)/"original_total.asc";landpath=Path(td)/"original_land.asc"
        dest.write_bytes(z.read(members[required[0]]))
        landpath.write_bytes(z.read(members[required[1]]))
        with rasterio.open(dest) as r,rasterio.open(landpath) as lr:
            t=r.transform
            meta={"source_rows":r.height,"source_columns":r.width,"source_bounds":list(r.bounds),
                  "source_transform":list(t),"source_nodata":r.nodata,"source_crs":str(r.crs),
                  "source_land_nodata":lr.nodata,"source_land_crs":str(lr.crs)}
            print(json.dumps(meta),flush=True)
            if r.count!=1 or r.width<4000 or r.height<1500 or r.width>4320 or r.height>2160:
                raise ValueError("Unexpected global 5 arcmin raster size")
            if lr.count!=1 or lr.shape!=r.shape or lr.transform!=t:
                raise ValueError("Original road and land rasters not aligned")
            for ds in (r,lr):
                if ds.crs is not None and ds.crs.to_epsg()!=4326:
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
            raw=r.read(1,masked=True);lraw=lr.read(1,masked=True)
            valid=~np.ma.getmaskarray(raw);lv=~np.ma.getmaskarray(lraw)
            v=np.asarray(raw.data,dtype=np.float64);l=np.asarray(lraw.data,dtype=np.float64)
            if np.any(~np.isfinite(v[valid])) or np.any(v[valid]<0):
                raise ValueError("Negative/nonfinite original valid road density")
            if np.any(~np.isfinite(l[lv])) or np.any(l[lv]<0) or valid.sum()<100000:
                raise ValueError("Invalid original land area or sparse road coverage")
            land,a,n=aggregate_land_density(v,valid,l,lv,x,y,xw,yw)
            qa={"source_valid_road_pixels":int(valid.sum()),"source_zero_valid_road_pixels":int((valid&(v==0)).sum()),
                "source_missing_road_pixels":int((~valid).sum()),"source_valid_land_pixels":int(lv.sum()),
                "source_positive_land_without_road_density_pixels":int((lv&(l>0)&~valid).sum()),
                "source_road_density_valid_range":[float(v[valid].min()),float(v[valid].max())],
                "source_land_km2_valid_range":[float(l[lv].min()),float(l[lv].max())]}
    l2,a2,n2=roll(land),roll(a),roll(n)
    for old,new in ((land,l2),(a,a2),(n,n2)):
        if not math.isclose(float(old.sum()),float(new.sum()),rel_tol=1e-12):
            raise ValueError("Original road area/length not conserved")
    return {"version":"v381b-publisher-original-GRIP4-2018-global-road-land-density",
        "source_doi":"10.5281/zenodo.6420961","source_zip_MD5":MD5,"source_zip_SHA256":sha,
        "source_density_units":"m of road per km2 original pixel LAND area","source_meta":meta,"source_QA":qa,
        "source_total_land_area_m2":float(land.sum()),"source_valid_road_land_area_m2":float(a.sum()),
        "source_implied_road_km_in_valid_land":float(n.sum()/1e9),
        "source_zero_density_NOT_imputed_for_nodata":True,
        "source_CRS_if_missing_assigned_from_publisher_5arcmin_WGS84_documentation":True,
        "source_historical_road_inventory_NOT_2026_full_road_network":True,
        "source_density_and_land_assumed_uniform_within_original_5arcmin_pixel":True,
        "source_global_GRIP4_road_density_by_grid":{"100000":grid_rows(land,a,n,W),"200000":grid_rows(l2,a2,n2,2*W)}}
def main():
    p=argparse.ArgumentParser();p.add_argument("--archive",type=Path,required=True);p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    o=run(a.archive);a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(o,allow_nan=False))
    print(json.dumps({k:v for k,v in o.items() if k!="source_global_GRIP4_road_density_by_grid"}),flush=True)
    print(json.dumps({k:len(v) for k,v in o["source_global_GRIP4_road_density_by_grid"].items()}),flush=True)
if __name__=="__main__":main()
