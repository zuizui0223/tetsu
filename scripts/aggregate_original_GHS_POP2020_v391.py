#!/usr/bin/env python3
"""v391: official GHSL 2020 population counts, frozen equal-area grid."""
import argparse,hashlib,json,math,shutil,zipfile,tempfile
from pathlib import Path
import numpy as np
import rasterio
from rasterio.warp import reproject,Resampling
from rasterio.transform import from_origin
from aggregate_original_GHSL2020_global_builtup_v376b import X_START,Y_END,WIDTH,ROWS,COLS,fingerprint
NAME="GHS_POP_E2020_GLOBE_R2023A_54009_1000_V1_0"
URL=("https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL/GHS_POP_GLOBE_R2023A/"
     "GHS_POP_E2020_GLOBE_R2023A_54009_1000/V1-0/"+NAME+".zip")
NODATA=-200.
def grid(src):
    dst=np.full((ROWS,COLS),-1.,dtype=np.float64)
    reproject(rasterio.band(src,1),dst,src_transform=src.transform,src_crs=src.crs,
              src_nodata=src.nodata,dst_transform=from_origin(X_START,Y_END,WIDTH,WIDTH),
              dst_crs="EPSG:6933",dst_nodata=-1.,resampling=Resampling.sum,
              num_threads=2,warp_mem_limit=256,init_dest_nodata=True)
    return dst
def roll(src,n):
    if src.shape!=(ROWS,COLS) or ROWS%n or COLS%n:raise ValueError("Wrong frozen grid")
    return src.reshape(ROWS//n,n,COLS//n,n).sum(axis=(1,3))
def output_rows(pop,area,n):
    p=roll(pop,n);a=roll(area,n);w=WIDTH*n;out=[]
    for y,x in zip(*np.where(a>0)):
        z=float(p[y,x]);s=float(a[y,x])
        if z<0 or not math.isfinite(z) or not 0<s<=w*w*1.03:
            raise ValueError("Invalid aggregate population or source-support area")
        out.append({"source_equal_area_grid_ix":int(x),"source_equal_area_grid_iy":int(a.shape[0]-1-y),
                    "source_grid_width_m":w,"source_GHS_POP2020_people":z,
                    "source_GHS_POP2020_approx_valid_pixel_area_m2":s,
                    "source_GHS_POP2020_approx_valid_fraction":s/(w*w)})
    return out
def run(archive,headfile):
    headers=headfile.read_text(errors="replace")
    lines=headers.splitlines()
    # curl -I -L can log redirects; parse final HTTP block only.
    starts=[i for i,l in enumerate(lines) if l.upper().startswith("HTTP/")]
    if not starts:raise ValueError("No publisher HEAD response")
    h={k.strip().lower():v.strip() for l in lines[starts[-1]+1:] if ":" in l for k,v in [l.split(":",1)]}
    if int(h.get("content-length","-1"))!=archive.stat().st_size or archive.stat().st_size<10000000:
        raise ValueError("Publisher HEAD Content-Length mismatch or implausible ZIP")
    sha=fingerprint(archive)
    with tempfile.TemporaryDirectory() as td,zipfile.ZipFile(archive) as z:
        if z.testzip() is not None:raise ValueError("Original publisher ZIP CRC failed")
        t=[a for a in z.infolist() if a.filename.lower().endswith(".tif") and not a.is_dir()]
        if len(t)!=1 or NAME not in Path(t[0].filename).name:raise ValueError("Unexpected source epoch/TIFF identity")
        p=Path(td)/"pop.tif";mask=Path(td)/"source_support.tif"
        with z.open(t[0]) as f,p.open("wb") as o:shutil.copyfileobj(f,o,4<<20)
        with rasterio.open(p) as src:
            m={"source_shape":list(src.shape),"source_crs":str(src.crs),
               "source_transform":list(src.transform),"source_dtype":src.dtypes[0],
               "source_nodata":src.nodata,"source_band_count":src.count}
            print(json.dumps({"source_population_raster":m}),flush=True)
            if (src.count!=1 or src.crs is None or
                not ("moll" in str(src.crs).lower() or "54009" in str(src.crs)) or
                abs(src.transform.a-1000)>1 or abs(src.transform.e+1000)>1 or
                src.nodata!=NODATA or src.dtypes[0] not in ("float32","float64")):
                raise ValueError("Original GHS-POP metadata differs from publisher docs")
            profile=src.profile.copy()
            profile.update(dtype="uint8",nodata=255,compress="deflate",predictor=1,count=1)
            total=0.;valid_n=0;zeros=0;maxval=0.
            with rasterio.open(mask,"w",**profile) as out:
                for _,window in src.block_windows(1):
                    a=src.read(1,window=window,masked=True)
                    valid=~np.ma.getmaskarray(a)
                    values=np.asarray(a.data,dtype=np.float64)[valid]
                    if np.any(~np.isfinite(values)) or np.any(values<0):
                        raise ValueError("Invalid valid source population")
                    if values.size:
                        total+=float(values.sum());valid_n+=int(values.size)
                        zeros+=int((values==0).sum());maxval=max(maxval,float(values.max()))
                    out.write(valid.astype(np.uint8),1,window=window)
            pop=grid(src)
        with rasterio.open(mask) as s:area=grid(s)*1e6
    if not 7e9<total<8.5e9 or valid_n<100000000:
        raise ValueError("Original full-world population census implausible")
    valid=np.isfinite(area)&(area>0)
    if np.any(valid&(~np.isfinite(pop)|(pop<0))) or np.any((~valid)&(pop>1e-6)):
        raise ValueError("Original projected count and support inconsistent")
    pop=np.where(valid,pop,0.);area=np.where(valid,area,0.)
    projected=float(pop.sum());relative=abs(projected-total)/total
    if relative>.01:raise ValueError("Global projected population sum changed >1%; inspect projection")
    rows={str(WIDTH*n):output_rows(pop,area,n) for n in (1,2)}
    for r in rows.values():
        if not math.isclose(sum(x["source_GHS_POP2020_people"] for x in r),projected,rel_tol=1e-11):
            raise ValueError("Original population not conserved in exact rollup")
    return {"version":"v391-original-JRC-GHS-POP-2020-1km-global-counts",
            "source_URL":URL,"source_archive_bytes":archive.stat().st_size,"source_archive_SHA256":sha,
            "source_final_HTTP_headers":{k:h[k] for k in ("content-length","etag","last-modified") if k in h},
            "source_meta":m,"source_full_original_valid_pixels":valid_n,"source_zero_pixels":zeros,
            "source_max_people_per_pixel":maxval,"source_full_original_population_sum":total,
            "source_projected_population_sum":projected,"relative_population_sum_difference":relative,
            "source_nodata_NOT_zero":True,"source_projected_valid_area_approximation":True,
            "source_global_GHS_POP2020_by_grid":rows}
def main():
    p=argparse.ArgumentParser()
    for k in ("archive","head","out"):p.add_argument("--"+k,type=Path,required=True)
    a=p.parse_args();o=run(a.archive,a.head);a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(o,allow_nan=False))
    print(json.dumps({k:v for k,v in o.items() if k!="source_global_GHS_POP2020_by_grid"}),flush=True)
    print(json.dumps({k:len(v) for k,v in o["source_global_GHS_POP2020_by_grid"].items()}),flush=True)
if __name__=="__main__":main()
