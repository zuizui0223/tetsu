#!/usr/bin/env python3
"""v380a independent original EarthEnv consensus tree cover, exact prior grid IDs."""
import argparse, hashlib, json, math
from pathlib import Path
from contextlib import ExitStack
import numpy as np
import rasterio
from rasterio.windows import Window
from pyproj import Transformer
from scan_fullworld_2025_iNat_plant_macro_grid_v369 import X_START,Y_START
X_END=17400000
Y_END=7400000
W=100000
NX=(X_END-X_START)//W
NY=(Y_END-Y_START)//W
CLASS_IDS=(1,2,3,4)
BASE_URL="https://data.earthenv.org/consensus_landcover/with_DISCover/"
def file_sha(p):
    with p.open("rb") as f:return hashlib.file_digest(f,"sha256").hexdigest()
def rollup(x):
    if x.shape!=(NY,NX) or NY%2 or NX%2:raise ValueError("Not fixed original global 100km grid")
    return x.reshape(NY//2,2,NX//2,2).sum(axis=(1,3))
def validate_layers(ds):
    ref=ds[0]
    if (ref.width!=43200 or ref.height!=17520 or
        ref.crs is None or ref.crs.to_epsg()!=4326 or ref.count!=1 or
        not np.allclose([ref.transform.a,ref.transform.e,ref.transform.b,ref.transform.d],
                        [1/120,-1/120,0,0],rtol=0,atol=1e-9) or
        not np.allclose(list(ref.bounds),[-180,-56,180,90],rtol=0,atol=1e-5)):
        raise ValueError("Original EarthEnv 30arcsec global 90N-56S raster metadata mismatch")
    for d in ds:
        if (d.shape!=ref.shape or d.crs!=ref.crs or d.transform!=ref.transform or
            d.count!=1 or d.dtypes[0]!="uint8"):
            raise ValueError("Original consensus tree layers not aligned unsigned8")
def source_rows(treearea,validarea,w):
    rows=[]
    for iy,ix in zip(*np.where(validarea>0)):
        va=float(validarea[iy,ix]);ta=float(treearea[iy,ix])
        if ta<0 or ta>va+1e-3 or va>w*w*1.05:
            raise ValueError("Original equal area approximate source-pixel support inconsistent")
        rows.append({"source_equal_area_grid_ix":int(ix),
            "source_equal_area_grid_iy":int(iy),
            "source_grid_width_m":w,"source_valid_pixel_area_m2":va,
            "source_tree_consensus_area_m2":ta,
            "source_tree_consensus_mean_fraction_of_valid_source_area":ta/va,
            "source_valid_area_fraction_of_grid_approx":va/(w*w)})
    return rows
def aggregate(paths):
    sha={str(c):file_sha(p) for c,p in zip(CLASS_IDS,paths)}
    area=np.zeros((NY,NX),dtype="float64")
    tree=np.zeros_like(area)
    qa={"joint_valid_source_pixels_in_grid":0,"joint_valid_sum_above100":0,
        "joint_valid_sum_above104":0,"source_mask_disagreement_pixels_in_grid":0}
    tf=Transformer.from_crs("EPSG:4326","EPSG:6933",always_xy=True)
    with ExitStack() as stack:
        ds=[stack.enter_context(rasterio.open(p)) for p in paths]
        validate_layers(ds)
        r=ds[0]
        xedge,_=tf.transform(r.transform.c+np.arange(r.width+1)*r.transform.a,
                             np.zeros(r.width+1))
        xm=(xedge[1:]+xedge[:-1])/2
        ix=np.floor((xm-X_START)/W).astype(int)
        dx=np.diff(xedge)
        if np.any(dx<=0) or ix.min()<0 or ix.max()>=NX or np.any(np.diff(ix)<0):
            raise ValueError("Unexpected projected source longitude axis")
        starts=np.r_[0,np.flatnonzero(np.diff(ix))+1]
        xbins=ix[starts]
        _,yedge=tf.transform(np.zeros(r.height+1),r.transform.f+np.arange(r.height+1)*r.transform.e)
        ym=(yedge[1:]+yedge[:-1])/2
        iy=np.floor((ym-Y_START)/W).astype(int)
        dy=-np.diff(yedge)
        if np.any(dy<=0):raise ValueError("Unexpected projected latitude axis")
        for a in range(0,r.height,64):
            n=min(64,r.height-a)
            y=iy[a:a+n]
            ok_y=(y>=0)&(y<NY)
            if not np.any(ok_y):continue
            arr=[d.read(1,window=Window(0,a,r.width,n),masked=True) for d in ds]
            masks=np.stack([~np.ma.getmaskarray(z) for z in arr])
            values=np.stack([np.asarray(z.data,dtype=np.int16) for z in arr])
            if np.any(masks & ((values<0)|(values>100))):
                raise ValueError("Valid original EarthEnv pixel outside documented 0-100; inspect NoData")
            good=masks.all(axis=0)
            mismatch=masks.any(axis=0)&~good
            sums=values.sum(axis=0)
            qa["joint_valid_source_pixels_in_grid"]+=int(good[ok_y].sum())
            qa["joint_valid_sum_above100"]+=int((good[ok_y]&(sums[ok_y]>100)).sum())
            qa["joint_valid_sum_above104"]+=int((good[ok_y]&(sums[ok_y]>104)).sum())
            qa["source_mask_disagreement_pixels_in_grid"]+=int(mismatch[ok_y].sum())
            # Four rounded percentage layers are from the same consensus
            # classification. Only tiny rounding excursions may be clipped.
            if qa["joint_valid_sum_above104"]:
                raise ValueError("Tree consensus sum exceeds rounding tolerance; no forest proxy estimated")
            pixelarea=dy[a:a+n,None]*dx[None,:]
            v=np.reduceat(np.where(good,pixelarea,0.),starts,axis=1)
            t=np.reduceat(np.where(good,pixelarea*np.minimum(sums,100)/100,0.),starts,axis=1)
            for j in np.flatnonzero(ok_y):
                area[y[j],xbins]+=v[j];tree[y[j],xbins]+=t[j]
            if a%2048==0:print(json.dumps({"source_rows_sofar":a+n,**qa}),flush=True)
        meta={"source_bounds":list(r.bounds),"source_pixel_width_degrees":r.transform.a,
              "source_raster_width":r.width,"source_raster_height":r.height,
              "source_nodata_by_class":[d.nodata for d in ds]}
    if qa["joint_valid_source_pixels_in_grid"]<100000000:
        raise ValueError("Original full world consensus valid land source unexpectedly sparse")
    rows100=source_rows(tree,area,W)
    a2,t2=rollup(area),rollup(tree)
    if not math.isclose(float(area.sum()),float(a2.sum()),rel_tol=1e-12) or not math.isclose(float(tree.sum()),float(t2.sum()),rel_tol=1e-12):
        raise ValueError("Source area not conserved across exact 2x2 rollup")
    return {"version":"v380a-EarthEnv-2014-original-full-consensus-tree-1km",
        "source_citation":"Tuanmu and Jetz 2014 GEB 23:1031-1045",
        "source_layer_urls":{str(c):BASE_URL+f"consensus_full_class_{c}.tif" for c in CLASS_IDS},
        "source_file_SHA256":sha,"source_raster_QA":meta,"source_pixel_QA":qa,
        "source_tree_classes":"1 needleleaf, 2 evergreen broadleaf, 3 deciduous broadleaf, 4 mixed/other trees",
        "source_total_valid_area_m2":float(area.sum()),"source_total_tree_consensus_area_m2":float(tree.sum()),
        "source_pixel_centres_assign_to_equalarea_grid_exact_projected_pixel_area":True,
        "source_tree_fraction_denominator_is_joint_valid_raster_area_NOT_all_cell_area":True,
        "source_tree_consensus_not_recent_2025_observed_canopy_or_forest_resistance":True,
        "source_global_100km_200km_tree_consensus":{str(W):rows100,str(2*W):source_rows(t2,a2,2*W)}}
def main():
    p=argparse.ArgumentParser();p.add_argument("--source-dir",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True);a=p.parse_args()
    o=aggregate([a.source_dir/f"consensus_full_class_{c}.tif" for c in CLASS_IDS])
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(o,allow_nan=False))
    print(json.dumps({k:v for k,v in o.items() if k!="source_global_100km_200km_tree_consensus"}),flush=True)
    print(json.dumps({k:len(v) for k,v in o["source_global_100km_200km_tree_consensus"].items()}),flush=True)
if __name__=="__main__":main()
