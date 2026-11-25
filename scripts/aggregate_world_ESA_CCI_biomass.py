#!/usr/bin/env python3
"""Official ESA CCI Biomass v7 0.1-degree woody AGB -> fixed EPSG:6933 grids.
Exact projected rectangle intersection; valid source pixels only.
"""
import argparse,hashlib,json,math
from pathlib import Path
import numpy as np
from scipy.sparse import csr_matrix
from pyproj import Transformer
from scan_fullworld_2025_iNat_plant_macro_grid_v369 import X_START,X_END,Y_START,Y_END
DOI="10.5285/6429d1aafe1e43b9b414e4a5a7f8b903"
URL="https://dap.ceda.ac.uk/neodc/esacci/biomass/data/agb/maps/v7.0/netcdf/ESACCI-BIOMASS-L4-AGB-MERGED-10000m-fv7.0.nc"
W=100000
NX=(X_END-X_START)//W
NY=(Y_END-Y_START)//W
YEARS=(2020,2024)

def edges(centers,bounds,lo,hi):
    v=np.asarray(centers,float);i=np.argsort(v);x=v[i]
    if x.ndim!=1 or len(x)<2 or not np.all(np.isfinite(x)) or not np.allclose(np.diff(x),.1,atol=5e-4,rtol=0):
        raise ValueError("Not documented 0.1-degree coordinate axis")
    if bounds is not None:
        b=np.sort(np.asarray(bounds,float)[i],axis=1)
        if b.shape!=(len(x),2) or not np.allclose(b[:-1,1],b[1:,0],atol=5e-4,rtol=0) or np.any((x<b[:,0])|(x>b[:,1])):raise ValueError("Invalid source coordinate bounds")
        e=np.r_[b[0,0],b[:,1]]
    else:e=np.r_[x[0]-.05,(x[1:]+x[:-1])/2,x[-1]+.05]
    if e[0]<lo-.001 or e[-1]>hi+.001:raise ValueError("Axis outside geographic domain")
    e=np.clip(e,lo,hi)
    if np.any(np.diff(e)<=0):raise ValueError("Nonpositive pixel extent")
    return i,e

def overlap_matrix(e,start,n,w=W):
    rows=[];cols=[];data=[]
    for i,(a,b) in enumerate(zip(e[:-1],e[1:])):
        if b<=a:raise ValueError("Nonmonotonic projected edges")
        for j in range(max(0,math.floor((a-start)/w)),min(n-1,math.floor((b-start)/w))+1):
            d=max(0.,min(b,start+(j+1)*w)-max(a,start+j*w))
            if d>0:rows.append(j);cols.append(i);data.append(d)
    z=csr_matrix((data,(rows,cols)),shape=(n,len(e)-1))
    expected=np.maximum(0.,np.minimum(e[1:],start+n*w)-np.maximum(e[:-1],start))
    if not np.allclose(np.asarray(z.sum(axis=0)).ravel(),expected,rtol=1e-10,atol=1e-5):raise ValueError("Source axis area lost")
    return z

def aggregate(values,Y,X):
    v=np.ma.asarray(values,dtype=float);good=~np.ma.getmaskarray(v)&np.isfinite(v.data)
    if np.any(good&((v.data<0)|(v.data>5000))):raise ValueError("Invalid biomass values or scale/mask")
    def red(z):return (X@(Y@z).T).T
    val=np.where(good,v.data,0.);a=red(good.astype(float));s=red(val);q=red(val**2)
    dx=np.asarray(X.sum(axis=0)).ravel();dy=np.asarray(Y.sum(axis=0)).ravel()
    area=float(dy@good.astype(float)@dx);mass=float(dy@val@dx)
    if not math.isclose(float(a.sum()),area,rel_tol=1e-10) or not math.isclose(float(s.sum()),mass,rel_tol=1e-10,abs_tol=1e-6):raise ValueError("Area or density integral not conserved")
    return a,s,q,{"valid_pixels":int(good.sum()),"valid_zero_pixels":int((good&(v.data==0)).sum()),"valid_area_m2":area,"density_integral":mass,"max_Mg_ha":float(np.max(v.data[good]))}

def rollup(v):return v.reshape(NY//2,2,NX//2,2).sum(axis=(1,3))
def records(a,s,q,w):
    out=[]
    for j,i in zip(*np.where(a>0)):
        n=float(a[j,i]);mu=float(s[j,i]/n);sd=math.sqrt(max(0.,float(q[j,i]/n-mu*mu)))
        if n>w*w*1.000001:raise ValueError("Impossible valid area")
        out.append({"source_equal_area_grid_ix":int(i),"source_equal_area_grid_iy":int(j),"valid_area_m2":n,"agb_Mg_ha":mu,"subpixel_sd_Mg_ha":sd})
    return out

def read_source(path):
    from netCDF4 import Dataset,num2date
    tf=Transformer.from_crs(4326,6933,always_xy=True);out={};qa={}
    with Dataset(path) as d:
        if not {"agb","lat","lon","time"}<=set(d.variables):raise ValueError("Expected source variables missing")
        v=d["agb"];dims=list(v.dimensions);unit=str(getattr(v,"units",""))
        if set(dims)!={"time","lat","lon"} or len(dims)!=3:raise ValueError("Unexpected AGB dimensions")
        u=unit.casefold().replace(" ","").replace("·","").replace("**","").replace("^","").replace("−","-")
        if u not in ("mg/ha","mgha-1","t/ha","tha-1","tonnes/ha","tonnesha-1"):raise ValueError("Unknown published biomass units: "+unit)
        def bounds(name):
            b=getattr(d[name],"bounds",name+"_bnds")
            if b not in d.variables:return None
            z=np.asarray(d[b][:],float);return z.T if z.shape==(2,len(d[name])) else z
        li,le=edges(d["lon"][:],bounds("lon"),-180,180);ai,ae=edges(d["lat"][:],bounds("lat"),-90,90)
        if len(li)<3500 or len(ai)<1400 or le[0]>-179.9 or le[-1]<179.9 or ae[0]>-55 or ae[-1]<80:raise ValueError("Incomplete global extent")
        tv=d["time"];years=[int(t.year) for t in num2date(tv[:],tv.units,calendar=getattr(tv,"calendar","standard"))]
        if len(set(years))!=len(years) or not set(YEARS)<=set(years):raise ValueError("Expected publisher epochs absent")
        xe,_=tf.transform(le,np.zeros(len(le)));_,ye=tf.transform(np.zeros(len(ae)),ae)
        X=overlap_matrix(xe,X_START,NX);Y=overlap_matrix(ye,Y_START,NY)
        meta={"shape":list(v.shape),"dimensions":dims,"units":unit,"years":years,"lon_extent":[float(le[0]),float(le[-1])],"lat_extent":[float(ae[0]),float(ae[-1])],"global_attributes":{k:str(d.getncattr(k))[:400] for k in ("title","product_version","version","source","institution") if k in d.ncattrs()},"agb_attributes":{k:str(v.getncattr(k))[:300] for k in ("long_name","_FillValue","missing_value","scale_factor","add_offset","valid_range") if k in v.ncattrs()}}
        print(json.dumps({"ESA_source_schema":meta}),flush=True)
        for year in YEARS:
            z=v[tuple(years.index(year) if k=="time" else slice(None) for k in dims)]
            if [k for k in dims if k!="time"]==["lon","lat"]:z=z.T
            a,s,q,check=aggregate(z[np.ix_(ai,li)],Y,X)
            if check["valid_pixels"]<500000 or check["max_Mg_ha"]<=0:raise ValueError("Unexpectedly sparse or zero world biomass")
            a2,s2,q2=map(rollup,(a,s,q))
            if any(not math.isclose(float(x.sum()),float(y.sum()),rel_tol=1e-12) for x,y in ((a,a2),(s,s2),(q,q2))):raise ValueError("2x2 rollup failed")
            out[str(year)]={"100000":records(a,s,q,W),"200000":records(a2,s2,q2,2*W)}
            qa[str(year)]={**check,"n100":len(out[str(year)]["100000"]),"n200":len(out[str(year)]["200000"])}
            print(json.dumps({"AGB_year":year,**qa[str(year)]}),flush=True)
    return out,meta,qa

def main():
    p=argparse.ArgumentParser();p.add_argument("--nc",type=Path,required=True);p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    if a.nc.stat().st_size<1000000:raise ValueError("Not full 10km source file")
    with a.nc.open("rb") as f:sha=hashlib.file_digest(f,"sha256").hexdigest()
    grids,meta,qa=read_source(a.nc)
    z={"version":"railstatus-ESA-CCI-v7-10km-woody-agb-2020-2024-v1","DOI":DOI,"source_url":URL,"source_sha256":sha,"source_bytes":a.nc.stat().st_size,"source_schema":meta,"qa":qa,"grids":grids,"limits":["woody aboveground dry biomass density, not all plant biomass or NPP","source pixels assumed piecewise constant","valid source pixels only; no missing to zero","spatial SD not retrieval uncertainty","no total mass claim"]}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(z,allow_nan=False))
if __name__=="__main__":main()
