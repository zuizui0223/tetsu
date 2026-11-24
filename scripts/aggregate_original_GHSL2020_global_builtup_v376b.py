#!/usr/bin/env python3
"""v376b genuine EC JRC GHSL 2020 GHS-BUILT-S global full 1km m2 → 100/200km.

Original full ZIP 152527564 bytes and v376a source ETag checked in workflow.
Resampling.sum aggregates source built-up m2 per Mollweide 1km pixel
to original EPSG6933 100km plant/rail equal-area model grid (v369/370),
NEVER a single built-up reading at grid centre.
200km built square metres are sums of original 100km children exactly.
No source original TIFF, raw observer GPS or location-specific records
rehosted; archive aggregated source-built fraction by coarse grid.
"""
import argparse,hashlib,json,math,shutil,zipfile
from pathlib import Path
import numpy as np

EXPECTED_ZIP_BYTES=152527564
SOURCE_URL=("https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL/"
    "GHS_BUILT_S_GLOBE_R2023A/GHS_BUILT_S_E2020_GLOBE_R2023A_54009_1000/"
    "V1-0/GHS_BUILT_S_E2020_GLOBE_R2023A_54009_1000_V1_0.zip")
SOURCE_DOI="10.2905/9F06F36F-4B11-47EC-ABB0-4F8B7B1D72EA"
X_START=-17_400_000
X_END=17_400_000
Y_START=-7_400_000
Y_END=7_400_000
WIDTH=100_000
ROWS=(Y_END-Y_START)//WIDTH
COLS=(X_END-X_START)//WIDTH
NODATA=-999.
BANDS=(100_000,200_000)
KNOWN_SOURCE_NAME="GHS_BUILT_S_E2020_GLOBE_R2023A_54009_1000_V1_0"

def fingerprint(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda:f.read(4<<20),b""):h.update(chunk)
    return h.hexdigest()

def source_tif(zip_file,workspace):
    zpath=Path(zip_file)
    if zpath.stat().st_size!=EXPECTED_ZIP_BYTES:
        raise ValueError("Not original JRC GHSL source ZIP Content-Length")
    sha=fingerprint(zpath)
    with zipfile.ZipFile(zpath) as z:
        if z.testzip() is not None:raise ValueError("Original JRC ZIP CRC failed")
        selected=[f for f in z.infolist()
                  if not f.is_dir() and f.filename.lower().endswith(".tif")]
        if len(selected)!=1:
            raise ValueError("Original JRC full globe 2020 1km ZIP must contain exactly one GeoTIFF")
        original=selected[0]
        if KNOWN_SOURCE_NAME not in original.filename:
            raise ValueError("Original JRC TIFF filename does not belong to source epoch/projection")
        destination=Path(workspace)/"source_original_GHSL2020_builtup_1km.tif"
        destination.parent.mkdir(parents=True,exist_ok=True)
        with z.open(original) as f, destination.open("wb") as out:
            shutil.copyfileobj(f,out,4<<20)
    return destination,{"global_GHSL_2020_publisher_full_ZIP_bytes":zpath.stat().st_size,
       "source_JRC_GHSL2020_original_fullZIP_SHA256":sha,
       "source_original_TIFF_name":original.filename,
       "source_original_TIFF_size_uncompressed_bytes":original.file_size}

def grid_sums(srcfile):
    import rasterio
    from rasterio.warp import reproject,Resampling
    from rasterio.transform import from_origin
    dest=np.full((ROWS,COLS),NODATA,dtype=np.float64)
    with rasterio.open(srcfile) as src:
        desc={"source_crs":str(src.crs),
              "source_original_raster_dimensions_rows_cols":list(src.shape),
              "source_x_resolution_m":abs(src.transform.a),
              "source_y_resolution_m":abs(src.transform.e),
              "source_nodata":src.nodata,"source_dtype":src.dtypes[0],
              "source_band_count":src.count}
        print(json.dumps({"source_original_GHSL_raster_schema":desc}),flush=True)
        if not src.crs or src.count!=1:
            raise ValueError("JRC global builtup 1km original raster unrecognized")
        if abs(src.transform.a-1000)>1 or abs(src.transform.e+1000)>1:
            raise ValueError("Original GHSL 1km Mollweide original resolution changed")
        if (not str(src.crs).lower().find("moll")>=0
            and "54009" not in str(src.crs)):
            raise ValueError("Original GHSL 2020 source is not Mollweide 54009")
        if src.nodata is None or src.nodata>=0:
            raise ValueError("Original GHSL 2020 raster expected negative nodata that is not built-up zero")
        reproject(
           source=rasterio.band(src,1),
           destination=dest,
           src_transform=src.transform,src_crs=src.crs,
           src_nodata=src.nodata,
           dst_transform=from_origin(X_START,Y_END,WIDTH,WIDTH),
           dst_crs="EPSG:6933",dst_nodata=NODATA,
           resampling=Resampling.sum,
           num_threads=2,
           warp_mem_limit=256,
           init_dest_nodata=True)
    if dest.shape!=(ROWS,COLS):raise AssertionError("Wrong full-world 100km EPSG6933 grid")
    valid=np.isfinite(dest)&(dest>=0)
    fractions=dest[valid]/(WIDTH*WIDTH)
    if len(fractions)<10000 or np.max(fractions)>1.005 or np.min(fractions)<0:
        raise ValueError("Original global builtup source pixel units or projected sum physically implausible")
    summary={"source_fullworld_100km_tile_dimensions":[ROWS,COLS],
             "source_100km_grid_cells_with_valid_GHSL_2020":int(valid.sum()),
             "source_100km_builtup_m2_total_all_valid_cells":float(np.sum(dest[valid])),
             "source_100km_fraction_min":float(np.min(fractions)),
             "source_100km_fraction_mean_unweighted_over_valid":float(np.mean(fractions)),
             "source_100km_fraction_max":float(np.max(fractions)),
             "source_100km_zero_builtup_cells":int(np.count_nonzero(fractions==0)),
             "original_source_GHSL_unit_m2_builtup_per_1km_pixel":True,
             "GHSL2020_source_nodata_not_equivalent_to_builtup_zero":True}
    return dest,valid,{**desc,**summary}

def rowwise_coarse_grid(dest,valid,source_m):
    if source_m not in BANDS:
        raise ValueError("Unsupported exact original rail-plant equal-area grid")
    factor=source_m//WIDTH
    if ROWS%factor or COLS%factor:
        raise ValueError("100km fixed global source grid cannot roll up to 200km")
    out=[]
    # Internal raster row zero is NORTH, while geographic grid iy=0 is SOUTH.
    h=ROWS//factor
    w=COLS//factor
    for new_iy in range(h):
        for new_ix in range(w):
            yy=(h-1-new_iy)*factor
            xx=new_ix*factor
            area=dest[yy:yy+factor,xx:xx+factor]
            ok=valid[yy:yy+factor,xx:xx+factor]
            if not np.any(ok):
                continue
            n=float(np.sum(area[ok]))
            frac=n/(source_m*source_m)
            if not 0<=frac<=1.005:
                raise ValueError("GHSL global 1km sum to large original equal area grid exceeded 100%")
            out.append({"source_equal_area_grid_ix":new_ix,
                        "source_equal_area_grid_iy":new_iy,
                        "source_grid_width_m":source_m,
                        "source_GHSL2020_total_builtup_m2_in_cell":round(n,1),
                        "source_GHSL2020_builtup_fraction_in_original_cell":float(frac),
                        "source_GHSL2020_valid_100km_children":int(ok.sum()),
                        "source_original_GHSL2020_builtup_m2_not_percent_observations":True})
    return out

def analyze(zipfile_source,workspace):
    original,provenance=source_tif(zipfile_source,workspace)
    try:
        dest,valid,qa=grid_sums(original)
    finally:
        original.unlink(missing_ok=True)
    a=rowwise_coarse_grid(dest,valid,100000)
    b=rowwise_coarse_grid(dest,valid,200000)
    sum100=sum(x["source_GHSL2020_total_builtup_m2_in_cell"] for x in a)
    sum200=sum(x["source_GHSL2020_total_builtup_m2_in_cell"] for x in b)
    if abs(sum100-sum200)/max(1,abs(sum100))>1e-10:
        raise AssertionError("100km and 200km source builtup m2 totals not conserved")
    return {"version":"v376b-publisher-original-GHSL2020-full-global-builtup-m2-100km-200km",
        "source_European_Commission_JRC_data_DOI":SOURCE_DOI,
        "source_original_GHSL2020_URL":SOURCE_URL,
        **provenance,
        "original_GHSL_raster_source_QA":qa,
        "original_global_GHSL2020_built_m2_by_100km_200km_cell_width":{
            "100000":a,"200000":b},
        "aggregate_source_QA":{
          "valid_GHSL_cells_100km":len(a),
          "valid_GHSL_cells_200km":len(b),
          "aggregate_total_built_m2_100km":sum100,
          "aggregate_total_built_m2_200km":sum200,
          "total_m2_conserved_in_grid_rollup":True},
        "scientific_boundaries":[
          "GHSL built-up land surface 2020 is source separate from 2025 plant observations and 2026 rail GIS",
          "One square km original source cell contains builtup m2, not total urban human footprint or railkm",
          "100km/200km builtup fraction derived by spatially weighted sum not center-only pixel",
          "No plant functional traits, railway effects, species niche or source rail seed dispersal evaluated here",
          "Published output contains 100km/200km coarse builtup area aggregate only, no source full original raster or private GPS"]}
def main():
    p=argparse.ArgumentParser()
    p.add_argument("--original-global-GHSL-ZIP",required=True)
    p.add_argument("--workspace",required=True)
    p.add_argument("--output",required=True)
    a=p.parse_args()
    report=analyze(a.original_global_GHSL_ZIP,a.workspace)
    out=Path(a.output);out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding="utf8")
    print(json.dumps({k:v for k,v in report.items()
        if k!="original_global_GHSL2020_built_m2_by_100km_200km_cell_width"},ensure_ascii=False),flush=True)
if __name__=="__main__":main()
