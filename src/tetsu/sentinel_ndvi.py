"""Source-controlled Sentinel-2 L2A reflectance and NDVI safety checks.

Explicitly preserve BOTH interpretations when Earth Search COG metadata
reports earthsearch:boa_offset_applied=true and raster:bands.offset=-0.1.
No ecological effect or multi-year harmonization can be asserted unless
that data-provider contradiction is independently resolved.
"""
from __future__ import annotations
import math
import statistics

CLEAR_SCL_PRIMARY={4,5}
CLEAR_SCL_SENSITIVITY={4,5,7}
BASELINE_CHANGE_YEAR=2022


def calibrated_reflectance_values(raw_dn, *, subtract_dn=0):
    """Return scaled reflectance from raw nonzero digital number.

    subtract_dn is 0 for 'provider COG already BOA-corrected' and
    1000 for 'raster:bands offset -0.1' interpretation. These are
    alternatives for validation, not equivalent physical corrections.
    """
    try:
        dn=float(raw_dn)
    except (TypeError,ValueError) as exc:
        raise ValueError("Not a numeric S2 DN") from exc
    if not math.isfinite(dn) or dn<=0 or dn>32767:
        return None
    if subtract_dn not in (0,1000):
        raise ValueError("Only documented 0 and -1000 calibration scenarios")
    scaled=(dn-subtract_dn)/10000.
    # Negative reflectances can occur after offset subtraction but they
    # do not support a stable NDVI where the pair sum is tiny.
    return scaled if -0.1 <= scaled <=1.5 else None


def calculate_ndvi(red_dn,nir_dn,subtract_dn=0):
    red=calibrated_reflectance_values(red_dn,subtract_dn=subtract_dn)
    nir=calibrated_reflectance_values(nir_dn,subtract_dn=subtract_dn)
    if red is None or nir is None:
        return None
    if min(red,nir)<0 or red+nir<=0.03:
        return None
    value=(nir-red)/(nir+red)
    if not math.isfinite(value) or not -1<=value<=1:
        return None
    return float(value)


def metadata_ambiguous(item):
    """Does provider say BOA already applied but band metadata says -0.1?"""
    props=item.get("properties") or {}
    flag=props.get("earthsearch:boa_offset_applied")
    offsets=[]
    for band in ("red","nir"):
        r=(item.get("assets") or {}).get(band) or {}
        entries=r.get("raster:bands") or []
        if not entries:
            return True
        offsets.append(entries[0].get("offset"))
    if flag not in (True,False):
        return True
    if abs(float(offsets[0] or 0)-float(offsets[1] or 0))>1e-7:
        return True
    return flag is True and any(abs(float(v or 0))>=0.09 for v in offsets)


def ndvi_branch_values(red_dn,nir_dn,item):
    props=item.get("properties") or {}
    baseline=props.get("s2:processing_baseline")
    flag=props.get("earthsearch:boa_offset_applied")
    if baseline is None or flag not in (True,False):
        return {"provider_flag_branch":None,
                "raster_offset_branch":None,
                "calibration_disagreement":True}
    raw=calculate_ndvi(red_dn,nir_dn,0)
    adjusted=calculate_ndvi(red_dn,nir_dn,1000)
    props_bands=(item.get("assets") or {}).get("red",{}).get("raster:bands",[])
    if not props_bands:
        return {"provider_flag_branch":None,
                "raster_offset_branch":None,
                "calibration_disagreement":True}
    offset=float(props_bands[0].get("offset") or 0)
    off=1000 if offset<-.09 else 0
    based_on_property=(raw if flag else
                      adjusted if float(baseline)>=4.0 else raw)
    based_on_band=adjusted if off==1000 else raw
    return {
        "provider_flag_branch":based_on_property,
        "raster_offset_branch":based_on_band,
        "calibration_disagreement":metadata_ambiguous(item)
    }


def summary(values):
    values=[float(v) for v in values if v is not None and math.isfinite(v)]
    if not values:
        return {"n":0,"median":None,"mean":None}
    return {
        "n":len(values),
        "median":round(statistics.median(values),6),
        "mean":round(statistics.mean(values),6)
    }


def paired_annual_ndvi(pairs, ndvi_by_year, year, branch):
    """Near-minus-far paired-site observation in ONE exposed rail corridor.

    A pair is two fixed projected points at the same along-rail transect
    with different offsets. Even many pairs are spatially correlated, not
    independent transport/railway treatments.
    """
    values=[]
    near=[]
    far=[]
    for near_id,far_id in pairs:
        a=ndvi_by_year.get(year,{}).get(near_id,{}).get(branch)
        b=ndvi_by_year.get(year,{}).get(far_id,{}).get(branch)
        if a is None or b is None:
            continue
        values.append(a-b);near.append(a);far.append(b)
    return {
        "near_minus_control_ndvi":summary(values),
        "near_ndvi":summary(near),
        "control_ndvi":summary(far)
    }


def prepost_paired_change(pairs,ndvi_by_year,after_year,branch,base_year=2017):
    changes=[]
    for a,b in pairs:
        base_near=ndvi_by_year.get(base_year,{}).get(a,{}).get(branch)
        base_far=ndvi_by_year.get(base_year,{}).get(b,{}).get(branch)
        new_near=ndvi_by_year.get(after_year,{}).get(a,{}).get(branch)
        new_far=ndvi_by_year.get(after_year,{}).get(b,{}).get(branch)
        if any(v is None for v in (base_near,base_far,new_near,new_far)):
            continue
        changes.append((new_near-base_near)-(new_far-base_far))
    return {
        "year":after_year,"baseline_year":base_year,
        "paired_near_minus_control_delta_NDVI":summary(changes),
        "not_causal":True,
    }
