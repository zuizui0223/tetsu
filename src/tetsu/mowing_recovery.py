"""Mowing-aligned projected vegetation/flower ROI from verified fixed plot photos.

Requires externally confirmed cutting events and human-reviewed masks.
This is a *descriptive* within-plot pre/post change, not a causal treatment
effect, not physical height or biomass, and not genetic adaptation.
"""
from __future__ import annotations

import csv
import datetime as dt
import math
import statistics
from collections import defaultdict


def _date(value):
    try:
        return dt.date.fromisoformat(str(value))
    except ValueError as e:
        raise ValueError(f"invalid ISO date {value!r}") from e


def _number(value, name):
    try:
        x=float(value)
    except (ValueError,TypeError) as e:
        raise ValueError(f"not numeric {name}={value!r}") from e
    if not math.isfinite(x) or x<0 or x>1:
        raise ValueError(f"{name} must be a fraction in [0,1]")
    return x


def summarize_mowing_roi(roi_rows, event_rows, min_before=2, min_after=2):
    """One record per plot/event, including failure reasons.

    Fixed windows: -28:-1 days and +7:+42 days relative to independently
    documented mowing event. Multiple events overlapping an assessment
    window invalidate the comparison. No interpolated/imaginary photographs.
    """
    if min_before<1 or min_after<1:
        raise ValueError("minimum repeated photos must be >= 1")
    eligible=defaultdict(list)
    rejected=defaultdict(int)
    seen=set()
    for r in roi_rows:
        uid=(r.get("image_id") or "").strip()
        if not uid or uid in seen:
            raise ValueError("missing/duplicate photo image_id")
        seen.add(uid)
        if r.get("quality_status")!="patch_pass" or str(r.get("valid_for_fixed_plot_flower_fraction")).lower() not in ("true","1"):
            rejected["not_verified_fixed_conspecific_patch"]+=1
            continue
        plot_id=(r.get("plot_id") or "").strip()
        site=(r.get("site_id") or "").strip()
        taxon=(r.get("taxon") or "").strip()
        if not all((plot_id,site,taxon)):
            rejected["incomplete_plot_identity"]+=1
            continue
        date=_date(r["photo_date"])
        cover=_number(r.get("plant_frame_fraction"),"plant_frame_fraction")
        floral=_number(r.get("flower_area_over_plant_area"),"flower_area_over_plant_area")
        eligible[(site,plot_id,taxon)].append((date,cover,floral))
    cuts=defaultdict(list)
    for e in event_rows:
        site=str(e.get("site_id") or "").strip()
        kind=str(e.get("event_type") or "").strip()
        evidence=str(e.get("evidence") or "").strip()
        if not site or not kind or not evidence:
            raise ValueError("management log missing site/event_type/evidence")
        if kind!="cutting":
            continue
        cuts[site].append(_date(e["event_date"]))
    for site in cuts:
        cuts[site]=sorted(set(cuts[site]))
    results=[]
    for (site,plot,taxon),obs in sorted(eligible.items()):
        for event in cuts.get(site,[]):
            # Other known cuts in either window make the single-event
            # calculation uninterpretable without treatment histories.
            overlap=[x for x in cuts[site] if x!=event and -28 <= (x-event).days <=42]
            before=[r for r in obs if -28 <= (r[0]-event).days <= -1]
            after=[r for r in obs if 7 <= (r[0]-event).days <= 42]
            reason=""
            if overlap:
                reason="overlapping_cutting_events"
            elif len(before)<min_before or len(after)<min_after:
                reason="insufficient_repeated_fixed_plot_photos"
            entry={
                "site_id":site,"plot_id":plot,"taxon":taxon,
                "cut_date":event.isoformat(),
                "n_pre":len(before),"n_recovery":len(after),
                "eligibility":"pass" if not reason else "restricted",
                "restriction_reason":reason
            }
            if not reason:
                med=lambda z,k: statistics.median(x[k] for x in z)
                for prefix,column in (("vegetation_frame_occupancy",1),
                                      ("floral_area_per_vegetation",2)):
                    b=med(before,column); a=med(after,column)
                    entry[prefix+"_pre"]=round(b,6)
                    entry[prefix+"_recovery"]=round(a,6)
                    entry["delta_"+prefix]=round(a-b,6)
            results.append(entry)
    return results,dict(rejected)


def load_csv(path):
    with open(path,newline="",encoding="utf-8") as f:
        return list(csv.DictReader(f))
