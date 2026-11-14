#!/usr/bin/env python3
"""v10.8 REAL Czech 2024 station area-mix vs alien flora, 39 source stations.

Original Kutlvašr et al. (2024) Preslia 96:247-266, DOI
10.23855/preslia.2024.247, publisher-corrected Appendix S2 (erratum).
Tests whether the proportion of station track-zone area is associated
with station alien-flora fraction beyond station total area, station
building-zone proportion and five original railway stretches.

Nonrandom one-country cross-section. An exploratory station habitat-mix
association, NOT a railway-vs-road or moving-seed causal estimate.
"""
from __future__ import annotations
import argparse, hashlib, json, math, random
from collections import Counter
from pathlib import Path
from statistics import mean, median

from audit_Czech_2024_real_railway_appendices_v107 import (
    audit_workbook, audit_presence_matrix_header_s1, SOURCE_DOI,
)

EXPECTED_S1_SHA="df6bc0b949e7f6a5cb57755dc8f9c86187223f66ddfbef959cf0e03ef7c46b06"
EXPECTED_S2_SHA="00c7cfdf5047854ba99d47dc8ec8c55b1b7df72731ce647c22b286ddef24a5a2"
RNG_SEED=20261114
PERM_DRAWS=9999
BOOT_DRAWS=4000


def solve_linear(a,b):
    """Solve full-rank small system with pivoting, without external numerics."""
    n=len(b)
    if n==0 or len(a)!=n or any(len(row)!=n for row in a):
        raise ValueError("OLS normal equation shape mismatch")
    m=[list(map(float,a[i]))+[float(b[i])] for i in range(n)]
    for j in range(n):
        k=max(range(j,n),key=lambda i:abs(m[i][j]))
        if abs(m[k][j])<1e-10:
            raise ValueError("Station model design singular")
        m[j],m[k]=m[k],m[j]
        v=m[j][j]
        for c in range(j,n+1): m[j][c]/=v
        for i in range(n):
            if i==j:continue
            scale=m[i][j]
            for c in range(j,n+1):m[i][c]-=scale*m[j][c]
    return [m[i][n] for i in range(n)]


def ols(y,mat):
    if len(y)!=len(mat) or not mat or len(mat)<=len(mat[0]):
        raise ValueError("Not enough independent stations for area-adjusted model")
    cols=len(mat[0])
    if any(len(row)!=cols for row in mat):
        raise ValueError("Design has mismatched columns")
    xtx=[[sum(r[i]*r[j] for r in mat) for j in range(cols)] for i in range(cols)]
    xty=[sum(r[j]*v for r,v in zip(mat,y)) for j in range(cols)]
    beta=solve_linear(xtx,xty)
    residual=[yi-sum(b*xi for b,xi in zip(beta,row)) for yi,row in zip(y,mat)]
    return beta,residual


def design(rows,extra=()):
    groups=sorted(set(r["stretch"] for r in rows))
    if len(groups)<2:raise ValueError("Only one independent original source railway stretch")
    return [[1.0,math.log(r["area_ha"]),r["building_share"]]
            +[float(r["stretch"]==group) for group in groups[1:]]
            +[r[k] for k in extra] for r in rows]


def partial_slope(rows,outcome,extra_controls=()):
    key="track_share"
    base=design(rows,extra_controls)
    target=[r[outcome] for r in rows]
    x=[r[key] for r in rows]
    _,ry=ols(target,base)
    _,rx=ols(x,base)
    ssx=sum(z*z for z in rx)
    if ssx<1e-10:raise ValueError("Original station track area fraction lacks independent variation")
    beta=sum(x*y for x,y in zip(rx,ry))/ssx
    return beta,rx,ry


def blocked_residual_permutation(rows,outcome,draws=PERM_DRAWS):
    """Exploratory within-stretch conditional residual permutation.

    Conditions on original published total station area (log) and
    building area fraction; approximate permutation p, no causal null.
    """
    slope,rx,ry=partial_slope(rows,outcome)
    rng=random.Random(RNG_SEED)
    groups={}
    for i,r in enumerate(rows):
        groups.setdefault(r["stretch"],[]).append(i)
    ssx=sum(z*z for z in rx)
    hit=0
    for _ in range(draws):
        vals=list(ry)
        for inds in groups.values():
            source=[ry[i] for i in inds]
            rng.shuffle(source)
            for idx,val in zip(inds,source):
                vals[idx]=val
        b=sum(x*y for x,y in zip(rx,vals))/ssx
        if abs(b)+1e-14>=abs(slope):hit+=1
    return (hit+1)/(draws+1)


def pct(x,p):
    s=sorted(x)
    if not s:return None
    z=(len(s)-1)*p
    i=int(math.floor(z))
    j=int(math.ceil(z))
    return s[i] if i==j else s[i]+(s[j]-s[i])*(z-i)


def stratified_station_bootstrap(rows,outcome,draws=BOOT_DRAWS):
    rng=random.Random(RNG_SEED)
    by={}
    for r in rows:by.setdefault(r["stretch"],[]).append(r)
    boot=[]
    for _ in range(draws):
        sample=[]
        for group in sorted(by):
            a=by[group]
            sample.extend(rng.choice(a) for _ in a)
        try:
            b,_,_=partial_slope(sample,outcome)
        except ValueError:
            continue
        if math.isfinite(b):boot.append(b)
    if len(boot)<draws*.90:
        raise ValueError("Insufficient bootstrap full rank within five original source railway stretches")
    return [pct(boot,.025),pct(boot,.975)],len(boot)


def effect_diag(rows,outcome):
    b,_,_=partial_slope(rows,outcome)
    p=blocked_residual_permutation(rows,outcome)
    ci,n=stratified_station_bootstrap(rows,outcome)
    groups=sorted(set(r["stretch"] for r in rows))
    loo=[]
    for section in groups:
        sub=[r for r in rows if r["stretch"]!=section]
        loo.append({"omitted_original_source_stretch":section,
                    "partial_beta_per_unit_track_share":round(
                        partial_slope(sub,outcome)[0],7)})
    return {"n_original_independent_station_locations":len(rows),
            "n_original_railway_stretches":len(groups),
            "partial_beta_per_ten_percentage_points_track_area_share":round(b*.1,7),
            "within_stretch_residual_permutation_abs_beta_p_exploratory":round(p,6),
            "source_stretch_stratified_station_bootstrap95_per_10pp":[
                round(v*.1,7) for v in ci],
            "n_valid_stratified_bootstrap_draws":n,
            "leave_one_stretch_out":[
                {**r,"partial_beta_per_10pp_track_share":round(
                    r["partial_beta_per_unit_track_share"]*.1,7)} for r in loo],
            "not_a_train_seed_dispersal_or_mowing_or_slipstream_causal_coefficient":True}


def recent_station_data(s1,s2):
    from openpyxl import load_workbook
    book=load_workbook(s2,read_only=True,data_only=True)
    sheet=next(iter(book.worksheets))
    records={}
    for data in list(sheet.values)[5:]:
        if not data or len(data)<11 or not data[1] or not str(data[1]).strip():
            continue
        code=str(data[1]).strip()
        if code in records:raise ValueError("Duplicated original station code")
        stretch=str(data[2]).strip()
        area=float(data[10])
        building=float(data[7])
        surround=float(data[8])
        track=float(data[9])
        if not (area>0 and building>=0 and surround>=0 and track>=0):
            raise ValueError("Original station has unknown/negative area")
        if abs(building+surround+track-area)>0.005:
            raise ValueError("Original source three station zones do not sum to S2 station total")
        recorded=int(data[3]);alien=int(data[4]);invasive=int(data[5]);native=int(data[6])
        if not (0<=invasive<=alien<=recorded and 0<=native<=recorded):
            raise ValueError("Source plant count invalid or S2 columns misaligned")
        records[code]={
            "station":code,"stretch":stretch,"area_ha":area,
            "track_share":track/area,"building_share":building/area,
            "surroundings_share":surround/area,"alien_n":alien,
            "invasive_n":invasive,"taxa_n":recorded,"native_n":native,
            "alien_fraction":alien/recorded,
            "invasive_fraction_of_aliens":invasive/alien if alien else 0,
        }
    book.close()
    if len(records)!=39 or len(set(r["stretch"] for r in records.values()))!=5:
        raise ValueError("Original erratum 39 stations/5 stretches not recovered")
    return records


def audit_and_run(s1,s2):
    check1=audit_workbook(s1);check2=audit_workbook(s2)
    if (check1["original_source_sha256"]!=EXPECTED_S1_SHA or
        check2["original_source_sha256"]!=EXPECTED_S2_SHA):
        raise ValueError("Publisher Czech railway original S1/S2 source SHA mismatch")
    m=audit_presence_matrix_header_s1(s1,s2)
    if (m["n_recent_station_columns"]!=39 or
        m["n_historical_station_columns"]!=33 or
        m["original_presence_marker_by_rec_hist"].get("rec__x")!=6941):
        raise ValueError("Original published station rec/hist taxon matrix changed")
    rec=recent_station_data(s1,s2)
    station_codes=set(rec)
    if station_codes!=set(m["station_s1_rec_any_marker_count_vs_s2_richness"]):
        raise ValueError("Original station flora vs station data IDs mismatch")
    mismatch_taxa=[s for s in rec
        if m["station_s1_rec_any_marker_count_vs_s2_richness"][s]["recent_marker_taxa"]
        !=rec[s]["taxa_n"]]
    mismatch_aliens=[s for s in rec
        if m["station_s1_rec_explicit_arch_neo_vs_s2_aliens"][s]["recent_arch_neo_tagged"]
        !=rec[s]["alien_n"]]
    rows=list(rec.values())
    primary=effect_diag(rows,"alien_fraction")
    secondary=effect_diag(rows,"invasive_fraction_of_aliens")
    robust=[r for r in rows if r["station"] not in set(mismatch_taxa+mismatch_aliens)]
    robustness={
      "n_original_stations_without_discrepant_species_or_alien_counts":len(robust),
      "alien_fraction_partial_track_share_beta_per_10pp":round(
        partial_slope(robust,"alien_fraction")[0]*.1,7)
        if len(robust)>=20 else None,
    }
    return {
        "complete":True,"version":"v10.8",
        "original_article_DOI":SOURCE_DOI,
        "original_source_S1_sha256":EXPECTED_S1_SHA,
        "original_corrected_S2_erratum_sha256":EXPECTED_S2_SHA,
        "n_recent_stations":len(rows),
        "n_source_rail_stretches":len(set(r["stretch"] for r in rows)),
        "n_recent_species_times_station_x_from_S1":6941,
        "n_recent_station_richness_exact_S1_to_S2":39-len(mismatch_taxa),
        "original_station_richness_disagreement_codes":sorted(mismatch_taxa),
        "n_recent_station_arch_neo_exact_S1_to_S2":39-len(mismatch_aliens),
        "original_station_alien_status_disagreement_codes":sorted(mismatch_aliens),
        "source_area_ha_range":[round(min(r["area_ha"] for r in rows),4),
                                round(max(r["area_ha"] for r in rows),4)],
        "source_track_zone_fraction_range":[round(min(r["track_share"] for r in rows),4),
                                           round(max(r["track_share"] for r in rows),4)],
        "source_aliens_fraction_range":[round(min(r["alien_fraction"] for r in rows),4),
                                      round(max(r["alien_fraction"] for r in rows),4)],
        "descriptive_source_alien_count_mean":round(mean(r["alien_n"] for r in rows),4),
        "primary_outcome_station_aliens_fraction_track_area_share_adjusted":primary,
        "secondary_outcome_station_invasive_share_of_aliens_track_area_share_adjusted":secondary,
        "original_disagreement_omit_robustness":robustness,
        "method":"OLS partial slope track-area share using fixed-effects for five original source stretches, log total area and building-area share, approximate source-stretch blocked residual permutation 9999 and within-stretch station bootstrap 4000; all exploratory noncausal.",
        "critical_limitations":"The study sampled nonrandom 39 freight-corridor stations across a single country with different survey-zone areas and source selection. Track area share is a station habitat composition proxy, NOT train traffic, mowing, direct seed arrival, thermal/wind stress or independent railway-vs-road treatment. S1 and corrected S2 3 station richness and 13 alien-count disagreements retained. Historical records cannot be used as equally sampled station-year absences. Source authors credited, source XLSX not redistributed.",
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--s1",required=True)
    p.add_argument("--s2",required=True)
    p.add_argument("--output",required=True)
    a=p.parse_args()
    out=audit_and_run(a.s1,a.s2)
    fp=Path(a.output)
    fp.parent.mkdir(parents=True,exist_ok=True)
    fp.write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps(out,ensure_ascii=False),flush=True)


if __name__=="__main__":
    main()
