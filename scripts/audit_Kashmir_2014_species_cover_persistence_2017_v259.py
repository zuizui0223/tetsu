#!/usr/bin/env python3
"""v25.9: Does original 2014 cover abundance rank predict whether a taxon
recorded at a Kashmir railway plot in 2014 is still recorded in 2017?

Author-original Dryad Matrix_plot.xlsx (SHA-gated by v25.4), same 27 true
spatial sites with 4 100m2 original zones in both years. Early cover is
an *ordinal numerical* within-plot predictor; units need not be percent.

AUC within site = P[2014 cover retained taxon > cover of taxon not
recorded 2017] + 1/2 P[tie]. Preserves the taxon risk set and source
site persistence/loss counts when labels are permuted. Does not treat
2014 unobserved taxa as true absences or actual seed colonizations.
Does not establish mortality vs imperfect detection, herbicide, mowing,
train-induced wind, site resource limitation or evolutionary adaptation.

No original species names, plot×taxon cover or site geography exported.
"""
from __future__ import annotations
import argparse
import json
import math
import random
from collections import defaultdict
from pathlib import Path
from statistics import mean,median

from analyze_Kashmir_same_sites_reference_year_v254 import (
    parse_xlsx_matrix, matching_sites
)

ZONES=("Track","Embankment","Interior Plot 1","Interior Plot 2")
GROUPS=("station","between")
SEED=259118
DRAWS=999


def early_numeric_cover(xlsx: Path, authentic_presence:dict) -> tuple[dict,dict]:
    """Re-read publisher source as numbers only after authentic SHA gate."""
    import openpyxl
    wb=openpyxl.load_workbook(xlsx,read_only=True,data_only=True)
    sheet=wb.worksheets[0]
    iterator=iter(sheet.values)
    header=next(iterator)
    if len(header)!=129 or len(header[6:])!=123:
        raise ValueError("Original Kashmir original 6 source columns + 123 taxa drift")
    numeric={}
    missing=0
    for row in iterator:
        if len(row)!=129:
            raise ValueError("Original source workbook row must preserve 129 columns")
        site=" ".join(str(row[0]).split())
        year=int(row[4])
        zone=" ".join(str(row[5]).split())
        key=(site,year,zone)
        if key in numeric:
            raise ValueError("Duplicate field site/year/plot")
        values={}
        complete=True
        for col,v in enumerate(row[6:]):
            if v is None:
                missing+=1
                complete=False
                continue
            if not isinstance(v,(float,int)) or isinstance(v,bool) or not math.isfinite(v) or v<0:
                raise ValueError("Malformed original nonnegative source cover code")
            if v>0:values[col]=float(v)
        if key not in authentic_presence:
            raise ValueError("Original numeric workbook extra site/year/zone")
        published=authentic_presence[key]
        if complete:
            if published is None or set(values)!=published:
                raise AssertionError("Numeric original source taxa do not reproduce strict source presence")
            numeric[key]=values
        else:
            if published is not None:
                raise ValueError("Previously flagged original missing cell discrepancy")
            numeric[key]=None
    if len(numeric)!=248 or missing!=4 or set(numeric)!=set(authentic_presence):
        raise ValueError("Kashmir source plot geometry/numeric cover missingness changed")
    return numeric, {"source_numeric_plots":len(numeric),"source_missing_original_cover_cells":missing,
                     "original_source_verified_presence_sets_exactly_reproduced":True,
                     "source_numeric_cover_interpreted_only_as_within_plot_rank":True}


def auc_of_early_covers(retained:list[float],lost:list[float]) -> float:
    if not retained or not lost:
        raise ValueError("Within-plot AUC requires both persistence and lost taxa")
    wins=0.
    for x in retained:
        for y in lost:
            wins+=1. if x>y else 0.5 if x==y else 0.
    return wins/(len(retained)*len(lost))


def site_early_persistence(early:dict[int,float],late:set[int]) -> dict:
    current=set(early)
    keep=current&late
    lost=current-late
    result={
        "n_original_2014_taxa":len(current),
        "n_persistent_recorded_2017":len(keep),
        "n_not_recorded_2017":len(lost),
        "n_newly_recorded_2017":len(late-current),
        "eligible_auc":bool(keep and lost),
    }
    if result["eligible_auc"]:
        retained=[early[t] for t in sorted(keep)]
        vanished=[early[t] for t in sorted(lost)]
        result["source_2014_cover_AUC"]=auc_of_early_covers(retained,vanished)
        result["mean_2014_cover_retained_numeric"]=mean(retained)
        result["mean_2014_cover_lost_numeric"]=mean(vanished)
        result["raw_ranks_all"]=sorted(early.values())
        result["n_positive_retained_labels"]=len(keep)
    return result


def quantile(a:list[float],p:float)->float:
    a=sorted(a)
    x=(len(a)-1)*p;i=int(x)
    return a[i] if i==len(a)-1 else (1-(x-i))*a[i]+(x-i)*a[i+1]


def source_stratum_analysis(sites:dict[str,dict],draws:int,seed:int)->dict:
    rows=list(sites.items())
    n=len(rows)
    if n<5: return {"identified":False,"n_sites":n}
    eligible=[(k,v) for k,v in rows if v["eligible_auc"]]
    out={
        "identified":len(eligible)>=4,
        "n_source_spatial_sites":n,
        "n_auc_eligible_original_sites":len(eligible),
        "n_2014_recorded_species_site_events":sum(v["n_original_2014_taxa"] for _,v in rows),
        "n_persisted_site_taxa":sum(v["n_persistent_recorded_2017"] for _,v in rows),
        "n_not_recorded_2017_site_taxa":sum(v["n_not_recorded_2017"] for _,v in rows),
        "n_newly_recorded_2017_site_taxa":sum(v["n_newly_recorded_2017"] for _,v in rows),
        "n_species_risk_sets_with_both_outcomes":len(eligible),
    }
    if not out["identified"]:
        out["no_mean_auc_reported_without_four_independent_sites"]=True
        return out
    vals=[v["source_2014_cover_AUC"] for _,v in eligible]
    observed=mean(vals)
    rng=random.Random(seed)
    null=[]
    for _ in range(draws):
        scores=[]
        for site,v in eligible:
            # Uniformly randomize persisted and lost labels among this
            # site's observed 2014 taxa, strictly preserving both counts.
            covers=v["raw_ranks_all"]
            n_ret=v["n_positive_retained_labels"]
            ids=rng.sample(range(len(covers)),n_ret)
            retained=[covers[i] for i in ids]
            other=[covers[i] for i in range(len(covers)) if i not in set(ids)]
            scores.append(auc_of_early_covers(retained,other))
        null.append(mean(scores))
    if abs(mean(null)-0.5)>0.07:
        raise AssertionError("Within-source-plot label-randomization expectation drift")
    loorange=[mean(vals[:i]+vals[i+1:]) for i in range(len(vals))]
    out.update({
        "site_equal_source_2014_cover_retained_vs_lost_AUC":observed,
        "fraction_sites_AUC_gt_0_5":sum(x>.5 for x in vals)/len(vals),
        "null_label_within_source_plot_mean":mean(null),
        "observed_minus_null_mean":observed-mean(null),
        "permutation_upper_rank_NOT_field_p":
          (1+sum(x>=observed-1e-12 for x in null))/(len(null)+1),
        "permutation_technical_q025":quantile(null,.025),
        "permutation_technical_q975":quantile(null,.975),
        "leave_one_original_source_site_out":{
            "min":min(loorange),"max":max(loorange)},
        "n_fixed_label_shuffle_draws":draws,
        "source_cover_used_only_as_within_plot_order_not_percent":True
    })
    return out


def run(path:Path,draws:int=DRAWS)->dict:
    if draws<99:raise ValueError("Need at least 99 within-plot label randomizations")
    species,source=parse_xlsx_matrix(path)
    raw,qc=early_numeric_cover(path,species)
    cohorts={g:matching_sites(species,g,ZONES) for g in GROUPS}
    if {g:len(v) for g,v in cohorts.items()}!={"station":15,"between":12}:
        raise ValueError("Original true common source cohort 15/12 drift")
    out={}
    for ix,(group,sites) in enumerate(cohorts.items()):
        byzone={}
        for zi,zone in enumerate(ZONES):
            rows={}
            for site in sites:
                before=raw[(site,2014,zone)]
                after=species[(site,2017,zone)]
                if before is None or after is None:
                    raise ValueError("Original site not in common 4-zone 2-year cohort")
                rows[site]=site_early_persistence(before,after)
            byzone[zone]=source_stratum_analysis(rows,draws,SEED+101*ix+43*zi)
        out[group]={
            "original_independent_spatial_sites":len(sites),
            "same_original_sites_all_2014_2017_zones":True,
            "source_2014_numeric_early_cover_to_2017_recorded_persistence":byzone
        }
    return {
        "version":"v25.9",
        "original_Dryad_source_sha256":source["original_file_sha256"],
        "original_source_provenance":source,
        "numeric_cover_reconciliation":qc,
        "metric":"Within-2014-plot rank AUC for persistence recorded again in 2017",
        "strata":out,
        "source_2014_cover_unit_not_assumed_percent":True,
        "cover_predicts_a_field_record_not_verified_sapling_survival":True,
        "loss_and_gain_not_certain_absence_or_colonization":True,
        "no_field_mowing_frequency_wind_gust_seed_transport_genetic_adaptation_measured":True,
        "not_separate_new_independent_railway_studies":True,
        "source_taxon_individual_covers_not_redistributed":True
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    p.add_argument("--draws",type=int,default=DRAWS)
    a=p.parse_args()
    out=run(a.source,a.draws)
    a.output.parent.mkdir(exist_ok=True,parents=True)
    a.output.write_text(json.dumps(out,ensure_ascii=False,indent=2)+"\n",encoding="utf8")
    for group,x in out["strata"].items():
        for zone,y in x["source_2014_numeric_early_cover_to_2017_recorded_persistence"].items():
            print(json.dumps({
                "version":"v25.9","group":group,"zone":zone,
                "n_early_species_site_records":y["n_2014_recorded_species_site_events"],
                "n_persisted":y["n_persisted_site_taxa"],
                "n_not_recorded_2017":y["n_not_recorded_2017_site_taxa"],
                "n_eligible_sites":y["n_auc_eligible_original_sites"],
                "source_original_within_plot_cover_auc":y.get("site_equal_source_2014_cover_retained_vs_lost_AUC"),
                "label_shuffle_mapping_rank":y.get("permutation_upper_rank_NOT_field_p")
            },ensure_ascii=False))


if __name__=="__main__":
    main()
