"""Repeated individual phenology logic for USA-NPN Open flowers observations.

NPN status 1=yes, 0=no, -1=uncertain. This module never imputes
unobserved days or interprets photo absence as biological no.
Event onsets are interval-censored; genuine site/rail treatment is unknown.
"""
from __future__ import annotations
from collections import Counter, defaultdict
from datetime import date

ALLOWED = {-1, 0, 1}


def _s(value):
    if value is None:
        return ""
    return str(value).strip()


def parse_status(value):
    try:
        out=int(str(value))
    except (ValueError,TypeError) as exc:
        raise ValueError("invalid observed NPN phenophase status") from exc
    if out not in ALLOWED:
        raise ValueError("unexpected NPN phase status")
    return out


def normalize_open_flowers(records):
    """Return one deconflicted record per site × plant × day.

    Conflicting simultaneous reported statuses are excluded. Duplicate
    agreements are deduplicated; unknown statuses excluded. All exclusions
    reported, no biologically unsupported zero imputation.
    """
    by_key=defaultdict(list)
    reasons=Counter()
    for row in records:
        if _s(row.get("phenophase_description"))!="Open flowers":
            reasons["other_phase_not_open_flowers"]+=1
            continue
        site=_s(row.get("site_id"))
        individual=_s(row.get("individual_id"))
        when=_s(row.get("observation_date"))[:10]
        if not site or not individual:
            reasons["missing_site_or_plant_id"]+=1
            continue
        try:
            d=date.fromisoformat(when)
            status=parse_status(row.get("phenophase_status"))
        except ValueError:
            reasons["invalid_date_or_status"]+=1
            continue
        by_key[(site,individual,d)].append(status)
    clean=[]
    for key, statuses in by_key.items():
        definite={s for s in statuses if s>=0}
        if len(definite)>1:
            reasons["conflicting_same_day"]+=len(statuses)
            continue
        if not definite:
            reasons["uncertain_no_imputation"]+=len(statuses)
            continue
        if len(statuses)>1:
            reasons["same_day_duplicate_rows_deduplicated"]+=len(statuses)-1
        site,individual,d=key
        clean.append({
            "site_id":site,"individual_id":individual,"date":d,
            "year":d.year,"status":next(iter(definite))
        })
    clean.sort(key=lambda x:(x["site_id"],x["individual_id"],x["date"]))
    return clean,dict(reasons)


def summarize_individuals(clean,prior_no_max_days=21,followup_max_days=21):
    """Descriptive incident flower intervals and potential return flowering.

    Require a qualifying prior observed NO before first observed YES in year.
    The 'first observed YES' is not necessarily true annual onset. A second
    flowering bout requires two explicit NO visits after a YES and subsequent
    YES within followup_max_days of last NO. Seasonal/detection explanations
    remain open even when the pattern is observed.
    """
    if prior_no_max_days<1 or followup_max_days<1:
        raise ValueError("day intervals must be positive")
    groups=defaultdict(list)
    for row in clean:
        groups[(row["site_id"],row["individual_id"],row["year"])].append(row)
    summaries=[]
    for (site,ind,year),rows in sorted(groups.items()):
        rows.sort(key=lambda r:r["date"])
        n_yes=sum(r["status"]==1 for r in rows)
        n_no=sum(r["status"]==0 for r in rows)
        # Interval-censored observed re-entry into flowering.
        all_gaps=[]
        last_no=None
        for row in rows:
            if row["status"]==0:
                last_no=row["date"]
            elif last_no is not None:
                gap=(row["date"]-last_no).days
                if gap <= prior_no_max_days:
                    all_gaps.append(gap)
                last_no=None
        # First yes of growing season preceded by explicit recently observed no
        first_yes=next((r["date"] for r in rows if r["status"]==1),None)
        first_yes_prior_no=None
        if first_yes:
            candidates=[r["date"] for r in rows if r["date"]<first_yes and r["status"]==0]
            if candidates:
                last=max(candidates)
                if (first_yes-last).days<=prior_no_max_days:
                    first_yes_prior_no=last
        # Distinct recurrent blooming is a strong pattern with 2 observed
        # NOs between two YES observations; gaps still need surveillance.
        seen_yes=False
        gap_nos=[]
        resumed=False
        for row in rows:
            if row["status"]==1:
                if seen_yes and len(gap_nos)>=2:
                    if gap_nos[-1]-gap_nos[0]>=__import__("datetime").timedelta(days=1) and (
                        row["date"]-gap_nos[-1]).days<=followup_max_days:
                        resumed=True
                seen_yes=True;gap_nos=[]
            elif seen_yes:
                gap_nos.append(row["date"])
        summaries.append({
            "site_id":site,"individual_id":ind,"year":year,
            "n_dates":len(rows),"n_yes":n_yes,"n_no":n_no,
            "both_yes_no_observed":bool(n_yes and n_no),
            "first_yes_observed_date":first_yes.isoformat() if first_yes else None,
            "first_yes_prior_no_date":first_yes_prior_no.isoformat() if first_yes_prior_no else None,
            "first_yes_bracket_days":(first_yes-first_yes_prior_no).days if first_yes_prior_no else None,
            "n_no_to_yes_pairs_leq_21_days":len(all_gaps),
            "has_two_no_interrupted_flowering":resumed
        })
    return summaries


def coverage_summary(rows,clean,excluded,summaries):
    return {
        "source_records_total":len(rows),
        "open_flower_rows_raw":sum(_s(r.get("phenophase_description"))=="Open flowers" for r in rows),
        "unique_valid_phase_visit_rows":len(clean),
        "yes_visit_rows":sum(r["status"]==1 for r in clean),
        "no_visit_rows":sum(r["status"]==0 for r in clean),
        "n_distinct_sites":len({r["site_id"] for r in clean}),
        "n_distinct_individuals":len({(r["site_id"],r["individual_id"]) for r in clean}),
        "n_individual_years":len(summaries),
        "individual_years_with_yes_and_no":sum(s["both_yes_no_observed"] for s in summaries),
        "individual_years_with_bracketed_first_yes":sum(s["first_yes_prior_no_date"] is not None for s in summaries),
        "individual_years_with_2no_interrupted_flowering":sum(s["has_two_no_interrupted_flowering"] for s in summaries),
        "excluded_or_deduplicated":excluded,
        "critical_caveat":"These are individual flower-status time series, not evidence of mowing, train exposure or adaptive flowering shift."
    }
