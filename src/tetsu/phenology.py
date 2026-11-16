"""Observed flowering-probability summaries with phenology quality gates.

Never convert first photographic flower detection to true first flowering date.
No inferential comparisons or claims about genetic adaptation are implemented.
"""
from __future__ import annotations

import csv
import datetime as dt
from collections import defaultdict
from pathlib import Path

ALLOWED_SITES = {"rail", "road", "background"}
ALLOWED_PHASES = {
    "unknown", "vegetative", "bud", "flowering", "fruiting",
    "flowering_fruiting"
}
FLOWER_PHASES = {"flowering", "flowering_fruiting"}
ALLOWED_EVENTS = {"cutting", "herbicide", "steam", "sheeting", "mulching", "unknown"}


def parse_date(value: str) -> dt.date:
    try:
        x = dt.date.fromisoformat(value)
    except ValueError as exc:
        raise ValueError(f"invalid ISO date: {value}") from exc
    return x


def load_events(path: str | Path) -> dict[str, list[dt.date]]:
    """Only independently documented cutting dates; other management stays excluded."""
    if not path:
        return {}
    times = defaultdict(list)
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if not {"site_id", "event_date", "event_type", "evidence"}.issubset(reader.fieldnames or []):
            raise ValueError("events require site_id,event_date,event_type,evidence")
        for row in reader:
            site = row["site_id"].strip()
            etype = row["event_type"].strip()
            evidence = row["evidence"].strip()
            date = parse_date(row["event_date"].strip())
            if not site or etype not in ALLOWED_EVENTS or not evidence:
                raise ValueError("events must include a site, accepted type, and evidence")
            if etype == "cutting":
                times[site].append(date)
    return {site: sorted(set(days)) for site, days in times.items()}


def normalize_visits(path: str | Path, events=None) -> tuple[list[dict], dict]:
    """Return reviewed nonunknown observations from fixed complete survey visits.

    Completeness means a survey had the opportunity to record nonflowering
    specimens; a missing row is NOT a nonflowering observation.
    """
    events = events or {}
    rows = []
    skipped = {"incomplete_survey": 0, "unknown_phenophase": 0}
    seen = set()
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {
            "visit_id", "site_id", "taxon", "obs_date", "site_class",
            "phenophase", "survey_complete", "observer_id", "effort_minutes"
        }
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"visits missing {sorted(missing)}")
        for raw in reader:
            visit_id = raw["visit_id"].strip()
            site = raw["site_id"].strip()
            taxon = raw["taxon"].strip()
            observer = raw["observer_id"].strip()
            site_class = raw["site_class"].strip()
            phase = raw["phenophase"].strip()
            date = parse_date(raw["obs_date"].strip())
            if not all((visit_id, site, taxon, observer)):
                raise ValueError("blank visit/site/taxon/observer ID")
            if site_class not in ALLOWED_SITES or phase not in ALLOWED_PHASES:
                raise ValueError(f"invalid site or phenophase {site_class}/{phase}")
            key = (site, visit_id, taxon)
            if key in seen:
                raise ValueError(f"duplicate site/visit/taxon: {key}")
            seen.add(key)
            flag = raw["survey_complete"].strip().lower()
            if flag not in ("true", "false", "1", "0"):
                raise ValueError(f"survey_complete not bool: {flag}")
            try:
                effort = float(raw["effort_minutes"])
            except ValueError as exc:
                raise ValueError("effort_minutes must be numeric") from exc
            import math
            if not math.isfinite(effort) or effort <= 0:
                raise ValueError("effort_minutes must be positive and finite")
            if flag in ("false", "0"):
                skipped["incomplete_survey"] += 1
                continue
            if phase == "unknown":
                skipped["unknown_phenophase"] += 1
                continue
            eligible_cut_dates = [d for d in events.get(site, []) if d <= date]
            last = max(eligible_cut_dates) if eligible_cut_dates else None
            iso = date.isocalendar()
            rows.append({
                "visit_id": visit_id, "site_id": site, "taxon": taxon,
                "obs_date": date.isoformat(), "year": date.year, "doy": date.timetuple().tm_yday,
                "week": iso.week, "iso_week_year": iso.year,
                "site_class": site_class, "phenophase": phase, "observer_id": observer,
                "effort_minutes": effort, "flower_observed": int(phase in FLOWER_PHASES),
                "last_verified_cut_date": last.isoformat() if last else "",
                "days_since_verified_cut": (date - last).days if last else None,
            })
    return rows, skipped


def weekly_summary(visits: list[dict]) -> list[dict]:
    """Descriptive flowering proportions only, not onset and not treatment effects."""
    acc = defaultdict(list)
    for row in visits:
        k = (row["taxon"], row["site_class"], row["iso_week_year"], row["week"])
        acc[k].append(row)
    result = []
    for (taxon, site_class, year, week), items in sorted(acc.items()):
        n = len(items)
        sites = {x["site_id"] for x in items}
        n_flower = sum(x["flower_observed"] for x in items)
        result.append({
            "taxon": taxon, "site_class": site_class, "iso_week_year": year,
            "iso_week": week, "n_complete_visits": n, "n_sites": len(sites),
            "n_flowering_visits": n_flower,
            "observed_flowering_proportion": round(n_flower/n, 6),
            "n_with_verified_prior_cut": sum(
                x["days_since_verified_cut"] is not None for x in items
            ),
            "screening_only": int(n < 5 or len(sites) < 3),
        })
    return result
