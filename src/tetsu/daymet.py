"""Strict daily Daymet source parser and independent pre-flowering climate summaries.

Daymet 1 km air temperature describes a gridded environmental covariate.
It is NOT ballast temperature, train slipstream, management or selection.
"""
from __future__ import annotations

import csv
import io
import math
from collections import defaultdict


YEARS = (2023, 2024, 2025)
FIXED_DAYS = (60, 90, 120, 150)
BASE_TEMP_C = 5.0


def parse_daymet(text: str, years=YEARS) -> dict[int, list[dict]]:
    lines = text.splitlines()
    header_position = next(
        (i for i, line in enumerate(lines)
         if line.lstrip().lower().startswith("year,yday,")),
        None
    )
    if header_position is None:
        raise ValueError("Daymet returned no year,yday header; login/error HTML or schema drift?")
    header_lines = "\n".join(lines[header_position:])
    reader = csv.DictReader(io.StringIO(header_lines))
    if reader.fieldnames is None:
        raise ValueError("empty Daymet table")
    field_map = {
        str(x).strip().lower().split(" (")[0]: x for x in reader.fieldnames
    }
    expected = {"year", "yday", "tmax", "tmin", "prcp"}
    if not expected.issubset(field_map):
        raise ValueError(f"Missing Daymet columns: {sorted(expected - set(field_map))}")
    rows = defaultdict(list)
    seen = set()
    for line in reader:
        if not any(value for value in line.values()):
            continue
        try:
            year = int(line[field_map["year"]])
            day = int(line[field_map["yday"]])
            tmax = float(line[field_map["tmax"]])
            tmin = float(line[field_map["tmin"]])
            rain = float(line[field_map["prcp"]])
        except (TypeError, ValueError, KeyError) as exc:
            raise ValueError("Invalid Daymet day or value") from exc
        if year not in years:
            raise ValueError(f"Daymet returned unexpected year {year}")
        if not (1 <= day <= 365) or (year, day) in seen:
            raise ValueError(f"Repeated/invalid Daymet day {year} day={day}")
        seen.add((year, day))
        if not all(math.isfinite(v) for v in (tmax, tmin, rain)):
            raise ValueError("non-finite Daymet weather values")
        if tmin > tmax or rain < 0 or tmax > 65 or tmin < -80:
            raise ValueError("Daymet physically invalid weather values")
        rows[year].append({
            "year": year, "doy": day, "tmax_c": tmax, "tmin_c": tmin,
            "precip_mm": rain
        })
    for year in years:
        if year not in rows or len(rows[year]) != 365:
            raise ValueError(
                f"Incomplete Daymet year {year}: {len(rows.get(year, []))}/365"
            )
        rows[year].sort(key=lambda x: x["doy"])
        if any(row["doy"] != i for i, row in enumerate(rows[year], 1)):
            raise ValueError(f"Missing Daymet day in year {year}")
    return dict(rows)


def summarize_annual_climate(rows, base_c=BASE_TEMP_C):
    if not isinstance(rows, list) or len(rows) != 365:
        raise ValueError("need full 365-day Daymet series")
    if not math.isfinite(base_c):
        raise ValueError("base temperature must be finite")
    cumulative_gdd = 0.
    cumulative_precip = 0.
    heat_sum = 0.
    metrics = {}
    for index, item in enumerate(rows, 1):
        if item["doy"] != index:
            raise ValueError("Daymet day order inconsistent")
        daily_mean = (item["tmax_c"] + item["tmin_c"]) / 2
        cumulative_gdd += max(0., daily_mean - base_c)
        cumulative_precip += item["precip_mm"]
        heat_sum += daily_mean
        if index in FIXED_DAYS:
            metrics[f"gdd5_to_doy{index}"] = round(cumulative_gdd, 2)
            metrics[f"precip_mm_to_doy{index}"] = round(cumulative_precip, 2)
            metrics[f"mean_air_c_to_doy{index}"] = round(heat_sum / index, 3)
    return metrics
