#!/usr/bin/env python3
"""Kashmir original Dryad 2014/2017: matched rail-versus-interior annual incidence.

Requires ORIGINAL, checksum-authenticated Matrix_plot.xlsx and an auditable
original-name life-span CSV. Does not infer life span from genus/family and never
publishes original occurrence rows. This is one study, not two independent years.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import random
import statistics
import sys
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from analyze_Kashmir_same_sites_reference_year_v254 import (  # noqa: E402
    DRYAD_SHA256, NS, PKG_NS, REL_NS, parse_xlsx_matrix,
)

ZONES = ("Track", "Embankment", "Interior Plot 1", "Interior Plot 2")
YEARS = (2014, 2017)
PAIRS = (("Track", "Interior Plot 1"), ("Track", "Interior Plot 2"),
         ("Track", "Embankment"), ("Embankment", "Interior Plot 1"),
         ("Embankment", "Interior Plot 2"))
KNOWN = {"annual", "perennial", "biennial"}


def canonical(text: str) -> str:
    return " ".join(text.split()).casefold()


def original_taxon_headers(path: Path) -> list[str]:
    """Read the 123 genuine column names, respecting XLSX shared strings."""
    if hashlib.sha256(path.read_bytes()).hexdigest() != DRYAD_SHA256:
        raise ValueError("Dryad Matrix_plot.xlsx SHA256 mismatch")
    with zipfile.ZipFile(path) as z:
        shared = []
        if "xl/sharedStrings.xml" in z.namelist():
            s = ET.fromstring(z.read("xl/sharedStrings.xml"))
            shared = ["".join(t.text or "" for t in si.iter(f"{{{NS}}}t"))
                      for si in s.findall(f"{{{NS}}}si")]
        book = ET.fromstring(z.read("xl/workbook.xml"))
        rid = book.find(f"{{{NS}}}sheets/{{{NS}}}sheet").get(f"{{{REL_NS}}}id")
        rels = ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))
        target = next((r.get("Target") for r in rels.findall(f"{{{PKG_NS}}}Relationship")
                       if r.get("Id") == rid and "worksheet" in (r.get("Type") or "")), None)
        if not target:
            raise ValueError("Worksheet location missing")
        target = "xl/" + target.lstrip("/").removeprefix("xl/")
        root = ET.fromstring(z.read(target))
        first = root.find(f"{{{NS}}}sheetData/{{{NS}}}row")
        if first is None:
            raise ValueError("No XLSX header")
        bycol = {}
        for cell in first.findall(f"{{{NS}}}c"):
            addr = cell.get("r", "")
            letters = "".join(x for x in addr if x.isalpha())
            pos = 0
            for ch in letters:
                pos = pos * 26 + ord(ch) - 64
            pos -= 1
            v = cell.find(f"{{{NS}}}v")
            raw = v.text if v is not None else None
            if cell.get("t") == "s":
                label = shared[int(raw)]
            elif cell.get("t") == "inlineStr":
                label = "".join(x.text or "" for x in cell.iter(f"{{{NS}}}t"))
            else:
                label = raw
            bycol[pos] = str(label or "")
    if sorted(bycol) != list(range(129)):
        raise ValueError("Source header not 129 columns")
    species = [bycol[i].strip() for i in range(6, 129)]
    if len(set(canonical(x) for x in species)) != 123 or any(not x for x in species):
        raise ValueError("Duplicate/empty original taxon header")
    return species


def life_span_crosswalk(path: Path, original_names: list[str]) -> tuple[dict[int, str], dict]:
    """Only explicitly sourced, exact original-name annotations are eligible."""
    index = {canonical(name): i for i, name in enumerate(original_names)}
    classified = {}
    seen = set()
    unmatched = []
    with path.open(newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        if not {"original_taxon", "life_span", "trait_source"}.issubset(reader.fieldnames or []):
            raise ValueError("Traits CSV needs original_taxon,life_span,trait_source")
        for row in reader:
            name = canonical(row["original_taxon"] or "")
            span = canonical(row["life_span"] or "")
            provenance = (row["trait_source"] or "").strip()
            if not name or name in seen:
                raise ValueError("Empty/duplicate source taxon in trait crosswalk")
            seen.add(name)
            if span not in KNOWN | {"unknown"}:
                raise ValueError(f"Unsupported life span for {name}: {span!r}")
            if span in KNOWN and not provenance:
                raise ValueError("Each classified life span needs an explicit trait_source")
            if name not in index:
                unmatched.append(row["original_taxon"])
                continue
            if span in KNOWN:
                classified[index[name]] = span
    if unmatched:
        raise ValueError(f"Trait crosswalk contains {len(unmatched)} unmatched original names; do not fuzzy-join")
    return classified, {"taxa_in_original": len(original_names),
                        "taxa_classified": len(classified),
                        "taxa_unclassified": len(original_names) - len(classified),
                        "crosswalk_sha256": hashlib.sha256(path.read_bytes()).hexdigest()}


def annual_stat(species: set[int], classified: dict[int, str]) -> dict:
    n = len(species)
    if n == 0:
        raise ValueError("Empty field zone has undefined annual fraction")
    annual = sum(classified.get(i) == "annual" for i in species)
    known = sum(i in classified for i in species)
    unknown = n - known
    return {"present": n, "annual": annual, "known": known, "unknown": unknown,
            "coverage": known/n,
            "annual_fraction": annual/known if known else None,
            "minimum_fraction": annual/n,
            "maximum_fraction": (annual + unknown)/n}


def q(sorted_data: list[float], p: float) -> float:
    loc = (len(sorted_data)-1)*p
    lo = int(loc)
    hi = min(lo+1,len(sorted_data)-1)
    return sorted_data[lo] + (loc-lo)*(sorted_data[hi]-sorted_data[lo])


def summarize(plots: dict, classified: dict[int,str],
              min_coverage: float = .7, bootstrap: int = 1999,
              seed: int = 285) -> dict:
    if not 0 <= min_coverage <= 1:
        raise ValueError("Coverage must be within [0,1]")
    source_sites = sorted({key[0] for key in plots})
    common = [site for site in source_sites if all(
        plots.get((site, year, zone)) is not None
        for year in YEARS for zone in ZONES)]
    result = {"original_sites": len(source_sites),
              "complete_four_zones_both_years_sites": len(common),
              "site_group": {"station": sum(not x.startswith("Between ") for x in common),
                             "between": sum(x.startswith("Between ") for x in common)},
              "minimum_known_trait_coverage": min_coverage,
              "source_years_not_independent": True,
              "effect_scale": "site-weighted delta strict annual fraction of species occurrences",
              "comparisons": {}}
    stats = {(site,year,zone): annual_stat(plots[site,year,zone], classified)
             for site in common for year in YEARS for zone in ZONES}
    for inner, outer in PAIRS:
        pair_name = inner + "__minus__" + outer
        result["comparisons"][pair_name] = {}
        for group in ("station", "between", "all"):
            candidate = [site for site in common if group == "all" or
                         ("between" if site.startswith("Between ") else "station") == group]
            per_site = {}
            for site in candidate:
                samples = []
                for year in YEARS:
                    a = stats[site,year,inner]
                    b = stats[site,year,outer]
                    if (a["coverage"] < min_coverage or b["coverage"] < min_coverage
                        or a["annual_fraction"] is None or b["annual_fraction"] is None):
                        break
                    samples.append({"delta":a["annual_fraction"] - b["annual_fraction"],
                                    "lower":a["minimum_fraction"] - b["maximum_fraction"],
                                    "upper":a["maximum_fraction"] - b["minimum_fraction"],
                                    "year":year})
                if len(samples) == len(YEARS):
                    per_site[site] = samples
            report = {"eligible_original_sites": len(per_site),
                      "ineligible_original_sites": len(candidate)-len(per_site),
                      "per_year": {}, "two_year_mean": None,
                      "coverage_extreme_bounds": None,
                      "site_cluster_bootstrap_technical_percentiles_not_field_CI": None}
            if per_site:
                report["per_year"] = {str(y):statistics.mean(
                    v[yi]["delta"] for v in per_site.values())
                    for yi,y in enumerate(YEARS)}
                effects = [statistics.mean(v[yi]["delta"] for yi in range(2))
                           for v in per_site.values()]
                report["two_year_mean"] = statistics.mean(effects)
                report["sites_positive"] = sum(x>0 for x in effects)
                report["sites_negative"] = sum(x<0 for x in effects)
                report["coverage_extreme_bounds"] = [
                    statistics.mean(statistics.mean(x[yi]["lower"] for yi in range(2))
                                    for x in per_site.values()),
                    statistics.mean(statistics.mean(x[yi]["upper"] for yi in range(2))
                                    for x in per_site.values())]
                if len(effects) >= 2 and bootstrap:
                    rng = random.Random(seed + len(inner) * 19 + len(outer) * 37 + len(group))
                    draws = sorted(statistics.mean(rng.choices(effects,k=len(effects)))
                                   for _ in range(bootstrap))
                    report["site_cluster_bootstrap_technical_percentiles_not_field_CI"] = [
                        q(draws,.025), q(draws,.975)]
            result["comparisons"][pair_name][group] = report
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True, help="Original Dryad version 99735 Matrix_plot.xlsx")
    parser.add_argument("--traits", type=Path, required=True, help="Original-name and evidence-labelled life-span CSV")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--min-coverage", type=float, default=.7)
    parser.add_argument("--bootstrap", type=int, default=1999)
    args = parser.parse_args()
    if args.bootstrap < 0:
        parser.error("--bootstrap must be nonnegative")
    plots, source_metadata = parse_xlsx_matrix(args.source)
    names = original_taxon_headers(args.source)
    classified, crosswalk = life_span_crosswalk(args.traits,names)
    result = {"version":"v28.5", "source": source_metadata,
              "trait_crosswalk": crosswalk,
              "life_span_is_not_genetically_measured": True,
              "validation":summarize(plots, classified, args.min_coverage, args.bootstrap)}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({"result_path":str(args.output), "trait_classified":len(classified),
                      "common_sites":result["validation"]["complete_four_zones_both_years_sites"]},ensure_ascii=False))


if __name__ == "__main__":
    main()
