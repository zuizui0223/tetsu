"""Usage: python kashmir_paired.py TETSU_REPO SOURCE_DIR > results.json

SOURCE_DIR contains Matrix_plot.xlsx from https://datadryad.org/api/v2/versions/99735/download
and canopy_height.txt and seed_mass.txt from the tetsu v24.5 LEDA sources.
Use tetsu commit 43884d11, Python with numpy, scipy, openpyxl.
"""
import hashlib
import json
import sys
from itertools import product
from pathlib import Path

import numpy as np
import openpyxl
from scipy.spatial.distance import pdist, squareform

repo, source = map(Path, sys.argv[1:])
sys.path.insert(0, str(repo / "scripts"))
from audit_SI_CZ_original_continuous_trait_coverage_v220 import binomial, read_leda

digests = {"Matrix_plot.xlsx": "529f3b97521933018e12623b61f531f012c67241cae2602a6675dbd534192474",
           "canopy_height.txt": "cbf5e761390d78a843512dfc1139acfc1154d46d6c1995a827354840a9f78a92",
           "seed_mass.txt": "34aa12c333fc1ba6de2e560a495ac5d693aa0722a09ecf68418ad94ac1c19cef"}
for filename, expected in digests.items():
    assert hashlib.sha256((source / filename).read_bytes()).hexdigest() == expected
height, _ = read_leda(source / "canopy_height.txt", "height_m")
mass, _ = read_leda(source / "seed_mass.txt", "seed_mass_mg")
with (source / "Matrix_plot.xlsx").open("rb") as f:
    book = openpyxl.load_workbook(f, data_only=True, read_only=True)
    header, *rows = list(book.worksheets[0].values)
names = [binomial(" ".join(x.split())) for x in header[6:]]
assert len(names) == 123
valid = [i for i, t in enumerate(names) if t in height and t in mass]
assert len({names[i] for i in valid}) == len(valid)
axes = (np.log([[height[names[i]], mass[names[i]]] for i in valid])
        - [-0.5802365754696401, 0.2936113186962433]) / [1.3801000545693773, 2.1657660559995824]
plots = {}
for location, _, _, _, year, zone, *values in rows:
    key = (" ".join(location.split()), year, zone)
    assert key not in plots and all(v is None or v >= 0 for v in values)
    plots[key] = np.asarray(values, dtype=float)
sites = sorted({k[0] for k in plots})
years = [2014, 2017]
zones = ["Track", "Embankment", "Interior Plot 1", "Interior Plot 2"]
assert set(plots) == set(product(sites, years, zones)) and len(sites) == 31

def eligible(site, year, zone, n, coverage):
    v = plots[site, year, zone]
    return (np.isfinite(v).all() and (v[valid] > 0).sum() >= n
            and (v[valid] > 0).sum() >= coverage * (v > 0).sum())

def analyze(cohort, year, rail, reference):
    if len(cohort) < 3:
        return {"n": len(cohort)}
    sides = [[plots[s, year, zone] for s in cohort] for zone in (rail, reference)]
    present = [[np.flatnonzero(v[valid] > 0) for v in side] for side in sides]
    profiles = [np.array([axes[p].mean(axis=0) for p in side]) for side in present]
    distances = [pdist(side) for side in profiles]
    delta = distances[0] - distances[1]
    matrix = squareform(delta)
    loo = [(matrix.sum()-2*r.sum())/((len(r)-1)*(len(r)-2)) for r in matrix]
    rng = np.random.default_rng(251117)
    sizes = np.minimum(*[[len(p) for p in side] for side in present])
    rare = []
    for _ in range(199):
        matched = [[axes[rng.choice(p, k, replace=False)].mean(axis=0)
                    for p, k in zip(side, sizes)] for side in present]
        rare.append(pdist(matched[0]).mean()-pdist(matched[1]).mean())
    tax = [pdist(np.array(side) > 0, metric="jaccard").mean() for side in sides]
    return {"n": len(cohort), "rail": distances[0].mean(), "reference": distances[1].mean(),
            "delta": delta.mean(), "leave_one_site_out": [min(loo), max(loo)],
            "richness_equalized_mean_delta": np.mean(rare),
            "species_Jaccard_dissimilarity_delta": tax[0]-tax[1]}

results = []
taxonomy = []
for group, population in [("station", [s for s in sites if not s.startswith("Between ")]),
                          ("between", [s for s in sites if s.startswith("Between ")])]:
    for rail, ref in [(zones[0], zones[3]), (zones[0], zones[2]), (zones[1], zones[3])]:
        complete = [s for s in population if all(np.isfinite(plots[s,y,z]).all()
                    and (plots[s,y,z] > 0).any() for y,z in product(years,[rail,ref]))]
        for year in years:
            pres = [np.array([plots[s,year,z] > 0 for s in complete]) for z in (rail,ref)]
            if len(complete) >= 3:
                observed = [pdist(a, "jaccard").mean() for a in pres]
                rng = np.random.default_rng(261117)
                sizes = np.minimum(pres[0].sum(axis=1), pres[1].sum(axis=1))
                draws = []
                for _ in range(199):
                    distances = []
                    for side in pres:
                        sampled = np.zeros_like(side)
                        for i, row in enumerate(side):
                            sampled[i,rng.choice(np.flatnonzero(row),sizes[i],replace=False)] = True
                        distances.append(pdist(sampled,"jaccard").mean())
                    draws.append(distances[0]-distances[1])
                taxonomy.append({"group":group,"zones":[rail,ref],"year":year,"n":len(complete),
                                 "all_species_distance_delta":observed[0]-observed[1],
                                 "richness_equalized_delta":np.mean(draws)})
        for n, coverage in [(5, .25), (10, .5)]:
            # Both years on exactly the same paired spatial cohort.
            cohort = [s for s in population if all(eligible(s,y,z,n,coverage)
                      for y,z in product(years,[rail,ref]))]
            results.append({"group": group, "zones": [rail, ref], "minimum": [n,coverage],
                            "single_year":{y:analyze([s for s in population if all(
                                eligible(s,y,z,n,coverage) for z in [rail,ref])],y,rail,ref) for y in years},
                            "years": {y: analyze(cohort,y,rail,ref) for y in years}})
print(json.dumps({"source_taxa":len(names),"both_traits":len(valid),"source_spatial_sites":len(sites),
                  "missing_cells":sum(int(np.isnan(v).sum()) for v in plots.values()),
                  "matched_both_years":results,"all_species":taxonomy}, indent=2))
