"""Run: python slovakia_rail.py <tetsu checkout> <sk.pdf> <canopy_height.txt> <seed_mass.txt>

Requires PyMuPDF and the dependencies of tetsu v24.5, commit 43884d11.
PDF: https://cbrb.sav.sk/usr/Janka/docs/BulletinSBS_2016_Majekova-Limanek.pdf
LEDA files: https://uol.de/f/5/inst/biologie/ag/landeco/download/LEDA/Data_files/
"""
import hashlib
import json
import math
import re
import sys
from pathlib import Path

import pymupdf

repo, pdf, height_file, seed_file = map(Path, sys.argv[1:])
sys.path.insert(0, str(repo / "scripts"))
import validate_Europe_3studies_82stations_original_LEDA_2D_trait_beta_v245 as model

for path, digest in [(pdf, "99465bdd7ed56985e65b1163f16eca4f702fc821619d84c89bf986d4d3736ae1"),
                     (height_file, model.LEDASHA["height"]), (seed_file, model.LEDASHA["seed"])]:
    assert hashlib.sha256(path.read_bytes()).hexdigest() == digest

with pymupdf.open(pdf) as doc:
    raw = "\n".join("\n".join(doc[p].get_text().splitlines()[2:]) for p in (4, 5, 6))
raw = raw[raw.index("Acer platanoides:"):raw.index("Zea mays: 17.") + len("Zea mays: 17")]
raw = re.sub(r"-\s*\n\s*(?=[a-z])", "", raw)
raw = " ".join(raw.split())
stations = [{"station": i, "source_names": set(), "region_group": "26b" if i <= 12 else "29"}
            for i in range(1, 20)]
names = set()
for entry in raw.split(";"):
    name, numbers = entry.strip().split(":")
    short, rest = name.split(" ", 1)
    if short.endswith("."):
        assert genus.startswith(short[:-1]), (genus, short)
    else:
        genus = short
    name = genus + " " + rest
    assert name not in names, name
    names.add(name)
    ids = []
    for item in numbers.split(","):
        bounds = re.fullmatch(r"(\d+)(?:[–-](\d+))?", item.strip())
        assert bounds, (name, item)
        first, last = bounds.groups()
        ids.extend(range(int(first), int(last or first) + 1))
    assert len(set(ids)) == len(ids) and all(1 <= i <= 19 for i in ids)
    for i in ids:
        stations[i-1]["source_names"].add(name)

counts = [len(s["source_names"]) for s in stations]
assert len(names) == 318
assert [counts[i-1] for i in (19, 15, 6, 12)] == [202, 171, 146, 134]
print("Source validation:", len(names), counts, file=sys.stderr)
height, _ = model.read_leda(height_file, "height_m")
seed_mass, _ = model.read_leda(seed_file, "seed_mass_mg")
# Frozen European v24.5 species-equal log scale, artifact 14296521971.
mu = [-0.5802365754696401, 0.2936113186962433]
sd = [1.3801000545693773, 2.1657660559995824]
traits = {t: ((math.log(height[t])-mu[0])/sd[0], (math.log(seed_mass[t])-mu[1])/sd[1])
          for t in height.keys() & seed_mass.keys()}
cohorts, coverage = model.source_cohort({"SK": stations}, traits)
sites = cohorts["SK"]
observed = model.get_metrics(sites, traits)
chains = [model.one_null(sites, traits, seed) for seed in (202611171, 202611172)]
draws = chains[0][0] + chains[1][0]
print(json.dumps({"source_taxa": len(names), "station_source_counts": counts, "coverage": coverage,
                  "station_matched_counts": [len(s["taxa"]) for s in sites],
                  "chains": [qa for _, qa in chains], "draws": len(draws),
                  "results": model.summarize(observed, draws)}, indent=2))
