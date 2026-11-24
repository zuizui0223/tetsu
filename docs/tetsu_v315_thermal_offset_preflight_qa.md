# v315 — Thermal-offset preflight QA (2026-11-22)

The earlier v314 preflight computes a community-weighted thermal indicator for matched railway/nonrail sites. This QA identifies important gates **before** fitting the proposed cold-region interaction.

## Required gates
- **Trait coverage:** v314 reports covered incidence fraction but does not exclude low-coverage comparisons. Predeclare a sensitivity grid at >=0.50, >=0.70 and >=0.90 coverage in *both* habitats; no filling 'x' with zero.
- **Source habitat:** within-railway ballast shoulder versus toe is not a railway/nonrail contrast. A national MIREN roadside reference distant from the rail site is not a local matched plot.
- **Climate replication:** the 19 Slovak stations are one geographically confounded transect; dependent station pairs do not create independent climatic regions. A slope fitted only within that line is not a cross-region climate-offset test.
- **Survey comparability:** record area, effort, year, and openness for each rail and nonrail plot; unmatched studies are descriptive only.
- **Taxon traits:** a thermal indicator is an ordinal preference index, not degrees Celsius. Species-level thermal niches from occurrences require sampling-bias and geographic-range controls.
- **Inference:** use study/region-level block resampling or leave-region-out checks; no causal claim of train dispersal, climate-limit release, or persistent colonization from observational CWM contrasts.

## Independent synthetic QA
A separate local Python implementation with an explicit coverage eligibility flag was run against three synthetic tests on 2026-11-22: missing indicator 'x' yields rail coverage 0.5 and fails a 0.7 threshold while retaining descriptive delta=2; complete coverage passes; mismatched within-site regional climate raises an error. **3/3 synthetic tests passed locally.** This is not a rerun of original ecological data, and no source-backed climate effect has been estimated.

## Next source-backed analysis
Recover the genuine study-by-site incidence and original thermal-indicator tables and audit whether paired *local* nonrail plots exist in at least two independent climate regions. If not, report the identifiability limit instead of a pooled climate interaction. Reproduce previously reported typical-taxon counts from original input files separately; do not select indicator species on the outcome data used to test the climate hypothesis.
