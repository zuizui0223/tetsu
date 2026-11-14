# tetsu v5.3 — Locked out-of-analysis-region US Southeast replication check

**Protocol frozen on 2026-11-10 JST, BEFORE querying its 50km RAIL/ROAD routes or associated species similarity outcome.**

## Why this specific region

The existing original v2.9 / v3.6 geographical source frame already fixed an independent 200×200km 50km-cell US Southeast block centered at longitude −78.60, latitude +35.80, separate from the US Midwest block. In original six-country ecology v4.4–v5.2, US Southeast was NEVER used to estimate the railway plant similarity coefficients or select the K80 effect. The original source-coverage audit v3.6 (obtained before examining its RAIL/ROAD/flora correlation) reported 12/16 GRIIS-listed introduced-plant cells passing both period K20 minimum, without testing K80.

We freeze this region to test an OUT-OF-ANALYSIS-REGION association in the same country but distinct climate and geographic location. Because we have already seen the six-country source outcomes and selected K80 posthoc, this is NOT an independent global country holdout or preregistered causal train effect, but a reasonable geographically untouched out-of-analysis-region sign/feasibility check.

## Frozen primary diagnostic

- Geographic unit: existing exact EPSG:6933 50km²? **No**: each source cell is a 50km×50km projected equal-area rectangle, 2,500km²; original 4×4 cells (16) in a 200×200km block. Do NOT reselect cell footprint to favor an effect.
- Taxa: country-GRIIS-listed introduced plant GBIF speciesKey positive indexed source counts, original country+source gating from v3.6; country status US only. Unlisted taxa UNKNOWN, not native.
- Source years: main period **2018–2025**; early 2010–2017 reported but not used to choose the main metric.
- Main recording depth K=**80** indexed introduced-plant RECORDS without replacement in each eligible cell, same set of cells also eligible in both periods. Species similarity response is expected rarefied intersection divided by expected union, NOT E[Jaccard], no observed ecological absences.
- Minimum observed introduced flora source gate: >=5 GRIIS-listed speciesKeys and >=20 recorded introduced observations per cell-period from v3.6, and >=80 indexed introduced source RECORDS in BOTH years blocks for this holdout, plus original >=80% country polygon and >=100 all-plant+3 datasets per period.
- Same source RAIL and ROAD physical paths from official 2026-10 Overture with 75km halo, **nearest connector inside each 50km cell** primary source anchor. Rail route over road route is the predictive exposure. This is not an actual operational freight service route.
- Partial rank: Spearman-like residual ranks for log(rail centroid route km / road centroid route km) and introduced-flora source similarity, controlling geographic great-circle route distance, physical ROAD path km, all-plant GBIF indexed record-count log ratio and log geometric mean, and GBIF all-PLANT **dataset-source overlap** on the SAME source cell-pairs. This cannot remove all source recording confounders.
- Hypothesis-compatible sign: **NEGATIVE**. More direct rail source paths relative to road route are associated with greater recorded introduced-flora similarity, conditional on the source controls.
- Source support: >=8 independent eligible 50km cells with both period source datasets and K80, and >=20 physically rail/road route-valid pairs. Pairwise observations share cells and are NOT treated as 20 independent plot replicates. A failed gate yields **NOT ESTIMABLE**, not a zero effect.
- Precommitted secondary sensitivities: K=20 and K=40 **on the same fixed K80 cell cohort**; original period 2010–2017; unadjusted source recording overlap; alternative v3.2 largest-rail-component source anchor. Do not choose the secondary result over the primary because its sign is favorable.
- Contradiction allowed: no effect, opposite sign or non-estimable US Southeast does NOT falsify train transport everywhere; it limits universal-source observational claims.
- No p-value: this one region, with dependent shared-cell pairs, cannot deliver independent geographic significance. A real causal seed-dispersal claim still needs temporal active freight and DOI-backed original repeated observations.

## Provenance and software recipe

Already fixed sources:
- [GBIF 10-block fixed geography and source recording opportunity](../scripts/audit_gbif_equalarea_macroblocks_v29.py).
- [Country-masked US GRIIS source species count per 50km cell](../scripts/harvest_griis_50km_10blocks_v36.py), original actual source run #42761140912.
- [Previous six-country fixed K80 test results and method](../scripts/analyze_joint_rarefaction_dataset_overlaps_v52.py), source run #42797571041, *already observed; not holdout*.

Independent new read-only source audits to run after this frozen document:
1. Overture US Southeast shortest physical RAIL route matrix via existing v3.2 implementation.
2. Overture same-cell shortest ROAD route matrix via existing v3.3 implementation, preserving <=2% missing-endpoint fail-closed source quality.
3. GBIF all-Plantae datasetKey full source marginal vs current DETECTED status in pre-existing country-GRIIS cells via v5.0 source protocol, including current API legacy-present incompatibility.
4. Joint results once these three are complete, report all outcomes independent of sign.

**Research interpretation:** This is a prospective geographically withheld source test of a posthoc-selected metric, NOT a preregistered randomized study. It advances macro replication only if honest source and geographic support gates are met; it cannot replace DOI-citable GBIF raw-event sampling correction or actual live train freight observations.
