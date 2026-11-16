# tetsu v8.9 — Directional freight–plant-record associations with origin/destination state fixed effects and constrained node-label nulls

**2026-11-12 independently verified actual original source workflow [SUCCESS #43160977392](https://github.com/zuizui0223/tetsu/actions/runs/43160977392)**, [original JSON artifact #13917882520](https://github.com/zuizui0223/tetsu/actions/runs/43160977392/artifacts/13917882520), [source analysis](../scripts/validate_US_FAF6_GRIIS_dyadic_FE_spatial_QAP_v89.py), [7 source/statistical tests](../tests/test_US_FAF6_GRIIS_dyadic_FE_spatial_QAP_v89.py). The full-repository descriptive smoke checks also passed.

## Main methodology: avoid treating 2,255 direction-pairs as independent biological replicates

Retain original v8.6 audited **48 U.S. Census state source units**, **2,255 directed origin→destination source comparisons**, verified 2018–21 origin COUNTRY-GRIIS-listed GBIF taxon recorded speciesKey categories, destination 2023–25 first-recorded candidates, and original 2022 BTS FAF6 state economic SCTG02+03 RAIL/TRUCK estimated freight tonnage. The v8.6 signed raw/K20 and both mode coefficient baselines are reproduced exactly before the v8.9 model.

K20 source-record expected first-recorded candidate fraction is a GBIF indexed-source outcome, **not real immigration/arrival**. A candidate taxon recorded earlier at origin but not destination is NOT proven biologically absent from destination. Country-level GRIIS status is not state-specific alien status. The 2022 freight economic OD is not observed train seeds or traffic along a physical railway.

**New model:** partially rank residualize recorded response and log1p RAIL/other-mode TRUCK economic OD against:
- 47 origin state source indicator fixed effects (first Census FIPS code reference), 47 destination state indicators, and intercept;
- other 2022 crop FREIGHT mode's same directed OD estimated tonnage;
- Census centroid log geographic separation, same Census region/division;
- log source taxon pool size in that pair.

The origin/destination state FE absorb many broad source differences (publisher/observer opportunities, prior plant population, agricultural/freight hub differences) without identifying their underlying causes. The source data remain dependent by state. FE can also remove real region-level freight/plant effects; the estimand becomes a *within-state-pair residual pattern across destinations/origins*, not the earlier pooled v8.6 rho.

**Constrained network-label structural null, NOT a calibrated causal p-value:** simultaneously relabel the SOURCE 2022 RAIL or TRUCK economic origin and destination state codes, preserving the original directed economic network's topology and exact freight weight multiset but reassigning the whole matrix to the geography/plant record sources. The primary labels are reassigned **within the same Census region (499 draws)**, followed by within Census division (299 draws) and unrestricted US-state code relabeling (199, less geographically credible). Controls and botanical recorded response remain unchanged; the source network is reranked and re-residualized for each. Fixed seed 20261112.

**Critical caveat:** exchangeability of historically determined state labels is NOT verified even under Census region constraints. Node-label QAP does not automatically solve network autocorrelation or known confounders. Therefore the null exceedance fractions below are **SOURCE REFERENCE DIAGNOSTICS**, not valid randomized-treatment p-values or 95% ecological confidence intervals. See Dekker, Krackhardt & Snijders (2007), *Psychometrika* 72:563–581, https://doi.org/10.1007/s11336-007-9016-1; and the 2022 network permutation nonindependence critique: https://pmc.ncbi.nlm.nih.gov/articles/PMC9617964/.

## Verified genuine source results

| Original mode, separately controlling opposite mode and state source FEs | K20 recorded taxon source FE residual rho | Source 499× Census-region constrained null 2.5–97.5% quantiles | Proportion region-constrained null absolute rho >= abs(observed) |
|---|---:|---:|---:|
| **RAIL** crop economic OD | **−0.014364** | **−0.048405..+0.026656** | **0.517034** |
| **TRUCK** crop economic OD | **+0.029635** | **−0.057779..+0.028981** | **0.244489** |

Additional source-definition robustness: raw *non-rarefied* source first-recorded taxon fraction fixed-effect partial rho **RAIL +0.011452** and **TRUCK +0.076036**, showing response metric dependency.

For **RAIL K20**:
- Census REGION conserved state labels: 499 trials, reference median **−0.008804**, source interval −0.048405..+0.026656.
- Census DIVISION conserved state labels: 299 trials, median **−0.026929**, source interval −0.061591..+0.011910.
- UNRESTRICTED state labels: 199 trials, median **+0.002444**, source interval −0.037147..+0.045687.

The actual rail result lies well within the central region- and division-stratified structural source null references. Trucks show a small positive source rho; its reference null is asymmetrically centered, so although the one-sided positive exceedance is ~0.022 under region-conserved labels, the two-sided absolute exceedance is ~0.244, and no significance or causal truck-specific dispersal interpretation is warranted.

**Decision:** there is no separately identifiable *rail* source-level positive association with destination 2023–25 first GBIF-recorded GRIIS-country-listed taxa after controlling origin/destination state source baselines. This reinforces earlier v8.3–v8.8 lack of robust positive railway seed-delivery signal. It does **not** establish that trains do not disperse seeds.

## What remains unresolved and next source audit

- 2022 estimated economic FAF state OD is not actual operating railway/track segment freight. State origin/destination direction cannot verify train-carried viable propagules or their deposition.
- 2018–21 vs 2023–25 GBIF recorded speciesKey categories are subject to nonrandom publisher, database backfilling, species synonym and taxonomy status, survey frequency and 4-year vs 3-year window biases; recorded plant *absence* is not verified.
- GBIF stateProvince is a text query. Prior v8.2–v8.8 analysis has NOT independently verified point records inside official Census state polygons. A separate ORIGINAL occurrence-lat/lon source audit v9.0 is therefore the correct next gate before stronger botanical interpretation.
- No real rail or road verge site plant plot monitoring, train/wind/mowing exposure, species regrowth/aerodynamics, phenotypic plasticity or heritable genetic adaptation is observed. All ecological mechanism hypotheses remain unverified.

Keep PR DRAFT and retain negative and non-estimable results. The supported macro result remains the 2022 FAF6 agricultural freight network's concentration and mode/directional structure, not a biological seed-distribution effect.
