# tetsu v8.4 — Does the US cargo–floristic relationship distinguish crop freight from resource freight?

**Actual source-verified run on November 12, 2026: [SUCCESS #43150760981](https://github.com/zuizui0223/tetsu/actions/runs/43150760981).** [Machine-readable real source output](https://github.com/zuizui0223/tetsu/actions/runs/43150760981/artifacts/13914206037), [reproducible Python model](../scripts/analyze_FAF6_crop_vs_resource_rail_floristic_v84.py), [7 source/mathematical tests](../tests/test_FAF6_crop_vs_resource_rail_floristic_v84.py). The full-repository smoke checks also passed.

## Question and fixed data universe

**Theory:** if commodity-associated propagule supply shapes recorded plant distribution, a 2022 **crop-linked** RAIL freight OD relationship should be distinguishable from a source-matched **non-crop resource** RAIL freight OD relationship, rather than simply detecting ANY rail cargo, rail hub or agricultural development gradient.

**Design:** same 49 U.S. state+DC geographic plant source units passing earlier objective GBIF record/status gates; same **1,176 dependent state pairs**; exact official 2022 FAF6.0 state-code economic OD commodity flow ZIP checksum as independent v7.9, and all 599,529 source rows, 51 source state codes, crop and resource source totals and crop OD links as independently verified v8.0; same 2025 Census TIGER state centroids/area checksum and v8.3 country-GRIIS-indexed GBIF 2018–25 plant source. Original v8.3 crop-only observed-plant K20 and raw coefficients reproduced exactly before comparing cargo groups.

The **crop-associated category** is SCTG02 cereals/grain and SCTG03 other agriculture. Comparison category SCTG10–12 resources is **not a verified seed-free or sterile freight control**; some resource freight might carry dust, soil, seeds or other propagules. All totals are 2022 **estimated economic tonnage**, not observed track segment traffic, cargo contamination or viable seeds.

**Primary response:** recorded U.S.-country-GRIIS-listed speciesKey flora ratio of expected 20-indexed-record intersection to expected union between pairs of states. Sensitivity: raw recorded-taxon Jaccard. These are positive **records**, not native–alien state-specific checklists, established presence, true absences or post-2022 introductions.

**Joint partial rank adjustment:** geographic Census state-centroid distance, within-Census region/division, Census state areas (mean and ratio), GBIF recording scale (mean and ratio), BOTH crop and resource TRUCK OD economic flows, BOTH cargo category RAIL/TRUCK node strengths, and **the other category's RAIL OD flow**. Numeric model coefficient is signed partial rho; neither OLS effect size nor causal gradient. Positive rho means greater estimated freight correlates with more similar recorded flora **under selected source conditions**.

## Actual 2022-vs-2018–25 results, same state pairs

| Economic OD comparison cohort | 2022 crop RAIL, partial rho | 2022 resource RAIL, partial rho | Crop minus resource descriptive rho | Source pair count |
|---|---:|---:|---:|---:|
| **ALL source-eligible 49-state pairs** | **−0.024108** | **−0.021684** | −0.002424 | 1,176 |
| Positive crop RAIL freight | +0.049589 | +0.052022 | −0.002433 | 141 (40 states) |
| Positive crop RAIL and TRUCK freight | +0.084335 | −0.031686 | +0.116021 | 133 (39 states) |
| **Positive RAIL cargo for BOTH crop and resource** | **+0.284174** | **−0.305962** | **+0.590136** | 41 (24 states) |

The **RAW Jaccard** corresponding crop/resource partial-rhos are:
- ALL 1,176: **−0.021295 / −0.014450**;
- crop RAIL+TRUCK positive 133: **+0.041736 / −0.056481**;
- BOTH crop+resource RAIL positive 41: **+0.212365 / −0.124166**.

The original economic crop RAIL source contains **155 positive directed interstate OD state links**; resource RAIL contains **142**, with only **24** common directed state OD links, versus 131 crop-only and 118 resource-only links. These are OD **economic** connections, NOT verified physical rail track lines, train journeys, or seed-dispersal events.

### Mandatory source-sensitivity stress tests of the apparent crop specificity

The subset of **41 shared-source economic state pairs** displays a notable cargo-category contrast under the full adjustment. Yet model specification and cohort selection matter. All 41-pair robustness cases were run, including unfavorable estimates:

| 41 matched source pairs, plant K20 response | Crop RAIL rho | Resource RAIL rho | Crop minus resource |
|---|---:|---:|---:|
| Full joint adjustment incl. rail/truck state-node strengths | **+0.284174** | **−0.305962** | +0.590136 |
| Drop all cargo node-strength proxies; retain distance, region/area, record effort and truck controls | +0.128035 | −0.176187 | +0.304222 |
| Simplified distance+GBIF record effort+both truck freight controls only | **−0.096442** | −0.139750 | **+0.043308** |
| Exclude largest 5 source pairs by combined cargo RAIL tonnage, full joint controls | +0.274569 | −0.335750 | +0.610319 |

**Delete-one-state from the 41 economic pair subset**, full adjustment only: 24 state codes were present, **23** refits remained estimable after deleting one state. In **23/23** estimable fits crop RAIL stayed positive (**+0.080881 to +0.512730**) and resource RAIL stayed negative (**−0.437328 to −0.064888**). The 24th deletion lost minimum source-pair/state support; it is NOT a failure or zero biological effect. There are not 23 independent ecological experiments.

**Crucial larger-cohort challenge:** In the **133 crop RAIL and TRUCK positive** source pairs, crop rho under full joint controls is only +0.084335, resource −0.031686. Remove cargo hub/node controls and it becomes crop **−0.011252**, resource **+0.030984**. With the simplified five-variable controls, crop **−0.048653**, resource **+0.047863**. A single observation-independent choice of which controls are appropriate cannot make the positive crop signal universal. In the FULL 1,176-pair nation-scale reference, both cargo categories are similarly weakly negative. In full-cohort leave-one-state-out, crop rho **49/49 negative** (−0.053754 to −0.009303), resource **49/49 negative** (−0.037271 to −0.007851). State deletion is not a valid dyadic p-value.

## Manuscript-safe interpretation

**Supported quantitative transport-source facts:** 2022 agricultural RAIL economic OD opportunities have very different geography from 2022 resource RAIL opportunities, and the 2022 RAIL vs TRUCK crop network is more interstate and destination-concentrated (v8.0). **Supported descriptive botanical result:** while crop and resource rail OD source associations are similarly small overall, a positive crop-versus-negative resource contrast appears in a **narrow, selected, cross-commodity active-cargo subset** under full/source-hub adjustment. It is robust to deleting any one eligible source state and top 5 shipment pairs, but is **not robust to all covariate choices and does not generalize to all original state pairs**.

**NOT established:** train-carried viable alien propagules, seed fall-off, transport-linked colonization; geographic arrival direction, freight service along actual physical railway lines; source-observer-independent GBIF flora; true local STATE-specific GRIIS alien status; causal mowing, gust, ballast heat, phenotypic plasticity or evolved stress adaptation. The 2022 economic cargo model is a single time point; 2018–25 plant source intersects it, so time direction and colonization are unidentified. SpeciesKey records share publishers and state source effort is highly uneven; partial rank adjustment does not fix this.

**Next stage:** prioritize independent dated commodity-seed carriage measurements or verified freight-track flow and independently surveyed rail versus road verges, with mowing/airflow exposure and seed/regrowth traits. A model with real SOURCE → DESTINATION new-colonization events, plus independent region/time replication, is needed before promoting this posthoc 41-pair contrast to evidence of railway seed dispersal.
