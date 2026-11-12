# tetsu v8.5–v8.8 — Genuine 2022 U.S. crop-freight directions versus earlier-origin and later-destination GRIIS-listed plant records

**Real, source-verified 2026-11-12; ALL original GitHub Actions passed.** This is ecological source-bias-aware MACRO exploration, NOT a causal demonstration of train-borne alien seed dispersal or trackside vegetation adaptation.

## Independently verified runs and inputs

1. [v8.5 temporal original GBIF indexed SOURCE #43153476962 SUCCESS](https://github.com/zuizui0223/tetsu/actions/runs/43153476962), [original 51-state source artifact #13915120892](https://github.com/zuizui0223/tetsu/actions/runs/43153476962/artifacts/13915120892). Current 2026 GBIF stateProvince speciesKey facets for 2018–2021 and 2023–2025, omitting ALL 2022 botanical records, using frozen Census 51-state/DC FIPS names, 2023-vintage U.S. country GRIIS-listed accepted species keys, ≥100 all-Plantae indexed records, ≥3 dataset UUID categories, ≥5 country-listed GRIIS taxa / ≥20 indexed GRIIS records, and ≤2% non-DETECTED/unknown status in each state/period. 51×2 original GBIF source period queries completed, no source errors. Early source quality 50/51, late 48/51, **48/51 BOTH qualifying**. Three were source-ineligible due to metadata quality, NOT botanical absence. **1,782 distinct USA-country-GRIIS speciesKey categories were recorded in 2023–25 but not in 2018–21 in at least one of those source states. NOT 1,782 true invasions.**
2. [v8.6 actual 2022 directed FAF RAIL and TRUCK crop OD versus source-year-separated GRIIS records #43154108559 SUCCESS](https://github.com/zuizui0223/tetsu/actions/runs/43154108559), [source result artifact #13914753197](https://github.com/zuizui0223/tetsu/actions/runs/43154108559/artifacts/13914753197), [code](../scripts/analyze_US_FAF6_directional_source_taxa_v86.py). Exactly **48 US states** in both GBIF bins, **2,255 distinct directional origin→destination pairs** (one eligible state pair below 10 candidate source-taxon categories), and **577,846 dependent origin-species category opportunity comparisons**. None are independent ecological replicates.
3. [v8.7 real SAME STATE/SAME PAIR indexed record rarefaction K10/20/40/80/160 and reversed freight direction #43154743438 SUCCESS](https://github.com/zuizui0223/tetsu/actions/runs/43154743438), [real JSON artifact #13914704313](https://github.com/zuizui0223/tetsu/actions/runs/43154743438/artifacts/13914704313), [code](../scripts/analyze_US_FAF6_temporal_rarefaction_sensitivity_v87.py). Original v8.6 K20 crop RAIL/TRUCK source coefficients reproduced **exactly** in every cohort before varying indexed-record opportunity.
4. [v8.8 real 2022 FAF6 crop SCTG02+03 RAIL versus resource SCTG10–12 RAIL same-time/direction origin→destination recorded taxa #43155204936 SUCCESS](https://github.com/zuizui0223/tetsu/actions/runs/43155204936), [actual result #13916185209](https://github.com/zuizui0223/tetsu/actions/runs/43155204936/artifacts/13916185209), [code](../scripts/analyze_US_FAF6_directional_crop_resource_temporal_v88.py). Official FAF6 published 599,529 state OD original rows, ZIP SHA256, 51-state FIPS, crop + resource commodity-mode marginal tonnages, original signed crop OD and earlier v8.6 source coefficients were independently verified, and **6/6 source integrity tests passed**. All framework smoke checks remained green.

### Temporal grain and observational meaning

Original U.S. BTS FAF6 2022 agricultural economic origin→destination shipment estimates, not observed seed loads or exact operated railway line routes, are used with source-year-separated GBIF indexed plant records:
- Origin state taxon is **recorded in 2018–2021** according to GBIF country-GRIIS speciesKey. Destination state taxon is **not recorded in 2018–2021** but **recorded in 2023–2025** in the same current 2026 GBIF source index.
- Every origin→destination pair first fixes the list of origin taxon keys absent from destination EARLY *positive records*. This is a **recording-defined candidate pool**, NOT a list of taxa verified absent from the destination.
- Response v8.6 K20: each candidate taxon's probability of appearing in 20 source-indexed post-2022 destination GBIF records, weighted by origin taxon's source-indexed probability of appearing in 20 earlier origin records; normalized over the original candidate pool. **This is not a seed movement probability, known colonization probability or E[ecological occupancy].**
- RAW fraction of candidate source taxa positive in destination's late GBIF indexed records is a mandatory outcome-definition sensitivity, and K10/40/80/160 vary only GBIF source-record inclusion depth on identical 48-state/2,255-pair geography.
- Target "2022" is an economic freight-base year strictly excluded from both plant record YEAR bins. This alone **does not create a quasi-experiment**: freight could have operated for decades and many older plant records may have been uploaded retrospectively. Pre period is four years, post three, and datasets/observers vary.

Controls in v8.6: same-direction 2022 TRUCK crop OD; Census geodesic distance and region/division; origin pre, destination pre and destination post GBIF all-plant source record counts; origin candidate recorded-taxon count; origin/destination land areas; four agricultural crop freight node outgoing/incoming strengths. These are source proxies only, NOT guaranteed climate, ecoregion, publisher, agriculture or land use adjustment. No dyadic p-values.

## Core v8.6 actual directional contrast: plant source records vs 2022 crop tonnage

Positive signed partial rho would mean larger **estimated RAIL crop shipment** coincides with a higher recorded origin-associated destination after-source fraction, conditional on the recorded source controls; it would STILL NOT prove actual train-carrying seeds.

| Original source cohort | Number of directional state OD pairs sharing states | K20 expected source first-record fraction vs crop RAIL partial rho | Same fraction vs crop TRUCK partial rho | RAW source first-record fraction vs crop RAIL partial rho |
|---|---:|---:|---:|---:|
| ALL 48 source-quality states | **2,255** | **−0.027935** | +0.037666 | +0.010019 |
| Crop RAIL original 2022 OD positive | **148** (30 origin, 35 destination states) | **−0.024113** | +0.101304 | +0.026230 |
| Crop RAIL and TRUCK both OD positive | **123** (28 origin, 33 destination states) | **−0.091896** | +0.114528 | −0.081802 |

The RAW source definition sometimes changes the sign, so even though the K20 result is reproducible, there is **no cross-source-definition positive railway transport effect**. Under v8.6, leave ONE source origin or destination state out of the ALL directed pairs and refit: RAIL K20 stayed negative in **48/48 origin removals** (rho **−0.035633..−0.018823**) and **48/48 destination removals** (rho **−0.036888..−0.017216**). These are within-US dependence sensitivity checks, NOT 96 country-level independent experiments or a hypothesis test.

## v8.7: mandatory fixed-source rarefaction depths and reversed OD freight direction

| GBIF source record rarefaction K | ALL 2,255 directed pairs: RAIL crop rho | 148 positive-Rail OD pairs: RAIL crop rho |
|---|---:|---:|
| 10 | −0.027849 | −0.015920 |
| **20** | **−0.027935** | **−0.024113** |
| 40 | −0.028532 | −0.035840 |
| 80 | −0.028944 | −0.052914 |
| 160 | −0.030029 | −0.061212 |

All fixed cohort K-depths preserve the **negative** RAIL correlation, though source count standardization is not observer/dataset bias correction. Reverse the **economic freight OD predictor direction** on identical plant source-response/covariate rows, without changing original positive-rail cohort: reverse-rail rho **−0.015982** for all pairs and **+0.113284** for the original positive-RAIL 148 pairs, versus true economic OD -0.027935/-0.024113. A reverse-direction transport predictor is a **direction diagnostic, not a verified placeboregression or independent train experiment**. It can correlate due to freight centrality/geography and nonrandom source records.

## v8.8: cargo-specific directional temporal model (crop vs resource railway)

Same original 48 source-quality states, 2,255 pairs, K20 taxon-record response; re-fetch identical official 2022 FAF6 source ZIP and retain separate directional estimated RAIL/TRUCK tonnages and origin/destination cargo-strength covariates. Each cargo RAIL coefficient is adjusted for the OTHER cargo's RAIL OD, BOTH TRUCK cargo OD, exact prior Census/GBIF source geography and BOTH commodity types' hub/source flow.

| Botanical source-recorded fraction and direction cohort | Original directional pairs | Crop RAIL rho | Resource RAIL rho |
|---|---:|---:|---:|
| **ALL 48-state directions** | **2,255** | **−0.017047** | −0.010768 |
| Crop RAIL OD positive | **148** | −0.046647 | −0.065249 |
| Crop RAIL + crop TRUCK both positive | **123** | −0.090519 | −0.018327 |
| **BOTH crop RAIL and resource RAIL positive** | **24** (9 origins, 14 destinations) | **NOT ESTIMABLE** | **NOT ESTIMABLE** |

The last source cohort fails prespecified **≥20 distinct origin states and destination states and ≥30 directional pair** gates. Do not relax this to reproduce v8.4's positive 41 **undirected** state-pair crop-specificity signal. The v8.4 selected-case positive crop–resource correlation was likewise not general on the full source data and depended on model controls.

Within the full directional v8.8 source, remove one state source: crop RAIL full coefficient remained negative in **48/48** origin-removal refits (−0.024813..−0.006680) and **48/48** destination refits (−0.027818..−0.006552). Resource RAIL full coefficient remained negative in 48/48 origin and 47/48 destination refits, one +0.000503. These numerical deletion checks are NOT p-values.

## Scientific conclusions — evidence ladder

**Supported robust MACRO TRANSPORT-STRUCTURE RESULT:** independently audited FAF6 economic source reports show domestic US agricultural RAIL freight strongly concentrates on state-interstate origin→destination paths (63.43% interstate among RAIL agricultural tonnage versus 17.76% TRUCK; 155 directed interstate economic rail OD links versus 1,347 truck) and major destination concentration. This is source-verified macro physical/economic movement opportunity, **not viable propagules**.

**Supported but negative/limited observational biology:** with recorded year separation and a directional source-pool proxy, RAIL agricultural cargo does not have a generally positive association with *source-indexed newly recorded* country-GRIIS taxa at destinations. Negative partial rank is stable for K10–K160 and deleted source states, but not a causal estimate of negative biological dispersal (or proof that freight never transports seeds). Weak positive TRUCK correlations are observation associations, NOT evidence trucks caused introductions. Noncrop RAIL gives similar small values nationwide, so **no durable crop-specific botanical signal**.

**Not established:** viable alien seeds riding freight cars or agricultural cargo, shed by trains, deposited/germinating along known operated track routes; state-specific invasive species colonization or arrival time; direct railway mowing, herbicide, mechanical slipstream gust, thermal/soil stress, species fitness plasticity or local genetic adaptation.

**Prioritize NEXT rather than another posthoc index search:** independently documented active freight corridors and commodity/cargo seed-contamination or railway-verge seed traps; source-authenticated state/rail-yard vegetation records with observer/provider-fixed dates; repeated rail vs road verge vegetation and measured mowing/gust regimes with true species regeneration/aerodynamics; multiple independent geographic and temporal replicates. The current state-year GBIF provenance/field-survey and FAF line-routing gap blocks identifiable transport × establishment/filtering effects.

### Original sources
- U.S. BTS FAF6.0: 2022 state/commodity/transport-mode economic OD estimates, [BTS](https://www.bts.gov/faf).
- U.S. Census 50 states + DC FIPS state code table, 2025 TIGER source centroid and land area metadata.
- GBIF 2026-11-indexed speciesKey recording-period source facets, original 2018–21 and 2023–25 record YEAR bins.
- Country GRIIS U.S. accepted speciesKey list; U.S.-country GRIIS status does not imply state introduction.
