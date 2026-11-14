# tetsu v9.5–v9.7 — Original iNaturalist observer/site continuity and source "first-recorded" botanical categories

**Verified original-source successful GitHub Actions, 2026-11-12:**
- [v9.5 original GBIF observer + 10/50km location source #43203802514](https://github.com/zuizui0223/tetsu/actions/runs/43203802514), [actual JSON artifact #13936058989](https://github.com/zuizui0223/tetsu/actions/runs/43203802514/artifacts/13936058989), [code](../scripts/audit_US48_iNat_observer_site_repeat_v95.py), 6 tests.
- [v9.6 necessary K20 same observer+site GRIIS source bound #43204765751](https://github.com/zuizui0223/tetsu/actions/runs/43204765751), [actual result #13935729819](https://github.com/zuizui0223/tetsu/actions/runs/43204765751/artifacts/13935729819), [code](../scripts/audit_US48_observer_site_strict_GRIIS_K20_gate_v96.py), 6 tests.
- [v9.7 ALL 7×2 observer/site and source first-record source correlations #43205672117](https://github.com/zuizui0223/tetsu/actions/runs/43205672117), [actual result #13936913219](https://github.com/zuizui0223/tetsu/actions/runs/43205672117/artifacts/13936913219), [code](../scripts/analyze_US48_GBIF_firstrecord_observer_site_v97.py), 6 tests. Project smoke validation also passed.

## Original data and limits

The GBIF iNaturalist Research-grade Observations dataset UUID is 50c9509d-22c7-4a22-a47d-8c48425ef4a7. Previous verified v9.3 showed that ONE dataset comprises **95.8247%** of all-Plantae indexed records in 2018–21 and **99.1155%** in 2023–25, median across the same 48 qualified Census-FIPS state source units. A SAME dataset UUID cannot imply SAME observer, survey plot or sampling effort.

For observer/site source QA, v9.5 fetched REAL GBIF ORIGINAL occurrence/search rows, 2 deterministic source result offsets (0 and 1,000) ×100 occurrence records each per U.S. state and year bin:
- **9,600 actual source rows, 2018–2021** (4 years).
- **9,600 actual source rows, 2023–2025** (3 years).
- **19,200 original occurrence rows across 48 states**. Original GBIF event year, scoped iNaturalist dataset UUID and GBIF record IDs are checked for source correctness. Individual observer names, recordedByID, gbifID, original longitude/latitude and personally identifying data are NOT archived.
- Observer identity is normalized from recordedByID if available, otherwise recordedBy. Coordinates are temporarily converted to projected EPSG:5070 and coarse 10km or 50km source-grid bins. Only state-level deidentified counts persist.

**FUNDAMENTAL SAMPLING LIMIT:** these are **deterministic NONRANDOM SEARCH RESULTS**, not probability samples from every GBIF/iNaturalist observation or observer. Missing a previously sampled observer in 200 late records does NOT mean that observer ceased participating. Same observer and same coarse 10km/50km grid do NOT verify a repeated standard botanical plot or equivalent plant detection.

### Actual original iNaturalist source recurrence, v9.5

| State-level original source statistic, 48 states | Verified value |
|---|---:|
| Median distinct sampled observer labels BEFORE | **73** |
| Median distinct sampled observer labels AFTER | **64.5** |
| States with at least one observer label in BOTH bounded samples | **47/48** |
| Median fraction of post source records whose observer label appeared in previous sample | **13%** |
| States with any repeated 10km original geographic cell | **48/48** |
| States with any repeated observer-label × 10km site | **43/48** |
| Median fraction of post records at 50km coarse tile already sampled previously | **88%** |
| Median post source-record count at same observer AND 10km coarse site | **4.5 / 200** |
| Median post source-record count at same observer AND 50km coarse site | **10 / 200** |

The 13% is a fraction within the selected **200-record** source sample, NOT true iNaturalist-wide observer retention. Exact observer labels, user IDs, observation GBIF IDs and private locations remain outside all archived artifacts.

### v9.6 K20 strict repeated observer + site source-support upper bound

For a LATE country-GRIIS-listed GRIIS plant source record to match the **same observer AND the same geographical tile** in EARLY source records, its count can never exceed either:
1. all LATE source records with repeated observer AND tile; OR
2. LATE GRIIS-country-listed source records from tiles also sampled EARLY.

Thus MIN(count1,count2) is a mathematical **necessary upper bound**, NOT an observed matched-observer GRIIS taxon count. At the already established K=20 plant-record rarefaction and >=40-state macro-eligibility threshold:

| Original GBIF coarse source grid | Median LATE source records at same observer×tile | States with >=20 all-Plantae late same-observer-tile records | States with >=20 GRIIS source records matching observer+tile EVEN POSSIBLE under upper bound |
|---|---:|---:|---:|
| 10km | **4.5** | 7/48 | **at most 3/48** |
| 50km | **10** | 16/48 | **at most 14/48** |

No observer-matched 40-state K20 FAF–plant model was fitted because the **necessary** source sample support failed. This does NOT show that all GBIF/iNaturalist source data lack repeated observers or that railway seed transport never happens. It only shows these original **nonrandom bounded source pages** cannot provide independently repeated observer+tile GRIIS samples.

## v9.7 observational-source association, explicitly POSTHOC

The complete 48-state original v9.5 aggregate was matched to independently verified v9.1 source records. Source outcomes:
1. Post-2022 country-GRIIS speciesKey category FRACTION that was not previously source-recorded in the same state (NOT new biological arrivals).
2. Alternative raw COUNT of post-only recorded country-GRIIS categories (NOT immigration count).

**All seven** coarse source observer/site continuity metrics were regressed descriptively against **both** responses; nothing was selected by significant p-values, and the exercise is labeled **POSTHOC exploratory**. Adjustment uses a partial rank regression controlling log 2023–25 all-Plantae indexed GBIF source record volume and the annualized 2023–25/2018–21 all-Plantae source record growth ratio. These controls do NOT correct all site/observer/provider biases. State-level source records have common provider/geography dependence, so no p-values or 48 independent country replicates.

| Observer and source site repeat metric | Raw rank rho vs post-only recorded taxon FRACTION | Adjusted partial rank rho vs post-only FRACTION | Adjusted partial rank rho vs absolute post-only COUNT |
|---|---:|---:|---:|
| Observer label repeated across source-period samples | **−0.319985** | **−0.330067** | +0.081869 |
| 10km geographic tile repeated, observer unrestricted | +0.287266 | −0.145850 | +0.107787 |
| 50km geographic tile repeated, observer unrestricted | +0.171654 | **+0.010609** | **+0.353371** |
| SAME observer and SAME 10km tile | −0.201460 | **−0.203583** | +0.049185 |
| SAME observer and SAME 50km tile | **−0.347158** | **−0.415382** | +0.072024 |
| Source GRIIS record occurred in previously sampled 10km tile | +0.093312 | −0.218233 | +0.081917 |
| Source GRIIS record occurred in previously sampled 50km tile | −0.106167 | −0.261386 | +0.176880 |

**Sensitivity**: after controlling source recording effort, SAME observer×50km tile v9.7 correlation with post-only recorded source speciesKey FRACTION remains NEGATIVE in **48/48** leave-one-Census-state-out descriptive refits (rho **−0.456612 to −0.347959**). SAME observer alone stays negative across all 48 refits (rho −0.386294 to −0.252556). The exact same observer×50km predictor has **rho +0.072024** against the alternative absolute LATE-ONLY taxon COUNT: the signed observation depends on denominator. No independent ecological causal inference.

### Scientific interpretation

**Supported descriptive recording-source finding:** despite the fixed GBIF iNaturalist dataset UUID, observer label/recording site composition in the original bounded sampled records differs across the two temporal bins; state-level observer×site recurrence tends to be higher where source-indexed later-only country-GRIIS taxon FRACTION is lower. This is consistent with observation-channel heterogeneity as an explanatory mechanism, but it does **not** prove that the observers or their turnover CAUSED apparent first records, let alone true biological immigration.

**What is not supported:** a biologically identified causal railway seed vector effect from FAF6 commodity OD, actual viable alien propagules carried by trains, active historical freight line flows, plot-confirmed rail-yard/trackside colonization, mowing/gust/stress adaptation, plant morphological plasticity or genetic evolution. Earlier v8.6–v9.4 national freight associations did not yield a durable positive railway first-recorded source signal, including when iNaturalist is excluded.

**Decision / next scientific identification gate:** avoid fitting new 40-state K20 biological transport models to this bounded sample. Obtain independently documented repeated same-observer/same-plot GRIIS plant source field surveys; true before-after site absence controls and source DOI/observer/provider continuity; actual freight car/seed trap contamination and dated RAIL traffic; matched rail/road verge field vegetation with recorded mowing/gust exposure and independent measured seed/regrowth traits. Maintain negative and unestimable results in Draft PR.
