# tetsu v9.1 — New recorded plant taxa coincide with a nationwide GBIF reporting shift

**Verified real 2026-11-12 source run:** [SUCCESS #43162581367](https://github.com/zuizui0223/tetsu/actions/runs/43162581367), [original 48-state result #13918307946](https://github.com/zuizui0223/tetsu/actions/runs/43162581367/artifacts/13918307946), [source code](../scripts/audit_US_GRIIS_firstrecord_observation_effort_v91.py), [7/7 unit tests](../tests/test_US_GRIIS_firstrecord_observation_effort_v91.py). Full repository smoke QA also PASSED.

## Key methodological problem

Previous 2022 freight × plant first-recorded analysis uses GBIF 2026 indexed occurrence event years 2018–2021 BEFORE (four years) versus 2023–2025 AFTER (three years), explicitly excluding 2022. A plant speciesKey not recorded in the early source window but indexed in the later source window is NOT necessarily a biological arrival or establishment; missing-source taxa are not confirmed absent.

This audit uses precisely the SAME 48 Census FIPS U.S. states with high-quality data in BOTH time periods, with the original valid country GRIIS checklist and original bounded GBIF source data; the three other state sources are excluded, not imputed zero. Annualization correctly divides source counts by four (early) versus three (late) event-year windows.

## Actual observed GBIF source process

| Original 48-state GBIF indexed source annual-record-rate post/early ratio | Minimum | **Median** | Maximum | States with ratio >1 |
|---|---:|---:|---:|---:|
| ALL Plantae speciesKey-indexed records/year | **1.166847** | **2.142561** | 3.842083 | **48/48** |
| USA-country-GRIIS-listed speciesKey-indexed records/year | **1.116180** | **2.063954** | 3.759035 | **48/48** |

**46/48** states show MORE late-window-only country-GRIIS-listed recorded speciesKey categories than early-only categories; just 2 show the reverse. The median source state has **122.5** late-only country-GRIIS indexed taxon categories, and the median fraction of the later-window recorded GRIIS taxa that were not in earlier records is **0.2455**. NONE of these numbers estimates true species immigrations.

Rank correlation of the late-only recorded source taxon-category count versus:
- Actual 2023–25 total ALL-Plantae indexed records: **rho +0.714282**.
- Actual 2023–25 recorded GRIIS country-listed taxon-category richness: **rho +0.847912**.
- Annualized ratio of late/early ALL-Plantae indexed records: **rho −0.035942**.
- Annualized ratio of late/early GRIIS-country-listed indexed records: **rho −0.094905**.

Late-only fraction divided by later reported GRIIS richness has different source associations: total later ALL-Plantae records rho **−0.751194**, annualized recording-rate ratio rho **+0.285714**. This is a denominator/source-composition warning, not a biological mechanism.

## Interpreting the results safely

The increase in indexed annual recording scale is **real as an observational-data process** and occurs in every source-qualified state, even accounting for different numbers of years. Absolute late recording volume relates strongly to the number of later-only source taxa, but their increase-RATIO does NOT: survey/publisher growth cannot yet be identified as the sole cause of the first-recorded categories. Species richness, climate, state area, provider/observer turnover, original GBIF taxonomic changes and retrospective upload all remain uncontrolled.

This reinforces prior v8.6–v8.9 evidence grading: source-recorded later-only taxon categories and small freight association coefficients are NOT actual 2022 train-carried propagules, true new state invasive plants, railway-verge seedlings, mowing or wind filtering, plasticity or local evolutionary adaptation. Only the macro 2022 FAF6 commodity network's geography and freight modes are directly supported.

## Next source integrity gate

GBIF stateProvince is an interpreted TEXT label, not necessarily coordinates in the named state. Independently compare ORIGINAL occurrence lat/lon to official Census 2025 state polygons (v9.0) before turning any state-indexed speciesKey source matrix into substantive biogeography. A bounded nonrandom API offset probe cannot estimate full-state mismatch prevalence without a probability-based original-record sample, and must not be described as such.

No original freight cargo seed contamination/track route use, repeated near-track vs roadside vegetation, actual grass-cutting/gust exposure or genetic adaptation has yet been measured.
