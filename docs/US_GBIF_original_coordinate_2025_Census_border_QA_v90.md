# tetsu v9.0 — Original GBIF stateProvince source geography QA against independent Census 2025 state boundaries

**Actual 2026-11-12 source run SUCCESS [#43162798791](https://github.com/zuizui0223/tetsu/actions/runs/43162798791)**, [aggregated source-verified JSON artifact #13918542955](https://github.com/zuizui0223/tetsu/actions/runs/43162798791/artifacts/13918542955), [reproducible original-coordinate quality audit](../scripts/audit_US_GBIF_stateProvince_original_coords_v90.py), [6 tests](../tests/test_US_GBIF_stateProvince_original_coords_v90.py). Full repository smoke QA also passed.

## Why it matters and exactly what was checked

Previous v8.2–v8.9 U.S. freight×plant source analyses use the GBIF 2026 search API with a stateProvince TEXT name filter. A source state name may not place the individual original occurrence inside that state, and invalid geography would compromise freight/pair alignment even before species biology is considered.

Use the **eight pre-frozen v8.1 Census FIPS states, selected independently of freight/plant effects**: California (06), Florida (12), Illinois (17), Iowa (19), New York (36), North Carolina (37), Texas (48), Washington (53). Two disjoint actual original GBIF event YEAR windows: 2018–2021 before the FAF6 2022 economic base-year, and 2023–2025 after, excluding 2022.

- Download original U.S. Census 2025 1:500k STATE polygon shapefile [cb_2025_us_state_500k.zip](https://www2.census.gov/geo/tiger/GENZ2025/shp/cb_2025_us_state_500k.zip), exact SHA256 in source output. Independently match shapefile FIPS/NAME/STUSPS to Census 51-state lookup before geographic audit. This is cartographic (generalized) boundary geometry, not a high-resolution land habitat map.
- Sample authentic GBIF occurrence/search **ORIGINAL coordinate-bearing Plantae species-rank positive record source** within each stateProvince and year bin. Source query is hasCoordinate=true / hasGeospatialIssue=false. Two deterministic nonrandom API offsets **0 and 1,000**, limit 100 each => 200 records × 8 states × 2 periods = **3,200 actual original GBIF occurrence rows**.
- Compare each original decimal latitude/longitude to the named state's WGS84 Census polygon, and use 5,000 m distance from projected Census shape for border/coastal generalization tolerance. Do NOT export precise coordinates, IDs, or individual observation records.
- Assess separately all Plantae and the authentic record subset whose GBIF speciesKey belongs to the U.S.-COUNTRY-level original GRIIS introduced checklist. The latter is NOT a confirmed state-alien classification.
- **This is not probability sampling of all GBIF state records.** In particular, only source result offsets 0 and 1,000 were inspected, so original records from beyond the first ~1,100 of each query population were unobserved. Provider and result-order bias may be strong; do not extrapolate percent matched to the whole United States.

## Original source verification outcomes

| Source-record category | Real GBIF original rows checked | Inside state's polygon or within 5 km boundary tolerance | More than 5 km outside named state | Far-outside fraction of PROBED rows only |
|---|---:|---:|---:|---:|
| **ALL Plantae records** | **3,200** | **3,176 (99.25%)** | 24 | **0.750%** |
| **US-COUNTRY-GRIIS-listed speciesKey records** | **716** | **714 (99.72%)** | 2 | **0.2793%** |

All **8/8** preselected Census states had source queries completed in both temporal windows, zero source request errors. All actual GBIF ID rows were deduplicated within source probes; no original coordinates or IDs are archived in the resulting artifact.

**Geographical heterogeneity in the original ALL-Plantae probes (200 GBIF records per state/year window):**
- **Illinois:** source coordinates more than 5 km outside Illinois polygon in **8/200 (4.0%) before** and **9/200 (4.5%) after**. Of the two original GRIIS-listed-probe records more than 5km outside any inspected state, both occurred in the Illinois BEFORE source subset (2 of its 32 GRIIS-listed probed original rows); not all Illinois source GRIIS records were misassigned.
- **Florida:** 0/200 BEFORE, 3/200 AFTER (>5km outside).
- **California:** 1/200 BEFORE and 1/200 AFTER.
- **Iowa:** 1/200 BEFORE, 0 AFTER.
- **New York:** 1/200 BEFORE, 0 AFTER.
- **North Carolina, Texas and Washington:** 0/200 far-outside in both windows.

### What this result supports

This supports a **bounded source-quality observation:** the named-state GBIF coordinate locations are highly concordant with official Census state boundaries *among the inspected original coordinate-bearing rows*, including GRIIS-checklist matches. StateProvince text is not automatically invalid as a geographic join. However the nonrandom, shallow-offset sample is insufficient to estimate the quality of all state source records, especially records beyond the first 1,100 API results, other U.S. states, records without coordinates and specimens concentrated among distinct publishing institutions.

**Limits:** a correct state location does NOT imply the taxon was truly INTRODUCED to that state, established there after 2022, absent beforehand, or moved by train. The GRIIS checklist is country-level, and year-bin GBIF indexed positive records are not independently designed field surveys. This source audit does not resolve provider/observer, taxonomy/backfill and opportunity confounding (see [v9.1 report](US_GBIF_2018_21_2023_25_reporting_effort_v91.md), which shows post-2022 annual recording effort >pre in 48/48 states), nor actual freight track use, viable seeds, rail-versus-road mowing/airflow, phenotypic plasticity or evolutionary adaptation.

**Priority after this audit:** source-provider-held-out original GBIF events with true geographical polygon validation across the geographic analysis universe, then actual cargo seed contamination/rail-yard surveys and independent wind/mowing vegetation plots. Do not treat a high bounded coordinate-matching fraction as proof that the causal railway plant-dispersal hypothesis is supported.
