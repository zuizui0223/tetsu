# tetsu v3.6 — 50km macro railway vs road *continuous* accessibility and observed introduced plants

**Status: actual 8-country source-backed work completed on 2026-11-10.** The central biological railway-dispersal hypothesis remains **unconfirmed**. No pair-level p-values are justified.

## Research question

Across genuinely broad, independent geographic regions, does a more DIRECT RAILWAY physical path relative to a ROAD path correspond to increased similarity of country-listed introduced plant assemblages, after matching/control for straight-line geography and GBIF recording opportunity?

This analysis moves beyond the earlier saturated "is this railway connected?" binary treatment. It uses continuously varying measured source physical network paths in geospatially identical 50km equal-area source cells. Physical route distance is **not** railway freight operation, active service, vehicle constraints, actual movement of seeds, observed seed establishment, or plant adaptation.

## Independent real-source baseline

- [v3.2 actual Overture RAIL physical shortest paths, CI run #42759357722](https://github.com/zuizui0223/tetsu/actions/runs/42759357722) — success, 5 independent unit tests, three predeclared source 200km blocks, 50km EPSG:6933 2500km2 nominal cells. Road comparison was not used to select rail source anchors. Two SOURCE-DEFINED anchor choices were both computed: nearest connector inside the cell to its centroid and nearest connector of the largest component. Mean traffic and actual active passenger/freight lines **unknown**.
- [v3.3 actual Overture ROAD physical shortest paths, CI run #42759664602](https://github.com/zuizui0223/tetsu/actions/runs/42759664602) — success, 6 tests, US Midwest + France Burgundy identical blocks and 75km physical graph halo, 846,408 and 1,028,680 source road segments. Both regions have all 120/120 50km cell pairs routed in the local mapped ROAD graph. Physical shortest route is approximate; driving legality, one-way/turn restrictions, bridges/closures and traffic are unmodeled.
- [v3.1 real introduced plant speciesKey source #42740320597](https://github.com/zuizui0223/tetsu/actions/runs/42740320597) — 15/16 US, 16/16 FR and 1/16 JP 50km cells met the BOTH-period source GRIIS-introduced species gate, given earlier plant source-record minimums. These are positively recorded GBIF species-key category counts, not certified absences, individual plants or original event row downloads.
- [v3.5 independent Natural Earth country polygon fraction #42760619186](https://github.com/zuizui0223/tetsu/actions/runs/42760619186) — 10 prior blocks, all 160 cells cross-checked against public-domain 10m admin-0 country polygons, exact source SHA in archived artifact. **132/160** cells have at least 80% area in their source country; together with two-period GBIF minimum source coverage, **121/160** remain (versus 137/160 record-only eligible). US Midwest retains **15/16**, FR Burgundy **16/16** such source cells. Natural Earth country area is NOT a validated terrestrial ecoregion/land-cover mask.
- [v3.4 real continuous-source GRIIS plant vs RAIL/ROAD association #42761289111](https://github.com/zuizui0223/tetsu/actions/runs/42761289111) — success, six strict synthetic tests, all route-anchor and response normalization sensitivities reported. Regionally independent replicate count is **two**, not 91/120 (these are pairwise comparisons sharing cells).

### Actual physical network geometry (not plant ecology)

| 50km study block | RAIL source: valid pairs | Median rail-centroid-route / geographical great-circle | ROAD source: valid pairs | Median road-centroid-route / geographical great-circle |
|---|---:|---:|---:|---:|
| Japan Tohoku | 105 | 1.719974 | not sourced in this pipeline | unknown |
| US Midwest | 105 | 1.950427 | 120 | 1.186433 |
| France Burgundy | 120 | 1.685462 | 120 | 1.196970 |

**Do not compare rail and road medians as if paired identically.** Within the species-qualified matched pairs, source **RAIL / ROAD route ratio** has the following directly paired quartiles: US Midwest **1.33606 / 1.71258 / 2.06445** and France Burgundy **1.16695 / 1.40289 / 1.71129** (Q25, median, Q75). These include distance from each source cell centroid to its mapped transport-network connector. A connected graph can still have widely differing continuous detour distances; binary connectivity was saturated in all three earlier v3.0 blocks.

### Observed plant association: *not* geographically replicated

Response is the ratio of expected species intersection to expected union under rarefying positive GRIIS-listed GBIF plant **indexed records to 20 per cell**, NOT the statistical expectation of Jaccard similarity. We residualize ranked similarity and ranked log(RAIL path / ROAD path) against ranked geodesic distance, ranked ROAD path, GBIF source record-count ratio and source geometric-mean count. **Within-block pairwise records are not independent**. This is a descriptive partial rank correlation without p-value or confidence interval.

| Region | 2010–17 partial-rank, rarefied source | 2018–25 partial-rank, rarefied source | 2018–25 raw observed speciesKey-Jaccard partial-rank | Distinct network/plant pair rows in later period |
|---|---:|---:|---:|---:|
| US Midwest | −0.044934 | **−0.264201** | −0.148354 | 91 |
| France Burgundy | +0.148512 | **+0.053809** | −0.099421 | 120 |

**Sign convention:** negative implies relatively shorter railway routes are associated with more similar observed introduced-plant species keys, holding the listed source covariates approximately constant. US late period points in this hypothesis-compatible direction; France late rarefied measure is near zero/positive, and France raw-Jaccard measure changes sign to negative. Early vs late US and France vary too. Thus **the result is not robust, general or causal**. No train-carried seed/spore conclusion is permitted. At this stage there is no independent count of five or more broad ecoregion blocks with convergent directions.

Both predeclared v3.2 rail anchor strategies gave identical correlations in these US and French regions, since their selected anchors belonged to the same relevant graph and did not alter pair distances. This does not repair missing biological controls.

## Newly verified v3.6: EIGHT countries, ten 50km geographic blocks

**[Actual eight-country source run #42761140912](https://github.com/zuizui0223/tetsu/actions/runs/42761140912) SUCCESS**, four offline source-join unit tests, 10/10 geographically frozen macro blocks, zero recorded source query errors. It re-used independently completed frozen country-polygon and GBIF counts and fetched **fully paginated GRIIS-listed introduced plant speciesKey counts** in 2010–17 and 2018–25, after requiring both-period source opportunity.

- 160 predefined 2500km² equal-area cells across 10 nonoverlapping regional blocks in eight countries;
- 132/160 have **>=80% surface of designated Natural Earth country polygon** (not a full terrestrial habitat mask);
- 121/160 pass this source-country mask and >=100 GBIF indexed Plantae records + >=3 GBIF source datasets in **both periods**;
- **121/121** country/source-eligible cells have complete two-period speciesKey facet retrieval;
- **105/160** cells have BOTH-period >=5 GRIIS listed introduced speciesKey categories AND >=20 introduced indexed records, plus all earlier eligibility gates.

| Country | Frozen macro blocks | All 50km cells | Both-period introduced-flora source-qualified |
|---|---:|---:|---:|
| JP Japan (Tohoku/Kyushu) | 2 | 32 | **3** (1 + 2) |
| US (Midwest/Southeast) | 2 | 32 | **27** (15 + 12) |
| CA Ontario | 1 | 16 | **10** |
| GB Yorkshire | 1 | 16 | **12** |
| FR Burgundy | 1 | 16 | **16** |
| AU Southeast | 1 | 16 | **16** |
| DE Saxony | 1 | 16 | **14** |
| NZ North | 1 | 16 | **7** |
| **Total** | **10** | **160** | **105** |

**This is source data availability, not 105 independent ecological plots or 105 rail treatments.** Eight-country taxon incidence now exists, but continuous 50km physical rail-vs-road network routes were only completed for US Midwest and France Burgundy in v3.4. Japan's low qualifying count reflects at least partly dataset/GRIIS source availability; no low biological alien diversity or train dispersal mechanism is inferred.

Current verified scripts: [v3.6 country-masked full-facet macro species acquisition](../scripts/harvest_griis_50km_10blocks_v36.py), [quality-mask tests](../tests/test_griis_50km_10blocks_v36.py), [real-source workflow](../.github/workflows/griis-50km-10block-countrymask-v36.yml). The original GRIIS country list, 2021/2026 source maps and historical observation periods remain different-year inputs, and original DOI-bearing GBIF event rows / observer effort remain unavailable.

## Evidence gaps / current inference boundary

1. Physical transport source is Overture **2026**, whereas plant records are 2010–17 and 2018–25. Historical route opening/traffic and actual active freight operations are missing.
2. 50km source cells do not yet carry authoritative RESOLVE terrestrial **ecoregion** IDs, climate/land cover, population density or local port exposure. Natural Earth country fraction filter alone is not biome matching.
3. GRIIS lists describe **country-reported introduced** status; unlisted flora cannot be assumed native, and lists have nonuniform vintage/coverage.
4. GBIF facet summaries cannot measure repeated plot-based absence, duplicated observers, coordinate uncertainty or single survey events. For true species occupancy/change estimates, use DOI-backed original occurrences and survey-aware availability.
5. One arbitrarily chosen representative rail/road connector per 2500km2 cell is not a real plant observation or station; 75km halo may truncate network-shortest paths.
6. Two broad blocks **cannot** support globally generalized positive or negative railway transport effects. Do not add significance tests based on treating >200 dependent cell pairs as independent samples.

## Next test (separate current workstream, no future result claimed)

The pre-frozen **10 blocks / 8 countries / 160 cells** now have a country-masked full speciesKey GRIIS census, with **105 two-period introduced-source-qualified cells**. The next ecological identification gate is to calculate actually comparable RAIL and ROAD continuous routes for more of these same broad blocks, verify actual ecoregions, land and event-level source effort, and seek independent geographic replication before any biological effect claim. Regional selection and all sensitivity settings must be frozen before examining ecological effect signs, not filtered toward positive effects.

Reproducibility:
- [Physical rail source paths](../scripts/audit_overture_50km_rail_paths_v32.py) and [unit tests](../tests/test_overture_50km_rail_paths_v32.py).
- [Physical road source paths](../scripts/audit_overture_50km_road_paths_v33.py) and [unit tests](../tests/test_overture_50km_road_paths_v33.py).
- [Real two-block introduced flora × matched route comparison](../scripts/analyze_macro_50km_rail_road_aliens_v34.py) and [unit tests](../tests/test_macro_50km_rail_road_aliens_v34.py).
- [Natural Earth country fraction](../scripts/audit_naturalearth_macro_country_v35.py) and [tests](../tests/test_naturalearth_macro_country_v35.py).
- [Country-masked ten-block GRIIS source census](../scripts/harvest_griis_50km_10blocks_v36.py).
