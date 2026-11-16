# v2.7 — Ten-country macro pilot: introduced plant composition, mapped railway connectivity, and actual road-density matching

**Verified real-source workflow:** https://github.com/zuizui0223/tetsu/actions/runs/42732247462 (2026-11-10 JST), SUCCESS; seven synthetic source/matching tests passed. All four independent-source artifacts were downloaded and the matching procedure ran successfully. This document records observed data diagnostics; **no train-mediated plant dispersal effect, p-value, or genetic adaptation** has been established.

## Scientific question and properly limited answer

Do mapped railway-connected small geographic cells share more **GRIIS-listed introduced plant species** than source-matched unconnected cells, and does this difference survive geographic distance, plant source-record count and **roadline length density** matching?

**Answer from the current selected urban windows: no convincing generalized positive effect.** The late 2018–25 all-cell median is only +0.000877 after road density matching, with 3 of 5 eligible city blocks positive. In the more focused both-have-rail comparison only two blocks are informative, the median is negative (−0.003231), and directions split one positive and one negative. In the earlier period the road-density-matched all-cell median is essentially zero.

These are descriptive source-record compositional measures, not ecological causal effects; the entire result must not be framed as proof that trains transport seeds or that railways have no effect.

## Real data input and provenance

- Country GRIIS introduced-plant backbone taxa: successful 9-country official GBIF index source audit https://github.com/zuizui0223/tetsu/actions/runs/42730506112. Accepted source country lists: JP 1,434; KR 32 (excluded by <100 checklist quality gate); US contiguous 3,837; CA 1,249; GB 1,575; DE 473; FR 1,618; AU 4,105; NZ 443. TW unassessed and excluded. Non-membership is **unknown**, never coded native.
- Full GBIF faceted speciesKey × 90 fixed small city-cell × two eight-year periods: https://github.com/zuizui0223/tetsu/actions/runs/42728479639. These are indexed RECORD counts per species key, not actual downloaded DOI-citable event rows and not reliable biological absence. 180/180 source cell-period bins retrieved; 65/90 cells passed an initial both-period recording-opportunity screen.
- Mapped Overture 2026-10-21 physical railway connector graph, clipped within fixed ten city boxes: https://github.com/zuizui0223/tetsu/actions/runs/42728551036. These are mapped physical rail connections and not current train service or traffic by route.
- **NEW** 2026 Overture ROAD and RAIL geodesic line lengths clipped to **each** of the same 90 cells: https://github.com/zuizui0223/tetsu/actions/runs/42731846908. All 10/10 source cities succeeded with zero errors. These are physical line km, *not* traffic or road-network journey times. Distinct area denominators are computed geodesically; no raw proprietary or third-party geometry is republished.

Country×city geographic windows: Sendai JP, Leipzig DE, Daejeon KR, Hsinchu TW, Madison US, Ottawa CA, York GB, Dijon FR, Canberra AU, Hamilton NZ. This is **one selected city per country**, not representative national or ecoregion-wide macroecological inference.

## Measured road and railway kilometer sources

| Fixed urban source window | Mapped rail km (clipped) | Mapped selected road km (clipped) | Cells with mapped road / 9 |
|---|---:|---:|---:|
| Sendai | 277.092 | 5,620.453 | 9 |
| Leipzig | 577.013 | 2,058.554 | 9 |
| Daejeon | 318.173 | 3,266.294 | 9 |
| Hsinchu | 302.102 | 3,340.501 | 8 |
| Madison | 135.479 | 3,514.041 | 9 |
| Ottawa | 189.810 | 4,547.449 | 9 |
| York | 200.154 | 1,147.662 | 9 |
| Dijon | 458.130 | 2,001.399 | 9 |
| Canberra | 20.679 | 4,153.254 | 9 |
| Hamilton | 106.695 | 1,486.050 | 9 |

Directly comparing total km across these unequal-sized city boxes is inappropriate. Per-cell km/km² is used as the covariate. Rail-line km in some cells will not match binary rail-connector presence when a mapped segment crosses a cell without containing a connector in it; no rail-presence imputation is allowed.

## Actual introduced-plant matched results

Similarity is the ratio of expected rarefied (20 indexed record) species-key intersection and expected union, **not E[Jaccard]**, based on GBIF records marked PRESENT and included in country-specific GRIIS introduced lists. The table is an **equal-city** summary of the median within-city matched pair difference connected minus disconnected. Matching uses centroid geographic distance, GBIF indexed record-count scale/balance, and optionally road line-length-density scale/balance. Cell pairs share source cells and are not independent observations.

| Source period & comparison set | Informative city blocks | Positive blocks / negative blocks | No road match median | Road-density match median |
|---|---:|---:|---:|---:|
| 2010–17 all source cells | 4 | 2/2 after road control | −0.013027 | **−0.000021** |
| 2018–25 all source cells | 5 | 3/2 after road control | +0.000328 | **+0.000877** |
| 2010–17 both cells have railway connector | 2 | 2/0 after road control | +0.008757 | **+0.018250** |
| 2018–25 both cells have railway connector | 2 | 1/1 after road control | −0.004249 | **−0.003231** |

In 2018–25 all-cell road-density adjusted comparisons (region-level median differences): Madison −0.001069; Ottawa +0.000877; York +0.001487; Dijon +0.002310; Hamilton −0.042680. Other cities fail treatment contrast support or matching gates, including Leipzig (all qualifying city cells source-connected) and Canberra (only one positive connected-cell pair). Japan's near-city alien cells are too few to identify connected vs disconnected pairs.

The late-period positive small median is **negligible**, mixed across cities, and disappears/reverses under the stricter comparison of mapped rail-present cells with one another. No statistical significance or railway propagule mechanism is identified. Negative test is also underpowered; cannot claim universal no railway effect.

## Additional quality findings

- Actual Overture road km exist in 9/9 cells for nine selected cities; Hsinchu has 8/9. Missing road km source in a cell is not treated as positive evidence of a road-free control.
- **Connector presence vs clipped physical rail length discrepancy** arose in some city cells (Dijon and Hamilton reported a single discrepancy in a qualified subset). This is not proof of a data error: a line may traverse a cell without a connector within it. Future exposure definition should consider **length-based rail passage** and **network-reachable access**, independently.
- Country-specific GRIIS checklist coverage varies, and route/station geocoding or real freight utilization was not verified.
- Source plant records have not been deduplicated at occurrence/event/observer level or assigned reliable coordinate uncertainty; record-count-based rarefaction cannot remove systematic GBIF source bias.
- Railway 2026 source maps postdate both 2010–17 and 2018–25 plant record windows; no treatment-onset temporal inference is possible.
- Road **density** is not road **connectivity**, traffic, travel time or port/freight logistic throughput.

## Priority next macro design

**Stop adding small near-identical urban source windows merely to chase a positive sign.** Instead freeze independent country × ecoregion 25/50 km equal-area blocks with comparable rail-connected and disconnected rail-present locations, and require explicit positive/negative source exposure support. Use full DOI-backed GBIF occurrence data and a country-introduced-plant status crosswalk to correct observation effort, native/alien status uncertainty and spatial precision. Introduce *actual ROAD network shortest-path accessibility* as a matched comparator plus climate, land cover, urban intensity and ports. Separate documented active freight routes from passenger, disused and freight-unknown rails.

Success criterion is a **replicated effect in independent regional blocks, robust to road and source controls**; if it is absent, report the null/uncertainty instead of manufacturing a railway-dispersal conclusion.

## Reproducibility

- Railway/road km source script: ../scripts/audit_overture_transport_cells_v26.py
- Six targeted road/rail source unit tests: ../tests/test_overture_transport_cells_v26.py
- Source workflow and real execution: ../.github/workflows/overture-cell-road-rail-density-v26.yml , https://github.com/zuizui0223/tetsu/actions/runs/42731846908
- Road-matched alien species community script: ../scripts/analyze_griis_rail_with_road_density_v27.py
- Seven source/matching tests: ../tests/test_griis_rail_with_road_density_v27.py
- Four-source joined analysis workflow: ../.github/workflows/road-adjusted-alien-rail-macro-v27.yml , https://github.com/zuizui0223/tetsu/actions/runs/42732247462
- Previous GRIIS source, no road control, v2.5: rail_griis_macro_results_v25.md
