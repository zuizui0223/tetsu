# tetsu v7.3 — real original GBIF events versus legacy facets: French GRIIS introduced flora × GIFT growth form × physical railway/road paths

**Measured Nov 11, 2026.** [GitHub Actions genuine source run #43031487832](https://github.com/zuizui0223/tetsu/actions/runs/43031487832) **SUCCESS**, using *previously independent successful source artifacts*:

- Original August 2026 public GBIF 81-object Parquet interpretive event audit [v6.9 #42978217481](https://github.com/zuizui0223/tetsu/actions/runs/42978217481) — source range **2018–2025**; 697,577 ALL-plant source bounding-envelope rows; 24,279 deduplicated source country-GRIIS strict-COL accepted introduced plant event–species counts across 16 frozen 50km cells; 18,991 events with known coordinate uncertainty ≤5km.
- Official May 2026 GBIF→COL26.5 strict accepted country key map [v6.3 #42844630098](https://github.com/zuizui0223/tetsu/actions/runs/42844630098). Species concepts absent or ambiguous in the map remain UNKNOWN.
- Independent public GIFT 1.2.1 growth-form accepted-name assignment [v7.1 #42980383943](https://github.com/zuizui0223/tetsu/actions/runs/42980383943), using legacy GBIF accepted species keys. Only unambiguous direct COL concepts are included; synonym-only original GRIIS alternative is excluded from the primary outcome.
- Original pre-frozen France Burgundy 4×4 EPSG:6933 50km cells from [v3.6 #42761140912](https://github.com/zuizui0223/tetsu/actions/runs/42761140912), plus physical 2026 rail [v3.2 #42759357722](https://github.com/zuizui0223/tetsu/actions/runs/42759357722) and physical road [v3.3 #42759664602](https://github.com/zuizui0223/tetsu/actions/runs/42759664602).

## What genuinely changed

Country-GRIIS introduced species concepts from original GBIF **event rows**, not only old GBIF indexed presence facets, now have GIFT growth form. The source match is:

| Strict accepted GRIIS events, 16 source cells | Herbaceous | Woody | Unmatched/ambiguous |
|---|---:|---:|---:|
| All original events | **13,052** | **9,462** | 1,765 |
| Original events, known uncertainty ≤5km | **10,215** | **7,339** | 1,437 |
| Distinct matched COL concepts among ALL original events | 186 | 65 | excluded from matched outcome |

There were **22,514/24,279 (92.73%) original event–species counts** assignable to strict unambiguous growth form, vs **17,554/18,991 (92.43%)** for the precise-location subset. These are **events, not 22,514 independently surveyed plants or species**.

Within the fixed 16-cell FR Burgundy geometry, all **120 cell-pairs** had physical RAIL and ROAD paths and both herbaceous/woody source outcomes in ALL original, precise-original and legacy facet cases. These are *dependent comparisons using 16 repeated cells*, NOT 120 independent sample blocks. Both physical network anchor options yielded numerically identical path rows for France, so they are NOT two independent robustness replications.

## Results: precisely the same physical pairs and covariates, three botanical source definitions

The predictor is log(RAIL shortest mapped physical route / ROAD shortest mapped route). Lower values imply relatively more direct rail access. We residualize ranked plant similarity and ranked rail/road route against source geodesic separation, physical ROAD path length, and legacy indexed ALL-plant observation-count ratio/geometric mean. Negative partial rho is hypothesis-compatible for rail-associated species similarity.

**2018–2025, nearest-in-cell physical paths, partial rank rho, no pairwise p-values:**

| Response source | Herbaceous K20 expected rarefied source overlap | Woody K20 | Herbaceous raw species Jaccard | Woody raw species Jaccard |
|---|---:|---:|---:|---:|
| 81-file original events, all source coordinate precision | **−0.011523** | **−0.149633** | −0.158238 | **+0.055146** |
| Same original events, known ≤5km precision only | **−0.057278** | **−0.208302** | −0.071956 | **+0.049118** |
| Original legacy 2026 GBIF source indexed species facets mapped to same strict COL-GIFT concepts | **−0.022836** | **−0.141642** | −0.107863 | **+0.047549** |

K20 is the ratio of expected cross-cell present-species intersection to expected union after hypergeometric 20-record rarefaction per growth group, NOT E[Jaccard]. Raw Jaccard only measures the unstandardized list of positive recorded taxon categories and is highly sensitive to effort/source records.

### Interpretation

1. In **this ONE French region**, a negative woody rail-relative-distance association persists when replacing old indexed facets with deduplicated original-event counts; the known ≤5km subset remains negative and is more negative descriptively. This is a **real source-backed single-region exploratory association**, not independent globe-wide validation or observed seed transport.
2. The **woody raw species Jaccard is POSITIVE** in all three data definitions, while the rarefied woody response is negative. Response choice strongly affects biological interpretation even with original event rows. Neither is a validated railway seed-propagule flux measure.
3. The earlier *multi-country* GIFT woody source pattern [v7.2 #43009955904](https://github.com/zuizui0223/tetsu/actions/runs/43009955904) **did not survive removing the dominant GBIF provider** (woody 4/4 negative BEFORE to 2/4 AFTER; median −0.091771→−0.024878 on identical source cells). These v7.3 original events may still be dominated by one provider; no source-independent confirmation has occurred.
4. The event-source counts are drawn from the **conservative 81-file A5 polygon+centroid union, not verified complete original GBIF**. The source detected a country A5 polygon-versus-centroid coverage disagreement. A genuine public source preview is not necessarily a full regional census. Only 2018–2025 was analyzed; **2010–2017 original events NOT scanned**.
5. Physical rail paths are based on Overture 2026 maps and contain neither actual historical active freight routes nor observed propagules. Cross-cell species records are not railway-verge plots. GIFT woody/herbaceous growth form is not directly measured cutting/wind tolerance; no plant adaptation, phenotypic plasticity or genetic selection inferred.

## Immediate next macroecological test

First partition original *event-grain* plant records by the actual largest GBIF dataset source before aggregation, retaining strict source-specific row/event dedup and known coordinate uncertainty, then compare FULL and WITHOUT the dominant provider on identical source cell pairs. This must be attempted in multiple geographically independent blocks; France alone cannot establish macroecological universality. Next, obtain independently dated **active freight connectivity** plus species-specific dispersal or cargo-pathway traits, and separate independently measured railway **mowing/wind/heat exposure** and road/verge controls. Only then can network-mediated delivery be statistically distinguished from environmental filtering and sampling artifacts.

Code: [real original event-vs-facet GIFT test](../scripts/analyze_fr_original_GBIF_GIFT_rail_macro_v73.py) · [strict synthetic QA](../tests/test_fr_original_GBIF_GIFT_rail_macro_v73.py) · [live Actions workflow](../.github/workflows/FR-original-events-GIFT-rail-macro-v73.yml). The raw observation IDs, precise locations and human names are never committed.


## v7.3b — 16-cell leave-one-cell-out sensitivity (source verified)

**[Actual extended source workflow SUCCESS #43033395401](https://github.com/zuizui0223/tetsu/actions/runs/43033395401).** Every one of the 16 France 50km cells was excluded in turn, then the same physical RAIL/ROAD route-ratio partial rank was recalculated on the remaining 15 cells (105 dependent cell pairs); all 16 refits were estimable. This is a robustness diagnostic, NOT 16 independent country/ecoregion replications or 16 significance tests.

| Plant source and group | Full K20 partial rho | K20 negative leaveouts | Raw Jaccard partial rho | Raw positive leaveouts |
|---|---:|---:|---:|---:|
| Original ALL events, **woody** | −0.149633 | 15/16 | +0.055146 | 13/16 |
| Original **known ≤5km**, **woody** | **−0.208302** | **16/16** | +0.049118 | **13/16** |
| Original ALL events, herbaceous | −0.011523 | 9/16 | −0.158238 | 0/16 |
| Original known ≤5km, herbaceous | −0.057278 | 15/16 | −0.071956 | 0/16 |
| Legacy indexed species facets, woody | −0.141642 | 15/16 | +0.047549 | 13/16 |

For the primary **known ≤5km woody K20**, all 16 excluded-cell refits had negative signed coefficients, range **−0.394052 to −0.036500**. But for the SAME exact source group using unstandardized species-list Jaccard, 13/16 refits were positive (range **−0.036542 to +0.171270**). **Response-definition dependence is NOT explained away by deleting any single 50km cell.**

This does not identify the true process: it could reflect source dataset/detection composition, a differing sensitivity of rarefied overlap to taxa with many records, taxonomic resolution or habitat effects. It remains only one France block with current (2026) map paths and observation records in 2018–25. Before promoting a railway dispersal effect, perform source datasetKey exclusion on ORIGINAL de-duplicated event rows and replicate among independent country–ecoregion blocks with independently verified freight service. Neither the apparent positive nor negative direction is claimed as universal.

