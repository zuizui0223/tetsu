# tetsu v2.5 — GRIIS-introduced plant composition vs mapped railway connectivity

**Real source result:** GitHub Actions run #42730951549, 2026-11-10 JST, SUCCESS. Six targeted synthetic validation tests passed, with three independently completed upstream source inputs reused as immutable Actions artifacts: GBIF geographic speciesKey facets 90 cells × 2 periods (#42728479639), Overture railway connector topology for 10 city windows (#42728551036), and country-specific GRIIS introduced-plant backbone keys (#42730506112).

**Research question:** Are geographically similar cells connected by mapped railway segments compositionally more similar in INTRODUCED plants than disconnected cells after limited geographic distance and recording-volume matching?

**Important scientific status:** No railway-mediated biological dispersal effect has been identified. Even the direct exploratory GRIIS-only association is unstable and tiny, and the stricter both-have-rails design has insufficient independent blocks.

## 1. Countries and record matching

Country-specific GRIIS (introduced and invasive) source list resolved to GBIF backbone species keys: JP=1,434; KR=32; US=3,837 (contiguous USA checklist, applied to Madison); CA=1,249; GB=1,575; DE=473; FR=1,618; AU=4,105; NZ=443. Taiwan unassessed and excluded, never coded as zero introduced taxa or as native. KR's 32-key source fails the >=100-key quality gate and is omitted from biological comparisons.

GBIF cells are nine geographic rectangles in each of ten city windows, with source speciesKey count facets fully paginated in two periods (2010–2017, 2018–2025). Cell inclusion requires >=100 all-plant indexed records and >=3 source dataset categories, plus >=5 actually observed GRIIS-listed introduced species and >=20 introduced-species record counts in the cell-period. Observed unlisted species are UNKNOWN, not assumed native. *Expected rarefied intersection / expected rarefied union* is computed at a 20-indexed-record draw; this is NOT mathematically the expectation of Jaccard and is NOT a standardized biological census. Cell pairs are dependent.

Railway connector source is Overture 2026-10-21.0, mapped physical connectivity within a city window. No proof of active services or current freight transport. A pair's disconnected status can be an artifact of clipping the graph at the window boundary.

Matched geography/source diagnostic: nearest available connected vs unconnected pair, without replacement, matched on approximate centroid geographic distance (max 5 km or 33% relative distance difference), log geometric mean all-plant indexed records (within 1) and log absolute record-effort ratio (within 1). Similarity outcome NOT used for matching. City-level medians are summarized as independent units; there are no pseudo-independent cell-pair p-values.

## 2. Real GRIIS-introduced-only city×period source results

| City/country | Intro listed observed species-key union, earlier / later | Early eligible alien cells | Late eligible alien cells | Early matched contrast | Late matched contrast |
|---|---:|---:|---:|---:|---:|
| Sendai JP | 44 / 63 | 1 | 3 | NA | NA |
| Leipzig DE | 80 / 187 | 4 | 9 | NA | NA (all source rail cells connected) |
| Madison US | 165 / 426 | 6 | 9 | −0.002407 | −0.001069 |
| Ottawa CA | 232 / 467 | 9 | 9 | +0.006310 | +0.000475 |
| York GB | 180 / 312 | 8 | 9 | −0.023647 | +0.000328 |
| Dijon FR | 84 / 183 | 7 | 9 | −0.069161 | +0.013344 |
| Canberra AU | 643 / 852 | 9 | 9 | NA (one connected source pair) | NA |
| Hamilton NZ | 27 / 123 | 1 | 7 | NA | −0.042680 |
| Daejeon KR | unassessed: country GRIIS only 32 keys | NA | NA | NA | NA |
| Hsinchu TW | unassessed: no matching country GRIIS list | NA | NA | NA | NA |

Note that "introduced species-key union" means taxon keys recorded in cells with adequate total-plant source coverage. It is not a fully observed site flora, and differences in early vs late union are NOT estimates of colonization or railway spread.

### Equal-weight city direction

| Period and treatment contrast | Informative city blocks | Positive / negative | City-level median matched similarity contrast | Adequacy |
|---|---:|---:|---:|---|
| 2010–2017 all source cells | 4 | 1 / 3 | **−0.013027** | below five-city design minimum |
| 2018–2025 all source cells | 5 | 3 / 2 | **+0.000328** | reaches screen but tiny and direction mixed |
| 2010–2017 only pairs with both rail cells | 2 | 2 / 0 | **+0.008757** | underpowered geographic blocks |
| 2018–2025 only pairs with both rail cells | 2 | 0 / 2 | **−0.004249** | underpowered and reverses sign |

The small positive late-period all-cell result does **not** support robust railway-induced non-native plant homogenization. It becomes negative when both sides contain mapped rail, leaving only two meaningful matched city blocks. No P-values, confidence intervals or train-seed transport assertions are warranted.

## 3. Why the present positive signal is weak

- All-cell comparisons conflate rails vs no rails, urban/green habitats, observer accessibility and conservation/survey pathways.
- The negative-control *ROAD physical connector graph and road intensity by cell* have NOT yet been joined, so rail association could be road association.
- Source countries are each one hand-chosen urban window, not random national or ecoregional replication, and remote window-boundary graph truncation creates false disconnected labels.
- Spatial cell polygons are WGS84 rectangular windows, not equal-area 25/50 km cells, and no matching environmental, port, management, land-cover or date-specific infrastructure overlays are included.
- Country GRIIS checklist recording standards and source vintages vary; absence from GRIIS is unknown status, not biological native status.
- Current Overture 2026 network exposure overlaps both GBIF historic periods, so early-vs-late contrasts are NOT a true railway opening-before-after design.
- The 20-record rarefaction corrects one component of record-count inequality, not dataset/observer geographic sampling differences.
- Main conditional mechanistic hypothesis (train-borne propagules) needs independent station/cargo/propagule observations and cannot be identified from spatial association alone.

## 4. Next substantive macro-analysis (not more source-count audits)

1. Add actual 2026 Overture ROAD connectedness/intensity and 25/50 km source road accessibility, with comparable cell-region buffers; compare rail association conditional on roads.
2. Broaden from ten urban rectangles to independent country × ecoregion blocks where both connected and disconnected rail pairs exist; use active rail / freight checks and spatial null models.
3. Use real DOI-citable occurrence records (not only aggregated source facets) to audit event-level sampling/coordinate uncertainty and deduplicate.
4. Where country alien taxonomy is insufficient (KR/TW), mark missing status or source a trusted alternative; never infer native species from GRIIS non-membership.

The study remains an ambitious globally replicated biogeography design, but the directly tested introduction-specific railway connectivity association is **not yet a positive robust ecological result**.

### Source files

- Actual v2.5 run: https://github.com/zuizui0223/tetsu/actions/runs/42730951549
- Source analysis: ../scripts/analyze_griis_alien_rail_source_v25.py
- Unit tests: ../tests/test_griis_alien_rail_source_v25.py
- Workflow: ../.github/workflows/griis-rail-alien-9country-v25.yml
- GRIIS 9 country source: https://github.com/zuizui0223/tetsu/actions/runs/42730506112
- GBIF source matrix: https://github.com/zuizui0223/tetsu/actions/runs/42728479639
- Overture source physical network: https://github.com/zuizui0223/tetsu/actions/runs/42728551036
