# tetsu v4.6 — Six-country 50km equal-area macro study: physical rail, road, circuits, and country-listed introduced plants

**Actual source-backed results: 2026-11-10 (JST)**. This file records both hypothesis-compatible and contradictory findings from the same preselected source data. It must NOT be used to claim a causal train-borne seed dispersal effect.

## Study question and independent geographic units

Do pairs of broad 50km geographic cells with comparatively shorter mapped RAIL routes versus ROAD routes, or greater alternative-path RAIL CIRCUIT redundancy, show higher **observed GRIIS-listed introduced plant species-key similarity** even after source/route-distance controls?

The original tetsu v3.4 result from US Midwest vs France Burgundy changed direction. We therefore extended to **four additional independent country blocks** from the pre-existing frozen 8-country/10-block EPSG:6933 50km source frame, selected by independent botanical RECORD coverage, not by biological rail correlation: Canada Ontario, Britain Yorkshire, Germany Saxony and Australia Southeast. The resulting six study blocks are one 200×200km regional window per country (US, FR, CA, GB, DE, AU), not randomized global ecoregions or country-level floristic census replicates.

**All SOURCE data are real** (not simulations): Overture Maps Oct 2026 rail and road junctions, shortest paths in the 75km spatial halo; Circuitscape-network-equivalent physical RAIL effective resistance (edge resistance proportional to geodesic mapped track length); GBIF indexed positive speciesKey counts in preselected 50km EPSG:6933 cells from two eight-year periods; GRIIS country-specific introduced plant species-key filters and country boundary >=80% source-area mask. GBIF faceted data are *not* original occurrence event rows / observer or date-specific sampling efforts / confirmed ecological absence.

The outcome controls and analysis were kept as v3.4 and v3.8: two independent centroid-to-network anchoring strategies (nearest in cell vs nearest of largest physical component in cell), two periods (2010–17, 2018–25), two plant similarity responses (observed species-key Jaccard and ratio of expected 20-record-rarefied intersection to expected union). Same ranked geography, ROAD shortest path and indexed record-count covariates; the source CIRCUIT partial rank additionally conditions on physical RAIL shortest path. A country is ONE independent summary block, not its 20–120 overlapping cell pair comparisons.

## Source acquisition and complete country network results (all successful)

| Region, fixed 200km wide | 2026 mapped RAIL pair routes among 16 cells, nearest source connector | Median rail-to-geodesic circuity | 2026 mapped ROAD pair routes among 16 cells | Median road-to-geodesic circuity | Physical mapped RAIL circuit Q median |
|---|---:|---:|---:|---:|---:|
| US Midwest | 105 | 1.950427 | 120 | 1.186433 | 0.525157 |
| FR Burgundy | 120 | 1.685462 | 120 | 1.196970 | 0.588462 |
| CA Ontario | 37 | 2.002217 | 120 | 1.360600 | 0.278382 |
| GB Yorkshire | 55 | 1.428759 | 120 | source audit artifact; all 120 route pairs | 0.757443 |
| DE Saxony | 78 | 1.448124 | 120 | 1.219154 | 0.792217 |
| AU Southeast | 21 | 1.487415 | 120 | 1.530415 | 0.331064 |

CA/GB/DE/AU actual rail source physical shortest-path workflow [#42780398147](https://github.com/zuizui0223/tetsu/actions/runs/42780398147): four matrix jobs SUCCESS. Actual matching road physical source [#42781107642](https://github.com/zuizui0223/tetsu/actions/runs/42781107642): four matrix jobs SUCCESS; in GB and AU the original source included some road line records without both explicit terminal connector IDs. Those records were **not** snapped to invented junctions: only exactly identified incomplete-endpoint source records were excluded under an explicit <=2% source-quality cap, with source exclusion counts retained in the result JSON. A missing mapped junction is not biological absence. The original FR/US sources are v3.2 rail [#42759357722](https://github.com/zuizui0223/tetsu/actions/runs/42759357722) and v3.3 road [#42759664602](https://github.com/zuizui0223/tetsu/actions/runs/42759664602).

Four additional physical mapped RAIL electrical resistance jobs [#42781597862](https://github.com/zuizui0223/tetsu/actions/runs/42781597862): four matrix jobs SUCCESS, 0 actual source errors. Q = 1 - (effective equivalent circuit track resistance)/(shortest-path equivalent track resistance). CA Q 0.278382, GB 0.757443, DE 0.792217, AU 0.331064; all positive, structurally varying, and mathematically consistent with their respective source physical rail networks. Do not interpret Q as number of train movements or seed flux. Official Julia Circuitscape software was previously crossvalidated **only on three synthetic networks**; all large actual Overture graphs use an independently tested SciPy grounded-Laplacian solver.

**Do not compare unpaired raw physical road and rail median circuities as a railway efficiency effect.** Networks have distinct source link density and endpoint topology; botanical models use only matched cell pairs for which both source network paths exist.

## Real macro ecological outcome v4.4: RAIL-vs-ROAD shortest path vs country GRIIS flora

[Successful actual 6-country plant × Overture RAIL/ROAD route source #42782246966](https://github.com/zuizui0223/tetsu/actions/runs/42782246966), with source/quality and leave-country-out tests.

Response: pairwise observed introduced species-key similarity, with the primary cautious version based on **20-indexed-record rarefaction**. Predictor: log(rail source-cell-centroid route / road source-cell-centroid route), additionally controlling ranked great-circle geographic separation, ranked road route distance and GBIF source record opportunity. A **negative coefficient sign** aligns with the hypothesis that relatively more direct rail paths coincide with higher recorded introduced flora similarity.

| Independent country block | Both-period source-eligible GRIIS cells | Paired plant-rail-road source rows (nonindependent) | Recent 2018–25 rarefied20 partial rho | Recent 2018–25 RAW Jaccard partial rho |
|---|---:|---:|---:|---:|
| US Midwest | 15 | 91 | −0.264201 | −0.148354 |
| FR Burgundy | 16 | 120 | +0.053809 | −0.099421 |
| CA Ontario | 10 | 22 | +0.201410 | −0.024219 |
| GB Yorkshire | 12 | 36 | −0.011291 | −0.186594 |
| DE Saxony | 14 | 66 | −0.069825 | −0.252495 |
| AU Southeast | 16 | 21 | +0.411986 | +0.151796 |

**Result:** The RAW non-rarefied source species-key Jaccard source-correlation sign aligns with the railway hypothesis in **5/6 countries** (median about −0.123888). However, the record-rarefied response aligns in **3/6** and opposes it in **3/6**, median **+0.021259** (near zero). Source country-one-at-a-time leave-one-out medians change sign; the same six regions' alternative source anchor strategy also changes recent rarefied sign in several blocks.

**Scientific interpretation:** a weak hypothesis-compatible source effect exists when species-key inventories are used as-is, but it does NOT survive a basic record-count rarefaction. These GBIF indexed occurrences are not survey-standardized plant absences, and source observer effort/dataset biases remain uncontrolled. Do not infer railway-mediated dispersal, floristic homogenization causality, or significance from 91/120 dependent within-country pair observations. The most recently tested six independent country windows are still only selected geography, not a random global sample.

## Real macro ecological outcome v4.6: RAIL electrical circuit redundancy versus flora (ROAD route held constant)

[Successful actual 6-country GRIIS species-key × railway circuit + source ROAD shortest path #42782796042](https://github.com/zuizui0223/tetsu/actions/runs/42782796042), with both original anchor strategies, two periods and raw/rarefied outcomes. Q = 1 - R_eff / Dijkstra distance on the SAME RAIL source connectors, not a biological propagule pathway probability. The network predictor tests *extra physical path redundancy at a fixed shortest railway route*, also adjusting for straight-line geographic separation, ROAD source shortest path, and simple GBIF record-count factors.

A **positive** Q partial rho is the *rail-alternative-route-facilitated plant homogenization* hypothesis-compatible direction.

| Independent country block | 2018–25 partial rho of rail Q vs RAREFIED20 plant similarity | 2018–25 partial rho vs RAW observed species-key Jaccard |
|---|---:|---:|
| US Midwest | −0.181667 | −0.296914 |
| FR Burgundy | +0.066599 | +0.259670 |
| CA Ontario | +0.020171 | +0.434321 |
| GB Yorkshire | −0.239513 | +0.354890 |
| DE Saxony | −0.138180 | +0.327265 |
| AU Southeast | −0.628758 | +0.284371 |

**Recent RAREFIED20 observations:** **2/6** independent countries have the hypothesis-compatible positive sign; **4/6 negative**, equal-country median **−0.159923**. Leave-one-country-out medians remain negative. **Recent RAW unstandardized Jaccard:** **5/6** countries positive, median ~**+0.305818**.

This is a striking **response-definition sign reversal**. Despite the robust *structural* finding that physical rail network circuits have alternative routes, the **ecological relation between that alternative-route measure and the observed alien plant species composition is not robust to recording-count standardization**. Do NOT celebrate the 5/6 raw-positive signs without this sensitivity qualification. The survey-source detection process may be driving the contrast; the current approach cannot identify what fraction.

Important: The six-country circuit test controls ROAD physical shortest path, **not** ROAD electrical circuit Q, because independent road circuit effective resistance has only been computed for the US and FR source regions. The separate US/FR road-Q-controlled v4.2 test [#42777378885](https://github.com/zuizui0223/tetsu/actions/runs/42777378885) already failed to yield a geographically replicated positive recent railway plant correlation. Do not claim six-country direct matched ROAD vs RAIL electrical resistance effects.

## Main limitations and decision

1. **Recording bias:** GBIF speciesKey facets are positive indexed record aggregates; no occurrence-level DOI download, duplicate-event/photographer adjustment, coordinate precision, standardized survey plots, biological absences or real abundance. RAREFIED20 is a count correction, not a full effort/detection adjustment. RAW and rarefied outcomes' opposite patterns are not independent evidence.
2. **Actual exposure:** 2026 Overture physical rail/road networks do not establish that specific rail lines were active or carried freight during 2010–17 and 2018–25; transit direction and volume are unobserved. A source rail line can be a habitat correlate rather than a seed-transport vector.
3. **Network source graph:** Overture source-boundary halo 75 km, source railroad vs road classes, centimeter/epsilon same-position connectors, truncated outside routes and source connector snapping affect the assumed electrical conductances. A circuit score is not biological dispersal probability.
4. **Macro sample:** 6 countries × 1 selected broad geography each. Country-block summaries are independent at a more reasonable level than the 50km pair comparisons, but 6 is NOT enough for continental/biome-universal causal conclusions; and 21–36 paired observations in AU/CA/GB are dependent and relatively limited.
5. **Ecoregions/land:** Natural Earth country-area fraction >=80% is not an actual RESOLVE/WWF terrestrial ecoregion or vegetation fraction mask. Port access, climate, urban intensity, grass cutting and observed plant traits are not yet included.
6. **Inference:** no valid pairwise hypothesis p-values, no properly evaluated conditional freight mechanism, and no evidence that high *physical* electrical conductance causes a shift in flora.
7. **Reproducibility:** Both anchoring choices and both response definitions must remain in the report even when they disagree. Do not cherry-pick source accessibility, metric or year to claim a universal positive outcome.

### Scientific decision

Use continuous railway shortest-path and Circuitscape-compatible electrical redundancy as **two distinct exploratory infrastructure descriptors** in a regional multi-source analysis. The main biological hypothesis **has not survived source-standardization and independent country contrasts**. Prioritize original occurrence-row and observer/dataset effort controls, historical active freight route documentation, true land/biome masks and independent country×ecoregion replication over extra pairwise correlations.

Source files and exact runs:
- [Four-country rail path workflow](../.github/workflows/rail-50km-expansion-v43.yml) and [actual run #42780398147](https://github.com/zuizui0223/tetsu/actions/runs/42780398147)
- [Four-country road path workflow](../.github/workflows/road-50km-expansion-v43.yml) and [actual run #42781107642](https://github.com/zuizui0223/tetsu/actions/runs/42781107642)
- [Six-country plant-shortest-route analysis](../scripts/analyze_6country_50km_network_plants_v44.py), [tests](../tests/test_6country_50km_network_plants_v44.py), [actual run #42782246966](https://github.com/zuizui0223/tetsu/actions/runs/42782246966)
- [Four-country Circuitscape-equivalent rail source workflow](../.github/workflows/rail-circuit-four-new-macro-v45.yml) and [actual run #42781597862](https://github.com/zuizui0223/tetsu/actions/runs/42781597862)
- [Six-country circuit/flora partial rank](../scripts/analyze_6country_circuit_macro_plants_v46.py), [tests](../tests/test_6country_circuit_macro_plants_v46.py), [actual run #42782796042](https://github.com/zuizui0223/tetsu/actions/runs/42782796042)
