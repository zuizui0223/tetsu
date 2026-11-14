# v3.0 — From ten small urban windows to independent 50km equal-area macro source regions

**Report dated 2026-11-10 (JST).** Real source success and negative/uncertain biological inference are distinguished. **No railway-caused plant dispersal, floral trait evolution, mowing, flowering or genetic adaptation** has been established.

## 1. Why change scale

Prior v2.7–v2.10 source science used one small 3×3 urban sample per country (9 small WGS84 cells), often with the same physical rail and road network in almost every cell. This causes a serious **positivity / discrimination failure**: even a perfect source cannot identify a railway-specific contrast if all cells are connected or if very few disconnected rail-present cell pairs exist.

Actual v2.8 Overture October-2026 ROAD graph [run #42737690925](https://github.com/zuizui0223/tetsu/actions/runs/42737690925): all ten geographical source regions successfully mapped from explicit ROAD connector IDs, with complete nine-cell centroid route distance matrix per city. **Nine cities have 36/36 physically road-connected cell pairs**; Hsinchu (TW) has 28/36 because only 8/9 cell centroids snap to the source road within the 5km gate. Road binary connectedness is largely saturated, unlike physical road route length / circuity. Road routes remain optimistic graph approximations excluding road access restrictions, turns, freight vehicles and map edges, not actual travel time or traffic.

The 5-source source-only GRIIS-introduced plant × mapped RAIL connector × observed GBIF speciesKey × Overture ROAD km density × ROAD physical path diagnosis is a **posthoc** sensitivity to earlier analysis, not preregistered independent validation. [Successful run #42738469033](https://github.com/zuizui0223/tetsu/actions/runs/42738469033), all 7 unit tests passed. Equal-city 2018–25 matched introduced-plant observed similarity difference (rail-connected minus rail-unconnected):
| Predeclared network road path/circuity matching sensitivity | Informative city windows | Positive city windows | Equal-city median source similarity difference |
|---|---:|---:|---:|
| Tight | 5 | 3/5 | +0.001421 |
| Main | 5 | 4/5 | +0.002749 |
| Loose | 5 | 3/5 | +0.000877 |
| Main restricted to both rail-present cells | **2 only** | 1/2 | +0.001545 |
| Tight restricted to both rail-present cells | 2 only | 0/2 | −0.001425 |
| Loose restricted to both rail-present cells | 2 only | 1/2 | −0.003231 |

The signal remains small, sensitive to matching conditions and biologically non-identifiable. Neither the positive all-cells value nor the 2-city negative/positive values warrant a railway dispersal claim, p-value or generalized no-effect conclusion. Spatially shared cells are non-independent pairwise records.

## 2. New source frame: fixed equal-area cells, actual real coverage

[Actual 50km GBIF source opportunity run #42738021553](https://github.com/zuizui0223/tetsu/actions/runs/42738021553): **SUCCESS**, seven offline tests, all **320/320** source cell×period queries, zero source errors.

Before seeing biological or railway outcomes, **24** 200km ×200km geographic blocks were frozen across **8 countries** (three geographically separated blocks per country; JP, US, CA, GB, FR, AU, DE, NZ). Each block contains 4×4 non-overlapping EPSG:6933 **50km nominal equal-area cells (2,500km² each)**, giving 384 cells / 960,000km² proposed complete *map footprint*, not 384 terrestrial valid ecological plots.

The first **ten geographically frozen blocks** (eight countries, with Japan and USA represented twice) actually returned:
- **160 unique equal-area cells**, nominal summed map footprint 400,000km²;
- **137/160 cells** pass both 2010–17 and 2018–25 source-opportunity gates of at least 100 indexed species-rank Plantae records and at least 3 source datasets per cell-period;
- 320/320 actual GBIF source cell-period count requests complete, 0 error;
- source counts are positive GBIF indexed records only, NOT unique species, non-native flora, plant abundance, independent field surveys or the species-by-cell matrix.

| Geographic fixed block / country | 50km cells | Both-period record + 3 datasets eligible | Source records 2010–17 | Source records 2018–25 |
|---|---:|---:|---:|---:|
| Tohoku Japan | 16 | 12 | 24,351 | 4,835 |
| Kyushu Japan | 16 | 9 | 8,103 | 5,038 |
| Midwest USA | 16 | 15 | 24,569 | 301,588 |
| Southeast USA | 16 | 13 | 19,080 | 244,426 |
| Ontario Canada | 16 | 10 | 21,183 | 269,077 |
| Yorkshire UK | 16 | 14 | 2,082,146 | 1,911,619 |
| Burgundy France | 16 | 16 | 2,033,433 | 698,318 |
| Southeast Australia | 16 | 16 | 201,146 | 406,685 |
| Saxony Germany | 16 | 16 | 118,709 | 908,055 |
| Waikato New Zealand | 16 | 16 | 26,664 | 160,086 |
| **TOTAL** | **160** | **137** | **4,559,384** | **4,909,727** |

The above country windows were selected by distinct geography, NOT verified terrestrial ecoregion IDs; no land/ocean fraction was yet masked. Count changes can reflect data submission/observer effort, not plants shifting. **Do not infer historical railway expansion from the 2026 rail map** or from GBIF first records.

## 3. Railway macro-level exposure, separate from ecological outcome

New [50km regional Overture rail physical connector source audit](../scripts/audit_overture_50km_macro_rail_v30.py), [unit tests](../tests/test_overture_50km_macro_rail_v30.py) and [three independent regional workflow](../.github/workflows/overture-rail-50km-macro-v30.yml) map the same EPSG:6933 equal-area source cells to actual 2026 Overture rail-line passage and connector components within a 75km halo around a frozen 200km study block.

This is deliberately **outcome-unseen network opportunity testing**: for each 16-cell block it reports numbers of rail-passage cells, connected rail–rail cell pairs and locally disconnected rail–rail cell pairs. A source block must have both connected and disconnected *rail-present* pairs (proposed >=3 of each) before a comparative rail effect is identifiable in that block. Geographic halo paths may still miss national connections. A source segment with rail passage but no connector in a cell remains rail-present (not coded rail-absent).

Do not present v3.0 macro railway ecological coefficients until that real source workflow succeeds and the *same 50km* cells have a country-introduced plant species-key matrix, road-network counterpart, and qualified source coverage.

## 3b. REAL source result: 50km physical rail connectivity is saturated

**Successful source-backed three-macro-block run:** [GitHub Actions #42739564813](https://github.com/zuizui0223/tetsu/actions/runs/42739564813), 2026-11-10, six unit tests, zero source errors. These macro rail blocks were geographically fixed **before seeing their train-plant outcomes**; 2026 source network topology only was queried.

| Actual predeclared 200km macro block | Overture rail segments in block+75km halo | 50km cells actually crossed by mapped rail | Both-rail cell pairs physically connected | Both-rail pairs disconnected in halo | Binary rail exposure positivity? |
|---|---:|---:|---:|---:|---|
| JP Tohoku | 3,252 | 16/16 | 120 | **0** | **No** |
| US Midwest | 13,697 | 15/16 | 105 | **0** | **No** |
| FR Burgundy | 18,253 | 16/16 | 120 | **0** | **No** |

**Major macroscopic structural conclusion:** Increasing city-scale source windows to **50km equal-area cells** does NOT rescue binary railway connectivity. Among rail-present cells, there are **zero** within-window mapped rail disconnected contrasts in all three independent broad geographic blocks. This is a positivity failure before plant composition is considered, NOT evidence that connected rail causes plant homogenization, and NOT proof every real train/freight network in each country is operationally connected.

**Consequent analysis reformulation:** Abandon binary connected/unconnected as the primary broad-region exposure in these blocks. Model *continuous railway-network accessibility/shortest-path excess or circuity*, ideally as a contrast against real ROAD-network accessibility, with active freight topology and country/ecoregion/climate/observer adjustment. Require source-based variability/positivity and external graph-halo sensitivity before a biological coefficient is estimated. A single connected component is still capable of varied shortest-path distances. Do not search only for less connected regions posthoc to manufacture a positive contrast.

## 3c. REAL macro botanical speciesKey coverage — GRIIS introduced-only (v3.1)

[Actual equal-area source GRIIS plant taxon run #42740320597](https://github.com/zuizui0223/tetsu/actions/runs/42740320597), 2026-11-10, SUCCESS; **five offline QA tests**, all three fixed source regions, no live source errors. Analysis source [full-facet species-key harvest](../scripts/harvest_griis_macro_50km_facets_v31.py), [offline tests](../tests/test_griis_macro_50km_facets_v31.py), [actual workflow](../.github/workflows/griis-50km-macro-species-v31.yml).

The original v2.9 pre-outcome source eligibility was used BEFORE fetching speciesKey facets. GBIF full speciesKey facets were paged to completeness and intersected with the independently sourced country GRIIS introduced plant key sets. A cell needs **>=5 GRIIS-observed introduced plant species keys and >=20 indexed introduced plant records in EACH period**, on top of the prior **>=100 all-Plantae records and >=3 source datasets per period**.

| Actual 200km macro block | Total nominal 50km source cells | Prior all-Plantae both-period source-qualified | Both periods fully speciesKey-paginated | BOTH periods sufficient introduced flora |
|---|---:|---:|---:|---:|
| JP Tohoku | 16 | 12 | 12 | **1** |
| US Midwest | 16 | 15 | 15 | **15** |
| FR Burgundy | 16 | 16 | 16 | **16** |
| **Pooled three fixed regions** | **48** | **43** | **43** | **32** |

**Interpretation:** The US and French blocks are comparatively source-rich in GRIIS-listed introduced plant observations, but the Japanese block FAILS most introduced-species opportunity thresholds despite passing many general plant opportunity thresholds. These are indexed record reporting opportunities, NOT a demonstration that Japan has few introduced plants, nor of actual rail corridor effects. Increasing record abundance or adding new climates cannot by itself rescue a **saturated binary railway-connection exposure**: all 16/16 JP, 15/16 US and 16/16 FR *rail-passage* cells share one mapped railway component in the first source network test. Consequently, these broad regions need a **continuous rail-network shortest path / circuitousness exposure** and a comparable source ROAD network path, not more binary connectedness tests.

There are **no rail–introduced plant biological effect sizes** for 50km equal-area cells yet. GBIF speciesKey records need actual DOI-cited occurrence files, duplicated-event QA, suitable coordinate uncertainty and country/ecoregion masks before high-confidence inference.

## 4. Macro analysis requirements, not completed results

- Keep eight-country geographically varied frozen blocks and later verify actual independent **ecoregion IDs**, land fraction, study-country overlap, ports, elevation/climate and land cover. Equal-area grid does not mean representative biogeographic sampling.
- Build full indexed GBIF speciesKey incidence after independent eligibility, with country GRIIS introduced-only status, and retrieve actual DOI-citable event rows for deduplication/coordinate uncertainty / observer effort if doing inference.
- Build physical rail and ROAD path networks in the same spatial block plus halo; without independent actual historical active/freight evidence, no train-mediated seed claim.
- Use actual independent source country×ecoregion regional block as unit; pairwise similarities within each region are dependent. Test true matching support, and report missingness and sign uncertainty. No p-hacking of calipers or naive 1st- vs 2nd-half species counts.

Sources:
- [v2.8 actual 10-city road network source](https://github.com/zuizui0223/tetsu/actions/runs/42737690925).
- [v2.10 actual 5-source posthoc ecological diagnostic](https://github.com/zuizui0223/tetsu/actions/runs/42738469033).
- [v2.9 actual 10-block 50km equal-area GBIF source census](https://github.com/zuizui0223/tetsu/actions/runs/42738021553).
- [Overture 2026-10-21 release](https://docs.overturemaps.org/blog/2026/10/21/release-notes/).
- [Overture road/rail segment connector topology](https://docs.overturemaps.org/guides/transportation/segments-and-connectors/).
