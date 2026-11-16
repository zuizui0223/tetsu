# tetsu v5.5 — Macro railway–plant analysis: source reliability and next data acquisition

**Verified: November 10, 2026 (JST).** This document unifies actual source-backed findings from v5.0–v5.5. It is not a railway-mediated seed dispersal result.

## Research question and provenance

Do relatively direct RAIL physical routes (compared with ROAD) or high all-path electrical resistance network redundancy (Circuitscape-equivalent Q) correlate with greater observed country-GRIIS-listed introduced-plant species-key similarity across 50km EPSG:6933 equal-area regions? Can recording effort, source dataset composition, and true event-level location uncertainty account for the apparent association?

The 2026 Overture mapped physical network is **not** verified active freight service or railway-side habitat for 2010–2025 GBIF records. GBIF indexed positive RECORD facets do not establish biological absence, field plant abundance, or unique independent observations.

## 1. Real GBIF dataset-source concentration (v5.0)

[Actual six-country full GBIF DATASET_KEY/DETECTED facets #42796286437](https://github.com/zuizui0223/tetsu/actions/runs/42796286437) **SUCCESS** in six parallel jobs. The source-eligible 50km GRIIS original study cells number **83** (US Midwest 15, FR Burgundy 16, CA Ontario 10, GB Yorkshire 12, DE Saxony 14, AU Southeast 16), or **166 historic cell-period full dataset-key profiles**, all complete under original record-count and status gates.

Late-period 2018–2025 median source cell top-one GBIF dataset share of **ALL PLANTAE indexed records**:

| Country block | Median share of SINGLE largest indexed plant dataset |
|---|---:|
| Canada Ontario | 94.37% |
| US Midwest | 91.95% |
| Australia Southeast | 54.37% |
| UK Yorkshire | 53.46% |
| Germany Saxony | 45.13% |
| France Burgundy | 28.13% |

Source dataset UUID concentration is a plausible survey-composition bias, not a recorder-level effort estimate and NOT introduced-only dataset concentration.

## 2. Real RAIL/ROAD source association before and after recording-source controls (v5.1)

[Actual same-pair source GRIIS-flora / dataset-overlap join #42797029122](https://github.com/zuizui0223/tetsu/actions/runs/42797029122) **SUCCESS**, six synthetic tests. Each original 50km source cell pair's all-Plantae dataset-Key overlap is the fraction of source dataset record proportions jointly shared between the cells; it is an imperfect **all-Plantae** proxy, not introduced-only effort correction. The SAME source cell pairs were used before and after adding the proxy.

For 2018–2025 introduced-flora K20-record-rarefied expected-intersection/expected-union similarity, physical RAIL/ROAD shortest route partial ranks (negative source rho = more-direct RAIL associated with higher observed GRIIS introduced similarity):

| Covariate controls | Country-equal rho median | Hypothesis-compatible negative countries |
|---|---:|---:|
| Geography + ROAD physical shortest path + indexed GBIF source record balance | +0.021259 | 3/6 |
| Add all-Plantae source dataset overlap | **+0.056600** | 2/6 |
| Add overlap + top-dataset concentration imbalance | **+0.039974** | 2/6 |

No robust general railway-specific positive ecological effect appears at K20.

## 3. Actual fixed K80 recording-depth plus source-dataset overlap test (v5.2)

[Real six-country K80-eligible identical-cell-cohort source sensitivity #42797571041](https://github.com/zuizui0223/tetsu/actions/runs/42797571041) **SUCCESS**, four tests. The same GRIIS source cells are retained at K20/K40/K80; K denotes indexed GBIF RECORDS per cell in both periods, not plants or detection probability.

| Recent 2018–2025 fixed geographic source cohort | Median rho, old source controls | Median rho plus source dataset overlap | Hypothesis negative directions after overlap |
|---|---:|---:|---:|
| K20 | +0.021259 | **+0.056500** | 2/6 |
| K40 | −0.021333 | **+0.011045** | 3/6 |
| K80 | −0.064249 | **−0.069033** | **5/6** |

The K80 hypothesis-compatible 5/6 result is **posthoc**: K80 was selected only after looking at previous raw vs K20 source discrepancies. Five of six countries have negative descriptive partial ranks (US −0.193542, FR −0.058652, GB −0.020952, DE −0.230200, AU −0.079415; CA +0.140425), but there are no pairwise-independent p-values, no recorded train operations, and no validated mechanism. The earlier 2010–2017 period at K80 had median source rho **+0.049069**, only 3/6 negative, so temporal transfer is weak.

## 4. Outcome-withheld US Southeast region: source support FAILURE, not a botanical null

The [US Southeast geographically untouched K80 protocol](holdout_us_southeast_v53_protocol.md) was committed BEFORE new railway/road physical source outcomes were queried. US Southeast was frozen within the original 10-block 8-country 50km registry and was not in the original six country-region railway/flora correlations.

[Three actual US Southeast source inputs #42798094566](https://github.com/zuizui0223/tetsu/actions/runs/42798094566) all **SUCCESS**. Physical rail 4,713 source segments and 120/120 50km railway cell pairs mapped; physical road 647,309 segments and 120/120 cell pairs mapped; GRIIS K20 original source minimum 12/16 cells; datasetKey source complete for all 12. Physical rail median route/geodesic path ratio 1.619842, road 1.171015; late maximum dataset source share median 0.918335.

[Actual postfreeze main outcome support #42798827089](https://github.com/zuizui0223/tetsu/actions/runs/42798827089) was executed under the locked eligibility minimum >=8 source cells and >=20 mapped source pair comparisons at K80. Actual K80 conditions yielded **4 source cells and 6 dependent pairs**. Thus the K80 source response is **NOT ESTIMABLE**, not positive, not negative and not a p-value. Do not change to K20 just to report a favorable held-out result. Details [here](us_southeast_heldout_K80_source_feasibility_v53.md).

## 5. Actual K80 capability across ALL ten originally frozen source regions (v5.4)

[Actual pre-outcome all-10-country-block source-depth screen #42799378672](https://github.com/zuizui0223/tetsu/actions/runs/42799378672) **SUCCESS**, four tests. All ten original 200km region boxes × 16 cells, spanning eight countries, were checked using immutable country-GRIIS speciesKey indexed source counts in both 2010–2017 and 2018–2025. No physical railway effect inspected.

| Per-cell minimum historic indexed introduced records | Qualified 50km source cells / 160 | Geographic country blocks with >=8 qualified source cells |
|---:|---:|---:|
| K20 | 105 | 7 |
| K40 | 95 | 7 |
| K80 | **85** | **6** |
| K160 | 69 | 5 |
| K320 | 62 | 3 |

At K80, the ONLY six blocks with >=8 eligible botanical source cells are US Midwest, Canada Ontario, UK Yorkshire, France Burgundy, Australia Southeast and Germany Saxony — exactly those used in the earlier posthoc six-country study. The other four original source regions have K80 eligible cells: US Southeast **4**, NZ North **3**, Japan Tohoku **0**, Japan Kyushu **0**. All are source recording-quality deficits, not confirmed ecological absence of nonnative plants.

**Therefore the existing ten-block source dataset cannot supply an unused K80-ready out-of-analysis geographic region.** More pairwise network calculation alone will not solve this.

## 6. OFFICIAL GBIF original DOI-bearing DWCA download request manifests: PREPARED, NOT SUBMITTED

[Real v5.5 manifest generator #42800153221](https://github.com/zuizui0223/tetsu/actions/runs/42800153221) **SUCCESS**, four offline source geometry and no-credential tests. [Generator](../scripts/prepare_GBIF_original_DWCA_download_manifests_v55.py) and [workflow](../.github/workflows/GBIF-DWCA-manifests-v55.yml) produced 10 official-format JSON predicate request files, one for every fixed geographical block, plus an SHA256 source index.

Each request is for Darwin Core Archive (DWCA), 2010–2025 inclusive, the original 16 anti-clockwise WGS84 50km source cell geometries with OR-within, country identity, Plantae kingdomKey=6, taxon rank species, geographic coordinates and no GBIF geospatial issue. **Legacy occurrenceStatus=PRESENT is omitted**, because previous current API tests show it returns zero despite original DETECTED indexed facet records. Detection/absence must be revalidated on actual full archive event rows.

**These files are prepared request PREDICATES ONLY:** no GBIF authenticated request has been posted, no original occurrence archive has been obtained, and NO DOI has been issued or verified. The user must authorize and authenticate a GBIF download with their registered GBIF username and password through a secure channel or GBIF's own website; do not put credentials in this repository or conversation. GBIF itself assigns a DOI when its asynchronous download succeeds, and cited datasets must be acknowledged.

Official GBIF reference:
https://techdocs.gbif.org/en/data-use/api-downloads

## 7. Decision

**Not ready for a causal railway-mediated plant dispersal publication.** The mapped rail/road shortest paths and Circuitscape-equivalent physical alternative-path redundancy are real measurable infrastructure quantities. The floristic observation source, however, varies strongly by dataset, record depth and geographic region. The 5/6 K80 source sign is posthoc, the independently chosen US Southeast K80 source lacks sample support, and no unused pre-frozen 10-block source region can provide K80 replication.

The next meaningful scientific bottleneck is **actual DOI-backed original GBIF event rows with status, independent observer/event/coordinate uncertainty, source country GRIIS taxonomy and REAL active historical freight routes**, followed by pre-frozen independent terrestrial ecoregion controls. Do not choose a rarefaction K, source region or electrical resistance scheme because it yields the desired sign.

No railway grass cutting adaptation, wind selection, flowering shifts, train seed contamination or actual seed delivery to new habitats has been measured in this source work.
