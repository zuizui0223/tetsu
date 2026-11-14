# tetsu v4.7–v4.9: rarefaction-depth sensitivity and GBIF true detection/status source audit

**Verified original results: November 10, 2026 (JST).** This methodological audit directly tests the primary source-observation sensitivity behind the 6-country tetsu GRIIS-introduced flora × physical RAIL-vs-ROAD network association. No actual propagule dispersal, rail freight or genetic adaptation has been verified.

## 1. Why 20 indexed GBIF RECORDS was a consequential choice

Earlier v4.6, same six independent 200km country macro blocks and 50km equal-area cells, found **5/6** railway-compatible partial rank directions with RAW observed indexed alien speciesKey Jaccard but only **3/6** railway-compatible directions when normalizing to expected species overlap under 20-record without-replacement rarefaction. Species observations are GBIF indexed positive RECORDS, not independent plants, systematic ecological plot occurrences or biological absence.

Actual original-source [v4.7 K20–K320 and fixed K80 source trial #42791684942](https://github.com/zuizui0223/tetsu/actions/runs/42791684942) succeeded, seven synthetic numerical and source quality tests. The source macro regions (US Midwest, FR Burgundy, CA Ontario, GB Yorkshire, DE Saxony, AU Southeast), both periods, physical rail/road route definitions and two outcome-blind source connector anchors were held fixed. Rarefaction probability was recomputed by numerically stable product sum of log1p terms rather than potentially cancellation-prone gamma differences; K20 reproduced prior aggregate rho **+0.021259**.

### 2018–2025 source-depth sensitivity, source-cell opportunity can change as K increases

| K indexed introduced records per source cell | Independent country blocks with geographic/network source support | Number of hypothesis-compatible (NEGATIVE) country rho values | Equal-country median descriptive partial rank |
|---:|---:|---:|---:|
| 20 | 6 | 3 | **+0.021259** |
| 40 | 6 | 4 | **−0.016212** |
| 80 | 6 | 4 | **−0.064249** |
| 160 | 5 | 3 | **−0.061336** |
| 320 | 3 | 2 | **−0.143032** |

The x-axis is K **records**, not equal independent species detections. For K160/K320, new source record thresholds remove geographic cells and country blocks, so apparent trend across K can be due to changing study coverage. K320 has only **3** independent country blocks and is not used to argue a general effect.

### Same-cell controlled K80 source cohort (the important quality control)

Keep the EXACT same cells with >=80 indexed GRIIS-introduced source records in *both periods* and evaluate them at K20/K40/K80. In 2018–25, all six country blocks can be compared:

| Same eligible 50km cells fixed from K80 source opportunity | Country blocks | Hypothesis-compatible country directions | Median partial rank |
|---|---:|---:|---:|
| K20 source rarefaction | 6 | 3/6 | **+0.021259** |
| K40 source rarefaction | 6 | 4/6 | **−0.021333** |
| K80 source rarefaction | 6 | 4/6 | **−0.064249** |

**Conclusion:** A sign change occurs even without changing the set of source cells. K20 by itself is arbitrary and under-informative. Going to K80 weakens the apparent counter-railway signal but does NOT support a strong general railway-dispersal effect. No country-weighted hypothesis p-value or causal seed transport result can be inferred from these pairwise dependent source associations.

A second genuine limitation is that RECORD rarefaction does not equalize datasets, observers, survey effort, timestamp density, spatial sampling routes or geographical coordinate uncertainty.

## 2. Critical source contract issue: old PRESENT vs current DETECTED / NOT_DETECTED

The original v4.8 GET-search raw record audit returned **0 original rows in 24 pre-fixed 50km country×period bins** despite previously indexed source speciesKey counts >0. This is a source semantic contradiction, not actual botanical absence.

Traceable real-source diagnostics:
- [v4.9 mode comparison #42792729516](https://github.com/zuizui0223/tetsu/actions/runs/42792729516) SUCCESS: first page, limit=1/300 and source facet limit=0 all report zero WITH legacy occurrenceStatus=PRESENT.
- [v4.9b pre-existing source-reconciliation run #42793141570](https://github.com/zuizui0223/tetsu/actions/runs/42793141570) intentionally FAILS strict historical revalidation, showing that dropping occurrenceStatus reproduces original GBIF country-cell source counts: US Midwest x1y1 old/new **68,254 / 68,254**, French Burgundy x1y1 **54,731 / 54,733**. Removing any other single tested filter did NOT restore results.
- [v4.9c explicit OCCURRENCE_STATUS facet #42793539814](https://github.com/zuizui0223/tetsu/actions/runs/42793539814) (program intentionally FAILS strict invalid API-filter contract): US source facet **DETECTED=68,254**; FR source facet **DETECTED=54,731, NOT_DETECTED=2**. Thus old count matches DETECTED exactly for the tested source pair.
- [v4.9d current GET status filter tests #42793849132](https://github.com/zuizui0223/tetsu/actions/runs/42793849132) SUCCESS of the *diagnostic*: GET occurrenceStatus=DETECTED or NOT_DETECTED yields HTTP **400** in both tested cells, whereas legacy PRESENT/ABSENT queries return 0 records. Direct original row GET data do **not** expose a normalized occurrenceStatus field in these tested first pages, despite status facet classification.

**Implication:** Do not use occurrenceStatus=PRESENT in new GBIF original search, and do not silently substitute unfiltered GET results as if individually DETECTED. Instead, audit exhaustive DETECTED/NOT_DETECTED *facet* counts in the same source cell-period and explicitly label raw first-page samples as unfiltered. Historical positive facet numbers are source-consistent with DETECTED in inspected original fixed cells, but not every country-cell has undergone the exhaustive status verification. Future DOI-backed complete GBIF occurrence download must resolve actual statuses at the row grain.

This is a **real API query/enum compatibility problem** discovered by crosschecking source cohort counts; it must be kept as a reproducibility warning in the paper. GBIF documented general occurrence statuses, and the official API supports source faceting, but recent returned indexed categories must be verified rather than assumed equal to older string filters.

## 3. Actual original GBIF row metadata diagnostic (reconciled v4.8b)

[Successful source v4.8b #42794536402](https://github.com/zuizui0223/tetsu/actions/runs/42794536402) sampled exactly 2 pre-fixed 50km cell *positions* (x1y1, x2y2) × two historical 8-year periods × six countries: **24 cell-period source strata**. Five synthetic source-schema tests passed. For EACH stratum, source DETECTED/NOT_DETECTED total counts were obtained exhaustively by STATUS FACET, then the first **300** unfiltered GBIF original occurrence rows were inspected *in memory*, with only metadata aggregates stored (no original observation IDs, human recorder names or precise coordinates published).

| Country, one 200km regional block | First-page original rows inspected | First-page rows missing coordinateUncertaintyInMeters | First-page rows missing recordedBy | Full 4 fixed-stratum source DETECTED indexed records | Full source NOT_DETECTED |
|---|---:|---:|---:|---:|---:|
| US Midwest | 1,200 | **468 (39.0%)** | 4 | 92,262 | 0 |
| FR Burgundy | 1,200 | **121 (10.1%)** | 119 | 381,287 | 81 |
| CA Ontario | 1,200 | **180 (15.0%)** | 3 | 108,175 | 0 |
| GB Yorkshire | 1,200 | **54 (4.5%)** | 74 | 406,038 | 86 |
| DE Saxony | 1,200 | **164 (13.7%)** | 130 | 46,596 | 0 |
| AU Southeast | 1,200 | **821 (68.4%)** | 65 | 178,719 | 17 |
| **TOTAL** | **7,200** | **1,808 (25.1%)** | **395 (5.5%)** | **1,213,077** | **184** |

**These are NOT global missing-coordinate rates.** Every one of the 24 source strata had >300 indexed records and GET returns the first 300 in an unspecified/nonrandom ordering. The full-status DETECTED/NOT_DETECTED counts, in contrast, cover ALL indexed matching records in those 24 precisely defined source-cell periods. Missing coordinateUncertaintyInMeters is a metadata limitation, not proof of inaccurate location. Duplicate rounded date/point/species identifiers are possible multiple legitimate observations, not confirmed false duplicated GBIF records. The source status facet reassessment reconciled prior GRIIS source record denominators within its fixed source tolerance for all 24 inspected slots.

Actual source quality script [v4.8b detected+raw-row QA](../scripts/audit_gbif_detected_original_rows_v48b.py), synthetic privacy/source tests [here](../tests/test_gbif_detected_original_rows_v48b.py), strict 24-bin workflow [here](../.github/workflows/gbif-detected-original-rows-v48b.yml).

## 4. Research implications

1. **The apparent 5/6 RAW Jaccard railway signals are method-sensitive.** Moving K20→K80 with exactly the same source-eligible country cells gives 3/6→4/6 negative route-rho directions and median +0.021→−0.064. Neither represents strong biogeographically universal railway-mediated plant homogenization; a signed partial rank is not a causal train effect.
2. **A large recording opportunity is not a robust observation process.** Full source DETECTED status is abundant in fixed cells, but coordinate uncertainty metadata omission varies strongly among initial source pages from 4.5% GB to 68.4% AU. GBIF all-Plantae indexed rows were NOT cleaned by actual coordinates, observer effort or dataset balancing in the previous ecological partial-rank models.
3. **Country GRIIS status is country-reported and incomplete.** Nonmembership remains unknown rather than native. Repeated event/source row filtering is needed before any true species occurrence/habitat/ecoregion comparison.
4. **Source method priority:** obtain original GBIF occurrence downloads with a **DOI and source dataset citations** (authenticated GBIF download); deduplicate by GBIF record IDs and publishable source event/observer IDs, check country GRIIS taxonomy, coordinate uncertainty and history, then build comparable country×ecoregion source-aware plant incidence. The bounded original-row GET API audit here is not a replacement.
5. **Rail ecological hypothesis remains unconfirmed.** Overture physical rail/road graphs are 2026, postdating 2010–25 GBIF observation records. True active historical freight connections, station transfers, repeated flowering/disturbance traits, propagule traps and independent ecoregion blocks are absent.

GBIF technical references:
- https://techdocs.gbif.org/en/openapi/v1/occurrence
- https://techdocs.gbif.org/en/data-use/api-downloads
- https://www.gbif.org/citation-guidelines
- https://gbif.github.io/gbif-api/apidocs/org/gbif/api/vocabulary/OccurrenceStatus.html

### Evidence readiness

**READY as a source-data-quality diagnostic with explicit nonrandom-sampling limitations. NOT READY as a railway-mediated seed dispersal study.** The original result direction varies with raw vs rarefied indexed recording definitions, the GBIF occurrenceStatus query API is not currently reproducible using the legacy PRESENT request, and true field-observed railway propagule movement remains unavailable.
