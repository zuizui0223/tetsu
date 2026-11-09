# v1.8 — Ten-country GBIF Plantae record source census (real API)

**Actual successful source run:** [Actions #42609761043](https://github.com/zuizui0223/tetsu/actions/runs/42609761043), 2026-11-09. All **20/20** source country × year windows retrieved, **six tests passed**, zero API errors. Aggregated source JSON artifact: gbif-macro-10country-2010-2025-source-counts-v18 (artifact ID 13668854141).

## Explicit population / caution

Frozen countries: JP KR TW US CA GB DE FR AU NZ. Two disjoint 8-year windows 2010–2017 and 2018–2025. Actual GBIF occurrence SEARCH counts, kingdomKey=6 (source-verified Plantae), species rank, latitude/longitude present, no GBIF spatial issue, country and year filters, all records vs human observations vs preserved specimens. Search request limit=0 retrieves *count only*, NOT original record rows and NOT any unique species count. The Plantae kingdom is **not a verified Tracheophyta-only set**, and no alien/native status, coordinate accuracy in meters or deduplication was checked.

| Country | 2010–17 all records | 2018–25 all records | Later / earlier | Human 2010–17 | Human 2018–25 | Preserved 2010–17 | Preserved 2018–25 |
|---|---:|---:|---:|---:|---:|---:|---:|
| JP Japan | 529,039 | 249,745 | 0.472 | 291,674 | 173,083 | 221,484 | 75,725 |
| KR South Korea | 327,723 | 525,137 | 1.602 | 260,639 | 471,775 | 65,248 | 53,183 |
| TW Taiwan | 140,914 | 1,170,453 | 8.306 | 95,907 | 1,153,520 | 25,452 | 13,559 |
| US United States | 3,469,078 | 23,798,357 | 6.860 | 2,343,471 | 23,147,219 | 933,684 | 559,144 |
| CA Canada | 308,199 | 4,752,318 | 15.420 | 181,887 | 4,692,903 | 103,403 | 19,243 |
| GB United Kingdom | 12,326,716 | 11,577,156 | 0.939 | 12,299,563 | 11,557,135 | 8,125 | 5,132 |
| DE Germany | 2,835,857 | 6,836,979 | 2.411 | 2,671,600 | 5,815,228 | 47,539 | 26,649 |
| FR France | 22,813,254 | 13,523,055 | 0.593 | 22,482,699 | 13,490,359 | 16,846 | 7,374 |
| AU Australia | 5,432,022 | 7,199,558 | 1.325 | 4,458,541 | 6,709,204 | 360,597 | 191,299 |
| NZ New Zealand | 264,295 | 1,254,081 | 4.745 | 185,342 | 1,070,131 | 64,632 | 30,972 |
| **TOTAL** | **48,447,097** | **70,886,839** | **1.463** | **45,271,323** | **68,280,557** | **1,847,010** | **982,280** |

**119,333,936 matching GBIF indexed occurrence records** total. Not 119M unique plant individuals, taxa, photos or downloaded rows. The public occurrence search record-retrieval API has a 100k overall cap; these count-only queries are not affected, but the underlying 119M occurrences have NOT been downloaded.

## Positive methodology-relevant result and limits

Source coverage succeeded in all ten countries and both time windows, with 2018–25 record counts ranging from 249,745 in Japan to 23,798,357 in the US (about 95-fold). The country changes range from 0.47x (Japan) to 15.42x (Canada), reflecting potentially enormous recording-coverage and source-composition differences. In particular, the millions of GB/FR HUMAN_OBSERVATION records cannot be assumed to represent equal individual observer or site efforts. **No claim that these changes represent railway-caused plant changes, non-native colonization, or true plant population growth.**

To enable actual railway macroecology:
1. Tracheophyta and country-specific alien-status crosswalk (regional status, not global species identity), accepted taxon keys, coordinate-uncertainty gates, site/dataset deduplication.
2. Obtain DOI-backed GBIF bulk downloads, equal-area 25/50km cell × taxon × period incidence, stratified target-group effort. The protocol's ≥100 valid occurrences and ≥3 datasets per cell-period remains a *threshold proposal*, not an achieved cell count.
3. Independently obtain active railway *graph topology*, roads, ports, environmental and observer controls. Road rail-distance pixels are not proof of train-carried propagules; freight-only is absent from some official rail layers.
4. Fit matched cell-pair community-similarity models for alien versus native assemblages with country/ecoregion spatial blocks and degree-aware rail-network nulls, report signed effects/uncertainty, not just predictive improvement.

**No new biological railway connectivity coefficient has yet been estimated.**

Source code: [auditor](../scripts/audit_gbif_macro_country_coverage_v18.py), [tests](../tests/test_gbif_macro_country_coverage_v18.py), [workflow](../.github/workflows/gbif-macro-country-coverage-v18.yml), [real Action #42609761043](https://github.com/zuizui0223/tetsu/actions/runs/42609761043).

GBIF source docs: https://techdocs.gbif.org/en/openapi/v1/occurrence and https://techdocs.gbif.org/en/data-use/api-downloads . Individual GBIF provider records have varying licenses and none was republished in this source-count-only pilot.
