# v7.0–v7.1 — Railway delivery × functional filtering macroecology: measured first-stage

**Status:** executable pipeline, not an affirmative seed-dispersal/adaptation result until the live GIFT source workflow successfully produces real source aggregates. The separate `v6.9` French original GBIF source workflow had a GitHub Actions numeric-leading job ID corrected on 2026-11-11; success must be confirmed independently.

## Novel, testable ecological question

> Across independently sampled countries, does mapped physical railway accessibility correspond to more similar introduced-flora assemblages **within different plant growth forms**, beyond physical roads and geographic/recording-opportunity covariates?

This is a first-stage **trait-stratified community-composition** question. It does *not* establish that any propagule was carried by a train, that any vegetation was actually cut, or that growth form arose through local railway adaptation. Subsequent documented mowing and air-flow exposure will require a separate habitat-specific field or high-resolution spatial data source.

## Exact non-posthoc existing source frame

- Six geographically preselected, **200 × 200 km** region windows: US Midwest, France Burgundy, Canada Ontario, UK Yorkshire, Germany Saxony, Australia Southeast; each 4×4 nominal 50km equal-area EPSG:6933 cells, country mask ≥80%, two indexed GBIF source periods 2010–17 and 2018–25.
- Existing verified country-**GRIIS-listed introduced** speciesKey × cell × period faceted-record source: run [#42761140912](https://github.com/zuizui0223/tetsu/actions/runs/42761140912). Countries with missing GRIIS designation are **unknown**, not native.
- Existing verified Overture Oct2026 physical RAIL and ROAD shortest-path graph source used for v4.4. Historical active freight/passenger operation is unknown; current graph is not automatically an exposure in the earlier years.
- **New independent source:** public [GIFT API](https://biogeomacro.github.io/GIFT/articles/web_only/GIFT_API.html), *growth form* metadata ID **1.2.1**, as documented in the [GIFT tutorial](https://biogeomacro.github.io/GIFT/articles/web_only/GIFT.html). The source driver retrieves GIFT taxon names and trait records, requires categorical source agreement ≥0.66, and excludes records derived solely from biased references/derivations (`biasref=0`, `biasderiv=0`). GBIF numeric speciesKey gets a strict source-current accepted Plantae species canonical binomial; only exact binomial joins to GIFT accepted work_species. Synonyms, changed taxonomy, ambiguous traits and name collisions are unknown/excluded and counted. No genus/family backfill.
- Original GBIF source record count minimums remain unchanged. For each group, a cell needs ≥3 separately observed trait-matched listed species categories and ≥20 positive indexed records; a country source trait coverage below 35% suppresses interpretation (reported as unestimable).

## What is actually calculated

1. For each plant growth group, **herbaceous vs woody**, apply the original source frequency-based 20-record rarefaction *within that group* to each qualifying cell. This means similar source records, not controlled presence/absence or independent plot surveys.
2. Compare pairwise **expected-intersection / expected-union** trait-group source similarity against `log(mapped RAIL shortest route / mapped ROAD shortest route)` using exactly the same two rail anchor definitions as v4.4. Control rank of geographic distance, mapped road distance and two GBIF all-plant record count proxies.
3. Report one **region-level** descriptive partial rho per country and group; all pairs within a country share cells and are **not independent n**. Never calculate pairwise p-values.
4. Only if **the same exact source cell pairs** satisfy both growth groups, compute a within-country paired group-difference diagnostic `rho_herbaceous − rho_woody`; otherwise show missing, not zero. Country-level summaries remain descriptive until enough independent blocks pass source gates.
5. Report the denominator for every source eligibility failure, taxon matching failure, unknown GIFT category and country outcome. All real aggregates are CI artifacts, not fictional local field adaptations.

**Interpretation:** A stronger negative rail/road partial rho among herbs, if robust across independent country blocks and to recording-source exclusion, is consistent with a **growth-form-dependent recorded assemblage association**. It is NOT sufficient to infer mowing selection or slipstream adaptation. Woody/herbaceous also does not identify resprouting strategy, flowering phenology, or dispersal mode.

## Next biological exposures (strictly separate)

| Biological component | Actual exposure necessary | Acceptable measure | False shortcut to reject |
|---|---|---|---|
| Propagule input | Dated active rail/freight, commodity movement or direct vehicle/station seed traps | Train movements/tonnage/freight terminals with start dates; species-specific seed trap assays | Mapped 2026 rail route = observed seed flux |
| Mowing filter | Actual cutting episodes at at-grade rights-of-way and road-verge controls | Geolocated maintenance logs, repeated same-plot vegetation cover/phenophase | 10m satellite vegetation class = mowing |
| Slipstream stress | Air-flow exposure at fixed distance from operating tracks | Safe operator-authorized gust sensors with controls and train passages | Train speed assumed from track density |
| Phenotypic plasticity | Same individual or conspecific population pre/post known event | Repeat photo masks in fixed frames and visit dates | Height from unscaled opportunistic photo |
| Adaptive evolution | Heritable variation or genotype contrasts across independent lineages | Common garden, genomic/environment association with demography control | Species mean growth form = local evolution |

The eventual **transport × stress × trait** interaction is only identifiable in spatial blocks where all three independently observed components align. It must *not* be estimated from source physical RAIL km alone.

## Reproducibility

- `scripts/harvest_gift_growth_form_v70.py`: public original trait source, strict taxonomy; output provenance and per-country matching denominator.
- `scripts/analyze_sixcountry_functional_rail_v71.py`: same-cell-pair group-specific 20-record-rarefied network contrast and independent-country summary.
- `tests/test_rail_functional_macro_v70_v71.py`: synthetic contract/fail-closed and geography tests; no tests are biological evidence.
- `.github/workflows/GIFT-sixcountry-functional-rail-v71.yml`: real frozen six-country GRIIS × 2026 Overture paths × newly sourced GIFT traits; CI artifacts retain results and missingness.

**Scientific stopping rule:** If the source GIFT API schema/taxonomic mapping fails, only 6-country source-map coverage was tested, not an ecological effect. If most countries fail ≥35% trait coverage or group paired-cell gate, do not reinterpret unavailable group as effect zero, relax threshold after inspecting signs, or claim a cross-country trait mechanism.