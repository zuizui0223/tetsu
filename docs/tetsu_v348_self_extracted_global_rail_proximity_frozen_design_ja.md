# v348 — 自前抽出する全球鉄道近傍植物（2026-11-24、結果確認前に凍結）

## Why this replaces the old literature-only shortlist
Previously available railway floras (CZ39, IT28, HRBA15, SI89, SK30, Kashmir, WA120) were selected *because original author plot-by-species datasets were found*, not from a geographically representative global sampling design. They cannot define the sampling frame for a 'global' ecological test. Instead, create an investigator-defined, **reproducible source-extraction frame** across all six permanently inhabited continents, independent of observed species composition or previously favourable effects.

### Pilot frame (FIXED before looking at any results)
A predeclared urban railway anchor in six broad biogeographic regions, 2025 calendar-year public wild vascular-plant observations:
- Tokyo, Japan (Asia) — approx anchor lat 35.6812, lon 139.7671
- Prague, Czechia (Europe) — approx anchor lat 50.0839, lon 14.4350
- Cape Town, South Africa (Africa) — approx anchor lat -33.9258, lon 18.4230
- Sydney, Australia (Oceania) — approx anchor lat -33.8830, lon 151.2070
- San Francisco, United States (North America) — approx anchor lat 37.7764, lon -122.3940
- São Paulo, Brazil (South America) — approx anchor lat -23.5357, lon -46.6340

Window fixed as center +/-0.012° longitude and latitude (~2–3 km each direction, latitude-varying actual area). One city/window per continent is **only an API/source/sampling opportunity pilot**, NOT independent replication of climate regions or representative of continents. No site was selected by plant temperatures or positive railway effect. Cities are heterogeneous in green space, urbanization and sampling practices.

### Two primary public source systems, never conflate
1. OSM actual way geometry, `railway=rail`, retrieved through Overpass GET (`out geom`) with an official OSM API map-XML fallback for *small* original bbox. Discard segments tagged `tunnel=yes`, `bridge=yes` and `service=yard/siding/spur`; record counts and route/source status; OSM rail geometry says a mapped railway, **not active train count, accessible rail ballast, precise trackbed boundary, or historically valid 2025 rail operation**. Persist no original detailed track geometry in GitHub outputs.
2. iNaturalist official `api.inaturalist.org/v1/observations`: 2025 whole-year, research-grade `iconic_taxa=Plantae`, wild `captive=false`, public georeferenced photos not required, per_page=200, bbox exact frozen; verify record `obscured != true`, `geoprivacy` open or null, taxon `rank=species` and exact binomial, `positional_accuracy` finite 0–30m, WGS84 valid. Report total_results; fetch all results **only if ≤1,200** (at most six sequential pages, API <=10k hard cap); otherwise source is *over cap*, do not treat first pages as a complete spatial sample. Respect rate guidance ≤1 request/sec and robots/HTTP 429. Keep no photographer handles, observer identifiers, original IDs, photos, exact coordinates, or raw taxon×location records in source artifacts.

### Frozen spatial design (site-level *proximity* NOT community plots)
Compute local metric distance to **all mapped OSM rail way segments** within each window using equirectangular projection centered on each study anchor. The inclusion decision uses *reported positional accuracy* `u`:
- strict `RAIL_NEAR`: distance to nearest rail centerline `d + u ≤ 80m`;
- strict `BACKGROUND_FAR`: `d - u ≥ 300m` **and** `d + u ≤ 1000m`;
- all other observations `uncertain buffer` excluded. If OSM source absent, no classification, NA.
- Pilot descriptive support: ≥10 near observations with ≥5 unique source taxa; ≥30 background with ≥10 distinct species; site must have real complete original API data; classify "insufficient" explicitly, not statistical zero.
- Matched local background should *eventually* constrain observation year/season, observer effort, land cover, distance to road and study plot area; the present urban 2025 API sample has none of those controlled. **Do NOT call the target a rail-biotope flora, valid species absence, habitat filter, dispersal rate, railway selection, adaptation, or causal/independent field test.** Citizen-science near rail often lies outside the right-of-way; geotags near station may describe platform gardens or street trees. The record pool is highly observation-effort biased.
- For positive evidence of an actual railway biological effect, require real matched on-track/near-track vs genuine outside-railway *standardized vegetation plots* with method, plot area, managed state, survey year, and direct original source species-by-plot matrix. Pilot observation proximity remains descriptive only.

### Source quality and scientific stops
- Output source counts by continent (HTTP retrieval success, mapped OSM actual rail ways/segments, iNat total count, full-pagination support, accuracy/cultivation/obscured drops, near and far observed unique species, exact support status); source URLs, fixed design and original public data citations. Raw coordinate geometry, source observers, IDs or rare species names are not republished.
- Do not extrapolate from 6 city samples to 'worldwide frequency of warmth'. A scientifically valid global study should predeclare multi-city *stratified probability sampling by climate × country × railway-use × urban/rural × protected areas*, then audit location, comparators and sampling effort. This first pilot checks whether observation-based study is even feasible.
- **Thermal affinity:** European EIVE does NOT apply to this six-continent pilot. The desired species global EIDC/GBIF climatic envelope trait remains blocked on original source ZIP HTTP, and no attempt to construct a trait from one city is valid. No global thermal effects or C4 adaptation asserted until an independent, complete global thermal species table is authenticated.
- Any source error/missing/429 is logged as **data unavailable**, never a negative biological effect. No opportunistic favorable city selection after results. OSM (ODbL) and iNaturalist source policies must be cited, city-level counts only.

### Implementation
`scripts/pilot_self_extracted_six_continent_rail_flora_v348.py`, tests `tests/test_self_extracted_rail_flora_v348.py`, CI `.github/workflows/tetsu-v348-six-continent-OSM-iNat-pilot.yml`. Research branch stacked on prior v347 PR #55, main not changed, raw publisher data not reproduced.
