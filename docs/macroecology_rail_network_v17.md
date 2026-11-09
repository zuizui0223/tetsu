# tetsu v1.7: macroecology of railway networks and plant filtering

Status (2026-11-09): analysis protocol and eight-new-region source-audit implementation. The 2011 JR East historical secondary reanalysis and 2021 satellite land-cover source results are real; railway-mediated dispersal, mowing effects, slipstream adaptation, evolutionary selection and flowering shifts remain UNTESTED.

## One central question

Across multiple independently sampled regions and plant species, do rail networks explain non-native plant redistribution beyond roads, geography, climate, land cover, trade/ports, urban development and botanical sampling?

Separate THREE processes instead of combining incompatible units:
1. Network delivery: spatial community composition and species incidence at 25–50 km macroecological grain.
2. Establishment filtering: independently verified active railway verge vs road verge and forest/open edges at 10–30 m or fixed quadrat scale.
3. Phenology / plasticity: repeated fixed-camera plots and complete monitored flower YES/NO visits with known cutting/maintenance dates and pre-flowering climate. Secondary, not proof of seed transport.

Transport does not imply establishment; establishment does not prove a seed or spore arrived on a train.

## Actually supported discovery baselines

- Abe et al. 2012 published JR East 2011 109-section questionnaire reanalysis: 1,334 taxon-by-habitat occurrence mentions. Equisetum arvense represents 24/70 (34.3%) trackbed mentions and 22/642 (3.4%) embankment slope mentions, a 10.01x relative compositional share. Both source warm/cold region aggregates show this direction. These are not independent plot records, and no causal effect is identified.
- Three fixed region GIS source normals (Sendai/Niigata/Nagano), 370 sites: 147 ESA 2021 tree-class, 24 grass-class, 117 built-up, 73 cropland, 8 water, one bare/sparse at 35m from railway line. Strict rail-road separation kept 9 tree and 1 water of 10; none mapped grass. Relaxed posthoc cutoff yielded two map-grass reconnaissance leads. **Strong source confound, not a railway habitat effect.**
- USA-NPN matched first-YES and Daymet sources have only two independent rail-near sites per species/year, changing temperature contrast signs; no flowering-onset effect.
- Single Utsunomiya LRT satellite NDVI sensitivity can change sign with clear-day/pixel selection; not a railway-caused effect.

## Primary macroanalysis: rail connectivity and alien plants

Spatial units: equal-area 25 km hexagonal (50 km sensitivity), nested in country x ecoregion blocks. Calendar windows: 2010–2017 and 2018–2025 for record effort, NOT a before/after treatment without independently verified rail/road historical dates.

Source ingestion contract:
- GBIF DOI-citable bulk occurrence downloads: country, date, species/accepted name, coordinate uncertainty, basisOfRecord, dataset key, taxon key, duplicate keys. Public search API caps total retrievable results at 100,000; do not confuse API-capped subsamples with entire GBIF flora. Large analyses must use GBIF authenticated bulk downloads and cite DOI.
- Regional alien-status checklists (Japan 1,753 source taxa previously audited; GRIIS or vetted national checklists for other countries). Never assign alien status from global species identity alone; unknown is missing, not native.
- Active railway topology and road network from OSM/Geofabrik or official licensed layers, preserving tunnel/bridge and railway lifecycle tags. OSM railway usage and service labels do not measure actual freight movements; Japan N02 lacks freight-only lines. Do not assume a passenger-only label or train count from absence of freight tagging.
- WorldCover 2021 background cover, climate variables, elevation, population, city/port exposure, coastal access and botanical reporting effort.

Pre-outcome effort eligibility: at least 100 de-duplicated eligible vascular-plant records per spatial cell-period across at least 3 independently sourced datasets; sensitivity at 50 and 200. Track exclusions by cell and rail exposure. Assess coordinate-uncertainty caps 100m / 1km / 5km separately. With 1km precision, do NOT classify a record as a 35m railway weed.

Primary comparisons:
- Non-native species composition: rarefied target-group-adjusted shared assemblage similarity / incidence between pairs of cells, and native species composition as negative-control sensitivity. A first GBIF record is NOT arrival or colonization date.
- Match cell pairs within country and ecoregion on 50–300km Euclidean distance, climate, urban exposure, port/coast access, road accessibility and target-group recording effort. Compare differences in rail **network** connectedness given similar ROAD connectivity.
- Fit ecological association with country/ecoregion effects and spatial blocking. Assess signed rail coefficients and uncertainty, not only prediction improvement.
- Leave-one-country/block-out reproducibility; blocked graph-rewiring/nulls preserving rail degree/length and spatial separation. Cell pairs sharing a site or route are not independent replicates. Use native/non-rail negative controls.
- Only after external freight evidence, test whether independently documented freight/port links amplify alien composition association. Association is NOT physical proof that seeds traveled on trains.

Macro launch target (proposed, NOT achieved): >=5 countries, >=20 country-by-ecoregion independent geographic blocks, >=500 sufficiently observed cells, >=200 regionally verified alien taxa, and adequate rail-connected/disconnected road-matched contrasts. Fail closed and report insufficient coverage otherwise.

## New eight-region replication of source selection bias (implemented)

Discovery set: Sendai, Niigata, Nagano 370 original normals; do not reuse in holdout replication.

Held-out geographic windows chosen BEFORE checking 2021 WorldCover classes: Sapporo, Akita, Tokyo-west, Nagoya, Osaka, Hiroshima, Fukuoka, Kumamoto. Different official N13 first-level meshes, validated map boundaries and projected UTM zones.

Code: scripts/audit_japan_rail_macro_holdout_worldcover_v17.py
Unit tests: tests/test_japan_rail_macro_holdout_worldcover_v17.py
Live source workflow: .github/workflows/japan-rail-macro-holdout-v17.yml

The 75m road separation threshold is the original v1.6 criterion, pre-fixed before new data. 25/50/100m are sensitivity screens. Read actual ESA WorldCover 2021 map source pixels and paired N02 passenger rail/N13 roads while suppressing exact coordinates and raw copyrighted asset uploads. At each region, report 2021 land cover classes among all geometry-frozen railside normals and among road-control-eligible normals. Only regions with >=20 fixed normals and >=3 strict eligible normals contribute to the descriptive independent-region sign/median summary. Zero eligible is missing comparison, never zero biological impact.

A repeated increase in tree-class representation would REPLICATE THE SAMPLING TRAP, not a railway ecological effect. Strict 35m centerline proximity cannot establish track-bed grass, actual at-grade operation, rail property, wind or mowing exposure. Publish region summaries and error audits only, not rail geometries or precise sites. No pseudoreplicated transect-level significance claims.

## How to make habitat traits and phenology macroecological

Stratify verified at-grade ground vegetation plots across ecoregions, rail traffic categories and maintenance types, with road verges / forest edges / unmanaged background comparison and corridor-cluster effects. Site-verified field species, consistent quadrat ROI and known treatment logs needed to test regrowth, shortness or flowering. Do not extract physical plant size or genetic adaptation from unscaled citizen photos. Equisetum spore cones are not angiosperm flowers, and spore/rhizome movement is different from seed transport.

For phenology, complete repeated non-flowering AND flowering observations, fixed species, thermal time, year and observer effort are required; first flowering photograph is not first flowering.

## Current evidence gate

**Real-source holdout completed on 2026-11-09**: GitHub Actions [#42608370719](https://github.com/zuizui0223/tetsu/actions/runs/42608370719), all 6 targeted tests PASS, full 8/8 source regions, zero source errors, aggregated downloadable workflow artifact. Across **2,169 new original rail-normal GIS positions** (not plant sampling plots), ESA 2021 tree cover is **453/2,169 = 20.9%**, mapped grass **63/2,169 = 2.9%**, built-up **1,374/2,169 = 63.3%**. Within the original strict road75 GIS condition, **28** source positions remain, **22/28 = 78.6%** tree-cover class and **1/28** grassland. These are **sampling selection proportions**, not comparable railway-vs-road ecological effects.

Original 8-new-region data:
| Holdout region | All GIS normals | ESA tree cover all | Strict road75 eligible | Strict tree class | Strict minus all tree fraction |
|---|---:|---:|---:|---:|---:|
| Sapporo | 116 | 40 | 4 | 2 | +0.1552 |
| Akita | 175 | 86 | 11 | 8 | +0.2358 |
| Tokyo-west | 682 | 100 | 5 | 5 | +0.8534 |
| Nagoya | 257 | 26 | 0 | 0 | undefined, not zero |
| Osaka | 358 | 58 | 2 | 2 | +0.8380 (below prespecified >=3 replication gate) |
| Hiroshima | 193 | 52 | 2 | 2 | +0.7306 (below gate) |
| Fukuoka | 301 | 68 | 2 | 1 | +0.2741 (below gate) |
| Kumamoto | 87 | 23 | 2 | 2 | +0.7356 (below gate) |
| **8-region source census** | **2,169** | **453** | **28** | **22** | not a pooled within-region causal estimate |

Under the **pre-outcome region support gate** (>=20 all / >=3 eligible strict), only **3 regions** are informative and **all 3/3 show enrichment**, with median region-level tree-fraction change **+0.235844**. Including two-eligible regions would give 7/7 directionally consistent, but that relaxes a fixed minimum-replicate gate and is **exploratory only**. Cannot claim statistical proof of universal forest enrichment, nor any botanical selection or railway adaptation. The new sites are 8 additional source regions (not random biomes), and the 370 v1.6 normals are a **discovery cohort**; together 11 regional boxes and 2,539 GIS site positions, never equate to 2,539 independent rail ecological experiments.

**Core network/plant macro analysis still awaits actual GBIF DOI-backed plant-by-cell data and independent rail/road topology, so no railway-dispersal effect size or p-value exists.**

Source/license reference:
- GBIF: https://techdocs.gbif.org/en/openapi/v1/occurrence and https://techdocs.gbif.org/en/data-use/api-downloads
- OSM: https://wiki.openstreetmap.org/wiki/Railways , ODbL https://www.openstreetmap.org/copyright
- ESA WorldCover 2021: https://esa-worldcover.org/en/about/about , DOI 10.5281/zenodo.7254221, CC BY4.0
- GRIIS country checklists: https://www.gbif.org/dataset/1288ee7d-d67c-4e23-8d95-409973067383
- Official Japanese N02 2025 and N13 2024 source URLs: see audit scripts. Historical N05 2025 is separately noncommercial-use restricted; do not redistribute raw N05 coordinates.
