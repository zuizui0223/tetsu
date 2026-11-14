# Tokyo real-data railway vs road photographic-occurrence audit — v0.5

**Execution date:** 2026-11-07 JST. **Source-backed results, but NOT a confirmed ecological effect, selection, plasticity, or phenological shift.**

## Provenance / reproducibility
- Corrected source + weight result: [GitHub Actions run #42376881325](https://github.com/zuizui0223/tetsu/actions/runs/42376881325), workflow artifact `tokyo-rail-road-plant-photo-observation-contrast`.
- Earlier audit [#42376261914](https://github.com/zuizui0223/tetsu/actions/runs/42376261914) used an erroneous `reference_population_count/selected_page_length` weight, resulting in **0.1159**. This is **superseded**, not the selected final value.
- The corrected page cluster weight is `n_available_pages` (page inclusion probability 1/`n_available_pages`, after a hash-chosen page), including a short last page; see `src/tetsu/transport_overlap.py` and the regression test in `tests/test_transport_overlap.py`.
- Code: `scripts/audit_rail_vs_road_inat.py`. Fixed first stage [protocol](rail_road_control_pilot.md).
- Public sources: 2025 MLIT N02 passenger railway, 2024 MLIT N13 5339 mesh road centrelines, iNaturalist research-grade plant photos observed in 2023–2025.

## Verified source scale
- N13 mesh 5339: **1,943,251** official road line features (complete file), of which **711,724** line-parts were indexed near the predeclared Tokyo metropolitan bounding box with a 0.03° halo.
- N02: **2,095** rail line-parts indexed near Tokyo (official passenger rail, excludes freight-only).
- Full source search for 4 focal herb species: **491** photographed observations (270 *Trifolium repens*, 76 *Solidago altissima*, 133 *Plantago lanceolata*, 12 *Erigeron annuus*).
- Spatial QC passed: **411 focal photos** (246 unique observers), **3,315 reference other-plant photos** (882 unique observers). Note: the reference is ONE fixed-hash-selected iNaturalist observation-ID page per year-month, not a complete other-plants census.
- Of the 491 original focal observations, 76 had missing or >100m spatial uncertainty and 4 obscured/private; 411 retained.

## Road-matched photographic composition (primary declared threshold)

Among observation points **≤100m from a public road** and rail **≤250m vs ≥1,000m**, exclude the intermediate rail band:

| Target-group photo taxon | Rail-near observation count | Rail-far count | Corrected relative photographic detection contrast |
|---|---:|---:|---:|
| All 4 focus herbs | **124** | **134** | **0.1333** |
| *Trifolium repens* | 69 | 77 | 0.1291 |
| *Solidago altissima* | 23 | 10 | 0.3222 |
| *Plantago lanceolata* | 27 | 46 | 0.0851 |
| *Erigeron annuus* | 5 | 1 | 0.5279 (too sparse) |

Reference other-plant records under the same primary condition have **weighted counts 3,837 near and 552 far**, with one entire ID page selected per month. Those are ***weighted photographic opportunities***, NOT the abundance of reference plants.

`relative contrast = (focal near / reference weighted near) / (focal far / reference weighted far)`, with 0.5 continuity correction. This contrast measures **composition in photographs**, not individual survival, plant fitness or a railway treatment. The selected page clusters are correlated and the indicated ratio has no independent-record sampling CI.

**Pre-fixed sensitivities:** rail ≤100m and road ≤100m: 0.1164. Rail ≤250m and road ≤150m: 0.1331. Month-stratified Mantel–Haenszel descriptive contrast: 0.1586. These directions agree but DO NOT solve confounding.

## Critical independent spatial matching gate — **failed**

In 2-km projected grid cells, after road-distance matching:
- **44** cells contain at least one focal and at least one reference photo (irrespective of rail distance).
- **10** cells contain at least one rail-near and at least one rail-far photo among either group.
- **Only 2 cells** contain both focus and reference photo types in BOTH rail-near and rail-far classes (full within-block 2×2 support).

This is far too weak to estimate a station/railway-specific effect while adjusting local land cover and urban environment. There is no independent external geographic validation; do not report inferential railway conclusions or causal p-values.

## Biological interpretation

**Positive confirmed engineering/data result:** The official road, rail, and observer/photo sources can be reconciled at consistent distances and audited with source SHA-256 and strict photo geoprivacy checks. A *quality problem* was quantitatively identified and a sampling-weight bug fixed.

**NOT demonstrated:** Railway-mediated species movement, railway-near alien enrichment, negative rail effect on any taxon, stronger wind-adapted growth form, accelerated or delayed flowering, mowing-driven plasticity or common-garden evolutionary adaptation.

A large concentration of photographic observations close to tracks and roads is likely bound up with accessibility, observation behavior, city planning and biological habitat differences. No photographed coordinate at 250m is independently validated as railway verge vegetation; source records can be city parks, public gardens or roadsides.

## Predefined next-phase upgrade

1. Independently verified **actual rail verge vs road verge vs background** plots (not nearest-line distance).
2. Match same taxa, land cover, urban context, road access and surveys within corridor pairs; at least several independent geographic blocks **with all four comparison cells**, not merely many photos from the same park.
3. Separate *plant community assembly* from *photo-taking behavior* with target-group full-floristic visit lists where possible, and repeat results across an independent region.
4. For cutting/flower ROI, implement verified fixed-quadrat complete photo masks and documented management dates; now available in `scripts/extract_photo_roi.py` (`patch_pass`) and `scripts/summarize_mowing_photo_roi.py`. Current naturalist photographs provide no validated fixed plot before-after data.
5. To infer train-induced wind or local adaptation: independent gust measurements, management, species turnover and common garden/genetic comparisons.

**Do not reinterpret this exploratory source audit as evidence against the original train seed dispersal mechanism.** The original mechanism tests remain open until field-verified habitat and freight transport data are available.
