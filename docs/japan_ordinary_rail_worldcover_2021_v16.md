# v1.6 — ESA WorldCover 2021 reveals tree-covered 2025-SCL4 railway candidates

**Source-backed result:** [Real COG run #42453959307](https://github.com/zuizui0223/tetsu/actions/runs/42453959307), executed 2026-11-08. This is **an independent land-cover source/product cross-check**, not a new biological railway effect. The code physically sampled the original 10 *fixed-by-geometry* Japanese railway near/road-control sample coordinates. It did not select or rank sites by WorldCover class, 2025 NDVI or the presence of grass.

## Authoritative land-cover source and correct chronology

- **ESA WorldCover 2021 v200**, 10m global landcover product, 11 categorical classes, based on Copernicus Sentinel-1 and Sentinel-2 remote sensing in 2021, DOI [10.5281/zenodo.7254221](https://doi.org/10.5281/zenodo.7254221). Product source [ESA WorldCover data access](https://esa-worldcover.org/en/data-access). Class **10 tree cover**, **30 grassland**, **40 cropland**, **50 built-up**, **60 bare/sparse**, **80 permanent water**, **90 herbaceous wetland**. The map is a machine-classified *land-cover interpretation*, not a species-level grassland survey.
- Source COG: `https://esa-worldcover.s3.eu-central-1.amazonaws.com/v200/2021/map/ESA_WorldCover_10m_2021_v200_N36E138_Map.tif`. The 3° tile covers **36–39°N, 138–141°E**, and therefore spans all three pre-fixed Sendai/Niigata/Nagano rail-road sample frames. Verified true **36000×36000, 1-band, EPSG:4326, 8-bit values, 0 nodata** source schema.
- Source construction: MLIT 2025 N02 passenger rail + MLIT 2024 N13 roads. Sample coordinates projected to UTM54 and back to EPSG4326 to read the ESA 2021 COG. **N02/N13 source selection is outcome-blind**; exactly the same road-isolated 10 sites from [v1.5](japan_ordinary_rail_road_separation_v15.md), not newly chosen forest locations.
- At each site, independently retrieve the **2021 WorldCover pixel at 35m from mapped rail**, the **1.6km road-control pixel** and the **3×3-pixel neighborhood** surrounding each near-rail point. Zero = unknown, never grass/bare/not-vegetated by assumption. No original GIS coordinates, WorldCover image data or source geometry are published.

## Actual independently retrieved COG classes

| Geographic cohort (fixed by GIS) | n rail / road matched candidates | 2021 WorldCover at near35m | 2021 WorldCover at far road-only controls |
|---|---:|---|---|
| Sendai / 仙山線 | **2** | **tree cover 2** | **tree cover 2** |
| Niigata / 磐越西線 | **4** | **tree cover 3, permanent water 1** | **tree cover 4** |
| Nagano / 吾妻線, しなの鉄道線, 篠ノ井線 | **4** | **tree cover 4** | **tree cover 3, cropland 1** |
| **TOTAL** | **10** | **9 tree, 1 mapped permanent water, 0 grassland or herbaceous wetland** | **9 tree, 1 cropland** |

At all nine near-rail 2021 tree-cover sites, **the 3×3 WorldCover window's plurality classification was tree cover**. For Niigata's remaining one source point, the 3×3 window's plurality was mapped permanent water. These are source **10m remotely classified dominant land-cover labels**, not independently measured actual ground flora or track grade.

**Critical source-data difference:** The same 10 source-selected near35m samples had *2025 Sentinel-2 SCL class4 = vegetation* on at least one summer source acquisition in [v1.5's real C1 audit #42448382497](https://github.com/zuizui0223/tetsu/actions/runs/42448382497), with 10/10 near and 10/10 road controls concurrently SCL4. This does **not** demonstrate 2021 tree-to-2025 grassland conversion. The sources have **different observation years and classification algorithms**; the 2021 water class may represent river/floodplain, bridge/tunnel effect, subpixel mixing, actual change or map classification error. No confident water→vegetation ecological transition is established. Tree canopy itself is naturally classed as 2025 SCL4 vegetation. Most importantly, **2025 SCL4 'vegetation' was never evidence of herbaceous railway verges**.

## New ecological interpretation and future quality gates

**Supported source result:** The ten *strictly road-isolated candidate 35m pixels* were mapped as trees or water in the independent 2021 categorical land-cover product, with 0 class30 grassland or class90 herbaceous wetland. This substantially **weakens their suitability as a direct herbaceous mowing/plasticity field-plot shortlist**, but it does NOT establish that grass, moss, spores, *Equisetum* rhizomes or any plant species are absent below trees or in local 2025 field verges.

The forest canopy finding raises a *distinct, falsifiable* landscape hypothesis: **rail corridors through forest may create narrow linear canopy gaps or disturbed edges even where 35m offsets remain mapped tree cover**. Demonstration would require cross-rail fine-scale canopy profiles, verified at-grade tracks and matched nonrail forest edges, not a high NDVI at a forest pixel. ESA 10m cannot resolve small grass patches or classify whether an adjacent area is railway-owned. 

For direct herbaceous ecology, reserve a new **exploratory (posthoc) candidate-discovery workflow** over all original **370 fixed geographical rail-road normals**, transparently retain the strict 10 and inspect class30/90 at other road-distance thresholds, before any attempt to fit rail effects. If there are no real mapped grassland sites, expand the geographic sampling frame to open/agricultural mainland or grassland rail contexts with an independently date-verified site-level aerial/field review. Any outcome-informed site selection must never become a hidden retrospective causal analysis.

## Attribution, script and repeatability

- Implementation: [`src/tetsu/worldcover_landcover.py`](../src/tetsu/worldcover_landcover.py), [`scripts/audit_japan_ordinary_rail_worldcover_2021.py`](../scripts/audit_japan_ordinary_rail_worldcover_2021.py).
- Source-sampling tests: [`tests/test_ordinary_rail_worldcover_2021.py`](../tests/test_ordinary_rail_worldcover_2021.py) (six checks for 3° tile codes, geographic coordinate conversion, categorical nosource, and true rasterio 3×3 source pixels).
- [Source run #42453959307](https://github.com/zuizui0223/tetsu/actions/runs/42453959307): **6/6 targeted tests successful, 3/3 source regions verified, 10/10 physical WorldCover pixels valid, source errors=0**.
- Original ESA product map credit: **© ESA WorldCover project 2021 / Contains modified Copernicus Sentinel data (2021) processed by ESA WorldCover consortium**. ESA 2021 and the 2025 Sentinel2 SCL share underlying satellite sensor classes (not independent ground-truth experiments), but the algorithms and reference years differ.
- The original railway **seed dispersal**, mowing and growth recovery, wind response, flowering/photo ROI, establishment and regional plant similarity hypotheses remain **untested** by this land-cover check. Their strongest confirmatory test requires actual site-level herbaceous flora and independently verified railway use/management exposure.
