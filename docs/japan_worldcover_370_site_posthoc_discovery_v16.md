# v1.6 — a real sampling confound: road-isolated railway pixels become overwhelmingly forest

**Real source run:** [GitHub Actions #42454338343](https://github.com/zuizui0223/tetsu/actions/runs/42454338343) executed 2026-11-08 JST. This is a **posthoc exploratory site-discovery and selection-bias audit**, NOT a biological railway treatment effect or a preregistered ecological test. Original 370 route-normal N02/N13 samples were selected by geography and road distances **before inspecting** ESA WorldCover vegetation classes. After seeing that all 10 previously road-isolated candidates had trees or water, we deliberately expanded the source inspection to all 370 fixed railway vs road normal offsets.

## Data and methods

- Fixed three source regions/ordinary passenger lines: Sendai/Miyagi (81 fixed rail-normal near/far sides), Niigata (154), Nagano (135), all sampled from 2025 official N02 ordinary passenger rail and 2024 N13 road centrelines. Near sample 35m normal to nominal rail source line, far candidate 1600m on the same side.
- 2021 ESA WorldCover v200 10m source tile **N36E138** from [official public COG](https://esa-worldcover.org/en/data-access), DOI [10.5281/zenodo.7254221](https://doi.org/10.5281/zenodo.7254221). Each site receives a 2021 map-center class and a 3×3 context read in raster EPSG4326, with projected N02/N13 sites transformed from UTM54. No restricted route coordinates or source imagery published.
- Historical-vintage limitation: 2021 WorldCover and 2024 current road geometry differ in time; no actual grade/tunnel/bridge, herbaceous flora, maintenance or train exposure observed.
- Source eligibility for a *paired modern road-only control* requires far<=50m from 2024 road, far>=1000m from N02 passenger railway, near>=300m from railway station, and a **near-road separation of 25m, 50m, 75m or 100m**. The original v1.5 fixed minimum is **75m**. All weaker thresholds (25m/50m) were explored only **after the original 10 sites failed mapped-ground-grass eligibility**; results must stay labeled posthoc and cannot be reclassified into original independent tests.

## Actual 2021 land-cover source for all 370 sites

| Mapped 2021 near35m class | Sendai n=81 | Niigata n=154 | Nagano n=135 | Pooled n=370 |
|---|---:|---:|---:|---:|
| Tree cover (10) | 30 | 60 | 57 | **147** |
| Grassland (30) | 3 | 5 | 16 | **24** |
| Cropland (40) | 6 | 47 | 20 | **73** |
| Built-up (50) | 42 | 34 | 41 | **117** |
| Permanent water (80) | 0 | 8 | 0 | **8** |
| Bare/sparse vegetation (60) | 0 | 0 | 1 | **1** |
| **Total original source sides** | **81** | **154** | **135** | **370** |

A grassland class at 35m from N02 is only a coarse remotely sensed 10m land-cover class, not a verified herbaceous railway embankment or authorized sampling plot. Bare/sparse 60 is not an observed flowering grass patch. Tree class10 can conceal ground-layer herbs.

## The actual road-separation/land-cover selection trap

| Eligibility cohort | n fixed sites | Near mapped 2021 tree cover | Near mapped grassland/herbaceous wetland | Tree fraction |
|---|---:|---:|---:|---:|
| Original 370 fixed GIS normals, no road-only control restriction | **370** | **147** | **24** | **39.7%** |
| Road-only control exists; rail-near >=25m from modern road | **51** | **26** | **2** | **51.0%** |
| Road-only control exists; rail-near >=50m from modern road | **19** | **14** | **2** | **73.7%** |
| Original strict control; rail-near >=75m from modern road | **10** | **9** | **0** | **90.0%** |

**Major supported descriptive result:** The stronger the source design's current-road isolation criterion, the more the eligible 35m railway-near samples lie in **2021 tree-covered landscapes**. This creates an ecological sampling confound: a comparison selected to avoid road adjacency can accidentally compare remote/forested railway surroundings with road controls in different landscape types, rather than actual disturbed railway verges. Do NOT call it selection caused by railway vegetation management; site land cover and road geometry are not randomized, and the original 370 normals/overlapping routes are not independent biological replicates.

The 10 original strict sites were WorldCover **9 tree + 1 permanent water, 0 grass**; they had all been classed vegetation SCL4 by 2025 Sentinel2 when cloud clear, a vivid example that optical SCL vegetation does NOT mean herbaceous rail-verge presence. The mapped water pixel could be bridge/river landscape or classification/year discrepancy, not evidence of new water disappearance.

## New posthoc grassland source opportunities at a weaker modern-road threshold

Among the original 370 fixed geometry samples, requiring near≥50m from current road rather than 75m yields **two 2021 grassland-class30** near-rail pixels **with otherwise valid source road-only controls** (far road≤50m, far rail≥1km, near station≥300m):

| Region | Official N02 route | Operator | 2021 mapped grassland near35m candidates | Threshold tradeoff |
|---|---|---|---:|---|
| Niigata | **信越線** | JR East | **1** | near 2024 road >=50m, but original >=75m criterion not met |
| Nagano | **小海線** | JR East | **1** | near 2024 road >=50m, but original >=75m criterion not met |

The 25m relaxation adds no further mapped grassland candidates with a valid road-only/rail-isolated/away-from-station control in these fixed regions, and the 75m original threshold retains zero. These **two** are real independently source-identified **2021 mapped grassland site reconnaissance opportunities**, not proof of actual trackside vegetation or new railway-linked species. They were discovered **after looking at** the 2021 WorldCover map. **Never** interpret their satellite NDVI or presence as an independent replication of a railway environmental filtering hypothesis.

## Hypothesis sharpened, not established

The most biologically useful question remains whether a railway **through** forest may create a narrow edge/linear open habitat for ruderal herbs, but 35m away can lie beneath forest canopy. That requires a cross-track canopy-gap/at-grade image analysis and independent nonrail forest edges; this analysis did NOT measure canopy gap width, train turbulence, track ballast, mowing, species or seed deposition.

**Next source check:** re-observe the exactly two posthoc 2021 grassland mapped candidates in 2025 Sentinel2 SCL and independently see whether official high-resolution GSI current aerial photography *is available*. The aerial mosaic is not a date-fixed 2021 map and availability does NOT prove a grass verge. Only approved detailed image review/field surveys can certify grassland on railway property. Keep all original 370 outcomes and rejection/selection rules auditable.

## Reproduction and licensing

- Source [`src/tetsu/worldcover_landcover.py`](../src/tetsu/worldcover_landcover.py), [`scripts/audit_japan_ordinary_rail_worldcover_expanded.py`](../scripts/audit_japan_ordinary_rail_worldcover_expanded.py), [`tests/test_ordinary_rail_worldcover_expanded.py`](../tests/test_ordinary_rail_worldcover_expanded.py).
- Original strict source selection: [`scripts/audit_japan_ordinary_rail_road_separation.py`](../scripts/audit_japan_ordinary_rail_road_separation.py), amended only to expose an outcome-blind **transient** complete fixed-geometry cohort, without changing its published original strict 10 GIS source eligibility.
- **All three regional N02/N13 WorldCover source joins succeeded, 370/370 source points and 10/10 strict original eligible verified, zero source errors, three targeted source tests passed.** [Real run #42454338343](https://github.com/zuizui0223/tetsu/actions/runs/42454338343).
- **© ESA WorldCover project 2021 / Contains modified Copernicus Sentinel data (2021) processed by ESA WorldCover consortium**; N02/N13 official Japanese government data acknowledged. Categorical land-cover labels have different 2021 vintage from roads (2024) and 2025 Sentinel2 C1 SCL.
- No raw copyrighted WorldCover COG, N02/rail geometry, N13 road geometry, precise sample coordinates, station positions or third-party aerial photographs redistributed.
