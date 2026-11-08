# v1.4 — Are any 10m Sentinel pixels actually suitable for railway-verge vegetation?

**New empirical rail-centreline proximity audit**: [successful GitHub Actions #42445463029](https://github.com/zuizui0223/tetsu/actions/runs/42445463029). Source: true 2023 Utsunomiya/Haga LRT MLIT N05 route geometry (fetched transiently; NONCOMMERCIAL source NOT redistributed), official 2024 MLIT N13 road centreline GIS, original 2017 Sentinel-2 SCL/red/NIR, and correctly calibrated 2024 and 2025 Sentinel-2 Collection1 L2A.

**Scientific status: direct source-measurement/quality feasibility, not a causal biological result.** This version corrects a serious defect in interpretation of earlier v1.2–v1.3: **140 metres from the official railway centreline is not trackside vegetation**. This is an urban surface tram mostly adjacent to paved roads. The new source samples the **actual 20m and 35m perpendicular offsets** from the same N05 rail centreline, with a fixed 1.6km comparator, and verifies road proximity and whether both physical pixels were *already vegetation in 2017*. Even at 20m from a rail centreline the 10m optical pixel and 20m SCL pixel may be carriageway, paved sidewalk, yard, street tree, roofs, or garden. No image has been independently verified as track-owned grass or gravel habitat.

## Frozen geographic design and data

- [`scripts/audit_lrt_true_centerline_20_35m.py`](../scripts/audit_lrt_true_centerline_20_35m.py): 48 evenly spaced along-rail anchors; both route sides when offset control is ≥1km from any route segment. 68 valid geographic side-pair definitions per offset.
- **20m and 35m** near-offset positions computed exactly by interpolating the existing paired **140m near/1600m far** points' common perpendicular vector; all sampling points frozen before inspecting any satellite greenness. Test: [`tests/test_lrt_rail_centerline_pixel_eligibility.py`](../tests/test_lrt_rail_centerline_pixel_eligibility.py).
- 2024 N13 roads: official line GIS meshes **5439 + 5440** and 137,286 local projected road parts, with 2024 distance to nearest **road centreline**. **Road within 25m of a pixel is not an image-derived road mask**, and road ≥25m is NOT proof that the entire 10m pixel is not pavement or roads.
- 2017-07-22 (first when both pair pixels SCL clear) / 2017-07-07 (fallback) exact same-date near and far SCL 4 vegetation / SCL 5 nonvegetated plus red/NIR NDVI; 2024-08-04 and 2025-07-30 reprocessed C1 L2A with -0.1 reflectance offset applied **exactly once**. Original legacy offset discrepancy was resolved on 49 exact same-date/source pixels in [2024/25 Collection1 cross-check #42441974063](https://github.com/zuizui0223/tetsu/actions/runs/42441974063).
- All points belong to **one correlated urban LRT** which began regular passenger service **2023-08-26**, but construction began in **2018**. Distances and satellite pixel classes do not constitute independently replicated management treatments.

## Actual suitability audit — key finding

| Eligibility diagnostic from real sources | 20m from actual rail centreline | 35m from actual rail centreline |
|---|---:|---:|
| Geographic rail-near versus 1.6km comparator pairs | **68** | **68** |
| Near point also within **25m of a 2024 road centreline** | **60 (88.2%)** | **45 (66.2%)** |
| Near-point median 2024 road-centreline distance | **13.20m** | **18.84m** |
| One 2017 date clear for BOTH pixels at fixed offsets | **17** | **18** |
| Both pixels **2017 SCL=4 vegetation** | **8** | **9** |
| Both 2017 vegetation & near 2024 road distance ≥25m | **3** | **4** |
| Stricter: both 2017 vegetation, near 2024 road ≥25m, 2024 road distance similar between near and comparator, 2017 near–far NDVI difference ≤0.15 | **1** | **2** |

This is a **positive source-based indication of the road/rail proximity confounding** in this specific recently constructed urban surface tram: even 20m from the line, 60/68 examined pixel centres are within 25m of a **current** road centreline. This is not the share of the tram line physically covered by road, nor the rate of grass presence. The 2017 SCL/road-independent eligibility leaves only 1–2 near/far pixel pairs, nowhere near independent biological replication.

## Measured NDVI and strict gating — no train effect

Define purely descriptive spectral contrast for a paired near/far sample, with the SAME historical near+far source acquisition and same physical pixel in 2024 or 2025:

`Δ_rel = (NDVI_near,post − NDVI_near,2017) − (NDVI_far,post − NDVI_far,2017)`.

| Year | Offset | Historical pair complete | Contrast median all baseline/post clear | BOTH 2017 vegetation only | Strict vegetation + current road context (n) |
|---|---|---:|---:|---:|---:|
| 2024 | 20m | 17 | **−0.148366** | **−0.294190** (8 pairs, all one block) | **−0.208449 (n=1)** |
| 2025 | 20m | 16 | **−0.279182** | **−0.328851** (8, one block) | **−0.255576 (n=1)** |
| 2024 | 35m | 18 | **−0.147604** | **−0.099457** (9, concentrated 8+1 in two blocks) | **+0.049537 (n=2)** |
| 2025 | 35m | 17 | **−0.194892** | **−0.155182** (9, concentrated 8+1) | **+0.021190 (n=2)** |

**Direction reverses at the stricter road- and vegetation-balanced gate between the 20m and 35m sampling offsets**; source units in this strict gate are just 1 vs 2, and remain all in **one** spatial block. Do NOT select the more attractive sign. These are mixed urban pixel greenness descriptions, **not railway-induced plant growth, suppressive track maintenance, seeds, flower/spore stage or local adaptation**. Current N13 roads could have been built during or after 2018 LRT construction; the only plausible historical-preconstruction road/habitat check requires independently dated 2017 ortho imagery or OSM historical records plus documented uncertainty.

## GSI historical high-resolution imagery source added, not yet a road classifier

- [`scripts/audit_lrt_gsi_2017_aerial_tile_coverage.py`](../scripts/audit_lrt_gsi_2017_aerial_tile_coverage.py): official GSI [year-specific orthophoto tiles 2007+](https://maps.gsi.go.jp/development/ichiran.html), `nendophoto{year}/{z}/{x}/{y}.png`. Z=17 sampling at equal distant LRT source anchors, 20m/35m and control points; tests for real PNG content versus HTTP404/transparent or solid blank. Samples **2016, 2017, 2018 and 2024** to check whether actual historical land-cover interpretation could be performed.
- **A tile dated 2017 should be reviewed visually for actual road carriageway and track/vegetation; the tile availability checker does not infer whether an unphotographed road existed and does not automatically classify it.** Even a nominal 2017 year layer may mosaic aerial photos with within-year/coverage limitations; inspect metadata/provenance before claiming a precise intervention timeline.
- The output reports **aggregate coverage only**; no GSI photos or restricted N05 route coordinates released. GSI source credits and permitted use conditions must be respected.

## Falsifiable ecological design still open

The main original questions remain (1) seed/spore/plant-material **transport** by railway versus ports and roads, and (2) **establishment / microhabitat filtering** under documented railway mowing, herbicide, wind stress and bare ballast. A useful biological result now requires *actually verified right-of-way vegetation sample polygons* and road-only comparators, preferably **multiple distinct non-road-adjacent rail corridors**. The 20m proximity metric is an intentionally strict **feasibility check** that rejects treating "rail near" as a biological treatment by mere mapped distance. The documented 2011 JR East horsetail trackbed prominence remains a separate observational historical result, not a modern controlled experiment.
