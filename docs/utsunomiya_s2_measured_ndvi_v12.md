# v1.2 — Utsunomiya LRT actual Sentinel-2 pixels: cloud-masked NDVI, not railway causal effect

**Actual pixel-level run:** [GitHub Actions #42441480708](https://github.com/zuizui0223/tetsu/actions/runs/42441480708). Artifact: `utsunomiya-actual-cog-ndvi-2017-2024-2025-two-offset-branches`. Source code [`scripts/measure_utsunomiya_s2_ndvi_pixels.py`](../scripts/measure_utsunomiya_s2_ndvi_pixels.py), strict safety/calibration math [`src/tetsu/sentinel_ndvi.py`](../src/tetsu/sentinel_ndvi.py). **These are real downloaded Sentinel-2 COG red/NIR/SCL pixel values**, not a catalogue search, synthetic test or claimed biological adaptation.

## Frozen measurement geometry and scene dates

- Independently verified 2023 Utsunomiya/Haga surface LRT source: MLIT historical N05-2025 (noncommercial-only GIS, not redistributed), official passenger service started **2023-08-26**, railway construction started **June 2018**, and testing took place before formal opening.
- Built **9 equally spaced anchors** along the verified 14.52km N05 rail geometry (not optimized for vegetation), ±140m **near offset** vs ±1,600m **far-offset control** on each possible corridor side, excluding far offsets that fall within 1km of another N05 line segment. Result **13 usable geographic near/far sides = 26 fixed physical pixel points**, all from **ONE** spatially correlated urban light-rail system.
- Matched **two neighboring Sentinel tiles** `54SVF` and `54SUF` within each year; same acquisition date per year, source STAC scene ids frozen before looking at pixel NDVI.
- Dates: **2017-07-07** (preconstruction), **2024-08-04**, **2025-07-30** (post service). The 25–28 day seasonal mismatch versus the preconstruction scene is an important phenological confound.
- Actual remote COG reads used red/NIR (10m) and scene classification SCL (20m), with only SCL **4=vegetation** and **5=nonvegetated** accepted; cloud/shadow/water/unclassified/unknown pixels excluded. A classified fixed sample may contain pavement, garden, street trees, roofs or construction. **No species or grassy railway strip is resolved at 10–20m**.
- User/rail geometry coordinates and raw 3-band pixels are not released, only small aggregated summaries with source traceability.

## What physically measured and how much was valid

| Year (real sensing date) | SCL/real NDVI provider usable pixel points / 26 fixed | Fixed within-year near/far point pairs both clear | Raw provider-calibration near−far median NDVI in those pairs | Reference caveat |
|---|---:|---:|---:|---|
| 2017-07-07 | **8** | **2** | **+0.063210** | Baseline/cloud coverage poor; valid points all SCL 4 vegetation |
| 2024-08-04 | **25** | **12** | **−0.105619** | BOA legacy metadata contradictory |
| 2025-07-30 | **24** | **11** | **−0.081330** | BOA legacy metadata contradictory |

**You cannot subtract the within-year medians directly**: the sets of clear matched points differ by date. With exactly the **same near and far points** SCL-clear on both baseline and post year, only **2 point pairs** qualify for 2017→2024 and only **2** for 2017→2025. Under the legacy item's **`earthsearch:boa_offset_applied` flag** calibration, descriptive *difference in near−far paired spectral contrast* is **−0.090955 NDVI** in 2024 and **−0.070353 NDVI** in 2025, respectively, but each is based on just two non-independent corridor sides. **No inferential or railway-specific vegetation-effect claim is justified.**

## A concrete source-calibration identification failure

The actual 2017 Earth Search STAC items have `s2:processing_baseline=00.01`, `earthsearch:boa_offset_applied=false` and band metadata `offset=0`: no source metadata conflict.

**2024 and 2025** have `baseline=05.11`, **`earthsearch:boa_offset_applied=true` AND each band says `raster:bands.offset=-0.1`**. These indicate opposite instructions about whether to subtract 1000 DN again from the current COG. This is a known legacy Earth Search source issue, not a tetsu result. [Element84 issue #71](https://github.com/Element84/earth-search/issues/71) and [Element84 discussion on Collection 1](https://github.com/Element84/earth-search/discussions/26). Element84 recommends **`sentinel-2-c1-l2a`** (reprocessed, more consistent radiometry) instead of the legacy `sentinel-2-l2a` where available.

The code computes both branches **without requiring the same pixels to be valid in both branches**. Requiring both would bias the source by discarding pixels with physically negative reflectance when subtracting 1000 DN twice:
- 2024 `provider_flag_branch`: **25** valid pixels and **12** near/far pairs, median near−far NDVI **−0.105619**.
- 2024 `raster_offset_branch`: **13** valid pixels and **5** near/far pairs, median **−0.320346**.
- 2025 `provider_flag_branch`: **24** valid pixels and **11** pairs, median **−0.081330**.
- 2025 `raster_offset_branch`: **11** valid pixels and **2** pairs, median **+0.143251**.

**The apparent contrast even changes sign across 2025 calibration interpretations.** Not the same retained pixel sample; cannot call the difference a measured biological effect. And 2017→2024/2025 **zero complete matched pairs** remain under the subtract1000 branch. A robust cross-era NDVI-change estimate requires a verified harmonized radiometric source and valid overlap.

## Independent scene classification (not plant species)

The same physical pixel sources also carry actual SCL land-class codes:

| Year | SCL near valid / class4 vegetation | SCL far valid / class4 vegetation |
|---|---:|---:|
| 2017 | 3 / 3 | 5 / 5 |
| 2024 | 13 / 4 | 12 / 4 |
| 2025 | 13 / 5 | 11 / 4 |

**These annual proportions are not a valid vegetation trend.** In 2017 18/26 pixels were excluded mostly for cloud, shadow or unclassified SCL, creating severe weather- and class-dependent selection. Not even the same valid points across years, and the SCL 4 class is an algorithmic vegetation classification that can vary by processing baseline.

## Next strict gate before any substantive vegetation change claim

1. Independently check the reprocessed **Collection 1 `sentinel-2-c1-l2a`** for the identical LRT place, dates and red/NIR/SCL bands, with correct scale/offset handling. 2017 source C1 historical reprocessing can be sparse; if so, use separately verified calibrated legacy 2017 with explicitly disclosed source.
2. Add one or more spatially and seasonally matched cloudy-scene alternatives and select source images **by independent local SCL clearance, never by NDVI direction**. Match the same spatial points and a comparable time-of-year before reporting a pre/post change.
3. Obtain real road-only and urban construction control pixels, stable 2017 land-cover/surface classification, and verify at-grade trackbed vs adjoining roads, cultivated greenspace and sealed urban surfaces. A 140m near-offset is **not a rail grassland plot**.
4. Replicate across independent railway corridors or plausible staggered rail interventions. Current 13 pairs in one LRT are pseudoreplicates for any railway-treatment inference.
5. Preserve scientific mechanisms: NDVI is *landscape greenness*; do not infer seed dispersal, rail-specific plant genetics, species flowering, post-cut recovery or train-slipstream adaptation from it.

**Actual accomplishment:** real cloud-filtered S2 pixel NDVI and source harmonization/selection bias successfully audited with reproducible CI and tests. **Unverified claim:** opening the LRT reduced railway vegetation, transported particular seed taxa, or selected wind/mowing-adapted plants.
