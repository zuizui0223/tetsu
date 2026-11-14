# v1.1 — New 2023 Utsunomiya LRT passenger service: botanical photo and satellite source checks

**Real-source code and complete audit runs:**
- [National N05 historical rail source](railway_history_n05_v11.md)
- [2023 LRT observed-source railway segment × iNaturalist photo opportunity](https://github.com/zuizui0223/tetsu/actions/runs/42438821083)
- [2017–2025 public Sentinel-2 source availability](https://github.com/zuizui0223/tetsu/actions/runs/42439134028)

### Independent event validity

An **entirely new Utsunomiya/Haga surface light-rail system**, operator-confirmed regular public service start **2023-08-26**, serving Utsunomiya and Haga town. The single matched official N05 2025 polyline for `宇都宮芳賀ライトレール線` / `宇都宮ライトレール`, year=2023, is approximately **14.52 km** long. These are **passenger operations**, not measured freight seed transport.

**Critical competing treatments:**
- **2017**: prior to railway *construction*; potentially usable satellite baseline, though resolution, observation calendar and regional land cover must be matched.
- **2018–2022**: construction started **June 2018**, test running November 2022; NOT an unexposed biological comparison period.
- **2023**: special construction/testing/opening transition year; exclude from simple pre/post.
- **2024–2025**: operation ongoing, vegetation management and surrounding road construction may also change.
- Road infrastructure is colocated with the LRT. Even near-rail vegetation may be street trees, planted gardens, or roadside grasses. A new passenger tram line cannot demonstrate freight seed transport.

Operator primary source: https://www.miyarail.co.jp/info/1143 . Construction background: https://u-movenext.net/information/164/ and JRTT/MLIT historical records.

## Observed iNaturalist plant photographs (2018–2025) — **coverage gate failed**

[Script](../scripts/audit_utsunomiya_lrt_photo_opportunity.py) uses the verified N05 LRT line in local UTM (not a generic city point), a fixed corridor envelope, iNaturalist plant research-grade photos, coordinate-quality ≤100m, no obscured locations, noncultivated status, and rail distances **≤250m** vs **1–3km controls**, with 250–1000m and >3km excluded. Observations were not selected for signs of flowering or post-opening enrichment; only aggregated counts and independent 1km grid/observer counts are exported.

| Year | Total photo records in source query | Strict geo-QC photos ≤250m | Strict geo-QC photos 1–3km control | Near 1km spatial grid cells | Independent near photographers |
|---|---:|---:|---:|---:|---:|
| 2018 construction start | 152 | **2** | 27 | 1 | 1 |
| 2019 | 16 | **0** | 3 | 0 | 0 |
| 2020 | 28 | **1** | 13 | 1 | 1 |
| 2021 | 15 | **0** | 0 | 0 | 0 |
| 2022 pre-operation/test | 26 | **0** | 10 | 0 | 0 |
| 2023 opening transition | 37 | **0** | 12 | 0 | 0 |
| 2024 operating | 31 | **1** | 7 | 1 | 1 |
| 2025 operating | 28 | **5** | 4 | 3 | 1 |
| **TOTAL** | **333** | **9** | **76** | *cannot add year-specific cells to get distinct across years* | *cannot add year-specific observers as independent people* |

Only **1 species** in **2025** had photographs in *both* distance bands; no year from 2018–2024 did. The 5 rail-near photos in 2025 came from just **one observer**, though spread across 3 spatial grid cells. Records alone therefore **fail the independent-site, balanced observer, and same-taxon comparison gates**.

It would be especially misleading to interpret the raw rise **3 near photographs (2018–2022, five years) → 6 (2024–2025, two years)** as a five-fold biological effect. Both timeframe length and camera/observer selection change. The 2018 "before" is also AFTER construction started. **No plant arrival, richness, flower timing, mowing, adaptation or railway effect is identified from these citizen photos.**

## Sentinel-2 source availability — stronger quantitative coverage opportunity, no vegetation measurement yet

[First public STAC seasonal audit](https://github.com/zuizui0223/tetsu/actions/runs/42439134028) independently queries Element84 Earth Search STAC (public Sentinel-2 L2A COG catalogue) for the same derived LRT study envelope at **June–September** of each year 2017–2025, source-scene cloud cover <20%. Both red, NIR and pixel-level scene-classification/cloud-mask **assets are present** for the returned items.

| Stage | Year | STAC scenes returned on first source page | Distinct acquisition dates on first source page |
|---|---:|---:|---:|
| Preconstruction | **2017** | **7** | **4** |
| Construction start | 2018 | 10 | 4 |
| Construction | 2019 | 14 | 5 |
| Construction | 2020 | 11 | 3 |
| Construction | 2021 | 21 | 8 |
| Construction/test | 2022 | 5 | 3 |
| Opening transition | 2023 | 10 | 7 |
| Operating | **2024** | **8** | **6** |
| Operating | **2025** | **11** | **7** |

These are **returned-page lower bounds, not exhaustive satellite scene censuses**: the provider indicates additional pages for all windows. Cloud fraction is **for the full satellite scene**, not proof a given 10m railway pixel is cloud free. To actually measure NDVI, separately retrieve equal-season cloud-masked red/NIR pixels in a spatially fixed grid, control for season and acquisition geometry and compare **independent road-only and city construction controls**. No pixel-valued NDVI, actual plant cover or ecological trend was computed in this audit. Code: [`scripts/audit_utsunomiya_sentinel_stac.py`](../scripts/audit_utsunomiya_sentinel_stac.py).

### Why this is worth doing

Fixed-location Sentinel-2 pixels are observable repeatedly even when no botanist photographs the rail corridor. A carefully masked long-term vegetation/green-cover curve could reveal a **surface-cover response to track construction and urban redevelopment**, whereas iNaturalist photographs near rail are too sparse for annual species/community tests. BUT **a 10m pixel is far broader than a narrow rail-side grass strip and often contains roads, pavement, planted trees, houses and bare construction**. It measures greenness only, not plant species or seeds, flowers, spore cones, train wind tolerance or adaptation. The effect of opening passenger operations **cannot** be isolated from construction and landscape reallocation without credible design and independent comparator corridors.

## Research decision

**Supported positive source-result:** a dated verified new surface rail event, complete N05 route-length/provenance verification, full annual Japan photo-source feasibility, and independently available preconstruction/construction/operating satellite scene metadata.

**Failed ecological inference:** direct species/photo before-after or moving-train propagule dispersal: 9 near observations across 8 years, 1 active photographer in each near-present year, and zero shared species/year in most seasons. The next appropriate quantitative measurement would be public cloud-masked multi-temporal 10m NDVI/change as a *land-cover* outcome with explicit road/construction controls, not speculative genetic wind adaptation or flowering advancement.

In parallel, keep the Japanese rail trackbed vs slope horsetail **historical secondary data** and actual operator vegetation-management (cutting vs herbicide vs growth retardant) hypotheses separate and correctly attributed.
