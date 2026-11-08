# v1.3 — Independent road and baseline land-class validity gates, one Utsunomiya LRT

**Scope:** *Source-quality and observed geographic-comparability audit*, NOT causal railway vegetation effect or replicated biological study. This protocol describes the data and eligibility rules fixed before seeing the present v1.3 source result; however the entire **v1.3 audit is exploratory/posthoc**, motivated by sign instability of the earlier v1.2 satellite NDVI contrasts. It is not an externally preregistered experiment.

## Background: what v1.2 identified and why it was not enough

[Actual Sentinel-2 pixel analysis](utsunomiya_s2_measured_ndvi_v12.md) independently validated 2024/2025 C1 reflectance calibration on 49 fixed ground pixels. A cloud-cleared 2017 July 7 vs 2024/25 comparison at 9 anchors changed its NDVI sign when an additional 2017 July 22 source date became available; a 48-anchor *posthoc sensitivity* returned a different direction. Source availability, acquisition timing, asphalt/trees and urban works create nonrandom comparison sets.

For one railway, a coherent descriptive NDVI contrast is not evidence of species dispersal, herbicide, mowing, train slipstream, or heritable plant adaptation. The candidate route is Utsunomiya/Haga LRT, regular operation since 2023-08-26, with **construction beginning in 2018**. The only true preconstruction remote-sensing source in this pilot is 2017; 2018–2022 is already exposed to construction changes.

## Public sources and original data restrictions

1. **Official MLIT 2024 N13 road centreline GIS**, two 1st-order map mesh archives **5439** and **5440**, needed because the corridor crosses 140°E. Data reference date **September 2024**, public CC BY 4.0. [N13 official metadata](https://nlftp.mlit.go.jp/ksj/gml/datalist/KsjTmplt-N13-2024.html), GSI land-survey acknowledgment. This is road **centreline geometry**, not a routing or actual traffic survey, and its coordinates are **seven years after** the preconstruction 2017 satellite scene. N13 **cannot** establish what road network existed before track construction; matching on post-road proximity can be post-treatment selection. All road distances must be labeled **2024 contemporaneous road proximity**.
2. **MLIT N05 historical passenger rail**, 2025 edition, **NONCOMMERCIAL-ONLY** due third-party historical source rights. Operator-verified Utsunomiya LRT surface line geometry ~14.52km, 2023 service start; no freight-only routes. Source is fetched temporarily and immediately discarded; **never redistribute restricted line geometry**.
3. **Sentinel-2 COG red/NIR and SCL classification**, 2017 legacy L2A (2017-07-22 or 2017-07-07, per-pair earliest valid same-day source) and reprocessed Collection1 2024-08-04 / 2025-07-30 (red/NIR offset corrected once). Source physical calibration cross-check [49/49 exact same pixel matches](https://github.com/zuizui0223/tetsu/actions/runs/42441974063). SCL=4 vegetation or SCL=5 bare/not-vegetated accepted; all other SCL classes considered unknown/unsupported for paired NDVI.

## Sampling and quality gates — frozen for this run

- **48 evenly spaced anchors** along the same surface LRT source line, ±140m near and ±1600m far on both sides, ignoring pairs where far is within 1km of any rail section. Every anchor pair is *correlated* with nearby anchors and is not an independently assigned rail treatment.
- **2024 road distance** from the nearest official N13 centreline. For both points, compare their distances. Primary road-availability QC: *both ≤100m of a 2024 road* and *absolute near/far road distance difference ≤50m*. Fixed sensitivity 100m difference. A station 140m from tram tracks is **not** verified railway right-of-way vegetation.
- **2017 same-acquisition SCL**: near/far must both be clear and of the same baseline classification (both vegetated 4 or both nonvegetated 5); **both vegetated 4** reported separately. SCL is an algorithmic land-class code, not a field survey. Both 2017 NDVI values must be valid on one single date, with Jul22 preferred over Jul7 independently of greenness.
- **Baseline spectral similarity:** absolute difference in the observed 2017 near vs far NDVI ≤0.15 (posthoc QA heuristic, not an approved regression-to-the-mean-free matched design).
- **Combined strict QC:** 2024 road comparability (≤50m distance difference, both road≤100m), same 2017 SCL class, and preconstruction 2017 NDVI similarity ≤0.15. If these yield too few comparable pairs, **report that failure and no railway/management effect**.
- Descriptive paired spectral contrast where available: `(NDVI near post – NDVI near 2017) – (NDVI far post – NDVI far 2017)`. Each pair uses exactly the same ground coordinates and a shared baseline acquisition date; different pairs may use different 2017 source dates. Summaries use **raw pair values** and also report occupied along-track blocks, but 48 anchors do not create 48 independent railway lines. **No confidence interval, p-value or causal interpretation**.

## Retention and scientific interpretation

Scripts:
- `src/tetsu/rail_road_landbalance.py` — safety bound on N13 downloads, streaming SHP/SHX to local temporary storage, geographically clipped segments, projected nearest-road geometry, quantitative quality gates, four along-corridor blocks.
- `scripts/audit_lrt_n13_road_land_balance.py` — fixed 2017 and 2024–25 remote-C1 time series integrated to N13 road distances.
- `tests/test_lrt_n13_road_land_balance.py` — synthetic source ZIP/CRS and distance fixture, guaranteed boundary classification, flat-NDVI storage and missingness, spatial-block safety.
- `.github/workflows/lrt-road-landcover-v13.yml` — run the actual source-backed complete audit and save **aggregate-only output**.

Only output **summaries, counts and original-source SHA256 hashes**, not raw N05 geometry, original N13 road lines, coordinates, individual source pixel DN or copyrighted images.

### Ecological effect is unidentifiable by these data alone

Even if a subset has matching 2017 SCL and similar 2024 road accessibility, the **2024 road network is not known preconstruction**, and the surrounding 10–20m pixels mix roads, rooftops, street trees, gardens, rails and bare construction. There is still only one physical tram corridor. No verified trackbed/slope/garden labels, train wind measurements, cutting/herbicide dates, repeated taxa/spore/flower ROI or propagule records. A high NDVI sample near tram lines does not establish railway grass adaptation. Modern original field plots and independently verified historical road-only controls are needed for railway-specific biological inference.


## Actual first run results: N13 roads + pre-2017 SCL + same ground pixel NDVI

**Successfully executed** [GitHub Actions #42444288968](https://github.com/zuizui0223/tetsu/actions/runs/42444288968) on 2026-11-08 JST, after fixing a scalar-vs-dict lookup error and adding regression coverage (7/7 targeted tests passed). The first two failed runs were an integration test import/lookup issue, not absence of road or satellite data. Full JSON result is retained as the `lrt-n13-2017-prelandcover-2024-roaddistance-48anchor-diagnostics` workflow artifact.

- N13 source mesh 5439: **1,267,863 national-mesh road feature lines**, with **95,426** projected line parts within the fixed study region bounding box. Download SHA256 `57b8849915784bb058b7e19ddc9b61c19936026f94587c9514bc302cd26eb72c`.
- N13 source mesh 5440: **613,014** road features, with **41,860** projected line parts indexed. Download SHA256 `7688c155d9f6ec355c1aac4e43f6dc8e28596bbbe46d386b6d687617561fe124`.
- **137,286** locally indexed projected road line parts total, measured at **68 fixed near/far transect pairs (136 physical points)**. 2024 official road-centreline distance median: **18.02m** at tram-near 140m offsets versus **16.95m** at far 1,600m offsets. This confirms that *modern road proximity is similar in aggregate*, not 2017 land use.

| Actual quality gate | Number of geometrically defined near/far pairs passing |
|---|---:|
| All fixed physical pairs | **68** |
| Both within 100m of 2024 official road | **62** |
| Both road≤100m, and absolute road-distance difference≤50m | **58** |
| 2017 baseline scene has both valid cloud-screened spectral values on the same date | **16** |
| 2017 baseline both in the **same** SCL class (both vegetation or both nonvegetation) | **12** |
| 2017 baseline both **vegetation** (SCL=4) | **8** |
| Road balance≤50m AND 2017 same SCL | **9** |
| Road balance≤50m, same SCL AND baseline |NDVI-near−NDVI-far|≤0.15 | **8** |

Only **two out of four** longitudinal road sections/blocks have eligible strict matched pixel pairs; *all observations remain from one LRT corridor*. A sample size of eight is not eight independent railway construction or wind treatments.

| Contrast among **same** fixed, cloud-valid spatial pairs | 2017→2024 median near−far NDVI change | n pairs | 2017→2025 median | n pairs |
|---|---:|---:|---:|---:|
| All baseline+post matched spectral pairs | **−0.073941** |16| **−0.116642** |15|
| Road proximity≤50m matched and 2017 same SCL land class | **−0.046250** |9| **−0.085644** |9|
| Road, 2017 same SCL and 2017 NDVI difference ≤0.15 | **−0.058956** |8| **−0.112438** |8|

The descriptive signs remain negative after these **posthoc comparability restrictions**, but there is **no causal railway effect identified**. In the strict eight pairs, the 2024 NDVI block medians were **−0.046250** (block with five pairs) and **−0.151576** (block with three); corresponding 2025 medians **−0.085644** and **−0.139232**. Two correlated along-track blocks in one city are insufficient to infer an LRT-wide rail operation response. Moreover the “same class” restriction also allows bare/bare 2017 pixels; a **separately labeled additional posthoc gate requiring both sides vegetated in 2017** has been implemented before paper-oriented biological interpretation. Never quietly treat bare/bare spectral differences as vegetation loss.

**Most important source-quality finding:** Road distance is not the principal immediate overlap obstacle (58/68 near/far pairs had similar 2024 road distance), whereas historical cloud/land-class/vegetation balance restricts the usable remote-sensing sample severely (16 valid baseline spectral pairs and 8 fully matched class/road/spectral pairs). Road geometry is **2024**, not an independent pre-2017 road control. The negative observation-time contrast, however repeated, could result from asphalt, housing, paving or LRT construction and **cannot be assigned to railway weed survival, flowering, propagule transport, wind, cutting or evolution**.
