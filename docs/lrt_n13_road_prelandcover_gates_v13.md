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
