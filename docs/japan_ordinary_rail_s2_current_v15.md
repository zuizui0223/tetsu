# v1.5 — real 2025 Sentinel-2 vegetation verification of Japan rail/road separated candidates

**Corrected actual source run:** [GitHub Actions #42448382497](https://github.com/zuizui0223/tetsu/actions/runs/42448382497), 2026-11-08 JST, including all 4 Sentinel C1 QA tests plus 3 neighboring-UTM-coordinate tests. Actual 2025 public reprocessed Sentinel-2 Collection1 L2A reflectance and scene-classification pixels were downloaded via COG HTTP ranges. Only aggregate route/site-source counts are retained. The exact source data are not experimental botany and no causal railway mechanism is estimated.

## Independent geographic site eligibility BEFORE looking at any satellite outcome

The companion [MLIT N02 conventional rail versus N13 road isolation audit](japan_ordinary_rail_road_separation_v15.md) selected exactly **10** samples from **5 route groups** across three geographically separate areas (Sendai/仙山線 two, Niigata/磐越西線 four, Nagano rectangle/吾妻線 two + しなの鉄道 one + 篠ノ井線 one). Rail-near pixels are nominally 35m perpendicular from N02 passenger railway centreline, >=75m from 2024 N13 road, >=300m from passenger station, with a paired 1600m road-only comparator <=50m from road and >=1000m from any N02 passenger rail. This is a *source-only sampling frame*, not inspected at-grade railway-owned vegetation.

For each region the source code queried 2025 summer public Element84 Earth Search Sentinel-2 **Collection1** with red 10m, near-infrared 10m and SCL 20m bands, and read exactly four independently STAC-ranked source acquisition days near August 1. Selection was based on calendar timing and scene-cloud metadata, **NOT the NDVI sign or plant cover**. Only **SCL class 4** is a positive vegetation classification; SCL class5 is not-vegetated and clouds/unclassified pixels are unknown, never actual plant absence. For matched phenological context, near and control must have usable SCL on the SAME acquisition date.

## Verified current vegetation SOURCE opportunity

| Region / ordinary route group | Preselected road-isolated near/control GIS pairs | Near SCL4 on >=1 selected summer date | Near & road-only control both SCL4 on same date | Near SCL4 on >=2 distinct selected dates |
|---|---:|---:|---:|---:|
| Sendai / 仙山線 | **2** | **2** | **2** | **1** |
| Niigata / 磐越西線 | **4** | **4** | **4** | **0** |
| Nagano / 3 routes | **4** | **4** | **4** | **4** |
| **Total** | **10** | **10** | **10** | **5** |

2025 fixed source dates were:
- Sendai: August 3, July 30, August 4 and July 27; only 2, 4, 2 and 2 usable matched pixel locations out of 4 at each date. Some scenes had no good pixels.
- Niigata: August 2, July 30, August 4, July 28; only the July30 scene had all 8 (four near+four control) strict clear SCL4 pixel readings. Therefore **0 repeated multi-day observations in Niigata is a CLOUD/SOURCE-COVERAGE result**, not evidence of short-lived near-rail vegetation.
- Nagano: August 2, July 30, August 4, July 28; 4, 7, 6 and 8 usable strict pixels at the 8 fixed measurement locations, including all near candidates observed SCL4 on multiple dates.

### Important independently discovered CRS failure corrected

The first 2025 source audit [#42448091802](https://github.com/zuizui0223/tetsu/actions/runs/42448091802) returned no Nagano pixel values because the C1 Sentinel tile used **EPSG:32653 (UTM53)** whereas our N02 sampling coordinates are **EPSG:32654 (UTM54)**. That was a source software-coordinate-mismatch, NOT a zero vegetation observation. The actual 2025 COG pixel reader was corrected to transform ground points into the raster's native source UTM CRS (53/54/55) before sampling. The corrected [run #42448382497](https://github.com/zuizui0223/tetsu/actions/runs/42448382497) physically retrieved all four Nagano candidate near-pixel vegetation states. Synthetic tests ensure sample site location is unchanged under cross-UTM conversion and refuse Web Mercator or unverified non-Japan projection as a satellite native raster.

## Limits and research next steps

- **A positive SCL4 vegetation class is NOT a confirmed botanical rail verge.** At 35m from source rail centreline, the 20m SCL pixel can be tree canopy, fields, surrounding natural forest, roadside growth or planted garden. Need actual above-ground/right-of-way imagery and operator-approved field plot verification; N02 does not mark bridges/tunnels/elevated track at sampling points.
- A road-only far control point that was SCL4 on the same date is a geographically eligible *current vegetation* comparison, not a standardized matched plant community or historical land-use control. The measurements are 2025 only and cannot establish **railway-induced change** or train-mediated transport.
- Satellite source sampling and road-source candidate selection are outcome independent, but **conditioning on observed 2025 vegetation** would bias any later estimate of flowering probability or ecological colonization. The source vegetation gate is for deciding where a field plot is possible, not an inclusion policy for treatment-effect modeling.
- The five lines cross three separate geographic regions, but 10 points are **correlated within corridors**, and neither service-frequency, freight cargo, operator cutting/herbicide histories, true microclimate, wind exposure, propagule arrival, local adaptation nor plant taxa have been measured. No significance test or species richness claim.
- Safe field research must not involve entering operating railway land without permission. Do not publish precise high-risk access points as a general field guide.

Next: use high-resolution local aerial/ground-source images to exclude tunnels/viaducts and identify genuinely track-owned herbaceous vegetation, then match vegetated **road verge** locations, control for land cover/temperature/observer, and use verified maintenance history or standardized fixed-plot photo ROI. The original Japanese freight-rail/port-alien plant assembly hypothesis remains independently testable via network data; N02 passenger railway lines alone cannot prove train-carried seed dispersal.
