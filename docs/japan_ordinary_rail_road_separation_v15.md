# v1.5 — Japan ordinary railway sites with independent road-only GIS controls

**Verified source run:** [GitHub Actions #42447780719](https://github.com/zuizui0223/tetsu/actions/runs/42447780719), 2026-11-08 JST. Actual MLIT N02 (2025) passenger-rail and N13 (September 2024) road-centreline geometries, each CC BY 4.0, were downloaded and joined in source-dependent local UTM. All GIS geometry remains transient; only aggregate counts and route labels are published.

**Stage: source sampling opportunity, NOT positive railway species/seed/mowing/wind/adaptation effects.**

## Prespecified independent source regions and eligibility

Study regions fixed before reading any biological outcomes: Sendai/Miyagi (N13 mesh 5740), Niigata (5639) and Nagano (5438), all within their corresponding national first-level map mesh. Conventional rail types 11=JR and 12=ordinary rail only. Exclude high-speed (operator code 1), municipal subway (operator code 3), funicular/monorail type13–17, tram type21–25, and explicit subway/high-speed route labels. **A first source run incorrectly let Sendai municipal subway Tozai Line in as an ordinary surface line; the operator-type3 bug was fixed and the counts below use the corrected code and a regression test.**

Choose railway source route × operator groups of at least 2.5km within the pre-fixed region, sample fixed side-normal near candidates at **35m** from centreline (about 1.2km along-rail step and 900m anchor spacing), and a geographically paired **1600m** far control on the same perpendicular side. A strict source opportunity requires:
- near rail point **≥75m** from any N13 2024 road centreline;
- near rail point **≥300m** from any N02 railway station line;
- far point **≥1km** from any N02 passenger rail;
- far road-control point **≤50m** from an N13 road centreline.

No photo, NDVI, flowering, species or plant presence is used to select these sites. All sites remain candidates until actual at-grade rail/right-of-way vegetation is independently verified.

## Actual GIS results

| Study region | N13 projected road line parts | N02 ordinary route groups >2.5km | Fixed geographic rail-side comparisons | Strict geometry-eligible rail/road pairs |
|---|---:|---:|---:|---:|
| Sendai/Miyagi | **183,706** | 4 | 81 | **2** |
| Niigata | **192,241** | 5 | 154 | **4** |
| Nagano study rectangle | **232,428** | 5 | 135 | **4** |
| **Total** | **608,375** | **14** | **370** | **10** |

| Region | N02 route label | Operator | Strict candidate pairs |
|---|---|---|---:|
| Sendai | 仙山線 | 東日本旅客鉄道 | **2** |
| Niigata | 磐越西線 | 東日本旅客鉄道 | **4** |
| Nagano rectangle | 吾妻線 | 東日本旅客鉄道 | **2** |
| Nagano rectangle | しなの鉄道線 | しなの鉄道 | **1** |
| Nagano rectangle | 篠ノ井線 | 東日本旅客鉄道 | **1** |

The five route groups across three regions are potential independent **rail corridors for future measurement**, not 10 random railway experiments. The Nagano region is a geometric rectangle rather than a prefecture-boundary assertion; 吾妻線 extends beyond Nagano into neighboring prefectures.

## What this does and does not show

The prior 2023 street-running Utsunomiya LRT pilot was highly road-co-located (60/68 samples just 20m from track centreline were <25m from roads), leaving 1–2 very strict vegetation comparison pairs. The **new positive geographical result** is that conventional railway alignments with a road-proximity contrast do exist, including more rural sections of 仙山線, 磐越西線 and 吾妻線.

However, N02 ordinary railway centreline data DO NOT establish whether a local section is in tunnel, on bridge, elevated, or at grade. Near-rail land may be forest, garden, roads, inaccessible railway property or ballast rather than grass. N13 from 2024 does not establish road surface width, road traffic or historical changes. N02 omits freight-only line and yard assets, and none of these geometries observes seed/spore transport, flowering, plant height, cutting/herbicide history, train winds, or genetic traits.

The next outcome-independent source QA stage recomputes the same strictly GIS-selected 10 candidates in memory and uses **2025 Sentinel-2 reprocessed Collection1 L2A** red, NIR and SCL. Select four 2025 summer acquisitions nearest August 1 by STAC calendar/date before reading actual vegetation values. Report source coverage and count of near-rail SCL4 vegetation; no biological effects or historical selection are allowed.

## Reproduction and attribution

- [Official MLIT N02 railway type fields](https://nlftp.mlit.go.jp/ksj/gml/datalist/KsjTmplt-N02-2025.html), [N02 railway classification code list](https://nlftp.mlit.go.jp/ksj/gmlold/codelist/RailroadSectionCd.html), [MLIT N13 2024 roads](https://nlftp.mlit.go.jp/ksj/gml/datalist/KsjTmplt-N13-2024.html).
- Source sampling: src/tetsu/ordinary_rail_corridors.py
- Main source audit: scripts/audit_japan_ordinary_rail_road_separation.py
- Regression QA: tests/test_japan_ordinary_rail_corridors.py
- GitHub Actions: .github/workflows/japan-ordinary-rail-road-isolation.yml

All source counts are derived from independently retrieved official public GIS; no original shapefile, image, observation ID or precise individual sampling coordinates were added to the repository. Field access to active rail corridors requires operator authorization and safe approved sampling.
