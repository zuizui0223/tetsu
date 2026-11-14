# v1.6 — Two mapped grassland rail candidates: actual S2 2025 and GSI aerial-source cross-check

**Source-backed successful execution:** [GitHub Actions #42454635103](https://github.com/zuizui0223/tetsu/actions/runs/42454635103), 2026-11-08 JST; four targeted tests passed and no source errors. This is posthoc sampling feasibility and data-quality investigation, NOT a railway biological response or ground-verified railway-owned grass cover.

## Candidate selection occurred BEFORE reviewing 2025 spectral data

The exploratory, source-audited [370 normal-rail transects](japan_worldcover_370_site_posthoc_discovery_v16.md) included just **two** nominal 35m-from-track GIS sites where ESA WorldCover 2021 assigned **grassland class 30**, the N02 rail-to-N13 road distance was at least 50m (relaxed from the ORIGINAL 75m cutoff), station distance at least 300m, the 1600m road-only comparator within 50m of a 2024 road and at least 1000m from any N02 passenger rail. These routes are **信越線, Niigata** and **小海線, Nagano**, both JR East. Selection is posthoc to ESA 2021 source classes, and cannot become confirmatory railway effect estimation.

## What the independent source check actually returned

| Real satellite/imagery indicator | 信越線 Niigata | 小海線 Nagano |
|---|---|---|
| 2021 ESA WorldCover near35m | Grassland (30) | Grassland (30) |
| 2025 S2 C1 near SCL4 on 1 or more of 4 summer dates | **Yes** | **No** |
| 2025 near SCL4 and 1.6km road-only control SCL4 on same date | **Yes (1 pair)** | **No (0 pairs)** |
| 2025-07-30 near NDVI (under SCL4 source gate) | **0.378209** | **Not eligible** |
| 2025-07-30 road control NDVI (under SCL4 gate) | **0.875961** | **No paired comparison** |
| 2025 near minus road-only NDVI same date | **−0.497752 (n=1)** | **Missing/null, not zero** |
| GSI current seamless aerial-photo Z17 | Real nonblank JPEG | Real nonblank JPEG |
| GSI current seamless aerial-photo Z18 | Real nonblank JPEG | Real nonblank JPEG |
| GSI Z19 tile | 404 unavailable | 404 unavailable |
| Verified actual surface railway verge grass | **No** | **No** |
| Verified cutting, herbicide, train wind or flower/spore response | **No** | **No** |

Dates were chosen by independent STAC metadata calendar timing before extracting near/far NDVI: Niigata returned both SCL4 on **2025-07-30**, while its other three selected summer dates had insufficient quality source pixels. Nagano had usable pixels on **2025-07-30** and **2025-07-28**, but did not register near-rail SCL4 vegetation on any of the four outcome-blind selected dates. A 2021 grassland label followed by a 2025 nonvegetation SCL status is NOT proof of grass loss, railway mowing, drought adaptation, construction or conversion: observation year, sensor/algorithm, land classification and seasonal vegetation may differ. Missing imagery never means biological absence.

A 2025 one-pair near-minus-control NDVI contrast of −0.497752 is merely **a descriptive current spectral difference**, not a rail effect: one corridor segment, potential differences in vegetation type, canopy, soil, mowing, observer selection and no actual sample plot validation.

**High-resolution source availability:** Both site areas returned nonblank **official GSI seamlessphoto JPEG** at zoom17 and 18, but zoom19 was unavailable. This confirms that detailed visual grade/ground-vegetation reconnaissance has a source. The mosaic has varying acquisition vintages and no image-based human classification was performed here; the availability result is not a verified at-grade railway grassland label or a 2017 historical road map.

## Next defensible ground-truth gate

A qualified botanist or independent visual review should determine at-grade vs tunnel vs bridge, tree canopy vs herbaceous ground, whether the nominated near-rail location is railway owned, and whether the far roadside is actually a vegetated road verge rather than another land-cover type. Require date/metadata for each review image and retain unknown rather than force yes/no labels. No entry onto operational railway land without permission.

Only with real plot evidence and dated railway cutting/herbicide records can the project move from remote source opportunity to **within-taxon plant/regrowth/flower-spore ROI**. The distinct main railway/freight/ports **propagule transport** mechanism remains untested by WorldCover, SCL or near-track NDVI.

## Reproducibility and attribution

- Actual source code: [two-site S2/GSI audit](../scripts/audit_japan_two_posthoc_grass_sites_S2_GSI.py); [four quality tests](../tests/test_japan_two_posthoc_grass_sites.py); [GitHub workflow](../.github/workflows/japan-posthoc-2-grass-S2-GSI.yml).
- ESA WorldCover 2021 v200: Zanaga et al. (2022), DOI https://doi.org/10.5281/zenodo.7254221. © ESA WorldCover project 2021 / Contains modified Copernicus Sentinel data (2021) processed by ESA WorldCover consortium.
- GSI aerial source: https://maps.gsi.go.jp/development/ichiran.html, copyright/attribution 国土地理院; no original source images or tile XY positions redistributed.
- MLIT N02/N13 and public C1 imagery source IDs/road eligibility are retained only as public aggregated geographic source count and route names, not original source geometries or exact 10m sample coordinates.
