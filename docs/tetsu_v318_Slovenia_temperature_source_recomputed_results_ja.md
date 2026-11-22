# v318 source-recomputed Slovenian BS–BT thermal indicator result (2026-11-22)

**Status: ORIGINAL-SOURCE WORKFLOW SUCCESS**, GitHub Actions run [45029366861](https://github.com/zuizui0223/tetsu/actions/runs/45029366861), exact analysis commit a7e866ee00e4718bfada0727648d83b5994f132d. Five synthetic tests passed. Original supplementary XLSX: Šabić et al. (2026) doi:10.1007/s44473-026-00155-z, SHA256 `26be7ca9513dd54136d4ec499d8413c1cb288045f87fe3cf8814731f2edac99f`. Source-verified v102 archive was reused; no original plot-by-species rows republished. Derived JSON artifact ID [14756709292](https://github.com/zuizui0223/tetsu/actions/runs/45029366861/artifacts/14756709292), SHA256 `9863ed2b3b26ad7dfe3b1fcf8c0a74162195735a95b449d7aaccd4aff586c010`.

## Actual result
Same railway survey location + same side pairs: 173; distinct field locations: 89; line sections: 15. Source temperature trait field: `EIV_T` (312 taxa: exactly 142 numeric, 135 `x`, 35 empty). Mean **observed incidence-based thermal trait coverage** among paired track sides: shoulder BS 0.480563, toe BT 0.503328. A zone pair is eligible only when both communities have at least one numeric trait. Each station/location's left/right differences are averaged **before** site-level inference.

| Minimum observed incidence EIV_T coverage in BOTH habitats | Eligible side pairs | Locations | Sections | Mean BS minus BT EIV_T | Shoulder warmer sites | Exploratory section-block bootstrap percentile 95% |
|---|---:|---:|---:|---:|---:|---|
| 0 | 171 | 89 | 15 | **+0.129208** | 59/89 | [+0.064880, +0.201071] |
| 0.5 | 60 | 47 | 11 | **+0.140544** | 28/47 | [+0.049958, +0.245322] |
| 0.7 | 3 | 3 | 3 | +0.017778 | 1/3 | [−0.206667, +0.380000] |
| 0.9 | 1 | 1 | 1 | −0.120000 | 0/1 | not estimable |

Bootstraps: 2500 fixed-seed draws resampling **entire original line sections**, averaging retained location-level contrasts. 95% percentiles are exploratory, conditional on source and eligibility; not country/study replication CIs. No standalone field p-value asserted. Full sample shoulder warmer in 59 locations, toe warmer in 30.

## Scientific inference
1. **Observed within-rail habitat filtering**: shoulder plants' recorded temperature indicator is higher on average than toe plants' under unrestricted and >=50% coverage comparisons. This is **consistent with** the original paper's inner-vs-outer thermophily description, not a newly independent country replication.
2. **Key obstacle**: requiring >=70% indicator coverage leaves just 3/89 locations. A high-coverage-only effect cannot be estimated reliably here. The +0.129 difference is an ordinal indicator-scale contrast, not +0.129 degrees Celsius or evidence of local microclimate magnitude.
3. **The climate-offset hypothesis (A) is NOT tested**: both BS and BT belong to the railway, there are no matched nonrail vegetation plots or exogenous regional thermal gradient in this source test. Study design and sampling widths are not proven equal. Neither temperature causality, train-carriage dispersal, persistence beyond climatic limits nor plant adaptation can be concluded.
4. New local/nonrail alternatives identified: Winkler et al. (2026) four Czech localities include real railway versus field margin/arable plots but publisher explicitly withholds original plot-level data; Dylewski et al. (2024/2025) 25 abandoned-railway plus 25 matched grassland reference plots are **on author request**. Rakotoson et al. (2025) RATP track/embankment original plot data also on request and both habitats are in railway environment. Details in [source-access audit](tetsu_v318_local_control_data_access_evidence.md).

## Reproducibility
`scripts/analyze_Slovenia_thermal_EIV_paired_v318.py`, `tests/test_Slovenia_thermal_EIV_paired_v318.py` and `.github/workflows/tetsu-v318-Slovenia-source-thermal-EIV-pairs.yml`, original source parser `scripts/analyze_Slovenia_rail_microhabitats_v103.py`. CI step **Recompute original observed same-site same-side temperature preference BS minus BT** succeeded. File names and source hashes are machine-checked. Source-derived aggregate only; no author original species-by-plot data redistributed.
