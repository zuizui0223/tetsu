# v0.9: Real 1-km Daymet weather at matched NPN rail-near vs rail-far flowering sites

**Source run:** [GitHub Actions #42430207432](https://github.com/zuizui0223/tetsu/actions/runs/42430207432), artifact `usa-npn-daymet-gdd-climate-control-v09`. Executed 2026-11-08 JST. **Independent weather source from ORNL DAAC**: https://daymet.ornl.gov/web_services.html ; tool citation Daymet Single Pixel Extraction Tool DOI **10.3334/ORNLDAAC/2361**.

**This is a successful SOURCE-BACKED CLIMATE-COVARIATE AUDIT, not proof of railway phenology shifts, local adaptation, mowing or wind exposure.**

## Reproducible inputs and quality control

- Plant observation source: USA-NPN Open flowers repeated status & site-phenometrics, 2023–2025, preselected *Cercis canadensis* (MD, ID 7) and *Prosopis velutina* (AZ, ID 317). Selection by positive species×site observation opportunity and geography, never flowering-time difference or weather results.
- Exposure geometry: station coordinate ≤250m of mapped FRA 2026 rail line vs station ≥1000m; far controls are outcome-blind geographically matched registered sites ≤75 km. **Station location does not establish individual plant rail-verge habitat or actual train traffic**.
- Weather: Daymet Single Pixel API at independently registered site coordinates, 1km pixel approximation, `tmax`, `tmin`, `prcp`, annual **365-day** Daymet calendar. All three years requested with one call per independent station.
- **All 16 unique site weather requests succeeded**; expected =16, returned=16; no source or coordinate errors.
- Prespecified GDD: `max(0,(tmax+tmin)/2 - 5°C)` summed from Daymet day 1. Cutoffs **DOY60, DOY90, DOY120, DOY150** were fixed before Daymet outcomes were inspected; below, DOY60 and DOY90 are the main **pre-bloom climate-sensitivity** screens. DOY120 can fall *after* the first observed YES in MD and AZ; **never treat it as a mediator at first flower**.
- Values = **near median across independent sites minus far median across independent sites**, with each site equal weight. No individual-visit pseudoreplication, no inferential confidence intervals or p-values. The far site set can differ across years.

## Real climate differences versus prior observed flowering midpoints

| Source organism / year | Prior pre-flowering GDD5 near−far to DOY60 (°C days) | Prior/early-season GDD5 near−far to DOY90 (°C days) | GDD5 near−far to DOY120 (°C days; post-onset in some cases) | Previously observed near−far first YES interval-midpoint (days; noncausal) | Site sample |
|---|---:|---:|---:|---:|---|
| *Cercis canadensis*, Maryland 2023 | **−10.165** | −18.145 | −28.570 | **−4.75** | 2 near, 4 far |
| *Cercis canadensis*, Maryland 2025 | **+15.400** | +19.280 | +22.340 | **−7.00** | 2 near, 2 far |
| *Prosopis velutina*, Arizona 2023 | −11.795 | **−31.310** | −28.615 | **+4.00** | 2 near, 2 far |
| *Prosopis velutina*, Arizona 2024 | −8.435 | **−7.885** | −0.725 | **+0.75** | 2 near, 2 far |
| *Prosopis velutina*, Arizona 2025 | +21.395 | **+20.310** | +23.060 | **+7.75** | 2 near, 2 far |

**The temperature contrast is not stable across source years**: the sign of near−far GDD changes in both candidate species. This is consistent with real environmental differences **and potentially changing choice of eligible far comparator sites**, not evidence that railway infrastructure generated the heat differences. Species flowering contrasts are themselves first observed YES after a previous NO, with uncertain onset dates (original widths often ~5–15 days) and extremely limited independent sites.

An overly simple warmer-site-earlier-flower story is **not uniformly concordant** with these observational signs:
- MD 2023 rail-near was slightly **cooler** through DOY60 yet photographed independent NPN first-YES midpoints were **earlier**;
- MD 2025 rail-near was **warmer** and first-YES appeared earlier;
- AZ 2023 rail-near was cooler and first-YES appeared later;
- AZ 2025 rail-near was warmer and first-YES nevertheless appeared later.

These combinations **do not isolate any railway effect**: temperature, precipitation, elevation, observer frequency, garden irrigation, phenophase recording uncertainty and site habitat remain confounded, and the sample is two rail-near sites per year.

## What this adds compared with v0.8

1. Establishes an independent real **1km daily macroclimate source** with date-specific GDD and precipitation for the same strictly matched site/year cohorts. This is a substantive quality gate beyond merely comparing DOY.
2. Demonstrates why using **GDD through day 120** to explain a first flower on day 80–100 would be temporal leakage (future weather). Pre-identified DOY60/90 contrasts are safer descriptive sources; proper interval-censored thermal-onset analysis would evaluate GDD at the **NO and YES bracket bounds** separately.
3. Shows that macroclimate alone cannot be assumed equal between geographically matched stations; 75km matching is NOT sufficient. A true climate-balanced subset and independent site-level environmental metadata are needed.
4. Leaves the **original Japan railway seed-redistribution**, **management selection**, **wind biomechanics** and **fixed-photo mowing/flower ROI** mechanisms as separate, untested questions, not contradicted by the climate audit.

## Scripts and sources

- Weather parser: [`src/tetsu/daymet.py`](../src/tetsu/daymet.py), hard checks for three complete 365-day years, expected field headers, no date duplicates, physically valid tmin/tmax/precip.
- Site selection and public climate API: [`scripts/audit_npn_daymet_climate.py`](../scripts/audit_npn_daymet_climate.py).
- Source execution: [`.github/workflows/npn-daymet-audit.yml`](../.github/workflows/npn-daymet-audit.yml).
- Synthetic schema and GDD validation: [`tests/test_daymet.py`](../tests/test_daymet.py).
- Prior real flowering audit: [`npn_observed_onset_midpoint_v08.md`](npn_observed_onset_midpoint_v08.md).
- Japan history of grass cutting versus herbicide weed composition: [`japan_rail_mowing_historical_baseline.md`](japan_rail_mowing_historical_baseline.md).

**Bottom line:** weather data and measured flowering statuses now join without source failures; the resulting seasonal-climate differences are real site covariates, **not proof that rail wind, railside microclimate or cutting shifted flower timing**.
