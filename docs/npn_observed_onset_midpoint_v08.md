# Observed Open flowers interval midpoint by site — real-data v0.8 diagnostic

**Source run:** [GitHub Actions #42391174556](https://github.com/zuizui0223/tetsu/actions/runs/42391174556), `npn-2023-2025-common-taxon-open-flowers-individual-status-check` artifact. Script: `scripts/audit_npn_candidate_flowering_status.py`.

**This is a DESCRIPTIVE photographic-independent NPN observation comparison, NOT a causal railway flowering, cutting or wind adaptation estimate.**

## Measure fixed before looking at directional differences

From each eligible monitored NPN individual-year (≥4 Open flowers status dates and at least one recorded YES and NO), retain:
- the latest **observed NO** date preceding the individual's first observed YES of year;
- the first observed YES date, provided the interval is **≤21 days**;
- observed-transition interval midpoint, measured as day of year; this is **not the actual first flower date**;
- per registered **site**, median of eligible individual-year midpoints;
- per exposure group (near/far) **median of site medians** with equal site weight; no record-/plant-level pseudoreplication.

Define **exploratory relative difference** = rail-near (station centre ≤250m FRA 2026 railway line) median minus matched rail-far (centre ≥1000m) median, with far station geographically within 75km of near station **and selected before observing flowering timing**. Negative is earlier recorded flowering toward rail, positive is later.

**Site eligibility conditional on positive Open flowers phenometric and repeated yes/no makes this an observed-event timing diagnostic only, NOT occurrence probability or population-level flowering probability.** No p-values or formal inference.

## Source-backed unadjusted interval-midpoint differences

| Species / state / year | Site-adjusted near midpoint DOY | Far midpoint DOY | Near − far, days | Sites with bracket near | Sites with bracket far | Typical bracket widths (near / far, days) |
|---|---:|---:|---:|---:|---:|---|
| *Cercis canadensis* MD 2023 | 91.50 | 96.25 | **−4.75** | 2 | 4 | 9 / 7 |
| *Cercis canadensis* MD 2025 | 80.00 | 87.00 | **−7.00** | 2 | 2 | 7 / 14 |
| *Prosopis velutina* AZ 2023 | 112.75 | 108.75 | **+4.00** | 2 | 2 | 4.75 / 11.5 |
| *Prosopis velutina* AZ 2024 | 115.75 | 115.00 | **+0.75** | 2 | 1 | 6.5 / 8 |
| *Prosopis velutina* AZ 2025 | 105.75 | 98.00 | **+7.75** | 2 | 1 | 14.5 / 8 |

**Observational description**: In two inspected years *C. canadensis* has earlier observed-onset bracket midpoints at the rail-near stations; in three inspected years *P. velutina* has later ones. The values are close to the bracket widths and sample of rail-near sites is only 2 per year. An average positive/negative sign is **not evidence for a railway effect**, nor proof of biological onset shift, plant stress adaptation or climatic selection.

Additional source limitations:
- Neither woody tree is a valid stand-in for repeatedly mown herbaceous trackside plants; site coordinates may reflect residential gardens or parks, not railway verge soil.
- Station rail adjacency measured against **2026 FRA rail**, not verified active passenger/freight status in 2023–2025.
- Latitude, road access, urban heat, climate/GDD, altitude, observer visit frequency, horticultural origin, watering and soil conditions are uncontrolled. In MD 2025, the closest geographically matched far site was about **73 km** away.
- We selected positive site-level flowering records before raw status lookup. Missing activity, no-flower site-years and site enrollment shifts can bias timing.
- Observed interval-midpoints are a convenience for **screening**; proper phenological modeling must use **interval-censored data** and incorporate days of site monitoring and observer/site clustering.
- A single active registered site may contain multiple plant individuals and years; do not use 7 individuals at 2 sites as 7 independent railway treatments.

## Next legitimate test

For phenology, find a larger herbaceous **actively monitored** rail/road/background fixed-plot cohort, validate true railway verge exposure and dates, and balance climate/GDD and management. For mowing-linked plasticity, use the fixed-frame patch ROI plus verified cutting-date analysis already implemented (no real before/after cutting data yet). For train-mediated seed redistribution, return to the independent Japan freight/port/rail road-network hypothesis.

**Supported statement:** source-backed repeated open-flower observations show measurable, variable site-group contrasts of a few days in two woody species. **Unsupported:** railway-specific environmental adaptation, genetic evolution, or causally shifted flowering.
