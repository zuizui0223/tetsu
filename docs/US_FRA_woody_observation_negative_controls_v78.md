# tetsu v7.8 — FRA train reports × plants: all-Plantae response comparator and closed-crossing infrastructure comparator

**Verified actual source workflow 2026-11-12:** [GitHub Actions SUCCESS #43126815842](https://github.com/zuizui0223/tetsu/actions/runs/43126815842), [original machine-readable output](https://github.com/zuizui0223/tetsu/actions/runs/43126815842/artifacts/13903114135). Seven unit tests passed, and source-tracked v7.6 woody coefficients reproduced exactly on both countries/blocks.

## Research hypothesis and evidence grade

Hypothesis: if actual railway movement causes introduced plants to appear or railway disturbance filters plant traits, the association should be more specific to railway-relevant trait responses than to general plant RECORDING intensity or unrelated map infrastructure.

What actually exists: **TWO original selected geographic US 200km blocks** (Midwest and Southeast), 14 and 13 qualified EPSG:6933 50km cells. FRA 2026 publicly mirrored OPEN road–rail crossing inventory with mixed reported train-count source years 2018–2026; original GBIF 2018–25 GRIIS-introduced indexed species facets; independently retrieved public GIFT 1.2.1 herb/woody growth forms. There are **no observed train-borne propagules, actual historic freight counts, verified trackside vegetation sampling, mowing frequency, wind gust measurements, or locally evolving genotypes**.

Same original eligible source cell cohort, exactly same log all-Plantae GBIF indexed record count + log count of currently open FRA road crossings covariates, and rank-regression gate from independently successful v7.6. No response-based geographical selection. No p-values from 14/13 spatially structured cell comparisons.

## REAL primary observational partial-rank rho, 2018–25 GBIF x mixed-year FRA reports

| Plant record outcome | US Midwest (n=14 cells) | US Southeast (n=13 cells) |
|---|---:|---:|
| GRIIS-listed **woody** expected richness after K20 indexed-record rarefaction | **−0.377168** | **−0.381404** |
| GRIIS-listed **herbaceous** expected K20 richness | −0.082780 | +0.229904 |
| GRIIS woody fraction among taxa with known GIFT growth form | −0.172056 | −0.534444 |
| **ALL Plantae** observed speciesKey category richness (source recording comparator) | **+0.353507** | −0.046017 |
| Plantae speciesKey categories NOT in source-country GRIIS list — **nativity UNKNOWN** | +0.146122 | +0.133851 |
| All GRIIS observed source speciesKey categories (raw) | −0.654119 | +0.118026 |

These values are correlations with **log1p(median reported through-trains/day at currently OPEN crossing inventory points)**, NOT a freight-train count. Conditional covariates (both source blocks): log all-Plantae indexed record count and log(1 + number of OPEN crossing rows).

The negative woody K20 sign is not a universal reduction of all recorded botanical categories: ALL Plantae does not have consistent negative rho. However this **does not validate the railway-specific biology**: GRIIS membership is incomplete, non-GRIIS taxa are not certified native; source datasets, observer effort, climate, urban/road gradients, and spatial record processes remain confounded.

## Comparison to infrastructure history (not a randomized negative control)

An additional **CLOSED crossing fraction** predictor, derived from all current crossing records in the exact same cell, is included only as an observational geography/history diagnostic. It is not active train traffic. Conditional on the same record/open-crossing covariates:

| Botanical response | Midwest closed-crossing-fraction partial rho | Southeast partial rho |
|---|---:|---:|
| Woody K20 expected GRIIS richness | **+0.476371** | **+0.501311** |
| GRIIS woody fraction among known growth-form taxa | +0.243218 | +0.494856 |
| ALL Plantae observed speciesKey categories | −0.300905 | −0.030142 |

Positive closed-crossing-fraction correlations alongside negative current reported traffic correlations are compatible with unmeasured transport geography, historical land use, vegetation and botanical reporting differences. They are **not** evidence that closed rail lines increase actual woody plants, nor experimental disproof of railway transport.

## Crucial positive-through-trains exposure boundary

Only cells whose median reported OPEN-crossing train count is **positive**:
- Midwest: **13** source cells, woody K20 rho **−0.225759**, woody known-growth-form fraction rho **−0.030670**, all-Plantae recorded species richness rho **+0.217279**.
- Southeast: **8** source cells, below original 10-cell support gate; all source-controlled active-train-only correlations are **NOT ESTIMABLE**. Do not relax the gate to force agreement with Midwest.

This distinction matters because some OPEN crossing inventory records carry zero reported through trains, and because an OPEN road crossing is not equivalent to measured freight activity on a line.

### Historical reporting limits

FRA is a CURRENT 2026 inventory; report years of underlying crossing traffic vary. Open crossings from the analyzed Midwest cells include 3 historical 2017-source reports, and Southeast includes 348 2017 reports (these older source reports are excluded from the v7.6/v7.8 train-count median). The date of latest crossing report does NOT validate past 2018–25 plant-time rail activity or operational network direction. The 2026 year/source categories are tracked separately rather than asserted equivalent.

## Scientific conclusion

1. **Technically supported:** independently joined 2026 US FRA crossing-level train count opportunities to frozen 50km GRIIS/GIFT botanical source data; robust exact baseline reproduction across n=14/13 cells; the apparent woody recorded-richness/train-intensity association differs from an all-Plantae recorded-species comparator.
2. **Macro ecological support insufficient:** only one Midwest block retains a minimum-sample active-positive-train gradient; Southeast and the a priori Pacific third US block fail necessary source support. Closed-crossing map history and GBIF source composition remain substantial rival explanations. No well-supported globally general woody railway association is established.
3. **Neither delivery nor environmental filtering is identified**. No actual freight seed carriage; no date-matched active service, railway verge colonization plots, measured cutting, slipstream gusts, plant stem-height plasticity or evolved adaptive change.
4. **Next gate:** independent true FRA Form 71 traffic-year verification or freight service-flow source linked to geographic track segments, independently measured rail-versus-road habitat/management controls, historical matched dates, and stronger geographically independent botanical replication. Avoid more posthoc trait/metric searches as a substitute for data that identify seed delivery and establishment.

Implementation: [v7.8 reproducible source code](../scripts/analyze_FRA_train_plant_negative_controls_v78.py), [tests](../tests/test_FRA_train_plant_negative_controls_v78.py), [workflow](../.github/workflows/US-FRA-traffic-plant-source-controls-v78.yml).
