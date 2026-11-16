# tetsu v10.0–v10.1: Roadside alien plants show higher gains and losses, plus species replacement rather than simple nested additions

## Provenance and actual verified analyses

The field survey belongs to the **Mountain Invasion Research Network (MIREN) consortium**, NOT to the tetsu project. Original global standardized MIREN Road Survey v3 2007–2023, [Zenodo 10.5281/zenodo.15749536](https://doi.org/10.5281/zenodo.15749536), MIREN data paper [Global Ecology and Biogeography 10.1111/geb.70284](https://doi.org/10.1111/geb.70284).

- [v9.8 raw CSV source SUCCESS #43212464976](https://github.com/zuizui0223/tetsu/actions/runs/43212464976): 17,899,241 bytes, 167,804 original species occurrence rows, confirmed MD5 b640195cc44f252ccb6ddb28eb6fff13.
- [v9.9 standardized paired road-verge field plots SUCCESS #43214312718](https://github.com/zuizui0223/tetsu/actions/runs/43214312718): 1,876 same-site/year road-edge Plot 1 versus interior Plot 3 survey pairs across 1,005 unique road-site sources, 24 usable region codes; first-visit region sign ROAD edge Alien richer 23/24 codes, 21/22 codes with >=20 distinct road sites.
- [v10.0 original road-plot revisit taxon gains/losses SUCCESS #43215926370](https://github.com/zuizui0223/tetsu/actions/runs/43215926370), [actual result JSON #13943811420](https://github.com/zuizui0223/tetsu/actions/runs/43215926370/artifacts/13943811420), [code](../scripts/analyze_MIREN_roadside_alien_resurvey_turnover_v100.py): 6 original-source tests PASSED.
- **[v10.1 real ROAD Alien/Native Jaccard beta species-turnover/nestedness SUCCESS #43538320965](https://github.com/zuizui0223/tetsu/actions/runs/43538320965)**, [actual complete JSON #14089493895](https://github.com/zuizui0223/tetsu/actions/runs/43538320965/artifacts/14089493895), [reproduction code](../scripts/analyze_MIREN_road_alien_beta_components_v101.py): 8 mathematical/source QA tests PASSED. Global repository smoke tests successful [#43538320830](https://github.com/zuizui0223/tetsu/actions/runs/43538320830) and [#43538324979](https://github.com/zuizui0223/tetsu/actions/runs/43538324979).

Original source protocol is direct vegetation inventory of equal-area 2m × 50m (100m²) mountain ROAD plots: Plot 1 immediately alongside road, Plot 2 intermediate 2–52m, Plot 3 interior 52–102m. Region-specific Alien/Native statuses are used, two original contradictory region×species nativity statuses conservatively set UNKNOWN; duplicate plot-year taxon lines deduplicated; 1A/2A/3A replacement plots are excluded and missing entire plot-years are NOT imputed as zero species. Original 26 region CODES are not necessarily 26 independent countries (data paper describes 25 mountain regions); 24 codes offer original strict Plot1/3 pairs.

## v10.0: Actual 484 ROAD original field plot revisits

484 distinct original ROAD sites across 11 region codes have repeated same-plot Plot1 and Plot3 at least twice, separated 5–16 years. Species first-to-last gain means LATE observed but not EARLY; loss means EARLY observed but not LATE. Neither is a true immigration or extinction event.

| Exact first-to-last MIREN Alien taxon count per original 100m² plot, 484 sites | ROAD edge Plot1 | Interior Plot3 |
|---|---:|---:|
| Mean recorded taxon gains | **3.640496** | 1.995868 |
| Mean recorded taxon losses | **2.807851** | 1.582645 |

- Mean roadside-minus-interior gain +1.644628; regional mean gain higher at ROAD edge in **10/11 codes** with >=10 repeated field sites.
- Loss is ALSO elevated at road edge, pooled Plot1-minus-Plot3 mean +1.225207.
- After dividing taxon gains by UNION of first and last taxa for each actual surveyed plot, and restricting to 337 sites with valid nonempty Alien taxa union at BOTH plots, mean roadside-minus-interior relative gain is **−0.036670**, region-level sign positive only **4/11**, negative 7/11. Absolute recorded gains being higher could be partly expected where richer source pools are present.
- Road edge Alien richness advantage CHANGE first-to-last: **5/11** region mean increases vs **6/11** decreases, median site-level change zero. There is no region-general temporal strengthening.
- Revisited field taxa are stronger local source observations than GBIF country-state first records, but imperfect plot detection and nonrandom roads still prevent true colonization and transport attribution.

## v10.1: Actual spatial ROAD species composition — replacement versus nestedness

Research question: Is near-road Alien flora primarily an enriched superset of interior flora (nestedness), or a different assembly (species replacement)?

For taxon groups observed at BOTH originally standardized roadside and interior plots of a same original ROAD site/year, define shared taxa a, roadside-only taxa b and interior-only taxa c. Jaccard species dissimilarity beta_jac = (b+c)/(a+b+c). Replacement/turnover beta_jtu = 2 min(b,c) / [a+2 min(b,c)]. Nestedness-resultant beta_jne = beta_jac − beta_jtu (Baselga-style identity; all mathematical tests pass).

**Eligibility is critical**: calculate components ONLY if both original Plot1 and Plot3 have at least one positively recorded taxon of that group. The original first-visit 1,005 unique ROAD sites comprise **583 Alien-present in both**, **266 Alien recorded roadside-only**, **17 interior-only**, and **139 neither plot recorded Alien taxa**. The remaining 422 cannot be silently assigned a zero turnover or a biologically meaningful Jaccard partition.

| Original first-year distinct ROAD site sources with Alien taxa on BOTH plots, n=583 | Mean Jaccard beta |
|---|---:|
| Total observed Alien assemblage dissimilarity | **0.695662** |
| Taxon species replacement/turnover component | **0.435966** |
| Nestedness-resultant dissimilarity component | **0.259697** |

Turnover accounts for **62.6692% of the SUM of observed Jaccard dissimilarities** in these both-positive pairs, rather than a mean of individual shares. At 379 of 583 both-positive pairs the turnover component exceeds nestedness; 184 show reverse (20 tied).

A region-balanced analysis using ONLY FIRST survey year at each original road site avoids counting repeated survey years as independent site replicates:
- Minimum **10 both-positive distinct original ROAD sites per region**: 17 qualifying region codes, **15/17** with greater region-mean turnover than nestedness; median across regions of regional mean turnover-minus-nestedness **+0.196753**.
- Minimum **20 both-positive sites per region**: 14 qualifying codes, **12/14** turnover dominant. ARS and AUN remain legitimate counterexamples favoring nestedness.
- Mean roadside-only recorded Alien taxa across ALL 1,005 first-year paired sources = **4.359204**, interior-only = 1.019900. This is an incidence description incorporating edge-richness difference, NOT a seed-source measurement.

**Native-plant methodological comparison**: Restrict unique first visits to actual years >=2017 because some older 2007/08 MIREN source protocols undercounted natives. Across 459 true original site-year pairs, 420 have Native taxa on BOTH plots; among these mean native beta_jac **0.763412**, beta_jtu **0.627570**, beta_jne **0.135842** (turnover constitutes **82.2059%** of summed dissimilarity). Native group sample differs from the Alien conditional sample, so this is not a formal randomized interaction test. But it demonstrates that species replacement is also dominant in NATIVE taxa — an Alien-specific roadside wind/mowing mechanism is not established.

## Scientific meaning and limits

**Supported real biological source pattern:** Road-edge Alien richness is higher in many MIREN mountain regions (v9.9) and both recorded gains and recorded losses are larger at road-edge repeated field sites (v10.0). Among sources with Aliens on both road-edge and interior, many paired assemblages exhibit substantial DIFFERENT-species replacement, not merely nested additions of roadside alien taxa (v10.1). This supports a descriptive road-associated spatial community difference.

**NOT established:** the mechanism (delivery of propagules vs local plant selection after germination). Roadside plots differ in light, openness, soil, maintenance, anthropogenic disturbance, existing pools and detection; the original data have no standardized mowing, vehicle gust/wind, cargo seed assays or genetic measurements. They are ROAD plots, not active railway plots. Native comparison shows turnover is not Alien-specific. No actual railway freight line, train-carried seeds, trackside establishment, phenotypic plasticity or evolved adaptation was measured.

This is a SOURCE REANALYSIS of an original MIREN consortium dataset, not a novel raw field dataset owned by tetsu. Cite consortium and discuss manuscript use with its original researchers. Never pool a Polish rail-flora 10km grid-centroid species list or a local metro construction survey as if it were same 100m² paired railway plots.

**Next experimental/source identification gate:** Acquire actual RAIL-versus-ROAD equal-area paired verge/interior vegetation plots in several independent ecoregions, dated active train frequency, measured rail/slipstream gust, documented mowing/herbicide and viable seeds from rail cargo/seed traps. Predefine same-region interaction rail_vs_road × distance_from_route × independently measured dispersal/regrowth trait; without this a direct train seed vector or rail-stress adaptation macro result cannot be claimed.
