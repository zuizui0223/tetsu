# tetsu v14.1–v14.2: REAL GPS-matched taxonomic identity versus functional profiles in Slovenian railway communities — geographic confounding changes the interpretation

**2026-11-15 verified real author-original source analysis.** Original field study by Šabić, Grašič & Rovanšek (2026), *Vegetation patterns, functional traits and invasion dynamics across railway microhabitats in Slovenia*, Plant Biosystems 160:154, DOI [10.1007/s44473-026-00155-z](https://doi.org/10.1007/s44473-026-00155-z). Authors performed all original botanical sampling in July–September 2022–2023; tetsu is a **secondary reproducible reanalysis**. Publisher open CC BY spreadsheet is source hash-authenticated SHA256 **26be7ca9513dd54136d4ec499d8413c1cb288045f87fe3cf8814731f2edac99f**, 401279 bytes. The original field sites were deliberately chosen with visible vegetation rather than randomly sampled entire network. The result is NOT independently replicated across countries and no original seed traps, train carriage, train aerodynamic exposure, or plant fitness are measured.

**Real original data executed, code/tests/Actions successful**:
- [**v14.1 actual author 89-rail-site plant identity versus four-guild functional similarity, first fixed-margin null SUCCESS #43667816461**](https://github.com/zuizui0223/tetsu/actions/runs/43667816461). [Source reader/analysis](../scripts/analyze_Slovenia_source_rail_functional_taxonomic_identity_v141.py), [correctness tests](../tests/test_Slovenia_functional_taxonomic_identity_v141.py), [source-provenance pipeline](../.github/workflows/Slovenia-original-rail-identity-functional-v141.yml). Original source workbook automatically downloaded as a previously SHA+CRC-verified GitHub Actions artifact from run 43538799032; **no invented vegetation observation or raw original publisher XLSX republished**.
- [**v14.2 source start/end GPS-based geographic and functional common-support analysis SUCCESS #43668179797**](https://github.com/zuizui0223/tetsu/actions/runs/43668179797). [Spatial matching script](../scripts/analyze_Slovenia_rail_spatial_matched_taxon_identity_v142.py), [tests](../tests/test_Slovenia_rail_spatial_matched_taxon_identity_v142.py), [original-data workflow](../.github/workflows/Slovenia-original-GPS-matched-identity-v142.yml).

## Objective and source unit

The main target for the prospective cross-study macro project is *species identity sharing above that expected from regional species pools, functional composition, environmental similarity and survey opportunities*, then its association with actual active railway movement/accessibility. **At this stage only original same-line-section membership and physical geographic separation are observed**; independent rail freight OD and road network route controls are **NOT** included. Thus current report is an in-source method/feasibility pilot, not the final network dispersal test.

Actual original published 92 Slovenian field locations on 15 rail line sections, four distinct railway microhabitat zones; 89 have legitimate SAME-SIDE BS (ballast shoulder)–BT (ballast toe) paired observations, **173** BS/BT matched sides pooled within each of 89 railway locations. No pseudoreplication of left/right side as two railway network locations.

The source records 312 original strategy/other taxon trait rows. **281** have actual BOTH author source diaspores and lifespan codes across the full trait table, but only **261** are observed in the real matched BS∪BT 89-site field cohort. The primary paired comparison retains only these actual **261** taxon codes, with **2729 original de-duplicated site×known-taxon incidences**; full-source site×taxon union **2814**, of which some traits missing, is separate sensitivity. Median original field site taxa-known fraction **0.97727**, minimum **0.86957**. No missing trait is assigned an invented wind/nonwind classification. Source sample original years 2022 **53** sites and 2023 **36**, single track **65**, double track **24**.

Functional representation is a deliberately simple **four-bin field species trait profile**, not a multidimensional continuous trait hypervolume or inherited genotype:
- source strict \`Lifespan.categories=annual\` vs other actual known lifespan source categories;
- source original published **\`Diaspore.strategy=Epilobium\` or \`Phragmites\`** as strict wind archetypes, versus other actual known original strategy archetypes.

The string \`Epilobium\` is a **strategy type** and **not** evidence a plant is of genus Epilobium; do not convert this trait proxy to measured wind dispersal by a train.

Taxonomic composition = source Jaccard similarity among same known taxa; four-guild proportional functional similarity = \(1-\frac12\sum_g|p_{ig}-p_{jg}|\). This is a total-variation overlap of coarse profiles, **not a functional beta phylogenetic hypervolume**. Source genus and regional plant status cannot be imputed from unavailable traits.

## v14.1 — genuine 89-site observation before adjusting geographic separation

All pairs: \(\binom{89}{2}=3916\). Pairs are dependent because they share original field sites: **340** pairs within the same author source railway section and **3576** pairs between sections; **15** source sections nested in **one** nation, not 3916 independent studies.

| Original source composition measure | Same author railway section | Different author railway section | Within minus between |
|---|---:|---:|---:|
| **Taxonomic Jaccard of 261 dual-known species** | **0.249413** | **0.196769** | **+0.052644** |
| Taxonomic Jaccard of ALL original source taxa (trait missingness sensitivity) | **0.244849** | **0.191331** | **+0.053518** |
| Four-guild functional profile compositional overlap | **0.843341** | **0.834281** | **+0.009060** |

Source guild-preserving null: 2×2 incidence checkerboard swaps preserve **every original site's exact number of species in each of four source guilds AND each original plant taxon's actual total source-site occupancy**. A stricter sensitivity additionally restricts swaps to original **survey year × SINGLE/DOUBLE track type**, thereby preserving species prevalence within each recorded year/type bin. Two random-walk chains, **180 draws each =360 correlated MC draws**, per chain **1,112,000 proposals** with source acceptance ~0.79 for full, ~0.75 for within-year/track. MCMC null values are correlated; they do not supply experimental p-values.

| Source 89-site identity within-minus-between | Observed | Expected under degree/guild-preserving null | Observed minus null mean |
|---|---:|---:|---:|
| Full species frequency/guild fixed-margin shuffle | **+0.052644** | +0.004499 | **+0.048146** |
| Same, but species frequency held within actual year×track type strata | **+0.052644** | +0.028787 | **+0.023858** |

All **15/15 leave-one-source-section-out** descriptive unadjusted gaps remain positive (+0.04854 to +0.06796). **But this cannot identify railroad dispersal**, since same-source-section sites can be substantially closer in space, sampled in the same year, and share environment and local geographic source species pools.

## v14.2 — ORIGINAL field GPS resolves the major geographic confounding

All **89 original field locations** have genuine published start/end **WGS84 GPS** fields; source start/end distances median **99.941m**, maximum **134.199m**, consistent with original 100m railway transect design. Source latitude 45.5323–46.67654, longitude 13.63413–16.22686; **89/89** recorded route-kilometrage locations read, but author source line names and kilometrage are **not** an externally validated train schedule, cargo direction or network topology.

| Source original site pair group | Direct great-circle distance median |
|---|---:|
| Same author railway section, 340 pairs | **17.196 km** |
| Different section, 3576 pairs | **75.239 km** |

**Critical source confounder**: same railway section corresponds to drastically shorter geographic separation. Without spatial control, the very large +0.052644 Jaccard identity gap is not interpretable as rail network movement. Geographic matching must be carried out before attributing the pattern to rail connection.

Predeclared geodesic distance bin boundaries in km: **0,5,10,20,40,80,160,500**; functional profile similarity bins: **0,0.7,0.8,0.9,1.000001**. Use only original bins containing at least **3 same-section** and **3 different-section** real site pairs; within-bin mean-gap contrast weighted by the smaller source pair count. The equal-site statistical dependence and original section/country grouping remain in place; results are descriptive (no pair-independent confidence interval).

| Source matching of ORIGINAL 89-site pair cohort | Source same-section pairs on support /340 | Different-section pairs on support /3576 | Number supported matched strata | Conditional same-minus-different source Jaccard |
|---|---:|---:|---:|---:|
| No matching | **340** | **3576** | — | **+0.052644** |
| **Same geographic distance bin** | **280 (82.4%)** | **3428** | 5 | **+0.005497** |
| **Same geographic AND four-guild functional overlap bin** | **235 (69.1%)** | **2951** | 14 | **+0.009855** |
| Same geographic, functional and same original survey-year/same original track-type match indicators | **191 (56.2%)** | **1307** | 17 | **+0.009500** |

These source supports are different subpopulations and their point estimates cannot be interpreted as a single unbiased progressive sequential effect decomposition. **Nonoverlapping source pairs are NOT ecological zeros and do not enter the matched estimator.** The original study retains all 15 source sections in leave-one-section-out v14.2: **15/15 estimable, 14/15 positive**, not proof of robust train-seed effects.

### Spatial-and-functional matched fixed-margins null (main falsification)

On EXACT same original source matched pair windows fixed in advance (geo+guild bins, 235 within /2951 between), recompute taxonomic Jaccard gap after checkerboard swapping exact guild site-richness and source taxon frequency margins; 2 source chains×160 draws=**320 CORRELATED samples**, each source chain **824,000 proposals**. The original matching strata & original GPS/functional margins remain fixed in all draws.

| Same original matched 89-site source cohort: main test | Observed conditional gap | Original null expectation | Excess above null |
|---|---:|---:|---:|
| Source guild counts / whole-study taxon occupancy preserved | **+0.009855** | +0.001347 | **+0.008507** |
| **Additional source taxon frequency preserved WITHIN 2022/23×SINGLE/DOUBLE strata** | **+0.009855** | **+0.018713** | **−0.008858** |

Illustrative MC source references **only**, not valid railway carriage p-values: whole-study null original 2.5–97.5 percentiles **[−0.008462,+0.010631]**, upper tail plus-one **0.05296**; stricter original year/track-type null percentiles **[+0.009346,+0.028009]**, upper tail **0.96885**. The interpretation does **not** hinge on crossing 0.05: the **SIGN OF EXCESS CHANGES** when actual survey year/type source recruitment prevalence is preserved. This defeats a stable generalizable same-section taxonomic species-identity surplus under strong observational confounder control. It does **NOT** demonstrate a true negative seed movement effect or nonexistence of train dispersal.

### Why the stricter null matters

When year and track type are allowed to shuffle, source uncommon species that occur preferentially in certain source years can create apparently section-specific identity. Restricting species swaps within year×track groups leaves that prevalence intact, while preserving source guild composition and all source site taxon counts. The source null expectation rises from +0.00135 to +0.01871 and exceeds original observed +0.009855. Year/type could themselves be associated with route selection and local conditions, so this is **a different conditional estimand and a sensitivity analysis**, not proof that survey year itself creates the actual local flora.

### Critical explanatory limits and data-readiness decision

1. **NO source RAIL-versus-NONRAIL comparison** in this Slovenian author workbook. Functional profile similarity 0.843 vs0.834 within vs between a rail section does **NOT** establish rail-specific convergence relative to a terrestrial habitat or ROAD flora. This result is an intra-RAIL section test only.
2. **NO actual source freight train OD, wind diaspore carriage, physical route train schedule, rail-to-rail connected-versus-unconnected geographic pair design** in this original workbook. Same section is only author-defined segment identity; cannot substitute for operational network distance, and site GPS distance is not true physical route distance.
3. **NO actual climate, geology, 100m transect canopy/substrate matching or species source regional plant pool** in present model. Matching on geography and broad four-guild type controls only coarse proxies; residual local species-pool confounding remains possible.
4. **NO pair-independent inferential p-values**: author 3916 site pairs reuse only 89 original locations nested in only 15 sections and 1 nation; all Monte Carlo fixed-margins null draws have chain autocorrelation. The 320-point ranges describe model-dependent diagnostic references, NOT cross-national confidence intervals.
5. **Original field sampling nonrandom**: railway staff screened site images and selected sufficient visible vegetation; this affects observed source taxa distributions. 2022/23 original survey followed 2021 herbicide policy change, and there is no pre-ban original matched repeated cohort or actual measured mowing, gusts, fitness, allele frequencies or plasticity.
6. **No site-level adaptation inference**: source species-level lifespan/diaspore codes are standardized traits, not same-species measured railway versus nonrail phenotypes, transplant or heritability.
7. **All source SPATIAL GPS was independently parsed from exact original author field data**; 15 source rail sections are NOT 15 independent study regions and do not supply a global effect.

**Core manuscript decision: V14.2 falsifies a naïve interpretation of positive same-section species sharing as train movement.** The empirical pattern is primarily sensitive to geography and baseline source-year/track-type occupancy. The research can still use the requested framework of *functional filtering versus extra species identity sharing*, but only for a better-designed, independent site-level railway species matrix paired with true historical active rail/road network accessibility, freight OD, habitat and regional pool controls. This study is an excellent development/falsification site, not confirmatory railway-dispersal evidence.

## Connection to the two-track v14.0 design

### B — primary prospective macro community analysis

The architecture now has a **real original-field data demonstration** of taxonomic similarity, coarse functional similarity, guild-preserving fixed row/column source-frequency nulls, true author GPS geographic matching, independent original section exclusion and no false pair-based degrees of freedom. It shows both statistical feasibility and critical falsification when stronger source controls are introduced. An across-country inference requires new source-level station/quadrat plant list sites, microhabitat metadata, source-native/non-native status, active rail vs road network and within-national geographic pools. New author-original Czech 2024 railway station S1/S2 published recent-historical station matrices and source Swedish 2024 trail/resurvey records are next potential data sources, **but not automatically equivalent rail-segment habitat observations**. The Czech 2026 Land field paper original 10m² habitat plots (83 taxa) states original per-plot raw data not publicly released, so cannot reconstruct actual 4 locality×5 habitat original incidence from marginal tables.

### A — separate propagule-stage evidence integration

Canada Gangadharan et al. 2017 actually caught freight-train spilled grain **mass** at 19 stations/sections (author whole-route 110t model extrapolation); Swiss Schulze et al. 2015 germinated **75/300** selected imported-wheat contaminated rape seeds in two laboratory trays (laboratory treatment, **not train-track colonization**); Austrian Pascher et al. 2016 detected sampled feral rape PLANTS across 60 sites/2014–15; Swiss Schoenenberger & D'Andrea 2012 observed 50 gene-marker-positive adult feral rape plants and Schulze et al 2014 detected GT73 in **2/48** non-GM maternal seed-pod samples at a port. These are DIFFERENT original dates, sites, response units and species/trait strata: do not simply multiply the marginal mean or claim complete seed-to-adult matching. Study papers may overlap a common Swiss cohort and must be deduplicated.

**No global railway seed carriage/germination/adaptation meta-analysis is claimed complete.** The primary process estimand remains source-specific, with structured evidence map while searching for independent viable on-track seed arrival and SAME cohort recruitment/reproductive trials. No stage-specific independent pairwise effect size is invented.

Original sources: Šabić et al. 2026 DOI [10.1007/s44473-026-00155-z](https://doi.org/10.1007/s44473-026-00155-z), Baselga (2012) DOI [10.1111/j.1466-8238.2011.00756.x](https://doi.org/10.1111/j.1466-8238.2011.00756.x), and [tetsu Draft PR #1](https://github.com/zuizui0223/tetsu/pull/1). Keep default main unchanged, do not merge publication-unsupported findings.
