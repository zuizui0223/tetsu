# tetsu v11.5–v11.6: source-verified side aggregation, exact section-flip and dispersal-class placebo sensitivity

## Original published source and reproducibility

This is **secondary analysis** of the original Šabić, Grašič & Rovanšek (2026) Slovenian railway field observations, *Plant Biosystems* DOI [10.1007/s44473-026-00155-z](https://doi.org/10.1007/s44473-026-00155-z), CC BY 4.0. Source original unaltered Supplementary 1 XLSX SHA256 `26be7ca9513dd54136d4ec499d8413c1cb288045f87fe3cf8814731f2edac99f`. These are **NOT newly collected field samples by tetsu**. Exactly **89 original field locations / 173 BS−BT same-side pairs / 15 original railway sections**; **84** original field locations have both left and right sampled sides; **261** source plant taxa have both known original lifespan and `Diaspore.strategy` codes and appear in an eligible BS/BT field plot.

**Real original-source workflows:**
- [v11.5 actual all 89 source locations and 32,768 section swaps SUCCESS #43570032263](https://github.com/zuizui0223/tetsu/actions/runs/43570032263); [script](../scripts/analyze_Slovenia_rail_side_aggregation_section_flip_v115.py), [source unit tests](../tests/test_Slovenia_rail_side_aggregation_section_flip_v115.py), [workflow](../.github/workflows/Slovenia-RAIL-side-section-exact-null-v115.yml).
- [v11.6 ALL original dispersal strategy classes and multiplicity-adjusted matched section flips SUCCESS #43570306000](https://github.com/zuizui0223/tetsu/actions/runs/43570306000); [script](../scripts/analyze_Slovenia_rail_dispersal_class_placebo_v116.py), [source unit tests](../tests/test_Slovenia_rail_dispersal_class_placebo_v116.py), [workflow](../.github/workflows/Slovenia-RAIL-diaspore-placebo-classes-v116.yml).

All source group labels, data corrections and opposite signs are retained. Never convert a source trait species count or source occupancy into original field wind-speed, actual train-borne viable seeds, plant biomass, phenotypic plasticity, adaptation, actual train traffic or mowing observations.

## v11.5: Side aggregation does NOT create the observed species-frequency weight reversal

Define each plant species' preference as **qₛ = (n actual field locations recorded only BS − n recorded only BT)/n actual field locations recorded in BS or BT**. At each original field location, left/right BS or BT plants are pooled in the v11.0 base. v11.5 separately computes side-matched species preference within each original location, then averages side preferences **within location** before the species/section statistics. No original side counted as an independently sampled railway section.

| Source original side design | Actual trait-known taxa | Equal-species strict-wind minus other BS preference | Same-taxa source site occupancy-weight preference | Exponent α where weight = occupancy^α gives zero |
|---|---:|---:|---:|---:|
| Union of original matched L/R sides within location, all 89 sites | **261** | **+0.189388** | **−0.098125** | **0.605497** |
| Average original side-matched BS−BT presence differences inside each site, all 89 sites | **261** | **+0.184566** | **−0.099183** | **0.600880** |
| Only original **84** dual-side field sites, L/R union | **258** | **+0.184750** | **−0.095037** | **0.610736** |
| Only original 84 dual-side field sites, side-matched within-site average | **258** | **+0.179952** | **−0.096148** | **0.605536** |

Weight exponent is **not** a biological adaptation threshold. It is a user-defined survey/estimation functional: α=0 weights every original species equally and α=1 weights every species by the number of original field sites where it was recorded at least once in matched BS/BT microhabitats. For same cohort it changes weight, not species identity.

## All 2^15 EXACT source railway SECTION BS/BT label-flip diagnostic

**No random Monte Carlo approximation**: For each of the original 15 geographical railway sections, interchange its BS and BT labels consistently for every source field location in that section, producing all **2^15=32,768** combinations. Keep exact source species site-occupancy opportunities, source strict annual vs other known lifespan coding and source microhabitat coexistence fixed. For every sign pattern calculate (i) species-equal wind/other preference, (ii) source-location opportunity-weighted wind/other preference, and (iii) **difference of these two estimands**.

| Original source statistic | Observed source contrast | Two-sided exhaustive section-flip tail fraction |
|---|---:|---:|
| Species-equal wind-other BS preference | **+0.189388** | **0.0740967** |
| Source-site opportunity-weighted wind-other BS preference | **−0.098125** | **0.1602783** |
| **Equal minus occurrence-weight estimand gap** | **+0.287514** | **0.0020752** |
| Same original taxa, **within-site original track-side average**: gap | **+0.283749** | **0.0020752** |

**Interpretation**: The strong observed **weight-scheme gap** persists as an unusual statistic under exhaustive grouped label swaps. Neither constituent's original direction passes the same 5% exploratory section-flip test. This stress test alone cannot validate that BS and BT zone labels are exchangeable or randomly allocated: **BS/BT areas/soil substrate/vegetation conditions are not experimentally exchangeable**. Do not describe p=0.00208 as a causal randomized field test of train wind, seed delivery, vegetation cutting, plasticity or adaptation. The one-country sampling may produce spatially confounded source preferences.

## v11.6: Source-diaspore strategy PLACEBO panel, no favorable-code cherry-pick

Source `Diaspore.strategy` has exactly eight originally published categories **Allium, Epilobium, Cornus, Bidens, Sparganium, Phragmites, Lycopodium, Zea**. Every source category plus strict wind (Epilobium + Phragmites) and optional broad wind (+ Lycopodium) was examined on the **identical source 261 species and 89 original field sites**, with original strict annual/other known lifespan stratification and unchanged occupation frequency.

Declare a category comparison **not estimable** if fewer than three source taxa belong to either category arm in **either** lifespan stratum. Five of the ten candidate contrasts are estimable, and Bonferroni for **all five** (including two overlapping pre-existing wind definitions) multiplies each original source section-flip gap tail fraction by 5. The other five must remain listed as non-estimable, **not ecological negative results**. These Bonferroni checks are still **posthoc exploratory**, based on the same source study and unverified section BS/BT label exchangeability.

| Original source strategy group vs its complement | Source classes in strict annual / other known lifespan | Equal-species contrast | Occupancy-weighted contrast | Exact 32,768 section-flip p for **gap**, unadjusted | Bonferroni over 5 identifiable comparisons |
|---|---|---:|---:|---:|---:|
| **Strict wind Epilobium+Phragmites** | **10 / 43** | **+0.189388** | **−0.098125** | 0.002075 | **0.010376** |
| **Epilobium type alone** | **10 / 38** | **+0.200928** | **−0.096289** | 0.003296 | **0.016479** |
| Wind plus Lycopodium spores | 10 / 46 | +0.167625 | −0.059401 | 0.011230 | 0.056152 |
| **Allium type** | 40 / 123 | −0.108436 | +0.038653 | 0.105286 | 0.526428 |
| **Bidens type** | 6 / 4 | −0.137289 | −0.114707 | 0.843140 | 1.000000 |

Original non-estimable category arms: **Cornus** 2 annual/21 other, **Sparganium** 2/5, **Phragmites** 0/5, **Lycopodium** 0/3, **Zea** 1/1. Their non-estimability is a *sample support problem*, not support for zero biological effect.

**Important constraint:** Of the 53 taxa in the original defined strict wind group, **48 (10 strict annual, 38 other known lifespan) are of original source 'Epilobium' dispersal-strategy type; only 5 other-known lifespan taxa are source 'Phragmites' type**. These labels are original **diaspore strategy archetype codes** and should **not** be mistaken for 48 species belonging to the botanical genus *Epilobium*. The apparent source wind-group behavior cannot be claimed to generalize across all wind dispersal strategies, especially Phragmites and Lycopodium with insufficient source lifespan-stratum support.

**What this DOES and does NOT add:** The source main result is not a mechanical outcome of any arbitrary dispersal label; the source Allium and Bidens archetypes have different patterns in this identical design. However, overlapping categories, limited taxon support, original field plot sampling geometry, species occupancy distributions and phylogenetic taxon correlations remain important alternative explanations.

## Independent ecological replication gate is still genuinely open

Only Slovenia has a verified original paired BS/BT plant×zone matrix linked with original dispersal strategy codes at this stage. This source-backed v11.5/v11.6 sensitivity work **does not create an independently replicated railway macro law**. The original global MIREN ROAD equal-area plots provide a useful reference for a *different trait/outcome (Alien vs Native)*; they cannot be pooled or subtracted to infer a RAIL-vs-ROAD wind-transport effect.

Independent source leads:
- **French RER 2021 394 phytosociological surveys (266 embankments, 128 track surveys) and 255 taxa**, Rakotoson et al. (2025), DOI [10.1111/wre.70025](https://doi.org/10.1111/wre.70025), original data available on request; do not fake the raw matrix.
- **Kashmir Himalaya original 2014/2017 transects**, Rashid et al. (2021), DOI [10.1111/1365-2664.13961](https://doi.org/10.1111/1365-2664.13961), actual Dryad metadata DOI [10.5061/dryad.0p2ngf20h](https://doi.org/10.5061/dryad.0p2ngf20h), raw `Matrix_plot.xlsx` publicly registered but original direct download returned HTTP 401 in tested workflow. No original independent biological model fitted.
- **German Stuttgart 2024/2025 railway station floristic study**, Winkler (2026) DOI [10.3390/ecologies7020052](https://doi.org/10.3390/ecologies7020052), explicitly states original detailed data not publicly available. A species list/management discussion in the text cannot be treated as a 89-site presence-absence matrix.
- **Czech railway/agricultural gradient, four sites × five habitats**, Winkler et al. (2026) DOI [10.3390/land15040523](https://doi.org/10.3390/land15040523), source detailed data not publicly available (supplement only species naming key); not a source harmonized same-zone wind-guild replication.
- **European Vegetation Archive (EVA)** and GIVD may have many original vegetation relevés but plot distribution, railway tags, plot sizes and governance/access need direct verification. Request under source governance and preserve plot-level sampling units.

Next aim is a genuinely independent **species×railway plot×site** dataset with *known* original area, source rail traffic and/or management exposure, plus harmonizable plant dispersal/life-history traits; only then test whether the source frequency-dependent wind pattern repeats across **countries**. Otherwise frame current work as rigorous **single-country descriptive functional sampling-estimand sensitivity** linked to an independent but **different** global ROADside alien/native benchmark.

Do not merge the Draft PR as a finalized causal paper. Original main branch remains unchanged.
