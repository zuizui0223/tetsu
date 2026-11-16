# tetsu v10.6: original RAIL wind-diaspore deficit persists WITHIN source lifespan strata

## Real original-source success

**[v10.6 SUCCESS #43561008890](https://github.com/zuizui0223/tetsu/actions/runs/43561008890)**; 6 decomposition/source tests PASS, original 2026 Slovenian railway field workbook checksum and v10.3/v10.4/v10.5 previous source outputs verified as unchanged. Reproducible [source decomposition](../scripts/analyze_Slovenia_rail_v106_annual_wind_decomposition.py), [tests](../tests/test_Slovenia_rail_v106_annual_wind_decomposition.py), [workflow](../.github/workflows/Slovenia-RAIL-annual-wind-composition-v106.yml). Original source authors Šabić, Grašič & Rovanšek (2026), DOI [10.1007/s44473-026-00155-z](https://doi.org/10.1007/s44473-026-00155-z), CC BY 4.0. **Original external authors' field data, NOT original field collection by tetsu.**

### Why this diagnostic matters

Original [v10.4](Slovenia_real_RAIL_diaspore_wind_vs_annual_v104.md) and [v10.5](Slovenia_rail_section_cluster_and_missing_trait_sensitivity_v105.md) document BS−BT **+0.054703** fraction strict annual taxa but **−0.050740** strict wind-diaspore (original Epilobium/Phragmites strategy). Such opposing signs need **not** imply separate ecological filters, because taxa bearing wind-diaspore traits may also be systematically different in lifespan. A two-group **exact identity** decomposes the fraction of wind taxa among species with BOTH lifespan and diaspore strategy coded:

wind fraction = P(strict annual) × P(wind|strict annual)
              + P(other known lifespan) × P(wind|other known lifespan).

BS−BT wind difference =
**annual-lifespan composition component**
+ **wind-within-lifespan component**.

This is the symmetric Kitagawa decomposition. IMPORTANT: “other known lifespan” includes mixed annual/biennial and other non-exact-annual original source labels; it must NOT be called strict perennial. Original traits are species-level and not measured changes in an individual phenotype.

### Real site/section inference

The published real railway source includes **173 true matched same-side BS/BT comparisons at 89 field locations in 15 sections**, but **36 source side pairs lack both lifespan categories in at least one BS or BT zone**, preventing their two-stratum decomposition. We did **not** code empty categories as zero. The eligible real subset comprises **137 matched zone-side pairs**, averaged within **82 distinct original field locations / all 15 railway sections**. Result is conditional on this subset, not an imputation for 89/173.

**6,000 deterministic SECTION-cluster bootstrap replicates** (not site/side as independent sections):

| Original source joint-known BS−BT contrast | Mean fraction difference | Section-cluster 95% interval |
|---|---:|---:|
| **Total strict wind fraction among jointly known lifespan × strategy taxa** | **−0.037128** | **[−0.072362,−0.004974]** |
| Annual/other-known lifespan **compositional shift component** | **+0.003296** | [−0.001877,+0.007022] |
| **Residual within-lifespan wind-strategy component** | **−0.040424** | **[−0.073793,−0.008518]** |
| Within strict-annual contribution | −0.021267 | [−0.037885,−0.003630] |
| Within other-known-lifespan contribution | −0.019157 | [−0.050466,+0.006331] |
| Original v10.4 known-wind fraction recomputed on SAME 82 eligible locations | −0.038679 | [−0.074255,−0.005739] |

All **15/15** leave-one-rail-section-out recalculations remain negative for the **within-lifespan** component; source range **[−0.048609,−0.032020]**. The 75 eligible sites inside source-eligible >=3-site rail sections show within-lifespan component interval **[−0.074008,−0.006211]** (does not cross zero). Some source sections and plant category terms differ in sign: do not ignore counterexamples.

**Original 312 trait-coded taxa inventory**: strict annual 76 species (10 strict wind, 62 other-known dispersal, 4 strategy-unknown); other-known lifespan 217 (45 strict wind, 164 other-known dispersal, 8 strategy-unknown); lifespan-unknown 19 (1 strict wind, 2 other-known dispersal, 16 strategy-unknown). These totals represent original source taxa, not transport observations.

### Scientific conclusion and noncausal firewall

**SUPPORTED original Slovenian RAIL field-species-composition result:** relative strict-wind-diaspore strategy depletion in BS is **not explained by a simple mathematical switch toward strict annual plants** among jointly trait-known source taxa. The life-history-conditioned component stays negative in section-cluster bootstrap. Wind and annual reflect distinguishable **descriptive trait axes**, although not independent mechanistic selection parameters.

**NOT supported / not tested:** train slipstream causing the wind trait gradient; wind-borne seed arrival, train-wheel/cargo dispersal; seeds' viability; mowing frequency, herbicide causal effect, seed mass, ballast microclimate exposure or phenotype plasticity; genetic local adaptation; RAIL-vs-ROAD causal contrast. The original one-country source was nonrandomly selected, BS/BT plot widths not equal-area verified, original 2022–23 after a 2021 herbicide policy change without treatment-before baseline. Original author findings and their own analyses must receive proper citation; tetsu is a secondary reanalysis.

**Next cross-country macro gate:** acquire original similarly interpretable Czech railway track/building/surroundings floristic matrices (Kutlvašr et al. 2024, DOI 10.23855/preslia.2024.247) and Slovak station plots (Májeková et al. 2025, DOI 10.1016/j.flora.2025.152780). Harmonize only data columns that truly correspond (within-region contrasts and source taxa); do not pool Czech stations and Slovenian BS/BT as identical sites. For a mechanistic separation of moving seeds and establishment, direct seed-trap / rail-carriage observations are still required.
