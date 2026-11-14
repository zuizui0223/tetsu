# tetsu v10.5 — original Slovenia railway wind-vs-annual trait result survives section clustering and unknown-trait worst-case bounds

**Verified real field source run: [v10.5 SUCCESS #43560671012](https://github.com/zuizui0223/tetsu/actions/runs/43560671012)**; synthetic methods checks 6/6 PASS, all three original v10.2/v10.3/v10.4 source artifacts downloaded and checksum-aligned, exact previous v10.4 result reproduced. Reproducible code: [v10.5 analysis](../scripts/analyze_Slovenia_rail_v105_cluster_missingness.py); [tests](../tests/test_Slovenia_rail_v105_cluster_missingness.py); [source-backed workflow](../.github/workflows/Slovenia-RAIL-section-bootstrap-trait-missingness-v105.yml). Source: Šabić, Grašič & Rovanšek (2026), *Plant Biosystems*, DOI [10.1007/s44473-026-00155-z](https://doi.org/10.1007/s44473-026-00155-z), CC BY 4.0. **This is independent reanalysis of authors' published Slovenian botanical field data, not a tetsu field survey.**

Original dataset: **89 true field locations across 15 source railway sections; 173 same-location, same-track-side ballast shoulder BS vs ballast toe BT zone pairs**. One original source side pair lacks an estimable strict wind-known-group fraction; 172 sides still support 89 location-level original wind estimates. Source provenance ZIP SHA256 **26be7ca9513dd54136d4ec499d8413c1cb288045f87fe3cf8814731f2edac99f**. Two track sides are **first averaged within their location**; 15 original railway sections, not 173 independently sampled transport networks.

## 15-section resampling: actual field observation, not causal confidence interval

For each contrast, resample **15 original railway sections**, preserving all nested field-location observations, **6,000** deterministic replicates (seed 20261114). Report equal-field-location mean, equal-section mean, both 2.5–97.5% intervals, and all-15 leave-one-section-out (LOSO) sign checks. Original sites were **nonrandomly selected in a single country**. These bootstrap intervals describe sample-section uncertainty, **not** probability sampling of all world railways, randomized treatment, train-wind or mowing effects.

| Original same-site BS − BT response | Equal-site mean | 95% section-cluster bootstrap interval | Equal-section mean | All 15 LOSO signs |
|---|---:|---:|---:|---|
| Strict annual share (original known-lifespan taxa) | **+0.054703** | **[+0.017619,+0.086547]** | +0.060223 | **15/15 positive** |
| Strict Therophyte share (original known-life-form taxa) | **+0.058871** | **[+0.020774,+0.092326]** | +0.064604 | **15/15 positive** |
| Strict wind-diaspore (Epilobium + Phragmites original strategy codes) | **−0.050740** | **[−0.088988,−0.012883]** | −0.048127 | **15/15 negative** |
| Broad wind + Lycopodium spores | −0.025074 | **[−0.059113,+0.009141]** | −0.023199 | 15/15 negative |
| Original reported cover units, NOT verified percentages | −35.179775 | [−41.782080,−27.363974] | −33.293953 | 15/15 negative |

Strict wind is lower at BS in the average of **10/15** original source sections, but 5 sections have mean wind BS-positive. Among the **11 sections with >=3 site locations**, strict-wind equal-site mean **−0.049690**, clustered interval **[−0.088598,−0.009919]**, equal-section mean **−0.042542**. Strict annual within the same eligible 11 sections: +0.050470 [0.008374,0.083200]. Do not turn 15 sections into independent countries, or suppress opposing section signs.

The optional addition of Lycopodium spores **widens clustered uncertainty across zero**: the strict-wind finding is not code-invariant across all broadly wind-movable propagules.

## Original strategy missingness — explicit worst-case sensitivity

A further diagnostic reconstructs **exact original wind-classified taxon counts** from each published species×zone row (not from rounded percentages), validates the matching v10.4 strategy-known count, and gives sharp all-taxa wind-fraction bounds when every originally missing Diaspore.strategy is allowed to be wind or nonwind in each zone:

- **All original sites, original same-side pairing:** sharp mean BS−BT bounds **[−0.074992,−0.021194]**. Even the most BS-favourable assignment to all unknown original zone records does **not** make this observed mean positive. This is a missing-trait bound, not train-induced seed transport proof.
- Original known-strategy taxon coverage BS−BT: mean **−0.006316**, 95% section bootstrap [−0.018222,+0.006984]. This is not a large directional coverage imbalance.
- Both sides source-strategy coverage >=90%: **83 locations/14 sections**, mean wind BS−BT **−0.043156**.
- Both sides complete 100% strategy labels: **69 locations/13 sections**, mean wind BS−BT **−0.057179**.
- Strategy-known BS & BT in at least one side pair: all **89** locations; **172/173** pairs at all thresholds up to 75% coverage.

Do not silently label original missing traits as native, nonwind or nonannual. These bounds specifically check original Diaspore.strategy data, not unmeasured biological seeds or undetected plant individuals.

## Scientific conclusion and remaining gates

**More robust, still purely observational:** within published 2022–23 Slovenia railway vegetation, original ballast shoulder species composition is relatively enriched in short-lived annuals/Therophytes and depleted in wind-diaspore strategy taxa, even under segment-level resampling and extreme missing-strategy classification. Lower cover is also robust in original source units. This is compatible with different **potential dispersal-versus-establishment life-history axes**, but their association could partly arise from **life-history composition and species-specific occupancy**, not independent process measurements.

**Not identified:** train-carried seeds, ambient versus train slipstream, measured mowing/cutting intensity, biomass recovery, physical seed traps, freight carriage, within-species phenotypic plasticity, genetic adaptation or causal road comparison. Observations were after the 2021 Slovenian herbicide ban and lack before-treatment baseline. BS and BT survey-zone widths differ/are unverified equal, so do not claim species density per 100m² or substitute for the independently sampled MIREN roadside plots.

**Next targeted diagnostic:** condition wind-strategy composition on strict annual/nonannual source lifespan within the original same-side BS/BT species×zone records. Does wind depletion persist *within lifespan strata*, or is it a mathematical consequence of more annual species in BS? Then expand to independently sampled countries with harmonized railway/road quadrats; do not treat Slovenia as a global macro replication.
