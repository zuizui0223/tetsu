# tetsu v12.7 — Exact original Czech vegetation CLASS-matched wind-diaspore plant comparison fails to reproduce a robust Slovenian-style frequency sign inversion

## Authenticated original author data, not new field observations

The analysis reuses real original 2026 published CzechVeg-Open vegetation data from Chytrý et al., full **117,739** author-original vegetation plot headers (76,386,498 bytes, official Zenodo MD5 `0bc920462d1f9eac473e7b2a231e6782`) and **2,510,340** original plant×plot records (382,939,947 bytes, original MD5 `b296e5e19f03fb3bad0a8347c317f8a3`), dataset [10.5281/zenodo.19231152](https://doi.org/10.5281/zenodo.19231152), plus actual Pladias original Sádlo et al. (2018) plant diaspore archetypes and source Raunkiaer lifeform features [10.23855/preslia.2021.001](https://doi.org/10.23855/preslia.2021.001). Original published plants belong to original researchers **not tetsu**.

**Actual original source GitHub verification: [v12.7 SUCCESS #43586719008](https://github.com/zuizui0223/tetsu/actions/runs/43586719008)**; original-source code [script](../scripts/analyze_CzechVeg_original_vegetation_CLASS_wind_source_v127.py), [tests](../tests/test_CzechVeg_original_vegetation_CLASS_wind_source_v127.py), [workflow](../.github/workflows/CzechVeg-original-vegetation-class-wind-v127.yml). The true original source original v12.5 species-by-exact-SYNTAXON and corrected v12.6 canonical class root receipt were reverified on download; no mock or missing botanical source observations treated as real. Raw author-original 76MB/383MB vegetation source and public Pladias XLSX were downloaded ephemerally in GitHub Actions, **not republished** as tetsu source data.

### Hierarchical original botanical vegetation coding

The official [Pladias vegetation hierarchy](https://pladias.cz/en/vegetation/) codes `XD` as class **Galio-Urticetea**, `XDA` as an alliance, `XDA01` as an association; the first two uppercase letters identify a class in this published classification. To avoid guessing, **we recognize a prefix class ONLY if the original full 117,739-source plot original `SYNTAXON` code matches a standard published pattern and the prefix is INDEPENDENTLY observed as a standalone exact two-letter original code in that original source**. The old sentinel `NA` is explicitly excluded. The original source has **35** valid standalone class roots, not 36, and **3,819/4,388** original labelled plot candidates with real study title/region/year/area were mapped to such classes. **569** could not be mapped (including 26 nonmissing but invalid/unanchored codes). Original source classes span multiple alliances/associations, so class-level coarsening increases comparability but reduces botanical specificity.

### Target estimand / sampling gate (NOT total railway effect)

For originally railway-associated vs road-associated site **TEXT** labels in the original `HABITAT` (STRICT) or `HABITAT+LOCALITY` (BROAD) fields, group original vegetation plots into common original **publication title × phytogeographic region × original decade × recorded source plot-area bin × actual same vegetation class**. Require source railway and road candidate plots in each block, with original source plot-count sensitivity min **1**, **2**, or **3** per group. Each source species' `rail observed prevalence − road observed prevalence` is standardized across blocks using the sum of block weights `min(n_source_rail,n_source_road)`. Compare actual Pladias **Epilobium+Phragmites dispersal strategy archetype** versus other original source dispersal archetypes within source `therophyte` vs other known life-form strata. The same source taxa compared under **equal-species** vs **observed original source plot-occupancy** weights. `therophyte` is NOT Slovenia's strict annual lifespan code. No fuzzy name matching or trait missingness imputation.

**This conditional community-composition estimator is NOT the source Slovenia original 89-site matched BS-vs-BT within-railway microhabitat effect. It is also not a causal total effect of RAIL: the vegetation class may itself be shifted by the railway environment.**

### Genuine v12.7 all original vegetation-class matching scenario results

| Original source site text | Source common study / region / class match | Min source plots each per block | Original shared studies | Source RAIL / ROAD plots | Pladias joint-trait source species (wind) | Equal species wind−other relative source prevalence | Occupancy-source-weighted prevalence | Opposite signs? |
|---|---|---:|---:|---:|---:|---:|---:|---|
| STRICT HABITAT | Same paper+region+decade+source area+source vegetation class | **2** | **0** | **0 / 0** | not estimable | **not estimable** | **not estimable** | not applicable |
| STRICT HABITAT | Same paper+region+decade+area+class | 1 | 1 | 4 / 1 | 29 (7) | +0.043103 | +0.139395 | No, but nrail/nroad and 1 paper wholly unsuitable |
| BROAD HABITAT+LOCALITY | Same paper+region+decade+source area+class | 1 | 48 | **189 / 325** | 667 (88) | **−0.000754** | **−0.013270** | No |
| **BROAD HABITAT+LOCALITY** | **Same paper+region+decade+source area+class** | **2** | **18** | **105 / 167** | **532 (76)** | **−0.002921** | **−0.023405** | **No** |
| BROAD HABITAT+LOCALITY | Same paper+region+decade+source area+class | 3 | 7 | 53 / 59 | 380 (54) | +0.000227 | −0.013460 | **Yes, but equal-species coefficient near zero and tiny n** |
| BROAD HABITAT+LOCALITY | Same paper+region+decade+class, area UNMATCHED | 2 | 20 | 134 / 206 | 572 (79) | −0.003340 | −0.025030 | No |
| BROAD HABITAT+LOCALITY | Same paper+region+class, area+year UNMATCHED | 2 | 21 | 142 / 243 | 588 (79) | −0.003244 | −0.023407 | No |

**Original source main BROAD class matched min2** has **26 genuine original study×region×decade×area×vegetation class blocks, 18 original published studies and 12 source class codes, effective block min-plot mass =78**. Genuine original exact Pladias dual trait plot×species coverage **66.79%**, the remaining recorded original source incidences did NOT have both original known traits and are not forced into "non-wind".

**No robust sign reversal** is reproduced when comparing species-equal vs occurrence-weighted estimates of the same wind archetype in this truly source CLASS-matched same-study region/year/area cohort. Original Slovenian positive equal and negative occurrence-weighted wind strategy **BS-vs-BT** contrast is not independently reproduced as an identical estimand. This is a genuine independent country's functional survey but **different rail–road and land-use/text comparison design**.

### Study-level sensitivity and uncertainty (not global causal CI)

Cluster the actual original **18 published study titles** in the main BROAD class+year+area+region cohort and sample original whole studies with replacement 1,200 times. Recomputed source plant guild contrasts **1,199/1,200** draws; **one** bootstrap sample lacked source minimum strategy/life-form cohort support. Conditional bootstrap percentile95:

| v12.7 BROAD original same-study class matching, source main min2 | Observed source estimate | Original source study-resample percentile95 |
|---|---:|---|
| Wind-other difference, **species-equal weight** | **−0.002921** | **[−0.013447,+0.007329]** |
| Wind-other difference, **source plot-occurrence weight** | **−0.023405** | **[−0.051187,+0.005788]** |
| **Difference between those two estimands** | **+0.020484** | **[−0.002863,+0.045664]** |

All **three** study-resample intervals **INCLUDE ZERO**. Only **19.02%** of identifiable source resampled source studies give opposite signs; the source 18 original study leave-one-out records and all null results are in the deposited JSON. The source-study bootstrap is **conditional on nonrandom paper sampling** and must not be called a global multiple-country confidence interval, evidence of treatment assignment, or causal test of rail mortality.

### Triangulation with previous source gates

- **v12.3 original author study region/decade/area** without vegetation class: BROAD 229 source RAIL vs 462 ROAD, 731 original dual-trait taxa, equal **−0.000581**, source-frequency **−0.013025**. STRICT 20/25 plots, 177 dual-trait taxa, **both positive +0.018702/+0.025242**. This strict-source sign is **not testable after class controls**.
- **v12.4 exact syntaxon** source full original study/date/area: BROAD **19/46 original plots**, 8 studies and 246 original strategy/lifeform source species, equal **−0.007371**, frequency **−0.046707** with 423/1200 original paper cluster resamples nonestimable. Strict source **0/0 original common class plot opportunities**.
- **v12.6 canonical CLASS** produces broader source standard overlap **105/167** across 18 studies, retaining more source class levels and reducing the magnitude of estimated frequency-weight difference. Yet its source study-resampled gap interval includes zero.
- **MIREN globally standardized ROADside Alien-vs-Native source** is an independently gathered source but tests a **different guild**; do not use as replicated wind-diaspore RAIL effect.
- **Slovenia Šabić et al. (2026) strict annual and wind-diaspore BS-vs-BT source gradient** remains conditional on real source same-site matched 89 railway locations and the nonrandom one-country original field survey; it does **NOT** demonstrate viable train-borne propagule transfer, actual cutting adaptation or slipstream selection.

## Research decision

There is now strong empirical evidence of a **sampling/estimand generalizability problem**: when historical railway-related botanical plots are pooled across incompatible land-use descriptions, original studies, original site-area/time differences and vegetation types, the apparent wind dispersal association (and especially a tiny species-equal coefficient's sign) changes. The original v12.7 class-matched independent Czech source does **not** provide robust field ecological validation of the Slovenian positive/negative frequency reversal. **A manuscript claiming universal railway-specific wind-mediated dispersal/adaptation is currently unsupported**.

A more defensible contribution for a methodology + macro study is: *Can plant functional railway-corridor contrasts be transported across source sampling regimes? Identifiability, community conditioning and species-occurrence weights.* Even that is a **descriptive original source comparison**, not a direct measured arrival vs established plant causal decomposition.

### What actually identifies the requested biological mechanisms

Direct railway **seed spillage/transport** requires original viable diaspores or field seed traps on railcars/cargo/vibration corridors and actual rail traffic/cargo data. A published real-world mechanistic candidate is the Austrian study of viable oilseed rape seeds spilled along transportation routes [10.3389/fevo.2017.00104](https://doi.org/10.3389/fevo.2017.00104), and independently mapped feral plants at actual railway freight sites. **This does NOT mean tetsu currently possesses those source sample matrices.** Train slipstream plant height selection, cutting/grass biomass and flower timing require measured train gusts, vegetation management events, phenotype/common-garden and repeated plot surveys. Do not infer them from original Sádlo source strategy code.

**Next high-value gate:** obtain truly **same railway section**, same area, same protocol and ideally **paired track/embankment** plot data in additional countries (2021 French RER data [10.1111/wre.70025](https://doi.org/10.1111/wre.70025) original **394** equal-area **2m²** surveys: 266 embankments, 128 track surveys; original data on request), and a direct recorded physical seed transfer dataset. Pladias substantial scholarly use should be discussed with original Governing Board as required by source database. Keep tetsu original [Draft PR #1](https://github.com/zuizui0223/tetsu/pull/1) DRAFT, main branch unchanged.
