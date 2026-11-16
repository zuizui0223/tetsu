# tetsu v12.3 — Same ORIGINAL source publication reveals that RAIL/ROAD wind-trait sign is study-selection sensitive

## Findings: real-data feasibility YES; causal or stable macro railway trait effect NO

**Original external botanical field sources, not original field data collected by tetsu:** Chytrý et al. (2026), *Czech Vegetation Database/CzechVeg-Open v1.1*, DOI [10.5281/zenodo.19231152](https://doi.org/10.5281/zenodo.19231152), 117,739 vegetation-plot headers and 2,510,340 original plant×plot records; Sádlo et al. (2018) genuine Pladias plant diaspores strategy archetypes and Czech source Raunkiaer life forms, database Chytrý et al. (2021) DOI [10.23855/preslia.2021.001](https://doi.org/10.23855/preslia.2021.001). Authors' public source files were fetched again, source Zenodo header MD5 `0bc920462d1f9eac473e7b2a231e6782` and species source MD5 `b296e5e19f03fb3bad0a8347c317f8a3` checked, Pladias species class inventory checked against full v11.9 original-source receipt. No raw source published in result ZIPs.

**[v12.3 REAL original source workflow SUCCESS #43584571818](https://github.com/zuizui0223/tetsu/actions/runs/43584571818)** · [script](../scripts/analyze_CzechVeg_rail_road_same_study_trait_v123.py) · [tests](../tests/test_CzechVeg_rail_road_same_study_trait_v123.py) · [workflow](../.github/workflows/CzechVeg-Pladias-same-original-study-traits-v123.yml).

## Exact source study design: NOT merely an artifact of rail/road broad text co-occurrence

The previous [v12.1 six source common-support scenario results](CzechVeg_Pladias_same_region_year_area_study_gates_v121_v122.md) selected observed RAIL-text and ROAD-text plot records within **original phyto-region × original source survey decade × original plot area class**, then compared original source Sádlo **Epilobium+Phragmites** wind-type strategies against the remaining taxon categories conditional on original Pladias therophyte/other known source life forms.

The real [v12.2 original bibliography audit](CzechVeg_Pladias_same_region_year_area_study_gates_v121_v122.md) revealed that many source rail and road text candidates are recorded in **different original source field studies**. Exact `REFERENCE_TITLE` presence is used to preserve original source original botanical field sampling provenance, rather than inferring a common source from a guessed publication year or broad taxon list. Missing `NA` original publication titles excluded, not collapsed into one bogus reference. Source original `STUDY` is a **public citation title**, not guaranteed the exact same field method or experimental randomization, and original RAIL site locations are still only habitat/locality TEXT candidates, not independently confirmed right-of-way tracks.

**Six predeclared same-publication comparison gates, both original RAIL labels retained:**

- S1 source study title × phytogeographic region × decade × original area class; min **1 RAIL and 1 ROAD** source plot per block.
- **S2 main** source study title × source region × decade × source area class; **min 2 original RAIL and 2 ROAD** source plots per block.
- S3 source study × region × decade (area NOT controlled); min 2/2.
- S4 source study × region × area class (year NOT controlled); min 2/2.
- S5 source study × region only (year and area NOT controlled); min 2/2.
- S6 source study only (region, area and date NOT controlled); min 3/3.

For every valid source study block, exact observed source species prevalence RAIL/nrail − ROAD/nroad weighted by `min(nrail,nroad)` to standardize original within-block rail vs road opportunity. Species with original source unknown Pladias strategy or life form excluded (not silently assumed nonwind). **Two different estimands from the same species** — species-equal and source detected plot-occurrence weighted — were compared. No biological p or causal treatment effect is claimed.

## Actual source data results — same-study filtering can reverse the previous sign

| Original source RAIL/ROAD site-label type | Source common study gate | Original shared studies | Matched original RAIL/ROAD candidate plot counts | Known both-trait source taxa (strict wind) | Species equal Sádlo wind-minus-other prevalence | Same original species, source plot occupancy weighted | Equal − occupancy-weighted |
|---|---|---:|---:|---:|---:|---:|---:|
| **STRICT HABITAT** | S1 original study-region-decade-area min1 | 7 | 26 / 44 | 250 (41) | **+0.013156** | **+0.016011** | −0.002855 |
| **STRICT HABITAT** | **S2 min2** | **5** | **20 / 25** | **177 (32)** | **+0.018702** | **+0.025242** | **−0.006540** |
| STRICT HABITAT | S3 same study+region+decade, no area | 5 | 24 / 46 | 240 (41) | +0.016767 | +0.030205 | −0.013439 |
| STRICT HABITAT | S4 same study+region+area, no decade | 5 | 20 / 25 | 177 (32) | +0.018702 | +0.025242 | −0.006540 |
| STRICT HABITAT | S5 same study+region, no time/area | 5 | 24 / 46 | 240 (41) | +0.016767 | +0.030205 | −0.013439 |
| STRICT HABITAT | S6 same study only | 3 | 22 / 48 | 231 (37) | −0.004871 | −0.018152 | +0.013281 |
| **BROAD HABITAT/LOCALITY** | S1 original study-region-decade-area min1 | 64 | 307 / 734 | 840 (112) | **−0.000172** | **−0.008089** | +0.007916 |
| **BROAD** | **S2 min2** | **36** | **229 / 462** | **731 (99)** | **−0.000581** | **−0.013025** | **+0.012443** |
| BROAD | S3 same study+region+decade, no area | 45 | 294 / 802 | 829 (112) | +0.000052 | −0.009315 | +0.009368 |
| BROAD | S4 same study+region+area, no decade | 37 | 236 / 544 | 742 (100) | −0.000939 | −0.015301 | +0.014362 |
| BROAD | S5 same study+region, no time/area | 45 | 304 / 867 | 838 (113) | −0.000251 | −0.009816 | +0.009565 |
| BROAD | S6 same study only | 45 | 469 / 1057 | 909 (119) | +0.000387 | −0.003868 | +0.004256 |

**Every source row is observational**; do not treat 12 overlapping sensitivity conditions as 12 independent geographical replications or posthoc p<0.05 tests. Original source 71,122 field plots have measured positive area, valid year and phyto-region **AND** meaningful original study title; 29,679 source area-eligible regional/year plot rows had no meaningful source study citation title and are conservatively excluded from v12.3 source within-publication comparisons.

### Direct contrast with the previously apparently supportive result (source restriction changes conclusion)

**STRICT HABITAT RAIL/ROAD candidates, common original region×decade×area class:**
- [v12.1](CzechVeg_Pladias_same_region_year_area_study_gates_v121_v122.md) across original studies: 104 RAIL / 206 ROAD source vegetation plots, 440 original joint-known-strategy/lifeform taxa, source **equal +0.000225 versus occurrence-weight −0.023503**, **opposite signs**.
- **v12.3** adding exact **SAME original publication** as a source study block: **20 RAIL / 25 ROAD**, 177 source joint-known trait taxa, **equal +0.018702 and occupancy-weight +0.025242**, **both positive**. The sign inversion **DOES NOT survive original same-study restriction** in the strict source label.

**BROAD HABITAT/LOCALITY RAIL/ROAD candidates, same source region×decade×area class:**
- v12.1 across studies: 756 RAIL / 1900 ROAD source plots, 1144 source joint-known plant taxa, **equal +0.000656 versus occurrence-weight −0.006559**, opposite signs.
- v12.3 same-publication S2: **229 RAIL / 462 ROAD**, 731 trait-known source taxa, **equal −0.000581, occupancy-weight −0.013025**, **both negative**. Again the original source *same-guild sign inversion changes when controlling study collection provenance*.

The dramatic sample reduction when controlling study title makes the strict cohort imprecise even if numerical directions seem positive. No clustered source study confidence intervals have been calculated; do not proclaim a biologically significant effect, source country-wide trait selection or a train-carriage seed vector.

### What this contributes and what it does NOT establish

**Observed and source-backed:** Study-of-origin selection and site-description precision can change direction of the original botanical Sádlo wind-group frequency-weighting contrast, even after original source locality/year/plot-area overlap standardization. The Czech independent country source allows actual same dispersal guild analysis on up to 2.51 million genuine source field records, but **strict comparability may leave only a few dozen vegetation plots**. Original bibliographic mixed-source flora is an important macro-data-quality bottleneck to mechanistic railway inference.

**Not supported by this source:** Universal opposite signs in Czech, independent reproduction of Slovenian BS–BT *within* railway source trait filter, random causal railway-versus-road direct seed dispersal, viability of propagules on railcar/cargo, wind-speed fitness selection, frequency of grass cutting, plant size recovery, flowering phenology, genotype local adaptation or measured biomass.

**Next rigorous identification test:** Original same-paper study blocks are not fully matched by vegetation CLASS or field methods. Source `SYNTAXON`/ `SYNTAXON_NAME` classifications and actual source coordinates can be examined for coherent vegetation formation, real railroad right-of-way and source land-use. For actual macro mechanism, independent **same railway section** BS versus BT plots across countries (or real standardized rail-vs-road study cohort) and observed physical seed/mowing/gust exposure are required. Pladias database large-scale scientific use requires original Governing Board consultation. Keep the tetsu Draft PR unmerged; main remains unchanged.
