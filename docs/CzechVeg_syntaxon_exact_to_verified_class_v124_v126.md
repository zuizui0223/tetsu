# tetsu v12.4–v12.6: original CzechVeg vegetation SYNTAXON matching reveals source selection bottleneck; class hierarchy permits wider diagnostic

## Original publication source and workflow receipts

**Study data belong to the authors of CzechVeg-Open v1.1 (2026), NOT original tetsu field collections.** Chytrý et al. (2026), original fully published 117,739 vegetation-plot source headers (source original 76,386,498 bytes, MD5 `0bc920462d1f9eac473e7b2a231e6782`) and 2,510,340 source plot×plant rows (382,939,947 bytes, MD5 `b296e5e19f03fb3bad0a8347c317f8a3`), source [10.5281/zenodo.19231152](https://doi.org/10.5281/zenodo.19231152). Original Pladias plant traits (Chytrý et al. 2021 [10.23855/preslia.2021.001](https://doi.org/10.23855/preslia.2021.001)) include Sádlo et al. 2018 original nine dispersal-strategy archetypes (Epilobium/Phragmites type = strict wind panel), and source Raunkiaer life forms (therophyte = **proxy**, not identical to strict annual lifespan). Names matched exactly, no genus-only trait guessing.

**Actual original-data workflows**, code, tests:
- **[v12.4 verified original 117739 plot SOURCE syntaxon inventory and study/area/year support SUCCESS #43585535257](https://github.com/zuizui0223/tetsu/actions/runs/43585535257)**; [script](../scripts/audit_CzechVeg_syntaxon_same_study_source_v124.py), [tests](../tests/test_CzechVeg_syntaxon_same_study_source_v124.py).
- **[v12.5 full genuine 2.51M plant source and same ORIGINAL syntaxon matched traits SUCCESS #43585966287](https://github.com/zuizui0223/tetsu/actions/runs/43585966287)**; [source script](../scripts/analyze_CzechVeg_syntaxon_matched_rail_road_traits_v125.py), [tests](../tests/test_CzechVeg_syntaxon_matched_rail_road_traits_v125.py); exact previously published original-study v12.3 source estimates also reproduced to 1e−12 before effect testing.
- **[v12.6 verified original vegetation CLASS coding hierarchy / source 117739 plot audit SUCCESS #43586479719](https://github.com/zuizui0223/tetsu/actions/runs/43586479719)**; [script](../scripts/audit_CzechVeg_original_vegetation_class_v126.py), [tests](../tests/test_CzechVeg_original_vegetation_class_v126.py). A regression test excludes original source `NA` sentinel from actual two-letter class codes: **35 genuine standalone source classes**, not 36 with a bogus `NA` class.

## v12.4 source classification completeness and EXACT SYNTAXON

The original source columns are `SYNTAXON`, `SYNTAXON_NAME`, `SYNTAX_OLD`, `SYNTAXON_OLD_NAME`. In the 4,388 original source rail-text or road-text plot candidates that also have a meaningful original publication title, an original phytoregion, source year and 1–1,000 m² plot area:

| Source original syntaxon coverage | Real plots |
|---|---:|
| Source candidates with valid same-study+region+year+area metadata | **4,388** |
| Of those, original nonmissing `SYNTAXON` | **3,845** |
| Original `SYNTAXON_NAME` | 3,845 |
| Original `SYNTAX_OLD` | 3,782 |
| Original `SYNTAXON_OLD_NAME` | 3,782 |
| Distinct source original modern `SYNTAXON` code values | **381** |
| Distinct source original legacy `SYNTAX_OLD` values | **503** |

**Source `SYNTAXON` is NOT uniform in rank**: e.g. source `XD` = *Galio-Urticetea* (vegetation class), `TDA` = *Arrhenatherion elatioris* (alliance), `MCA04` = *Phragmitetum australis* (association). Matching an exact source code is valid only as **same original vegetation code**, not evidence that all plotted codes represent equivalent hierarchy ranks.

Original RAIL-text vs ROAD-text groups are not original randomized treatments; these are historically selected source vegetation patches and site text indications, NOT independently checked active railway track-boundary locations. Original source comparison gate: common original source **publication title + phytogeographic region + decade + original source plot area class + identical `SYNTAXON`**, and at least 2 source observed plots in EACH original rail-text and road-text arm within a block.

| Source original rail/road text | Common-study/region/decade/area without vegetation restriction | Same block PLUS exact source modern syntaxon with ≥2 source plots each | Distinct source study labels remaining |
|---|---:|---:|---:|
| STRICT HABITAT-only rail or road terms | **20 RAIL / 25 ROAD** | **0 / 0** | **0** |
| BROAD HABITAT or LOCALITY site terms | **229 RAIL / 462 ROAD** | **19 RAIL / 46 ROAD**, **8** source common blocks | **8** |

Source strict HABITAT category also has 0 exact-code source blocks at even 1 per arm with all original paper, region, decade and area matching. The **0 is non-identifiability from lack of overlap, NOT a measured zero wind-diaspore biological effect**.

## v12.5 full species × original habitat code conditional source result

The full 2,510,340 source plant×plot records were read again from author-original verified Zenodo. All originally recorded missing plant traits remain missing; source strict wind archetypes are original *Epilobium* + *Phragmites* **strategy type labels, NOT botanical genus identities**. For each exact common study/region/date/area/syntaxon source block, compare per-species **rail-plot source prevalence − road-plot source prevalence**, weighted by source minimum rail/road plot opportunities per block. Then compare all source strategy-known taxa species-equal vs original plot-occupancy weighted within source Pladias therophyte/other life-form strata. This is a **conditional vegetation type comparison**, not a total causal rail land-use effect.

| Original comparison | Actual source RAIL / ROAD relevés | Source jointly trait-known taxa (wind-type) | Equal source taxa | Original plot-occupancy weighted |
|---|---:|---:|---:|---:|
| STRICT exact `SYNTAXON` min2 per group | **0/0** | **not identified** | **not estimable** | **not estimable** |
| BROAD exact modern `SYNTAXON` min1 | **117/179** | 542 (72) | **−0.000894** | **−0.015396** |
| **BROAD exact modern `SYNTAXON` min2** | **19/46** | **246 (36)** | **−0.007371** | **−0.046707** |
| BROAD original `SYNTAX_OLD` min2 (separate historical classification) | **14/36** | 194 (28) | **+0.004553** | **−0.053119** |

Original modern vs legacy syntaxon matching produces **different selected source cohorts and even a sign change of the species-equal coefficient**. Do NOT cherry-pick historical class assignments to claim a statistically robust positive rail effect. Source modern source known joint-trait coverage in the broad min2 panel is **68.84%** of original matched plot×plant incidences, not all species observed.

**Study-cluster source-bootstrap, broad exact modern syntaxon min2:** original source **8 original study titles**, 1,200 bootstrap reweight draws, only **777 estimable**, **423 non-estimable** because study resampling removed required life-form/wind archetype support. Conditional on the estimable draws:
- source species-equal 95% study resample interval **[−0.036284,+0.014012]**, includes 0;
- source plot-frequency-weighted source interval **[−0.076717,−0.009828]**, negative in that selected conditional subset;
- **difference (species-equal − occupancy weighted) [0.010554,0.063490]**. This interval is conditional on only estimable resamples and 8 historically selected original source publications, **not a valid independently replicated country-level causal p value**.
- source opposite signs in 30.63% of estimable draws; leave-one-original-study-out had sign reversal in **1/8**. Thus the original source population of study selection matters.

## v12.6: verified actual vegetation CLASS coarsening is scientifically interpretable, with real source caveats

Official [Pladias Czech vegetation classification](https://pladias.cz/en/vegetation/) reports a hierarchy where original standardized codes `XD` (Galio-Urticetea) designate **vegetation class**, `XDA` an alliance, `XDA01` an association. The original hierarchy uses the first two uppercase letters for class. We therefore construct an original-class identifier **ONLY** when the original code follows the actual Czech source pattern, and its two-letter prefix exists independently as a **standalone source class code** elsewhere in the SAME official 117,739 original plot headers. We leave all missing or unrecognized codes UNKNOWN.

- **35** valid standalone original source class roots (source `NA` correctly EXCLUDED).
- Of the 4,388 Czech original study/date/region/area labelled original plot candidates: **3,819** map to source-anchored classes, **569** have unknown/unmappable classes, of which 26 had source nonmissing but unrecognized codes. **No missing classes imputed by intuitive habitat names.**
- Exact original class match with **strict** source HABITAT rail-text vs road-text leaves **0** blocks at >=2 source plots each, even after scientifically justified rank coarsening.
- **BROAD** original HABITAT/LOCALITY candidate rail vs road, same source study, phytoregion, decade, area class, shared original class:
  - min **1** source rail and road per block: 86 original blocks / 48 source publications / 18 original vegetation classes / **189 RAIL + 325 ROAD** plots.
  - min **2** source rail and road per block: **26 blocks, 18 original studies, 12 source original classes, 105 RAIL + 167 ROAD** plots, effective minimum plot-pair block mass **78**.
  - min **3** per block: 9 blocks, 7 original studies, 6 original vegetation classes, **53 RAIL + 59 ROAD** plots, minimum block mass 44.
- Thus source class-level matching improves data retention relative to source exact individual syntaxon, but can obscure heterogeneity between associations within the same class. This is a **statistical sensitivity sample support gain**, not evidence that source rail habitat is genuinely operational railway land or comparable to Slovenian BS−BT original same-track sides.

**Critical causal note:** Vegetation class may be **mediated by railway disturbance**. Conditioning on that vegetation class can remove true habitat-sorting effects or introduce selection/collider artifacts. It is useful as a *within-broad-community descriptive comparison* and common-support gate but not a causal adjusted treatment effect. Original source railway references are site TEXT, not physical train seed movements or train wind field data.

### Next actual source guild test

The genuine source-backed v12.7 workflow compares strict vs broad source site descriptions and source verified classes using original Pladias diaspores on the SAME original author source full vegetation matrix; cluster tests at source-publication level. Do not claim a v12.7 biological result before its real original source Actions run completes. Full source code, source-verified proof and all identified or nonidentifiable source groups belong in [tetsu Draft PR #1](https://github.com/zuizui0223/tetsu/pull/1). Keep original default main branch unchanged.

Publication attribution and data use rules: cite original CzechVeg and Pladias data source papers; Pladias requests consultation with the Governing Board before substantial publication based on its database. No raw 76MB/383MB field-source CSVs nor raw original Pladias XLSX exported in result artifacts.
