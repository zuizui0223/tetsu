# tetsu v12.1–v12.2 — CzechVeg railway-vs-road source comparison after original region, date, area and original publication controls

## Executive evidence statement

**Actual original published CzechVeg-Open v1.1 2026 vegetation data and public Sádlo et al. (2018) Pladias dispersal traits were analyzed in full, with original file checksums.** Independent original field data: Chytrý et al. (2026), CzechVeg [10.5281/zenodo.19231152](https://doi.org/10.5281/zenodo.19231152) (117,739 original vegetation plot headers; 2,510,340 species-by-plot records), and Pladias [10.23855/preslia.2021.001](https://doi.org/10.23855/preslia.2021.001), nine real Sádlo 2018 dispersal-strategy archetypes and Raunkiaer lifeforms. These are **other researchers' field data, NOT novel tetsu original observations**.

- **[v12.1 complete original-source data+biological join SUCCESS #43583982369](https://github.com/zuizui0223/tetsu/actions/runs/43583982369)** — six prespecified **regional / year / plot-area** overlap designs for strict HABITAT-only RAIL-vs-ROAD source text vs broader HABITAT/LOCALITY RAIL-vs-ROAD source text, original full taxon×plot source checksum reverified, matched original Pladias species by exact accepted Latin names, all unknown strategy or life-form taxa omitted rather than coded nonwind. Real [script](../scripts/analyze_CzechVeg_rail_road_region_year_area_overlap_v121.py), [tests](../tests/test_CzechVeg_rail_road_region_year_area_overlap_v121.py), [workflow](../.github/workflows/CzechVeg-Pladias-original-region-year-area-v121.yml).
- **[v12.2 complete original-publication study-reference source provenance SUCCESS #43584215131](https://github.com/zuizui0223/tetsu/actions/runs/43584215131)** — original `REFERENCE_TITLE` provenance across 117,739 official source field plots; study ID treated as **valid original publication title only if nonblank, non-NA and meaningful**; original species and trait data not required for this gate. Real [script](../scripts/audit_CzechVeg_shared_railroad_reference_studies_v122.py), [tests](../tests/test_CzechVeg_shared_railroad_reference_studies_v122.py), [workflow](../.github/workflows/CzechVeg-original-same-study-railroad-v122.yml).

## v12.1: Six ORIGINAL source common-support definitions (NOT favorable-results subset choice)

Always require recorded 1–1000m² source plot area and source phyto-region and valid original YEAR (1850–2026). Exclude missing source plot-year and source regional code, and do not treat source bibliographic-only railway mentions as site labels. For all regions/area/period blocks, retain a block only when BOTH source railway-text and road-text candidate vegetation plots exist with >=specified original measured plot counts, then **within that block** compute original observed strict wind-typed taxon source prevalence rail/nrail − road/nroad. Standardize over blocks with weight `min(nrail,nroad)`, avoiding unmatched source years/area/regions. Compare strict wind Epilobium+Phragmites archetypes vs other source original Sádlo strategy categories, adjusting for source Pladias `therophyte` vs other known lifeform (NOT equivalent to Slovenia exact `annual`). Species-equal and source matched-plot occupancy weighted estimands use identical original source taxa.

Predeclared matching definitions:
- **A**: original phyto-region × 10-year bin × **exact original plot area**; >=2 source plots in each arm
- **B**: original region × 10-year bin × plot area class (1–4, >4–16, >16–64, >64–256, >256–1000 m²); >=2 each
- **C**: original region × 20-year bin × plot-area class; >=2 each
- **D**: original region × 10-year bin, NO area matching; >=3 each
- **E**: original region × plot-area class, NO year matching; >=3 each
- **F**: original region only, NO year or area matching; >=5 each

**All six original comparisons retained; NO multiplicity-controlled hypothesis test or independent source-country sampling claim.**

| Original botanical rail/road text cohort | Matched rule | Original matched RAIL plots | Original matched ROAD plots | Source original joint-known diaspore×lifeform taxa | Species-equal wind-minus-other RAIL relative prevalence | Source observed-plot occupancy weighted | Source estimand signs differ? |
|---|---|---:|---:|---:|---:|---:|---|
| Strict original HABITAT | **A exact area** | 58 | 78 | 266 | +0.002417 | **−0.046645** | Yes |
| Strict original HABITAT | **B decade + area class** | **104** | **206** | **440** | **+0.000225** | **−0.023503** | Yes |
| Strict original HABITAT | C 20yr + area class | 108 | 227 | 452 | +0.000182 | −0.023922 | Yes |
| Strict original HABITAT | D region + decade, no area | 119 | 312 | 552 | −0.004388 | −0.036088 | No |
| Strict original HABITAT | E region + area, no year | 107 | 218 | 456 | +0.001337 | −0.023833 | Yes |
| Strict original HABITAT | F region only | 127 | 322 | 583 | −0.001746 | −0.025654 | No |
| Broad HABITAT+LOCALITY | A exact area | **468** | **999** | 985 | **+0.001109** | **−0.008744** | Yes |
| Broad HABITAT+LOCALITY | **B decade + area class** | **756** | **1,900** | **1,144** | **+0.000656** | **−0.006559** | Yes |
| Broad HABITAT+LOCALITY | C 20yr + area class | 911 | 2,217 | 1,194 | +0.000470 | −0.005020 | Yes |
| Broad HABITAT+LOCALITY | D region + decade, no area | 1,001 | 2,410 | 1,217 | +0.001602 | −0.000694 | Yes |
| Broad HABITAT+LOCALITY | E region + area, no year | 1,058 | 2,523 | 1,233 | +0.001170 | −0.003165 | Yes |
| Broad HABITAT+LOCALITY | F region only | 1,249 | 3,181 | 1,302 | +0.001347 | −0.000988 | Yes |

**Most important finding:** A **small** weight-dependent sign switch can appear in Czech railway-vs-road candidate assemblages after matching area/decade/region. However, **all Czech species-equal estimates are extremely close to zero** (+0.0002 to +0.0024 in matched reversal designs), and the original observed-plot occupancy estimand is modestly negative. The raw v11.9 broad COUNTRY-POOLED comparator had *both positive* (+0.001131/+0.003015) — stronger source-geometry restriction materially changes its biological target. Never call tiny +0.000225 or +0.000656 a strong wind guild preference or proof of the Slovenian effect. Original sample exclusion is substantial. No v12.1 cluster bootstrap CI was estimated: sign switches are exploratory descriptive observations ONLY.

**Uncertainty/quality:** The strict-source B rule retains 20 eligible original source region×decade×area blocks across **12** phyto-regions (total 104/206 actual plots, 78 effective min-plot-pair mass). The broad-source B rule retains **126** original blocks across **50** phytogeographical regions (756/1,900 plots, 602 effective min-plot-pair mass). Source exact jointly known strategy+lifeform plot–taxon coverage in these matched cohorts is **about 66–68%**. Missing traits/taxa can bias source frequency and guild preference, so absent trait codes must remain unknown. The exact source accepted species names and generated source Pladias XLSX biological aggregate/source provenance remain v11.9 verified; original CzechVeg exact official published header MD5 `0bc920462d1f9eac473e7b2a231e6782`, species CSV MD5 `b296e5e19f03fb3bad0a8347c317f8a3`.

**Structural limitations:** These original Czech sample sets are NOT same-location BS versus BT, their vegetation communities and original study methods differ, and the railway-text candidates have not been independently confirmed to lie within an active railway right-of-way. Original field original source columns `SYNTAXON, SYNTAXON_NAME, SYNTAX_OLD, SYNTAXON_OLD_NAME` may contain vegetation classification, but have **NOT yet been harmonized in a defensible way**. Regional/year/area adjustment alone therefore cannot isolate environment, management or dispersal mechanism.

## v12.2: Is the original sample even observed in the SAME BOTANICAL STUDY?

Use `REFERENCE_TITLE` in the original 117,739 source plot header records, rather than globally pooling heterogeneous historic author's field studies. Blank, `NA` and dubious original reference labels cannot identify shared actual study; exclude rather than placing them into a single fake "NA source" cohort. Original association with publication title does **not** establish matched vegetation types, standardized source observations or original railway management.

| Source metadata gate | Strict original HABITAT rail/road plot candidates | Broader original HABITAT/LOCALITY candidate plots |
|---|---|---|
| Eligible original plots with source area | **195 RAIL, 681 ROAD** | **1,580 RAIL, 5,672 ROAD** |
| Source plots missing a meaningful `REFERENCE_TITLE` | **128 RAIL, 391 ROAD** | **486 RAIL, 2,137 ROAD** |
| Original publication with ≥2 RAIL and ≥2 ROAD plots in SAME original source title | **7 source studies; 30 RAIL, 80 ROAD** | **75 source studies; 626 RAIL, 1,381 ROAD** |
| SAME original publication AND phytogeographical region AND survey decade AND plot area class, ≥2 each source arm | **6 original strata; 20 RAIL, 25 ROAD** | **44 original strata; 229 RAIL, 462 ROAD** |
| Same original publication, region, decade, area class, ≥5 each source arm | **0 original strata, 0 RAIL/ROAD** | **11 original strata; 91 RAIL, 200 ROAD** |

**Meaning:** Much of the original Czech rail-vs-road signal from matched *regions* may still be comparing vegetation surveyed by completely different publications, observers and floristic missions. The stricter original bibliography matching collapses the strict original HABITAT subset to **only 20 RAIL and 25 ROAD plots** before checking whether it has enough independent *wind-strategy taxa* for a real biological guild model. Historical data **can be used for macro descriptive studies**, but this original cross-paper source selection cannot identify railway adaptation or physical seed delivery. Broad-text candidates provide greater observed sample but may include railway-adjacent—not genuinely on-track—vegetation, requiring source validation.

Next method gate: combine original study title × real source year × original phytoregion × area class **AND original observed plant diaspore/lifeform on the SAME plots**, then report if actual Sádlo wind vs nonwind groups are even estimable, instead of claiming an unsupported causal result. A registered same-study source analysis is implemented in [v12.3 script](../scripts/analyze_CzechVeg_rail_road_same_study_trait_v123.py) and separate reproducible source workflow; do not claim its biological result until real author-source workflow is complete.

**Any final paper must cite original CzechVeg and Pladias teams and consult Pladias data management for substantial reuse. No observed train seed carriage, maintenance/mowing dates, actual gusts, plant biomass, flowering plasticity or genetic adaptation were measured.** Original tetsu Draft PR remains unmerged, main unchanged.
