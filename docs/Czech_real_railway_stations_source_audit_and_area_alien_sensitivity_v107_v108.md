# tetsu v10.7–v10.8: Original independent Czech railway stations – documented non-significant station-area proxy result

**External original field study:** Kutlvašr J., Turková S., Hejda M., Vojík M., Kadlecová M., Berchová Bímová K., Pyšek P. & Pergl J. (2024), *Railways as a source of alien plants*, *Preslia* 96: 247–266, [DOI 10.23855/preslia.2024.247](https://doi.org/10.23855/preslia.2024.247). We use authors' original source XLSX appendices and **corrected Supplementary Table S2 erratum**, not any newly collected field observations from tetsu.

**Original-source actual successful GitHub Actions:** [v10.7 original publisher S1/S2 schema and 39-station presence audit #43561722534](https://github.com/zuizui0223/tetsu/actions/runs/43561722534), [v10.8 actual station area/alien proportion analysis #43562007783](https://github.com/zuizui0223/tetsu/actions/runs/43562007783). Source code [v10.7](../scripts/audit_Czech_2024_real_railway_appendices_v107.py), [v10.8](../scripts/analyze_Czech_2024_rail_station_area_aliens_v108.py) and tests [v10.7](../tests/test_Czech_2024_real_railway_appendices_v107.py), [v10.8](../tests/test_Czech_2024_rail_station_area_aliens_v108.py). All read-only publisher XLSX, no copyrighted raw publication tables redistributed.

## Direct provenance and source interpretation

- Original S1 plant list/rec-vs-hist presence matrix publisher source **79,261 bytes**, SHA256 **df6bc0b949e7f6a5cb57755dc8f9c86187223f66ddfbef959cf0e03ef7c46b06**; 1,060 rows including headers, 74 columns, 1,053 unique taxon labels across all recent AND historical source records, **NOT** 1,053 recent species.
- Original CORRECTED S2 station characteristics publisher source **13,533 bytes**, SHA256 **00c7cfdf5047854ba99d47dc8ec8c55b1b7df72731ce647c22b286ddef24a5a2**; original 39 railway station identifiers, all matching S1, with station/railway stretch, total taxa, alien/invasive/native counts and 3 distinct zone areas.
- Original S1 has 72 station/period columns: **39 recent and 33 historical**, and exactly **6,941 'x' markers in REC columns**, historical letters a,b,c,d (alone/combined) record historical period categories. All 39 S1 station codes exactly match the corrected S2 station codes. Historical Pladias records have uneven sampling and were grouped by the ORIGINAL authors for description; **no pseudo-absence or temporal turnover inference** from S1 hist columns.
- Recent S1 presence taxon count vs S2 total taxon count agrees exactly at **36/39 stations** (3 mismatches DV/SV/TI, source differences +4/+1/+1). Exact original appendices are preserved; no silent code changes or row/label rewriting. Recent alien source label (original explicit archaeophyte/neophyte labels) matches S2 counts at **26/39** stations; there are small additional source-specific discrepancies and not all redlist categories are alien labels. Use **S2 corrected published counts** as the primary station area comparison and carry S1-S2 disagreements as sensitivity.
- Original corrected S2 station **PRM** has original missing \u201cNA\u201d zone area; no unknown area set to zero. Excluded from numeric area model: **38 complete-area stations** across **all 5 original railway stretches**. 35/38 source plant counts exactly match S1, 25/38 source alien counts match simple explicit archaeophyte/neophyte tag parse. Sensitivity excludes all disagreement stations, retaining 23.

**Not the same design as Slovenia:** Czech data have 39 large station complexes with three zones (track/embankment, buildings, wider surroundings) of unequal area and a station-level source matrix; Slovenia original v10.3–v10.6 use 89 same-site ballast shoulder BS versus ballast toe BT matched zone-side differences and original wind/lifespan trait metadata. Do not pool sites, percentages, sampling areas or ecological effects across these two papers as if equivalent randomized RAIL-vs-ROAD quadrats.

## New v10.8 original 38-station ecological association: positive estimates, no clear evidence

Independent 2020–2022 source observation, NOT an identified physical railway delivery effect. Model source station alien fraction and invasive share of aliens as separate outcomes; predictors track-and-embankment share of total station area, adjusted for \u006clog(total station area), station building-area share and five original source stretch fixed effects. 9,999 within-stretch **approximate residual permutations** (two-sided) and 4,000 within-stretch station bootstrap replicates (exploratory).

| Outcome | Partial slope per +10 percentage points of station track-area share | Source-stretch-stratified 95% bootstrap interval | Approximate source-stretch-blocked permutation p |
|---|---:|---:|---:|
| Original S2 aliens / total recorded station taxa | **+0.0088911 (+0.89 pp)** | **[−0.0166077,+0.028394]** | **0.3224** |
| Original S2 invasive aliens / all alien taxa | **+0.0045737 (+0.46 pp)** | **[−0.0140756,+0.0203179]** | **0.5630** |

Both intervals include 0. **Do not claim the source station area composition is a significant or robust positive predictor.** Leave-one-of-five original railway stretches out for the *total alien fraction* beta remains numerically positive in **5/5**, from +0.0037939 to +0.0146804 per 10 pp, but it lacks support from the 5-stretch sample's broader uncertainty and the within-stretch conditional permutation.

The same track-area share coefficient among only 23 source stations with completely matching original S1/S2 species/alien counts drops to **+0.0039955 (+0.40 pp)** per +10pp (no inference significance claimed). Thus the small area proxy association is neither mechanistically interpretable nor robust to all source-data quality restrictions.

Original recent complete-area source station areas range **2.4585–87.6885 ha**; track-area fraction **0.2946–0.8522**; alien fraction **0.2895–0.5630**. The published authors' original independent result **track zones richer in taxa/aliens than surroundings** remains original work and is NOT replaced by our new station-area share analysis; zone richness and station alien *proportions* are different estimands. Original paper's more detailed abundance/management and phyto-geography interpretations must be attributed to their authors.

## Ecological inference – controlled separation of scales

1. **Slovenia observational fine-scale RAIL trait axis:** annual/strict Therophyte higher toward ballast shoulder; strict wind-diaspore strategy lower, persists when conditioning on lifespan category (v10.6). This is one-country source functional filtering, NOT proven wind or mowing cause.
2. **Czech independent station-area axis:** station-relative track-area size alone does NOT consistently predict alien composition fraction once total area/other zonal area/rail stretch are conditioned (v10.8 uncertain small estimate). This is not a refutation of the Slovenian microhabitat pattern.
3. **MIREN true worldwide ROAD field benchmark:** road edge enriched in alien richness, with high turnover and no widespread trend towards increasing edge advantage over repeated observation; **not** a same-ecoregion RAIL-vs-ROAD cause estimate.
4. **Still not proven:** direct loaded train propagules, real train aerodynamic seed displacement, mowing schedule/recovery, within-species plasticity, heritable adaptation or rail-mediated long-distance plant arrival.

**Next genuine multi-region identification gate:** Obtain comparable source *railway zone-level* raw species incidence matrices and survey area, directly paired road habitat quadrats, and independent freight/rail track use and actual mowing or gust exposure. Replicate the annual vs wind (strict categories + lifespan-adjusted component) on multiple railway sections across several countries, then hierarchical section/country effect; consider seed-traps/seed contamination data to separate delivery from establishment. Without direct movement/management or suitable natural experiments, call any findings **observational functional sorting**, not local adaptation or train-borne dispersal.

## Attribution and integrity

DO NOT claim these papers are tetsu original raw field observations. Original field paper and DOIs above must be cited; the reanalysis is reproducible at the level of published raw data obtained from the publishers and durable GitHub Actions audit receipts. No PR merge; original main remains untouched.
