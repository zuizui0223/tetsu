# tetsu v13.1–v13.3 — genetically traced feral oilseed rape recruitment along Swiss freight railway and port interfaces, with full source-unit and 8-year monitoring audit

## Executive research synthesis (source cutoff 2026-11-14)

**Observed with independent original government + peer-reviewed author field sources:** (1) Original 2011–12 Swiss railway-area survey found glyphosate-tolerant genetically engineered *Brassica napus* feral plants in **4 of 79 surveyed railway areas**, heavily concentrated at rail stations and port transfer sites; (2) original 2013 Basel freight station and Rhine port follow-up molecular analyses identified **GM plant genotypes, additional glufosinate events at the active port, and GT73 DNA in seeds collected from non-GM maternal plants**; (3) original Basel cantonal monitoring of two localities from **2013–2020** observed different dynamics in the discontinued-cargo freight station and the actively importing port under ongoing control. **All are genuine, publicly verified original authors' field observations reanalysed by tetsu, not original tetsu data.**

**Crucial noncausal limit:** Original studies did not catch and genotype the *same actual cargo grain* and then track its seed-to-plant development to a mature adult and next-generation viable seedling. The original Rhine port actively transfers **SHIP cargo into railway cars or trucks** and wheat can be contaminated with oilseed rape; not every genetically traced feral plant can be attributed uniquely to a train. Different original sampling definitions, non-random source site selection and herbicide removal prevent estimation of a causal freight-cargo pressure effect on GM relative fitness or spontaneous railway adaptation.

## Original-source provenance / reproducible author data

- **Schoenenberger N. & D'Andrea L. (2012).** *Surveying the occurrence of subspontaneous glyphosate-tolerant genetically engineered Brassica napus L. along Swiss railways*, *Environmental Sciences Europe* 24:23, [DOI 10.1186/2190-4715-24-23](https://doi.org/10.1186/2190-4715-24-23). Original peer-reviewed **8-page CC BY 2.0 PDF** actually obtained and 59 Table1 source rows natively parsed, file SHA256 **b514908d43b99ac59b6cd557ab471a283d22434089c557d64e982df4c0cc6c8a**, 386,075 bytes. **[v13.2 exact primary source 79 RAIL-area field audit SUCCESS #43592628808](https://github.com/zuizui0223/tetsu/actions/runs/43592628808)**, source [original-PDF machine analysis](../scripts/analyze_Swiss_2012_full_railway_GM_sites_v132.py), [tests](../tests/test_Swiss_2012_full_railway_GM_sites_v132.py), [workflow](../.github/workflows/Swiss-original-2012-79-railway-GM-ferals-v132.yml).
- **Schulze J., Frauenknecht T., Brodmann P. & Bagutti C. (2014).** *Unexpected Diversity of Feral Genetically Modified Oilseed Rape Despite a Cultivation and Import Ban in Switzerland*, *PLOS ONE* 9:e114477, original [DOI 10.1371/journal.pone.0114477](https://doi.org/10.1371/journal.pone.0114477). Original full **2013 field survey from the same freight station and Rhine port** with 7 sampled dates port and 6 dates station, PCR-confirmed GT73 event, glufosinate bar marker and source maternal seed-pod DNA. Original publisher live 206,136-byte HTML downloaded with SHA256 **771747599192d1bfd2ff64e9f580ea0ebedd9c0f5b16ec93723cfae498af6155**, actual author numerical sentences independently validated. **[v13.3 original 2013 genotypes/seed pollen field source SUCCESS #43593118581](https://github.com/zuizui0223/tetsu/actions/runs/43593118581)**, [real author source site numeric transcription](../data/Swiss_Schulze2014_PLOS_2013_rail_port_GM_OSR_leaf_samples_v133.csv), [source script](../scripts/analyze_Swiss_2013_original_PLOS_GM_OSR_pollen_v133.py), [tests](../tests/test_Swiss_2013_original_PLOS_GM_OSR_pollen_v133.py), [workflow](../.github/workflows/Swiss-2013-original-PLOS-GM-pollen-seeds-v133.yml).
- **Bagutti C. & Hamburger D. (2020-12-17).** *Gentechnisch veränderter Raps: Kantonale Überwachung 2014–2020*, Kantonales Laboratorium Basel-Stadt and Swiss Federal Office for the Environment (BAFU), [official original government PDF](https://media.bs.ch/original_file/9ff1074caebc95c112c02d023afc713a69e28f42/2020-gvraps.pdf). Original public-source 969,514-byte PDF SHA256 **a74b7ea5013962a6af64acc56454a1bec6949b7f53439a18afcca15f2ac68174**. Real printed original govt Tables1–2 **2 sites × 8 years (2013–2020)** original annual estimated plant population and separate genetically tested plant samples; source native PDF text checked against both 8-year sequences, original official whole study header **7,908 plant tests, 459 GM-positive tests**. **[v13.1 true source governmental 8-year panel SUCCESS #43592152232](https://github.com/zuizui0223/tetsu/actions/runs/43592152232)**, [source table transcription](../data/Swiss_Basel_CantonalLab_2020_GM_rapeseed_station_port_2013_2020_v131.csv), [source script](../scripts/analyze_Swiss_Basel_rail_freight_port_GM_rapeseed_2013_2020_v131.py), [tests](../tests/test_Swiss_Basel_rail_freight_port_GM_rapeseed_2013_2020_v131.py), [workflow](../.github/workflows/Swiss-Basel-original-government-feral-OSR-v131.yml).

Each observation comes from original peer-reviewed author/government sources, not generated RAIL data or fictitious source original plant samples. Original government 2020 and PLOS 2014 site cohorts overlap; DO NOT count overlapping site/years/plant lineages as independent country replications.

## 1 — Original Swiss national railway survey (2011–2012), **79 surveyed RAIL areas / 58 feral sites / 4 genetic marker hotspots**

| Actual original field biological denominator | Verified original source estimate |
|---|---:|
| Original **rail-related areas surveyed** | **79** |
| Original areas where feral *B. napus* was observed | **58 / 79 (73.4%)** |
| Original rail-related areas **without** recorded feral rape, absent from published species-genetics table | **21** |
| Original Table1 rows with sampled feral plants and genotypes | **58 actual railway areas PLUS 1 off-rail field negative control** |
| Total original paper 2011–12 plants in Table1 **including field negative control** | **2,403** |
| Actual feral plant genetic marker tests at railway areas, **excluding 5 negative control field plants** | **2,398** |
| Confirmed original Swiss CP4 EPSPS-positive plant samples in actual railway areas | **50 / 2,398 = 2.085%** |
| Railway areas with at least 1 confirmed CP4-positive plant | **4 / 79 = 5.06%** of all surveyed RAIL areas; **4/58 = 6.90%** of areas with feral plants |
| Mean of **58** source railway-area within-sample CP4 marker fractions (one equal-weight area each) | **4.735%** |
| CP4-positive plant source samples occurring in the **three largest positive rail-port sites** | **49 / 50 (98%)** |

Original author 2012 50 site-positive plant specimens are **not 50 distinct railway arrivals**, may derive from imported contaminated crop seed, local previously spilt seed bank, or locally reproduced seed generations. The 79 railway areas were not randomly selected from the country's railway network. Original source larger feral populations were subsampled by sampling at regular spatial distances, which affects the sample-weighted denominator.

### All FOUR actual primary-source hotspots confirmed

| Original surveyed site | Original source plant samples | CP4-positive plant samples | Site-specific positive fraction |
|---|---:|---:|---:|
| **Lugano railway station (Ticino)** | 23 | **21** | **91.30%** |
| **Kleinhüningen port unloading structure (Basel)** | 18 | **16** | **88.89%** |
| **Basel St. Johann railway station** | 13 | **12** | **92.31%** |
| **Muttenz Auhafen port-side railway area** | 47 | **1** | **2.13%** |

Original 2012 source CH1903 map coordinates confirmed by PDF source and converted via standard EPSG:21781→EPSG:4326 for display only (source east/north exact in source JSON). Those are **two railway stations and two port–rail handling sites**, not proof that all GM feral plants entered via train rather than ship transfer.

**Site-size-only diagnostic (noncausal):** Under hypothetical random uniform assignment of the **50 GM-positive labels across the 2,398 actual plant tests**, preserving every one of the 58 original site sample sizes, 25,000 conditional randomizations yielded a mean **24.09 source sites** with any positive marker, compared with the actual **4**. No simulated null sample produced ≤4; plus-one Monte Carlo tail estimate **1/25,001 ≈ 0.000040**. This is **NOT an experimental p-value for rail transport or herbicide selection**. It only shows the spatial clustering of confirmed GM positives is far beyond what an artificial homogeneous per-sampled-individual null would predict; true seed lots, cargo logistics, site-specific glyphosate exposure and historical propagule inputs invalidate uniform biological exchangeability.

**Key mathematical point:** the source national individual-sample-weighted CP4 fraction **2.085%** differs from the equal-weight mean of 58 railway-area CP4 fractions **4.735%** (factor **2.27**) even with the exact same observational samples. This is another REAL weight-estimand distinction, but uses **one plant genotype phenotype** and must NOT be pooled with Slovenian **wind dispersal archetype × BS/BT preference** numerical effects.

## 2 — Original 2013 Basel 2-site study shows **local reproductive transgene transfer to seeds**, not inferred seed viability

Original *PLOS ONE* 2014 authors directly observed real feral plants and performed event-specific genotyping from real plant leaves and maternal pods.

| 2013 original PLOS site / plant sample grain | St. Johann freight RAIL station | Kleinhüningen Rhine PORT transfer |
|---|---:|---:|
| **B. napus LEAF samples** used for marker genotype analysis | **279** | **2,291** |
| Event-specific **GT73** positive OSR leaf samples | **64 / 279 = 22.94%** | **106 / 2291 = 4.63%** |
| Additional source glufosinate resistance **bar** marker positive OSR leaf samples, reporting known | not explicitly source counted | **22** |
| Recognized bar-related source event-specific sample counts | not reported | **MS8 4, RF3 13, MS8+RF3 3 (20 total)** |

**Pollen-related direct genetic transmission:** Among original Rhine-port samples collected from **48 NON-GM maternal oilseed rape plants with seed pods**, original source found transgenic **GT73 DNA in TWO maternal seed-pod sample pools** (2/48 = **4.17% of MATERNAL SAMPLE UNITS**, **NOT a per-seed fertilization rate**, not viability and not seedlings). One maternal sample was collected near GT73 plants, the second approx **200 m** from the nearest detected GT73 plant. This is **direct molecular evidence of genetic material in seeds from non-GM plants**, consistent with successful pollen movement/outcrossing; no offspring germination or seedling survival assay from those specific seeds was published.

- Source original overall **2,787 total collected plants of OSR and relatives** differs from 279+2291=**2,570** genotype-sampled *B. napus* LEAF specimens. Residual 217 cannot be treated as another 217 independently sampled *B. napus* GT73 tests; it includes relatives/other plant sampling.
- **CRITICAL overlapping 2013 sampling definition contradiction**: original Swiss Cantonal **2020 retrospective** panel counts original 2013 genetically-tested *B. napus* **St. Johann 100**, **Rhine port 2274**, whereas *PLOS ONE* 2014 original 2013 source **279** St Johann OSR leaves and **2291** Rhine port. These are **not the same cohort or denominator**, despite overlapping field sites and calendar year. The difference may reflect site area/sampling definition/study inclusion and publication methods, not a fabricated 100→279 population increase. They are deliberately preserved in separate source data models and never added for site-year pseudoreplication.
- 2013 PLOS original recorded glyphosate spraying in port on **April25, July25, October8** and in station on **April29**; true dates do not establish randomized comparable site herbicide exposure or exposure-dependent selection coefficients. Port source shipments of seed goods (mainly wheat) were continuing, sometimes with contaminated rapeseed, whereas St. Johann had stopped seed-cargo handling around **2009**.
- Source paper genetic event classes do not prove double herbicide adaptation in situ vs source imported stacked cultivars without matched original seed-cargo genotypes.

## 3 — Official Swiss 2013–2020 monitoring (2 high-risk site series, not random macro population)

The original 2020 Canton government source reports one original annual estimated plant abundance \`Nhat\` and separate plant DNA-testing sample count per site-year. Repeated controls and changes in seasonal sampling cannot identify survival of the same wild plants.

| Original monitoring year | St. Johann station: estimated feral *B. napus* | St. Johann plants genotyped | Rhine port: estimated feral *B. napus* | Rhine port plants genotyped |
|---|---:|---:|---:|---:|
| **2013** | 639 | 100 | 5,438 | 2,274 |
| 2014 | 84 | 49 | 9,247 | 1,333 |
| 2015 | 7 | 7 | 6,756 | 1,260 |
| 2016 | 23 | 23 | 5,461 | 1,216 |
| 2017 | 0 | 0 | 1,038 | 417 |
| 2018 | 1 | 1 | 459 | 338 |
| 2019 | 16 | 16 | 3,259 | 528 |
| **2020** | 2 | 2 | 752 | 344 |

These author original separately annual estimates vs DNA tested sample counts EXACTLY match full government source pdf; 16 annual population cells +16 tested cells verified, source original total **7,908 genotyping tests**, **459 GM marker-positive tests (5.80% of genetic tests, not a population-wide prevalence estimate)**. Since 2013 monthly and 2014+ twice annually, **raw year-to-year changes are not uniform surveillance effort**.

- St. Johann station was previously a seed-cargo freight transfer point, but source government reports **no handling germinable cargo there since the 2009 closure**. Original report records **no confirmed GM OSR from 2016 onward** at the station, even though non-GM feral plants were still occasionally observed (e.g. source 2019 estimated 16 OSR plants). This is **NOT proof no crop seeds continued to arrive by other channels**.
- The Rhine port continues seed goods shipments and has GM OSR and multiple event types over monitoring years. The original government report attributes possible repeated introductions to contaminating cargo and notes that glyphosate control may favor herbicide resistant lines **conditional on population survival and new arrivals**. Actual 2013–2020 population estimates fell at both source sites, and port spikes 2014 and 2019 should prevent a simplistic constant exponential selection trend.
- The source original 459 positives over 7908 tests cannot be decomposed by source site/year without independent exact original government source Fig1/Fig2 genotype counts. Do NOT invent per-year GT73/GM counts from annual population or total source positive fraction.
- **Two sites, both managed, different source types, different cargo/seed rain, source repeated observations and nonrandom sampling:** station-vs-port comparison is a strong **mechanistic NATURAL HISTORY CONTRAST**, **NOT causal difference-in-differences** or independently replicated population-level freight intervention experiment.

## Connection to multi-country RAIL seed source / established plant literature

| Original source study (distinct original field grain) | Observation in biological process | Supports | Cannot estimate |
|---|---|---|---|
| **Canada Gangadharan et al. (2017)**, original RAIL track cargo kernels 19 sites / 1,979 repeated traps / 134km, published original [10.1111/acv.12336](https://doi.org/10.1111/acv.12336) | **Actual moving TRAIN cargo grain physical deposition**: reported mean warm rates 1.12/.87/2.75 g/m²/day in 2013/14/15; source original annual approx 110 tonne whole-route model extrapolation, canola 25% of spilled grain mass | actual **cargo physical propagule input** | viability of collected seeds, same-seed to plant establishment, wild noncrop anemochory |
| **Switzerland Schoenenberger & D'Andrea (2012)**, 79 original railway-area survey | **Actual adult feral plant presence and GM protein phenotype** with marked import-genotype at rail/port gateways | crop escape and site concentration, some herbicide survivors | what specific train seed arrived or seeded sampled individual |
| **Switzerland Schulze et al. (2014)**, two Swiss source sites 2013 | **GM genotype-confirmed feral adults and GT73 DNA in seeds from non-GM maternal plants** | reproduction/outcrossing of trait at a port railway interface | daughter seedling viability and the exact train source |
| **Switzerland Basel Canton (2020)**, annual 2013–2020 same two sites | estimated feral plant population trajectories and GM molecular tests under managed continuing-import port vs ceased-freight station | dynamic field persistence/decline under cargo/management context | causal freight stoppage, independent seed banks/management selection coefficients |
| **Austria Pascher et al. (2016)** 60 source transport sites×2 years | sampled feral *B. napus* adult presence at ordinary stations/road sections, 16/20 stations vs12/22 road sites positive; source stratified exact p≈0.0825 | second independent *COUNTRY* adult feral crop presence | cargo origin or *train* freight flux, direct viability and plant fitness |
| **Slovenia Šabić et al. (2026)** 89 original BS/BT railway plant field sites | population-/community-level functional **plant occupation** and source wind archetypes | descriptive source plant-trait filtering | train seed input, germination, mowing and adaptation |
| **CzechVeg original 2026** source 117739 vegetation plots, matched true Pladias community class cohort 105 rail-text/167 road-text, 532 taxa | source species-pool occupation, original strict-wind trait differences; both group effects and study-cluster gap intervals include zero | source-sampling sensitivity and limits to generalizable rail trait filtering | genetically tracing train-spilled seeds to distinct established wild plants |

**Conclusion:** These are many truly empirical observations at different links of the ecological chain, but not one single-site, same-year, same-original-species dataset measuring all links. Even Swiss genetically traced GM feral plants cannot uniquely distinguish **ship-to-train cargo**, repeated imports, persistent soil seed-bank, and pollen-mediated local reproduction. The next causal mechanism test requires matched cargo grain genotype/seed germination/seedling genotyping at the same rail route/time/station, management schedules, and multiple independent countries or networks. Mowing/strong wind adaptation of native railway plants remains biologically unknown.

**Recommended accurate manuscript framing:** *Anthropogenic propagule pressure and habitat management jointly constrain feral crop persistence at transport gateways: limits and opportunities for railway macroecology*. Empirical evidence supports the *existence* of each separate field mechanism, but a universal **railway-specific causal selection coefficient** or global meta-analysis of linked arrival→establishment is **not** identified.

Original copyright: Swiss 2012 CC BY 2.0, PLOS 2014 CC BY 4.0, Swiss government 2020 original public source, original Canadian and Austrian papers retained under respective publisher author citation. Pladias original substantial reuse would require its Governing Board consultation. Keep tetsu original [PR #1](https://github.com/zuizui0223/tetsu/pull/1) DRAFT, main original tree content unchanged.
