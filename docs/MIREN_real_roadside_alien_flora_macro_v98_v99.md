# tetsu v9.8–v9.9 — Independent standardized road-verge vegetation survey reveals near-universal roadside alien-species enrichment (not yet a RAIL effect)

**Original dataset and attribution:** Seipel, T., Haider, S. and MIREN consortium, [MIREN survey of plant species in mountains, v3](https://doi.org/10.5281/zenodo.15749536), observed 2007–2023. Data paper: Haider, S. et al. (2026), *Global Ecology and Biogeography* 35:e70284, [doi:10.1111/geb.70284](https://doi.org/10.1111/geb.70284). All actual underlying plant data are the consortium's, NOT original field data collected by tetsu. The 2026 data paper asks users to notify the MIREN Steering Committee of planned substantial analyses and discuss data contributor involvement for manuscripts.

**Independently executed source-backed GitHub Actions:**
- [v9.8 exact original 17.9MB MIREN 2025 source CSV and checksum/schema audit SUCCESS #43212464976](https://github.com/zuizui0223/tetsu/actions/runs/43212464976); [original raw and schema artifact #13941483126](https://github.com/zuizui0223/tetsu/actions/runs/43212464976/artifacts/13941483126); [source reader](../scripts/audit_MIREN_road_vegetation_real_source_v98.py).
- **[v9.9 full 2007–23 24-region paired road-vs-interior plant richness analysis SUCCESS #43214312718](https://github.com/zuizui0223/tetsu/actions/runs/43214312718)**; [complete original-source analysis JSON #13943505434](https://github.com/zuizui0223/tetsu/actions/runs/43214312718/artifacts/13943505434); [reproduction script](../scripts/analyze_MIREN_standardized_road_edge_vegetation_v99.py); [six targeted source-identity/matching tests](../tests/test_MIREN_standardized_road_edge_vegetation_v99.py). Six targeted tests and latest full repository descriptive smoke checks pass.

## New independent biological data: a genuine sampled local plant assemblage rather than a citizen-science source-indexed presence

MIREN researchers directly sampled vascular plant species composition using standard 2m×50m (=100m²) vegetation plots at mountain ROAD sites, with up to 3 plots per site on a gradient:

- **Plot 1:** road-adjacent strip, aligned parallel to road.
- **Plot 2:** intermediate zone approximately **2–52m** from road.
- **Plot 3:** interior vegetation approximately **52–102m** from road.

All plots in a same Region/Road/Site/Year group were observed by the source field protocol in the SAME survey year; comparison does not depend on GBIF stateProvince text labels or arbitrary indexed-record rarefaction. Original taxon field is Accepted.Name.MIREN, and **Status** is Alien/Native/Unknown **according to source regional native-status classifications** (not necessarily known invasive in every locality).

**Exact source v3 MD5**: b640195cc44f252ccb6ddb28eb6fff13. CSV 17,899,241 bytes; **167,804 original species-occurrence rows** and **6,854 distinct name categories**, matching the published order of magnitude. Source contains **26 distinct region CODES**; paper describes **25 mountain regions** and we do not invent an equivalence between the code count and independent eco-geographical sampling units. This analysis has **24 codes with actual same-site/year Plot 1 and Plot 3 recorded**. The excluded region codes **AUV** and **HWI** have no matching original STRICT 1/3 site-year cohort; they are **not** negative biological results.

### Original quality control before ecological comparisons

- **7,223 repeated identical-taxon-plot-year records removed** in strictly coded Plot1/Plot2/Plot3 source. A real duplicated species observation in one plot-year counts as ONE biological species category, not additional richness.
- **Two region×species taxon nativity classifications conflict** in original source (CNN: *Tephroseris kirilowii*; ITA: *Patzkea paniculata*); both are conservatively reclassified "Unknown" rather than arbitrarily considered Alien or Native.
- **915 source records with modified Plot suffix** 1A/2A/3A (replacements after plot destruction) excluded from strict SAME-site original comparison; these replacements cannot automatically be treated as identical previous vegetation plots. These valid distinct replacement records remain in the raw retained source but are not mixed into main biological comparison.
- Observed Plot1+Plot3 same site-year pair must have data in BOTH plots; missing plot-year in the long species-table is never imputed as confirmed 0 species. Distinct plot-year source observations in strict plots: **6,419**.
- Early 2007/2008 protocol often concentrated on aliens rather than a complete native flora census; therefore primary outcome is **observed Alien richness**, while native fractions/total richness and 2007/08 temporal analyses are secondary/limited. In MIREN, disturbance intensity and forest cover are not consistently measured; their absence prevents direct inference on mowing or wind.

### Primary region-level source-balanced spatial result

**First surveyed year of each geographically unique mountain road-site** (Region + Road + Site; site revisits excluded) is the primary independence/quality sensitivity. **1,005 unique road-sites** in 24 code regions were sampled with both roadside Plot1 and interior Plot3 in at least one year.

- **23/24** source region codes have positive mean per-site-year Alien species richness (Plot1 - Plot3).
- Under **>=20 geographically unique matched ROAD sites per code** (originally defined nontrivial macro sampling threshold): **21/22 source codes positive**, one negative, zero exact ties. Across these 22 source codes, **median of code-average Alien richness difference = +2.155953 species per 100m² roadside plot**, mean code-average +3.129523.
- The independent geography remains *region code*, not 1,005 independent countries. Source sites are spatially structured and road choice was not randomized; do NOT interpret these sign counts as proof of direct roadside traffic/wind/herbicide causal mechanisms.

**Full site-year observational cohort** (including repeated years, as sensitivity only): **1,876 Plot1-vs-Plot3 paired road-site/year source comparisons** across **64 road identifiers nested in 24 code regions**, 1,005 unique road-sites. Source roadside minus interior Alien richness difference:

| Genuine same-site/year pair statistic | Real v9.9 finding |
|---|---:|
| Plot1 Alien species richness minus Plot3 richness (pooled mean) | **+4.086354 spp / 100m²** |
| Pooled median source pair richness difference | **+2 species / 100m²** |
| Paired site-years with more roadside Alien spp | **1,346 / 1,876** |
| Paired site-years with same Alien richness | 354 / 1,876 |
| Paired site-years with more interior Alien spp | 176 / 1,876 |
| Region codes positive among ALL paired site-years | **23/24** |
| Region codes positive among >=20 paired source site-years | **21/22** (code-equal median regional mean difference +2.232143) |
| Native spp Plot1 minus Plot3 pooled mean | −0.417377, **heterogeneous**, and older native surveys may be incomplete |

**Complete three-plot same-site/year check**: **1,561** source triplet observations with all 1, 2, 3 plot vegetation records. Alien richness median per 100m² Plot1 **5 species** > Plot2 **3** > Plot3 **1**. Mean signed paired Alien richness differences Plot1-Plot2 **+2.575913**, Plot2-Plot3 **+1.175529**, Plot1-Plot3 **+3.751441**. The median-distance gradient supports a real roadside-associated plant community pattern, not a railway traffic seed-carriage measurement.

**Mandatory temporal/protocol source sensitivity**, including unfavorable results:
- Original **first** and **last** measured source year of each road-site: source region signs **23/24** positive in BOTH. This avoids summing repeat five-year surveys as independent source sites.
- Only **2017 or later** actual survey plot pairs: **20/21** positive source region codes (>=10 matched site/year per code).
- Only **2022 or later**: **11/12** source region codes positive (>=10 pairs).
- **484** source road-sites have original Plot1-vs-Plot3 pairs in >=2 actual survey years. Across the **11** region codes with >=10 repeated road-sites, mean change in roadside Alien richness advantage from each site's earliest to latest survey year is **positive in 5 regions, negative in 6**. Median site-level earliest-to-latest ROAD-edge advantage change is **0**; **NO consistent temporal strengthening** observed. The static roadside alien enrichment is much more regionally general than the change through time.

## What this contributes to the actual railway-seed-dispersal project

The previous v7.4–v9.7 train-seed macro work relied heavily on GBIF *positive indexed species records*, in which original 2018–25 "new records" are dominated by iNaturalist observational sources and do NOT establish new invasions. Here an INDEPENDENT source—a standardized, field-surveyed mountain-road vegetation matrix with region-specific Alien status—produces a much stronger and more consistent **near-road vs interior Alien richness concentration** across geographic regions.

**This is an actual observational biological ROAD disturbance/filtering pattern** in repeatedly surveyed equal-area vegetation plots, and can anchor a future matched railway-verge vs road-verge comparison. It is NOT yet a newly discovered global railway-effect mechanism:
1. Plot differences may be caused by road-associated propagule delivery, roadside light, soil, substrate, habitat opening, human management, vehicle disturbance or existing plant source pools. **Mowing frequency and measured wind gusts are unavailable** in this original MIREN dataset.
2. Actual railroad verges, rail freight cars, operated railway traffic and live dispersed seed counts are absent. **Cannot compare RAIL/ROAD biological effect sizes without equivalently sampled railway plots.**
3. No common-garden, within-species growth/reproduction after cutting, genetic markers or plasticity experiment; the regional Alien richness gradient does NOT prove evolutionary adaptation.
4. The source survey flora gradient has already been a motivation of the MIREN initiative; it must be credited as a **reanalysis** of consortium measurements, not portrayed as uniquely discovered original plot data.

## Next decisive source acquisition

Recruit/locate *original equivalent railway vegetation relevés* with region-coded Alien status, replicated equal-area plots at verified 0/2–52/52–102m distances from operating tracks, **independent management (mowing/herbicide)** records and measured traffic/slipstream wind as separate predictors. Seed sources can be tested with freight car / railway-yard seed traps and verified cargo/route exposure. Use geographic region as independent macro-level unit and original paired plot-year as within-site contrast, not GBIF record-pair p-values.

**Known direct railway-data leads** (MUST NOT be pooled with MIREN as though equivalent):
- [Poland railway-tracks vascular plant dataset, Nobis & Nobis (2026), GBIF DOI 10.15468/qu2242](https://www.gbif.org/dataset/63289930-6581-43d7-a587-c980f4ae5bca): only ~45 source species across southeast Poland, coordinates refer to **10km ATPOL grid centers**, NOT original track plots or matched road/interior distances, source 2004–05. It documents surveyed railway flora, but is not a region-matched experimental control.
- [India Surat Metro Rail Phase I GBIF floral survey](https://www.gbif.org/dataset/b0c10ea1-bff9-4148-80e9-35b8e1fec6e2): source metadata ~222 plant species and strip/transect vegetation within/near project corridor; no assured equivalence to MIREN 100m² matched 1–3 plots, and actual operational rail seed dispersal is NOT measured.
- [USGS Canaan Valley historic railroad embankment wetland vegetation/environmental plot data, DOI 10.5066/P9KUATM7](https://doi.org/10.5066/P9KUATM7): independent measured hydrology/vegetation near abandoned/historic embankment, NOT active train wind/carried seeds.

**Study identification gate:** Do NOT call the Poland/Surat/USGS datasets a multi-country standardized RAIL sample or merge their species lists into MIREN, since design, habitat, survey effort and historical train operation differ. The independent roadside benchmark remains strong, but future rail verification is a separate required data-collection phase.
