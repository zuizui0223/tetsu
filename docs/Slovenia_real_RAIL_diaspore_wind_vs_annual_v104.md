# tetsu v10.4: Real railway microhabitat annual-life-history enrichment does NOT imply wind-dispersed propagule enrichment

**Independent original published railway field data reanalysis**, not new raw tetsu observations and not inferred from GBIF occurrence counts.

Study original authors: Šabić A, Grašič M & Rovanšek A (2026), *Vegetation patterns, functional traits and invasion dynamics across railway microhabitats in Slovenia*, [DOI 10.1007/s44473-026-00155-z](https://doi.org/10.1007/s44473-026-00155-z), CC BY 4.0 Supplementary file 1. Provenance source checksum SHA256 26be7ca9513dd54136d4ec499d8413c1cb288045f87fe3cf8814731f2edac99f, original XLSX 401279 bytes.

**Actual verified source analysis:** [v10.2 authenticated supplementary original source Action SUCCESS #43538799032](https://github.com/zuizui0223/tetsu/actions/runs/43538799032), [v10.3 89 railway sites' annual and therophyte source trait analysis SUCCESS #43539536809](https://github.com/zuizui0223/tetsu/actions/runs/43539536809), **[v10.4 wind-diaspore/annual same-location source result SUCCESS #43540382128](https://github.com/zuizui0223/tetsu/actions/runs/43540382128)**, [actual v10.4 all source sensitivity JSON artifact #14091901171](https://github.com/zuizui0223/tetsu/actions/runs/43540382128/artifacts/14091901171), [reproducible source code](../scripts/analyze_Slovenia_rail_wind_diaspore_vs_annual_v104.py), [5 integrity/algebra tests](../tests/test_Slovenia_rail_wind_diaspore_vs_annual_v104.py).

## Question: biological dispersal mode and environmental life history may differ

Original 2026 Slovenian railway surveys include **92 actual railway field locations, 482 habitat-zone records, 6331 species-by-zone rows (6295 unique EPPO-by-zone incidence), and 312 species with a Diaspore.strategy / lifespan / life-form source record**. Fifteen actual railway sections; source sites were chosen near vegetated tracks rather than randomly within a country. Original dates 2022–2023 are after a 2021 herbicide ban, with no treatment-before baseline.

Real same-location same-track-side matched railway microhabitat pairs: **173 BS (ballast shoulder) vs BT (ballast toe) side pairs**, averaged within the **89 distinct field locations**, then summarized across **11 of the 15 railway sections with >=3 source field locations each**. BS and BT physical transect area is NOT known to be strictly equal, so study species-level composition shares among trait-classified taxa, NOT richness per identical 100m². Neither transport volume nor true local train slipstream or mowing pressure is measured.

The original author's **wind-dispersal diaspore guild** is exactly source strategy values **Epilobium** and **Phragmites**. Our separate OPTIONAL broad sensitivity includes **Lycopodium** (wind-borne spores possible, not the strict original strategy). The dataset's Allium (autochory), Bidens (epizoochory) and Cornus (endozoochory) source codes are kept as additional nonwind comparison strategies; all missing Diaspore.strategy entries have UNKNOWN source value and are excluded from denominators, NOT treated as nonwind.

### Actual joint functional source contrasts

| Same 89 actual railway locations, BS minus BT | Mean compositional share difference | BS greater locations | BT greater locations | ≥3-sampled-section sign BS/BT |
|---|---:|---:|---:|---:|
| Strict **annual** lifespan species (v10.3) | **+0.054703 (+5.47 pp)** | 60/89 | 29/89 | **8/11 BS** |
| Strict **Therophyte** life form (v10.3) | **+0.058871 (+5.89 pp)** | 62/89 | 27/89 | **8/11 BS** |
| **Strict wind-dispersed diaspore strategy Epilobium + Phragmites** (v10.4) | **−0.050740 (−5.07 pp)** | 31/89 | **58/89** | **7/11 BT** |
| Broad wind+Lycopodium spores sensitivity (v10.4) | **−0.025074 (−2.51 pp)** | 37/89 | 52/89 | 6/11 BT |
| Allium/source autochory | +0.017541 | 48/89 | 40/89 | 7/11 BS |
| Bidens/source epizoochory | +0.000047 | 47/89 | 35/89 | 6/11 BT |
| Cornus/source endozoochory | +0.012113 | 44/89 | 37/89 | 8/11 BS |
| Original vegetation cover source units (v10.3; NOT necessarily %) | **−35.179775 units** | 7/89 | **81/89** | **11/11 BT** |

**Direct same-location orthogonal trait contrasts:** Of all **89** original site-averaged BS−BT results with BOTH annual and wind-diaspore source fractions defined, **43** have annual species fraction BS higher but wind-diaspore fraction BS lower; **17** have both higher; **15** both lower; **14** annual lower but wind higher. This descriptive distinction in source life-history strategies is compatible with two different functional constraints rather than "wind-dispersed plants necessarily preferentially occur closer to tracks."

### Robustness and unfavorable regional cases

- Strict wind BS−BT by survey year: **2022 (53 locations) −0.040161**, **2023 (36) −0.066314**. Thus a lower wind-strategy share at BS appears in both observed years, but no repeated same-site before/after time series.
- By track type: **SINGLE (65 locations) −0.036203**, **DOUBLE (24) −0.090111**; DOUBLE includes only 2 sections with >=3 locations, so section-wide inference is especially limited there.
- Optional broad wind + Lycopodium source group: BS−BT still negative **−0.025074** overall, however among single-track locations effect only **−0.007757**, suggesting source taxonomy coding matters.
- At segment level 11 eligible source sections, wind-BS higher in **4** and wind-BT higher in **7**. Independent source oppositely signed sections must not be suppressed: source segment IDs 20, 34, 64, 80 are BS wind-positive; other sections with sufficient source site observations BT wind-positive.
- Exact source EPPO taxonomy issue LOTHE (source traits/actual 4 plant records) vs DOYPH (species codebook) is retained and annotated, not silently repaired; 36 repeated same-zone EPPO incidence records deduplicated, all 482 published zone richness totals reconciled. One cover value 200 means source cover cannot automatically be described as %; excluding its location leaves **88** for cover sensitivity.

## Biogeographic/mechanistic interpretation

**SUPPORTED observational biological finding**: In these actual Slovenian RAIL field plots, near-track ballast shoulder has less measured source vegetation cover and greater relative frequency of short-lived annual / Therophyte plants, whereas **wind-dispersal-strategy plants are not enriched at that same near-track shoulder**. A negative BS−BT wind-diaspore share is present across the aggregate and most but NOT all rail sections, suggesting that wind dispersal potential (source delivery-related trait) and short-lived growth form (establishment-related trait) are NOT an identical macro filter.

**NO causal railway mechanism established**:
1. A plant coded "Epilobium/Phragmites" wind-dispersal may be moved by **natural ambient winds** rather than train-generated slipstream. Source data measure neither moving seed trajectories nor actual train gust speed or freight cargo release. A lower wind trait fraction in BS does NOT disprove physical seed carriage by trains.
2. Short-lived plant enrichment may be due to ballast substrate, temperature, drainage, disturbance from maintenance, sampling width or development history, rather than measured mechanical wind/mowing selection. No true manipulation or site-year treatment repeated data are present.
3. Actual one-country 15 rail segments and 89 clustered locations are NOT 15 independent countries nor a global empirical law. Comparisons need additional independent rail ecoregions and equivalent equal-area road controls.
4. The true MIREN worldwide ROAD plot dataset supplies a separate standardized roadside Alien-richness enrichment/turnover baseline v9.9–v10.1, but **do NOT subtract percentages from Slovenian RAIL microhabitats** (different design, sampled widths, geography, nativeness definitions, protocols).
5. No common-garden/genomics or within-species repeated true phenotype/fitness measures: these traits are source taxon-level and do NOT prove adaptive genetic change or plasticity.

The next identification gate is matched equally sized RAIL vs ROAD botanically surveyed site-zone pairs across multiple independent ecoregions with original traffic/seed traps and measured mowing/gust exposure; stratify freight delivery, ambient wind and microhabitat disturbance separately. More GBIF state-level correlations are not a substitute for direct propagule/management measurement.
