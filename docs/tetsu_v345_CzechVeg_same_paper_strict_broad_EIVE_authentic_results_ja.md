# v345: Same-publication and same-block source label audit (authentic source, 2026-11-24)

## Result
The previous v343 CzechVeg railway text vs road text EIVE indicator difference was +0.049381 in 32 original papers. Strict publisher HABITAT-only labels in v344 retained 5 original papers and produced -0.052983. v345 re-fetches the authentic source and separates original-paper composition from more stringent field-label and actual original-plot selection.

| Original cohort and comparison | Source papers | S2 matched blocks | Study-equal rail minus road EIVE |
|---|---:|---:|---:|
| v343 broad label, all eligible original studies | 32 | 39 | +0.0493809 |
| v345 broad label restricted to the exact 5 strict-eligible source papers | 5 | 6 | -0.0085057 |
| v345 strict HABITAT label on exactly these 5 papers | 5 | 6 | -0.0529830 |
| v345 broad label restricted to 5 identical original S2 design blocks | 5 | 5 | -0.0388482 |
| v345 strict HABITAT label on identical 5 S2 design blocks | 5 | 5 | -0.0589434 |
| v345 identical original plot IDs (5 blocks, 17 rail and 16 road original plots) | 5 | 5 | -0.0557362 in BOTH label estimators |

**Finding:** Selecting only the original source papers that qualify for the strict HABITAT-only design, while retaining BROAD label classification, changes the reported thermal signal from +0.049381 to -0.008506. Thus the positive broadly defined Czech source effect is strongly sensitive to which original papers are represented, independently of the next label-sharpening step. Restricting to exactly the same five original papers but changing which original vegetation plots are included by HABITAT label changes the study-equal mean again from -0.008506 to -0.052983. These differences are **not** causal decompositions, because source plot selection and field habitat differ. The original field plots held identical in both calculations yield exactly equal computed EIVE temperature values (software source identity test).

Main EIVE source gate: >=3 independent indicator systems, 2014/source plot source taxa >=5, valid matched EIVE species >=3 and >=50% coverage, source original study × phytoregion × decade × plot area class, >=2 original plots in both RAIL/ROAD arms per S2 block. 70% coverage sensitivity gives broad same-five-paper +0.0321633 versus strict -0.0232831; on the exact common S2 blocks broad +0.0018207 vs strict -0.0292435. In the identical original plot-ID subset, both label modes yield -0.0150672. This sensitivity is inconclusive with just five papers.

## Authentic sources and validation
[Full authentic GitHub Actions run #45246412311 (SUCCESS)](https://github.com/zuizui0223/tetsu/actions/runs/45246412311): source 117,739 CzechVeg original vegetation plot headers and 2,510,340 original plant plot records re-downloaded and verified against official Zenodo 19231152; MD5 headers 0bc920462d1f9eac473e7b2a231e6782, species b296e5e19f03fb3bad0a8347c317f8a3. Independent EIVE1.0 source Zenodo 7534792 MD5 0f0ce19cd2a781eeed027d8b0acadbcb. Five synthetic original site/paper overlap and null-support invariant tests passed. No original plot species rows, geographic coordinates, or source publication identifiers exported.

v345 design was frozen following v343/v344 result inspection, therefore exploratory. Read protocol docs/tetsu_v345_CzechVeg_common_paper_label_source_gate_frozen_ja.md, code scripts/analyze_CzechVeg_common_paper_label_thermal_v345.py, test tests/test_CzechVeg_common_paper_label_thermal_v345.py.

## Scientific decision
Do not claim independently replicated, modern active-railway thermal selection from broad historical Czech RAIL/ROAD textual associations. The original data labels HABITAT and LOCALITY are not verified operational ballast. There are only five source studies meeting strict rail/road classification and six study-region-decade-area S2 blocks. Avoid treating them as independent climatic replication or randomized railway/nonrail nature controls. The European EIVE T preference is an ordinal plant preference indicator, not measured temperature in C.

A more rigorous global railway comparison still requires at least two geographically distinct railway and genuinely adjacent nonrail natural-habitat paired botanical field datasets with original plot species lists, matched area, survey effort/season, and independent trait calibration.

### External eligibility candidates
1. Dylewski, Mackowiak & Dyderski, first online 2024, Land Degradation & Development 2025 36:614-629, DOI 10.1002/ldr.5383. 25 abandoned rail sites and 25 adjacent semi-natural grassland reference sites. Authors expressly state the original data are available **upon reasonable request**, not as a freely downloadable source matrix. Abandoned tracks do not identify active train dispersal.
2. Everingham, Hemmings & Moles (2019), DOI 10.1111/aec.12719. The journal lists actual field species Appendix S4 (XLSX), but earlier official Wiley requests returned HTTP 403, and the original rail embankment vs fire-trail/artificial habitat comparator needs original plot-level pairing verification.
3. Winkler et al (2026), DOI 10.3390/land15040523. Five railway-agriculture habitat classes, four localities, 2021, 10 m² plots; article states complete raw data withheld, its published supplemental Table S1 contains taxon names but no full plot × species matrix.
4. Historical Tikka et al (2000), DOI 10.2307/1478915. 92 semi-natural grasslands and 90 railway/road verges studied in Finland, but railway and roads are mixed in an original aggregated design, and accessible exact source plot × species data not yet verified.

Both v343 broad source EIVE and v344/v345 strict-source counterexample remain as observational sensitivity examples, not additional independent climate regions. Main GitHub branch unchanged. Research Draft PR #54 remains unmerged.
