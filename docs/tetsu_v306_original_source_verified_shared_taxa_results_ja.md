# tetsu v306 — SOURCE VERIFIED: actual Slovenian freight and species sharing

2026-11-20 JST. Source-authentic Actions SUCCESS #44746851895:
https://github.com/zuizui0223/tetsu/actions/runs/44746851895
Canonical aggregated original-source artifact #14633721347:
https://github.com/zuizui0223/tetsu/actions/runs/44746851895/artifacts/14633721347

First-party author original 401279-byte species/trait XLSX SHA256 matched, original operator station PDF and 2020 Eurostat real two-direction network 9 segment/27 field site/5 source railway line matching succeeded. 89 original sites, 36 dependent between-segment dyads. Five source numerical tests and independent aggregate QA PASS. Three no-urban freight + geographic distance original v305 betas reproduced to 1e-9; site-pair Jaccard and total beta exactly sum over per-taxon contributions.

## New data: no-urban freight coefficient (the SAME source v305 no-urban baseline)

| Habitat | Real 2020 freight beta | Global 89-site richness-and-taxon-frequency null mean | Real minus null | Technical upper rank |
|---|---:|---:|---:|---:|
| BS ballast shoulder | +0.006309 | +0.003668 | +0.002642 | 0.3835 |
| BT ballast toe | +0.008956 | -0.000260 | +0.009216 | 0.0602 |
| BS union BT | +0.015538 | +0.006035 | +0.009504 | 0.0752 |

The earlier original v305 1km urban-adjusted coefficients are BS +0.004509, BT +0.007504, union +0.011860, **not the same coefficients as the new no-urban taxon accounting**. The original 89 site-level WorldCover covariates were NOT available in the archived aggregate, so this study does not claim new urban-adjusted per-species attribution.

## The crucial falsification: constrain year AND original single/double-track survey layout

| Habitat | Year x single/double layout frequency-fixed null beta | Real minus this null | Technical upper rank |
|---|---:|---:|---:|
| BS | +0.004522 | +0.001787 | 0.4110 |
| BT | +0.011291 | **-0.002335** | 0.6617 |
| BS union BT | +0.015557 | **-0.000019** | 0.5038 |

Main union positive freight association is **not exceptional** after preserving each source site richness and each source species occupancy WITHIN original survey-year and track-layout strata. This directly weakens the interpretation that a distinct network-transport force explains its coefficient. It is also NOT proof that train dispersal does not occur.

Among all 48 possible true source freight label arrangements WITHIN original national lines only, observed no-urban geometric adjusted beta has absolute size at least as large as the rearrangements BS 6/48, BT **1/48**, union 9/48. This is a conditional/technical comparison, NOT a randomized traffic treatment p-value; only 2 of 5 source lines have variation in multiple matched traffic intervals. BT's relative rank and its non-extreme source year-layout fixed plant null refer to different hypothetical models and must both be reported.

## Species driving the positive coefficient

In BT, 28 explicitly source-annotated alien taxa contribute +0.010675 to observed total +0.008956; nonexplicit aliens including unclassified taxa jointly contribute -0.001718. Strict annuals 60 taxa: +0.009809. Original source wind-archetype Epilobium/Phragmites 47 taxa: +0.011012. Source strategy archetype Bidens 9 taxa: -0.000616. Guilds overlap and their positive and negative coefficients should not be summed across traits.

| Original species source code | Taxon | Additive BT freight slope contribution | Native-status annotation |
|---|---|---:|---|
| ERIAN | Erigeron annuus | +0.00607 | explicitly invasive alien |
| ERICA | Erigeron canadensis | +0.00377 | explicitly invasive alien |
| PICHI | Picris hieracioides | +0.00286 | missing |
| MEUOF | Melilotus officinalis | +0.00269 | missing |
| EUOEU | Euonymus europaea | +0.00238 | missing |
| AMBEL | Ambrosia artemisiifolia | +0.00229 | explicitly invasive alien |
| SOOGI | Solidago gigantea | +0.00168 | explicitly invasive alien |

These individual coefficients are *accounting components*, not independent effect significances or biological dispersal effect estimates. Source Nativeness.status missing does not equal native. The source Diaspore.strategy field contains archetypes, not observations of train-borne attachment or seed viability. Source Bidens archetype is **not** a positive primary driver in these three aggregate models.

In BS: total +0.006309, explicit alien +0.004335, strict annual +0.000228, source wind-archetype +0.004937.
In union: total +0.015538, explicit alien +0.010394, strict annual +0.008487, source wind-archetype +0.009250, source Bidens -0.001762. Positive coefficients may be larger than total because other species make negative contributions.

## Scientific decision

The previously attractive two-scale fingerprint is still an exploratory observation, **not an identified mechanism or universal railway phenomenon**. The original local BS/BT negative vs freight and the distant union positive vs freight are separate conditional regressions on the SAME one publication; a natural non-rail local flora contrast, freight source-to-destination, field-year traffic, management and seed arrival are absent. No railway-wide transport/recolonization causal claim.

First-priority next audit is keeping the real same-site WorldCover covariates, and testing whether the few explicit alien source categories have a *frequency-conditioned* unexpected contribution, rather than merely a large accounting component. Independent 2022/23 railway site×taxon and traffic count datasets are required before comparative meta-analysis.

Original source details and methods: docs/tetsu_v306_shared_species_traffic_falsification_prereg_ja.md; scripts/analyze_Slovenia_traffic_shared_species_v306.py. This is a DRAFT research PR stacked on PR #34. Main unchanged.
