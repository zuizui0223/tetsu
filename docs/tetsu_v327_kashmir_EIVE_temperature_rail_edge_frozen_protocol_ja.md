# v327 frozen source-analysis protocol — Kashmir railway vs interior thermal EIVE (2026-11-23)

## Primary scientific question
In the independent Kashmir railway field study (Rashid et al. 2021, DOI 10.1111/1365-2664.13961), do **track-zone communities** have a different EIVE 1.0 species temperature-affinity centre or width from **paired Interior Plot 2** at the **same source site, same survey year and same original 2×50m field plot area**?

**This is a railway-edge vs local interior observational contrast**, not randomized climate treatment, globally replicated climate-offset effect, train seed carriage, or the claim that the interior is pristine unaffected habitat. The ~130km Kashmir corridor has elevation 1591–1741m; it is *one* climatic corridor.

## Source-gated design, frozen before source result
- Original vegetation: Dryad version 99735 official ZIP, exactly one `Matrix_plot.xlsx` (SHA256 `529f3b97521933018e12623b61f531f012c67241cae2602a6675dbd534192474`). Preserve 31 field sites, 16 stations and 15 between-station sites, 2014/2017 repeated, four original source zones, 123 source taxon columns, four originally missing species cells. Do not assign these missing cells biological absence. Use source-authenticated `parse_xlsx_matrix` from v254 and independent accepted binomial parsing.
- Independent European *trait* source: EIVE 1.0 SM08 original XLSX from Zenodo `10.5281/zenodo.7534792`, MD5 `0f0ce19cd2a781eeed027d8b0acadbcb`; exact binomial match **only**, no automatic synonym guessing or geographic niche extrapolation. EIVE niche **centre** `EIVEres-T`, **width proxy** `EIVEres-T.nw3`, contribution count `EIVEres-T.n` **>=3**. Count all unmatched original observed species in per-zone trait denominators. Do not infer Asia-wide physiological tolerances from European consensus traits.
- **Primary**: Track minus **Interior Plot 2** (farthest), all original spatial sites if source nonmissing in that year, each community **>=5** EIVE matched species and **>=0.50** of original source observed taxon indices covered. Assess 2014 and 2017 separately, never treat repeated years as independent new studies. Report stations vs between as site strata.
- **Sensitivity**: (a) nearer Interior Plot 1, (b) lower 0.30/3 coverage and EIVE n>=1, (c) same spatial sites eligible in both years for within-site temporal change, (d) leave-one-site-out sign/range, (e) no matched sites / trait coverage -> explicit *nonestimable* instead of defaulting zeros.
- Each field site contributes one BS/Interior CWM delta. Report number of eligible sites, station/between counts, mean and positive signs, mean trait coverage for *all* site-years, and the distribution of trait matched species. No standalone source p values unless block/stratum dependence is addressed. Comparisons are exploratory hypotheses chosen using prior European data.
- Prior v288 experience with European LEDA trait missingness is NOT a priori proof of EIVE missingness; fail closed if EIVE species coverage is poor. This is a pre-specified hypothesis transfer check to non-European railway observations.

## Interpretation gate
If the direction is positive in India with strict coverage, that is an **external observational railway-edge result across source regions**, but still not proof of global railway microclimate effects or climatic range-limit release. If no strict eligible pairs survive, report the limitation rather than lowering thresholds and claiming the result confirmatory. A corridor-specific negative effect would directly challenge simple globally uniform warm filtering. Do not compare numbers directly across different EIVE scales, source plots and habitat grains without describing their contrasts.

No original biological source matrix or GPS should be republished in the repository.
