# tetsu v28.5 — Kashmir annual microhabitat replication and study eligibility gate

Verified on 2026-11-18 JST. This is a SOURCE-AWARE PREFLIGHT + six synthetic mathematical tests. **No new Indian annual proportion effect is claimed** because the study's original 123-taxon life-span annotations have not yet been authenticated. Main is unchanged.

## Independent field studies (not independent revisions or years)

|Study|Actual original design|Current annual-lifeform result|Use|
|---|---|---|---|
|Slovenia Šabić et al. 2026, DOI 10.1007/s44473-026-00155-z|89 unique field locations, 173 same-side ballast-shoulder/ballast-toe pairs, 15 railway route sections, sampled 2022/23|Author-original verified strict annual BS−BT = **+0.054703**; 60/89 sites positive|One original source study; BS/BT differs from Kashmir Track/Interior|
|India Rashid et al. 2021, DOI 10.1111/1365-2664.13961; Dryad 10.5061/dryad.0p2ngf20h|31 spatial sites, four original 100m² zones, 2014 and 2017, 123 taxa. All-four-zone complete two-year sites: **27** (15 station, 12 between)|Original source XLSX SHA verified; **123 taxon names not yet mapped to original life span**|A second independent country source, but zero new quantitative annual effect from v28.5|
|France Rakotoson et al. 2025 DOI 10.1111/wre.70025|394 genuine 2m² vegetation plots (128 track /266 embankment), 255 taxa|Original article reports annual/ruderal enrichment on tracks and biennial/perennial competitive flora on embankments|Independent published sign; source raw matrix and variance on request|
|China–Laos Li et al. 2026 DOI 10.1016/j.jenvman.2026.130380|14 rail stations, 84 vegetation/soil/fungi plots|No identical annual effect verified|Candidate only|
|Brazil Mendes, Sebastiani & Melo 2026 DOI 10.1590/2175-7860202677029|**Published Nov 13, 2026**; four railway stretches, 293 morphospecies (277 identified, 16 unresolved), **400 one-m² herb quadrats** near railway-platform edge|No matched nonrail plot data shared; original published checklist via DOI 10.6084/m9.figshare.34015290|New South American independent rail-flora survey, **not** yet an annual-gradient meta estimate. The 2023 Mendes thesis may share source sites, not another independent study|
|Czech Kutlvasr et al. 2024 DOI 10.23855/preslia.2024.247|39 stations, station-by-taxon source data|Railway microhabitats mixed|Not an authenticated ballast-to-verge contrast|
|Italy Adorni & Ghillani 2020 Parma|28 stations, 438 taxa|Railway microhabitats mixed|Not a matched zonal annual contrast|

Brazil original article: https://new.scielo.br/j/rod/a/6HK8dshbFGTdRQ3NhTDFFtC/?format=html&ilang=pt_BR&lang=en

## Three independent estimands

1. H1: site-paired direction of a measured annual/perennial filter. BS vs BT, track vs embankment, track vs interior and other definitions are **different physical gradients**, not automatically interchangeable standardized effects.
2. H2: taxonomic versus functional beta convergence against an independently sampled *local non-rail species pool*. Within-rail frequency-fixed null distance is NOT the same measure.
3. H3: traffic/network-mediated delivery beyond geographic distance, road connectivity and site flora; no species list or annual filter directly demonstrates train-specific seed movement or within-species genetic adaptation.

## Kashmir v28.5 source-authenticated analysis gate

Script: scripts/analyze_Kashmir_annual_microsites_v285.py.
Tests: tests/test_Kashmir_annual_microsites_v285.py.

This imports the existing actual Dryad source parser scripts/analyze_Kashmir_same_sites_reference_year_v254.py. Source original SHA256 is 529f3b97521933018e12623b61f531f012c67241cae2602a6675dbd534192474 for the author-published Matrix_plot.xlsx. There are 4 original missing species cells; do not impute 0. The original source 31 sites reduce to 27 spatial sites with **all four original zones in both years**; the study has 123 separate original species columns.

The strict independent trait crosswalk CSV must have columns **original_taxon,life_span,trait_source**. Allowed span values: annual, biennial, perennial, unknown. Every non-unknown classification must have evidence/provenance. Names must match the ORIGINAL taxon headers; an unmatched name fails rather than fuzzy linking to a congener. Unknown taxa are not coded as perennial. The source paper analysed 81 frequently occurring taxa using annual/perennial and herb/grass factors, but the per-taxon correspondence is **not authenticated** yet, even though published Supporting Appendices S1/S2 are listed.

Per site/year/zone, strict annual species fraction = n annual / n with known life spans. Per matched site/year, primary contrast = Track annual fraction MINUS Interior Plot 1 annual fraction. Predeclared sensitivity: Track−Interior Plot 2; Track−Embankment; Embankment−Interior1/2. All use the same source-complete 27-site cohort. Main results retain sites passing known-life-span coverage at least 70% in both zones and **both years**. Group years within source spatial site; unit of inference is the original site, not a plot-year and not a pair of taxa. Group station/between is stratification, not traffic intensity. Classification unknowns yield bounds if all are annual vs all nonannual. Bootstrap 1999 resamples original sites and its 2.5%/97.5% interval is technical, **not randomized-field causal CI**. No original incidence matrix or site coordinates are written out.

Run after acquiring the genuine source workbook and evidence-tracked crosswalk:

    python -m unittest discover -s tests -p 'test_Kashmir_annual_microsites_v285.py' -v
    python scripts/analyze_Kashmir_annual_microsites_v285.py --source source_original/v254/India/Matrix_plot.xlsx --traits source_original/v285/Kashmir_verified_lifespan.csv --output results/audit/v285_Kashmir_annual_summary_only.json --min-coverage 0.7 --bootstrap 1999

## Global universality stop conditions

The target >=8 truly independent field studies on >=3 continents is still **unmet**. A random-effects 95% prediction interval that does not intersect zero means that a *comparable exchangeable future study* is predicted to have the same sign within that model; it is NOT a guarantee that all countries have the effect, and small-study variance/heterogeneity estimates are unstable. Multiple years of a single original source study never count as independent k. Leave-country/continent-out requires truly independent source studies with comparable effect and sampling uncertainty; failing to beat a naive region prior is failure of transport.

Meta-regression aridity, local flora, ballast exposure, mowing, cargo/passenger and network topology only after sufficient independent observations and authentic covariates; do not manufacture p-values, SE or Hedges g from insufficient published summaries. GBIF/iNaturalist railway 30m versus 100–1000m buffers are volunteer-observation corroboration, subject to observer/provider/geospatial bias, not independent rail quadrat experiments. MIREN surveys are ROADS, not RAILWAYS.

Next source priority: (a) authenticate the Kashmir author taxon-by-life-span table for 123 names; (b) run paired-source contrast and missing-trait bounds; (c) obtain original French 394-plot matrix and its real field/block sampling covariance; (d) seek additional physical rail-versus-nonrail plots beyond Europe before making any universal meta statement.