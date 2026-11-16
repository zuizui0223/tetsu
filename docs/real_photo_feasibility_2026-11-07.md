# First real iNaturalist phenophase × rail proximity pilot — 2026-11-07

**Data evidence, not evidence of evolutionary adaptation, mowing response, or rail-related flowering change.**

## Fixed query and code
- [Pilot protocol](inat_rail_phenophase_pilot.md)
- [Reproducible implementation](../scripts/audit_inat_rail_phenophases.py)
- [GitHub Actions success and original aggregated JSON artifact](https://github.com/zuizui0223/tetsu/actions/runs/42373021202)
- N02 source: official 2025 Japanese *passenger railway* line geometries only.

## Actual observations

Eight predeclared taxa × three locked Tokyo, Osaka and Sendai bounding boxes × two iNaturalist Flowers and Fruits annotation values, 2018–2025:
- API-reported explicitly annotated event counts: **189** (185 with visible open flowers, **only 4 without visible flowers or fruit**).
- Photos and observed coordinates satisfying un-obscured, noncaptive and <=100 m uncertainty: **132**.
- Of those, **37** fall within 250 m of a mapped railway line, 50 in the 250–1,000 m bin and 45 beyond 1,000 m.
- **0/48** annotation-by-species-by-city strata required capping at 100 observations for the selected query.
- Since a photograph can have more than one annotation value, these sums are not guaranteed distinct observation IDs; the JSON consists of separate stage-level aggregations, not a complete individual-level taxon sample.

| Predefined metro bounding box | Annotated counts | Geolocation pass | <=250m of passenger railway | No flowers or fruits visibly annotated |
|---|---:|---:|---:|---:|
| Tokyo | 136 | 104 | 28 | 1 |
| Osaka | 43 | 23 | 9 | 3 |
| Sendai | 10 | 5 | 0 | 0 |
| **Total** | **189** | **132** | **37** | **4** |

The species distribution was extremely unbalanced:
- *Trifolium repens*: **114 annotated flowering**, **2 no flowers or fruits** across the three bounding boxes.
- *Solidago altissima*: **61 flowering**, **1 no flowers or fruits**.
- *Plantago lanceolata*: **8 flowering**, **1 no flowers or fruits**.
- Other five taxa combined: **2 flowering**, **0 no flowers or fruits**.

**Consequent decision: explicitly reject rail-vs-background flowering onset or flowering probability inference from this particular labelled pilot.** The four "no flowers or fruits" observations do not support the contrast, and absence of any annotation is not nonflowering. Even the labelled nonflower state refers to photo *visibility*, not necessarily individual biological state.

The presence of multiple labelled, precise-location flowering photographs is useful as a **source of candidates** for independent photograph quality auditing; the apparent near-track photographs might be from station parks and sidewalks rather than rail vegetation.

## Next dataset strategy

1. Use the **complete species photo records rather than only labelled flowering records** to select licensed, geo-quality-checked photographic examples near and far from rail. Freeze an identifier-hash-based sample, not newest or most attractive flower photos.
2. Produce an attributed contact sheet and manually review target individual identity, whole-plant view, whether the camera cuts off the plant, and whether flower and entire plant can *both* be delineated.
3. For sufficiently good images, derive **pixel geometry only** from inspected plant and flower masks. If no defensible masks can be produced for enough photographs, do not claim morphology or ROI differences. Strict manual QC trumps apparent abundance of photos.
4. For actual flowering response and within-genotype plasticity, use repeated fixed-view plot/individual photos around *documented mowing events* with matching roadsides and uncut sites, temperature/GDD, site, species and observer.
5. For train wind, gather independent measured exposure. For heritable adaptation, use common-garden / genetics; do not label observational habitat filtering as evolutionary adaptation.

## Literature and ethics

- iNaturalist terms and annotation definitions: https://help.inaturalist.org/en/support/solutions/articles/151000191830
- Rate-limited API and bulk-data guidelines: https://www.inaturalist.org/pages/api+recommended+practices
- Photo bias and incomplete views: https://pmc.ncbi.nlm.nih.gov/articles/PMC6976896/
- iNaturalist annotation quality in large datasets and taxonomic differences: https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.70081
- Verify CC licenses, attribute original photo authors, do not reproduce restrictive photographs, and never enter railway operating land to gather data without authorization.
