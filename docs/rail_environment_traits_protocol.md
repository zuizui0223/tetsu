# Railway life-history and phenotypic response workstream (v0.3)

**Scope:** Supplement, do not replace, the existing seed-transport / network workstream.
**Status:** hypotheses, photo annotation specification, and quality-controlled measurements only; **no empirical railway effects have been measured.**

## Main question
Do repeated railway disturbances filter or induce changes in plant architecture, reproduction and flowering schedules, and are those differences stronger under cutting and slipstream regimes than in matched roadside plant populations?

The ecological systems must be kept separate:
1. **transport** (seed movement between sites);
2. **environmental filtering** (only resilient taxa persist);
3. **phenotypic plasticity** (same genotype/individual changes its phenotype with exposure);
4. **evolutionary adaptation** (heritable response to selection).

A roadside-vs-rail phenotype difference alone identifies none of 2–4 unambiguously.

## Testable predictions and competing explanations

### P1 — mowing escape
Frequent/early cutting is associated with low flowering positions, basal/rosette growth, more detectable resprouting after cutting, or changed timing of fruit maturation. *But*: species turnover, mowing by roadsides, mowing season, selective herbicide, irrigation and observer framing could recreate the association.

### P2 — train airflow
Train-induced turbulent wind/slipstream exposure is associated with low and flexible growth forms or lower canopy stature, after matching cutting, ambient wind, heat, ballast substrate and vegetation context. *But*: a railway line, distance to tracks or a train's scheduled speed is NOT a measured air-flow exposure. Barrier shielding, distance, train morphology and terrain alter flow. If no safe measured peak gusts and exposure, classify P2 as a proxy-only exploratory test and avoid "wind adaptation".

### P3 — flowering schedule
At matched GDD, climate, species and sampling effort, repeated photo-based **flowering probability across DOY** differs near rail vs road/away and may be interrupted/restarted after documented cutting. Distinguish temporary suppression of visible flowers from genuine seasonal earlier/later flowering; if no pre- and post-cutting repeated observations, call this a **distribution of observed phenophases**, not onset shift.

### P4 — persistence / regrowth
When the **same individual or fixed quadrat** is imaged with unchanged geometry before, immediately after, and later after cutting, relative recovered projected area and flowering are quantifiable. This is evidence of *within-season response*, not automatically genetically based plasticity or adaptation.

## Photo measurement contract

Each image requires: unique image and site IDs, taxon (validated), date and local date confidence, provenance/license, observer ID where available, site class, view class, whole-plant/occlusion/truncation flags, phenophase, plant mask, and **separately annotated** reproductive mask when available. Multi-species scenes must have an unambiguous target individual.

**Core derived metrics (from valid masks):**
- Plant foreground area fraction: plant-mask pixels / frame pixels — **photo composition/scale-sensitive**, use primarily as QC, not as biological size in crowdsourced images.
- Plant bounding-box aspect ratio width/height; silhouette solidity: mask area / pixel-corner convex hull area. **Both view-sensitive**; compare within taxon and matched full-plant view strata.
- Floral ROI fraction: flower pixels within plant mask / plant pixels. This is *visible projected flower fraction*, NOT reproductive biomass investment and NOT phenological timing.
- Flower vertical position: (bottom of plant bbox - flower-mask centroid y) / bbox height, 0 near base and 1 near apex. Only for complete, upright, comparable camera views; it is not an absolute flower height.
- Flower present / bud / fruit / none / unknown: phase labels and their image-confidence separately. A missing reproductive mask is **unknown**, not biological absence.

**Quality gates:** masks share dimensions, flower mask is almost entirely inside plant mask, target plant covers at least 1% of frame, no foreground within 1 pixel of the image boundary for whole-plant inference, camera view comparable, photo source and observer splits. Retain rejected observations and rejection reasons. Optional PIV/slipstream sensors must be separate from citizen photo inference.

The `scripts/extract_photo_roi.py` pipeline produces image-level **descriptive features only**, with reasons. Images and masks are inputs; no neural segmentation model is asserted or trained.

## Phenology sampling design

Sample repeated images at fixed safe-access plots or leverage existing records with **flowering AND nonflowering** observations. Collect species, site, DOY, year, local GDD, recent mowing days, train/road classifications, photographer and plot survey effort. Use a seasonal binomial model:

`Pr(flower present) = logit^{-1}(s(DOY) + GDD + rail_exposure × s(DOY) + documented days_since_cut + (1|species) + (1|site-year) + (1|observer))`

For a general across-species analysis, allow species-specific seasonal effects/slopes. Separately test **flower/fruit succession** and **secondary flowering peak after cutting**. Do not use the first photographed flowering day as true onset when nonflowered visits and effort are unavailable. With presence-only flower photos, report observation-time distributions only.

Negative controls: road verges with similar mowing, nearby nonmanaged grass, same plant species from similar temperature/urban strata, old/unused rail (after confirming maintenance), and publicly accessible areas outside strong slipstream. Avoid assigning zero cutting intensity to unobserved dates or management policy.

## Small field design, independent of railway rights-of-way

Repeated quadrat photos from **publicly accessible and legally safe** locations, or with operator authorization: ≥10 independent corridors (start with ~10 matched groups, not a power claim), each rail-adjacent, road-adjacent and background plots; repeat at 1–2-week intervals plus any documented cut events. GPS uncertainty must not put plants falsely inside rail land; retain offset/measurement uncertainty. *No entering tracks or railway maintenance property without permission.* For slipstream instruments, work only with operator approval and safety planning.

### Minimum metadata for a management event
`site_id,event_date,event_type,evidence,cut_height_cm(optional),herbicide(optional),confidence`.
Types: `cutting`, `herbicide`, `steam`, `mulching`, `sheeting`, `unknown`. Never equate any vegetation suppression with cutting. Railway operator cutting calendars often are **not public**.

### Contrasts and causal strengthening
- Difference-in-differences / event-study on **verified cuts** and repeated *fixed* plot photos, matched uncut roadside plots; pretrend diagnostics.
- Individual plant repeated measures / standardized camera angle is stronger plasticity evidence than unrelated observer photos; common-garden/transplant + genetic differentiation is needed to claim heritable local adaptation.
- Factorial `cutting history × measured slipstream × soil/moisture/heat` with site/species random effects when enough independent treatments exist; otherwise report each descriptive contrast separately.
- Out-of-site and out-of-region validation; exclude duplicate burst photos, account for observer and sampling-season bias.

## Immediate analysis priority

1. **Photo quality and feasible traits**: build annotated sample with plant and flower masks and classify each photo as full plant vs close-up and known vs unknown flowering status. Evaluate manual labeling agreement and within-species consistency.
2. **Flowering time**: assemble observation date + phenophase for same species from rail, road and matched background, obtaining nonflowering denominator or repeated plot visits; audit GDD and observation effort.
3. **Cutting chronology**: operator notices + repeated same-plot imagery, record `unknown` when not documented.
4. **Slipstream**: physical or defensible exposure data last; exploratory rail-distance proxy never called "strong wind adaptation".
5. Connect trait analysis to the existing rail transport hypothesis only after both workstreams provide independent evidence.

## Literature and source pointers (checked 2026-11-07)

- Railway Technical Research Institute, vegetation control and fast summer regrowth: https://www.rtri.or.jp/rd/seika/2022/2-21.html
- Tokyu railway weed management includes cutting, herbicides and weed-suppressing sheets: https://www.tokyu.co.jp/railway/safety/weeding/
- JR Central notes actual cutting calendars are not public: https://faq.jr-central.co.jp/61ffebc2-b485-4402-809c-dcec8052e4a2
- US Federal Railroad Administration on turbulent train slipstreams/peak gusts: https://railroads.dot.gov/sites/fra.dot.gov/files/2022-07/Arup%20Aero-A.pdf
- Road verge mowing differences and flowering success: https://doi.org/10.1111/j.1654-109X.2007.tb00528.x
- Photo phenology bias and detection best practices: https://pmc.ncbi.nlm.nih.gov/articles/PMC10457241/
- Yucca photo phenology whole-plant and observer bias: https://pmc.ncbi.nlm.nih.gov/articles/PMC6976896/
- CNN flower/nonflower classification from iNaturalist: https://pmc.ncbi.nlm.nih.gov/articles/PMC8801702/
- PhenoVision flowering/fruiting annotator on community images: https://doi.org/10.1111/2041-210X.70081
