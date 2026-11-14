# tetsu v13.0 — DIRECT measured grain/seed-kernel deposition from moving freight TRAINS (not just observed feral plants)

## Main finding: a genuine physical RAILWAY cargo propagule input is measured, but original crop-seed VIABILITY and realized establishment remain unmeasured

**Original authors Gangadharan A., Pollock S., Gilhooly P., Friesen A., Dorsey B. & St. Clair C. C. (2017), *Animal Conservation* 20:391–400, original peer reviewed article [DOI:10.1111/acv.12336](https://doi.org/10.1111/acv.12336), [original publisher full text](https://onlinelibrary.wiley.com/doi/10.1111/acv.12336/full).**

This is not an additional tetsu field study; published original Canadian Pacific original samples were actually collected by the university group. They physically **counted kernels** from spilled freight train cargo along the original **134km** rail stretch in **Banff and Yoho National Parks**, Canada. There are **19** independently monitored original rail field sites, sampled approximately every **14.1±2.7 days**, repeated 2013,2014,2015 April–October and 2014–15 winter. The source sample count is **1,979 collection-device field samples** (652 screens, 614 baskets, 614 boxes, 99 lids) nested in 19 sites. **1979 is NOT 1979 independent railway stretches or a seed viability sample.**

**The source authors measured deposition of kernels/physical grain** (wheat+barley, Brassica canola seed, peas+soybean, lentils, canary, flax), not the phenology or fitness of plants later developing from those kernels. Subsequent grizzly bear attraction was the original paper's ecological focus, not a germination bioassay. Unlike Austria Pascher original 2016's only feral-PLANT observations or von der Lippe et al 2013's direct CAR seed movement, this Canadian study is **physical sampled cargo grain deposited by actual moving TRAINS**. However source **viability/germination of deposited kernels** was **NOT measured**, and original seed-to-established flora chain is unlinked.

### Source-provenance and actual code

**[v13.0 successful authentic published numerical source reanalysis: GitHub Actions #43590341083](https://github.com/zuizui0223/tetsu/actions/runs/43590341083)**. [Exact manually author-original published source summary CSV](../data/Canada_Gangadharan2017_original_train_grain_seasonal_summary_v130.csv), [reproducible script](../scripts/analyze_Canada_true_train_spilled_grain_2017_v130.py), [source arithmetic and false-viability safeguard tests](../tests/test_Canada_true_train_spilled_grain_2017_v130.py), [original-source workflow](../.github/workflows/Canada-direct-TRAIN-spilled-grain-source-v130.yml).

Original article full text and its text were independently read from the publisher via scholarly web retrieval; the **GitHub Actions job itself could NOT download the Wiley publisher HTML, and the report's source-access receipt explicitly says so**. No site×visit raw matrix has been acquired or fabricated; secondary reanalysis is of exactly original author printed summary numbers only. Source CSV author numbers are manually transcribed and stringently checked against published printed values, not represented as byte-identical original scientific instrument output.

Original article provides several supplementary methodology/model tables (Appendix S1-S5), but those are NOT a 1979-row raw plot-by-grain sample matrix. No raw site×visit model was fitted. Original paper's predictors and R² are quoted as originally published, NOT re-estimated from reconstructed fake points.

## Actual source physical train deposition data by original year and collection method

| Original measurement window | Actual collector | Original reported all-grain mean | Original author SD | Original author **134km extrapolated** spilled tonnes |
|---|---|---:|---:|---:|
| April–October **2013** | Basket | **1.12 g m⁻² day⁻¹** | 1.23 | **46.0 t** |
| April–October **2014** | Basket | **0.87 g m⁻² day⁻¹** | 0.71 | **35.6 t** |
| April–October **2015** | Basket | **2.75 g m⁻² day⁻¹** | 5.52 | **112.8 t** |
| November 2014–March **2015** | Screen (winter-only deployed method) | **1.52 g m⁻² day⁻¹** | 2.37 | **44.6 t** |

Original publication pooled warm-season true author reported rate **1.64 g m⁻² day⁻¹ (sd 3.60)** across the three warm-season observation years; **not equal** to a simple mean of separate yearwise displayed values (1.58), since site/visit methods and pooling weights differ.

- Original author warm-season model projections mean (46.0+35.6+112.8)/3 = **64.8 tons** (source SD of yearly projected mass 41.9 tons).
- One *observed* winter model projection = **44.6 tons** across 134km.
- Adding one winter projection to the mean warm-season projection = **109.4 tons**, consistent with source original manuscript's **rounded "≈110 tons annually"**. It assumes the single winter monitoring interval is representative of other winters. **110 tons are an author's model extrapolation, NOT material actually physically weighed or recovered along 134 km**.
- Warm-season author reported 2015 mean deposition rate / 2014 = 2.75/.87 = **3.16×**, while annual model-extrapolated mass 112.8/35.6 = **3.17×**. This source year variability is large, with original 2015 SD5.52; **three years alone do not test time trend, climate effect, crop genetic adaptation or global growth**.
- Original source grain mass composition rounded percentages: wheat+barley **46%**, canola *Brassica* spp. **25%**, peas+soybean **18%**, lentils **10%**, canary **1%**, flax **1%**. The printed original class percents sum to **101% due rounding**; do not normalize, interpret 25% as exact canola kernel survival probability, or multiply annual tonne extrapolation into a measured viable canola seed delivery by location. Crop grains are often seeds but physical presence doesn't establish germination capability.
- Original collection method differences are large: compared with source basked collection rates, box/lid/screen source captured **77%/37%/24%** respectively in the paper's method comparison. Original 2014–15 winter was **screen-only**, so raw winter/summer differences are not automatically a controlled phenological effect.
- Original source study seasonal grain shipment/month were more consistently selected as predictors than local train speed/topography, and the author reported models explaining only **9–31%** of deposition variation. Thus cannot predict precise site grain flux from generic corridor class or national freight tonnage without new matched instrument data.

## Critical ecological connectivity insight

True train cargo grain input directly demonstrates a physical anthropogenic resource subsidy **at the source railway track**; this creates possible ecological linkages not captured by plant richness-only data:

1. **Transport delivery** — original Canadian freight trains directly dropped physical grain kernels (captured and counted); grain type/season varied.
2. **Plant recruitment** — source did **NOT** test whether the captured kernels were viable or emerged as seedlings at the track. Do NOT count the 2017 original grain samples as feral Brassica establishment.
3. **Community shift** — actual original Slovenian 89 within-rail matched BS/BT plant plot observations and Czech national railway-related survey species lists describe **already present plants**, but do NOT share the Canadian seed deposition sites/taxon-years or determine source seed arrivals.
4. **Wildlife food-web response** — original Canadian paper considered bears attracted to grain, and estimated **4.77×10⁸ kcal gross energy** from the 110-tonne annual resource projection, an energy amount comparable to 42–54 grizzly bears' annual caloric needs (the authors expressly warn that much of the grain is inaccessible, unavailable or consumed by other wildlife); this is **NOT evidence that 42–54 bears actually fed from or were supported by train grain**. Any bear-train collision causal hypothesis requires behavioral/collision matched data.
5. **Railway wind/slipstream, cut stress, size and timing** — NONE are experimentally measured in the Canadian author grain-trap study. Direct physical grain loss from hopper cars via leaks/trickles is a different source mechanism from passenger train-induced wind transport of native diaspores or adaptive flowering.

### Source evidence ledger now has a genuinely RAIL arrival anchor

| Original study/country | Source observation | Reanalysis availability | What it establishes |
|---|---|---|---|
| **Gangadharan et al. (2017), CANADA** | **Direct train-spilled physical grain kernels from 19 sites, 1979 repeated field samples, six cargo categories** | v13.0 published season summary real, original raw visit/source taxon kernel counts unavailable | **Physical train cargo grain delivery on railway tracks confirmed**; deposited crop seed viability/colonization NOT known |
| **Pascher et al. (2016), AUSTRIA** | **Sampled feral Brassica napus PLANTS at original 60 sites×2 years** | v12.8 source 60×2 original numeric table independently reverified against actual government publisher PDF 120/120 cells, 0 errors | Observed feral plants at transport-related sites, not cargo seeds or quantified post-arrival establishment |
| **von der Lippe et al. (2013), ROAD by car** | **Marked physically moved four plant species seeds under 48 km/h car airflow** | v12.9 original 4 taxa×five repeated source pass levels, 20 published summary cells | Physical source seed morphology×ROAD vehicle slipstream, NOT direct train wind |
| **von der Lippe & Kowarik (2007), ROAD tunnels** | **Viable vehicle-deposited seeds, motorway tunnel seed traps** | Only original peer-reviewed paper source published result (no raw 2007 data modeled by tetsu) | Road vehicle viable seed input by field traps, not rail |
| **Šabić et al. (2026), SLOVENIA** | **89 same-rail BS/BT actual plant community observations** | v11.0–v11.6 actual original flora trait/source analysis | Conditional species-frequency/rail microhabitat functional filtering, NOT seed traps |
| **Chytrý et al. (2026), CZECH Republic** | **117,739 vegetation original source relevés, selected RAIL-text/ROAD-text cohorts** | v12.7 exact Pladias plant traits, 18-study source community-class matched plots 105/167, 532 known plant taxa, all study-cluster CI zero crossing | Independent-national flora study, **not same railway within-track plot design** and no direct dispersal/colonization rates |

**Resulting credible macro science question:** When does transport-corridor cargo/vehicle-induced propagule INPUT correspond to actual colonization and to changes in functional traits, after maintaining site-year-species sampling equivalence? Published data demonstrate parts of the chain, but their source domains, vegetation communities, years and outcomes do not overlap sufficiently to estimate a joint causal effect.

## Next substantive data gate

Find and obtain original Canadian **19×visit×grain type seed count/weight** records (authors/source data repository on request if not publicly available), and link to genuinely matched original plant vegetation surveys and **germination experiments on the very same collected grain**. Separately, France RER original equal-area track/embankment 394 plot vegetation raw source available upon author request [10.1111/wre.70025](https://doi.org/10.1111/wre.70025), and Slovenia original 89 paired BS/BT plots can serve as established-community filters. For macro-multi-country generalization, demand at least two independent actual **train seed-trap** studies or explicitly frame one Canadian rail-delivery study as a source case rather than a globally representative effect.

No physical propagated seed train-sourced and established same-site outcome is yet independently identified. No train aerodynamic stress, grass-cutting adaptation or evolved timing mechanism is established. Respect original study authors' source governance; Pladias requests Governing Board discussion before substantial reuse. Keep **[tetsu PR #1](https://github.com/zuizui0223/tetsu/pull/1) as DRAFT**, do not merge new work into original default main.
