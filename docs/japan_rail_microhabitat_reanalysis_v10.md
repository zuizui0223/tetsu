# v1.0 — Japanese railway vegetation: original reanalysis of published microhabitat counts

**Created 2026-11-08.** [Verifiable GitHub Actions reanalysis](https://github.com/zuizui0223/tetsu/actions/runs/42436745865). [Reconstructed source-count CSV](../data/literature/abe_2012_railway_microhabitat_occurrences.csv). [Auditing and calculation code](../scripts/reanalyze_abe_2012_rail_microhabitat.py). **All reported source cells and regional margins were cross-checked** against the paper.

## Citation / actual input design

Abe, Kuribayashi, Tsuyuki & Ogasawara (2012), **"Vegetation of Railroad Site in East Japan"**, *Journal of Weed Science and Technology* 57:159–163. DOI [10.3719/weed.57.159](https://doi.org/10.3719/weed.57.159), **Table 4, printed p.161**.

- JR East railway 2011 management questionnaire, **109 surveyed railway sections** across 12 branches, reported **24 plant families and 62 taxa**. Geographic split: **41 warm-region sections** in Tokyo/Yokohama/Hachioji/Chiba/Omiya/Mito/Takasaki; **68 cold-region sections** in Nagano/Niigata/Sendai/Morioka/Akita.
- Table 4 is **top 10 reported taxa per region plus a residual "other" class**, with mentions in **six rail microhabitats**: embankment slope (法面), formation base (施工基面), trackbed (軌道), outside track, boundary and other.
- Our curated CSV retains **11 × 2 = 22** published taxa-or-residual records; the list of top 10 taxa **differs between regions**. An absent row is never imputed as a biological absence.
- Cell counts are **reported weed occurrences** by microhabitat, potentially from the same railway section/manager/species; they are **not independent plants, a random treatment, or a census of botanical survey quadrats**.

## Verified original-data reconstruction

Every row's six microhabitat counts adds to its published total, and all regional field totals reconcile exactly:

| Published source count | Warm region | Cold region | Pooled |
|---|---:|---:|---:|
| All species×habitat reported occurrences | **592** | **742** | **1,334** |
| Embankment slope occurrences | 247 | 395 | **642** |
| Formation base occurrences | 229 | 162 | **391** |
| Trackbed occurrences | 16 | 54 | **70** |
| Off-track | 30 | 16 | 46 |
| Boundary | 10 | 5 | 15 |
| Other | 60 | 110 | 170 |

The row and column margins form a useful reproducibility checksum; automated tests refuse a damaged source CSV.

## New descriptive ecological reanalysis (not a new field experiment)

The clearest source-backed microhabitat association is **`Equisetum arvense` / スギナ**:

| Metric | Warm | Cold | Pooled |
|---|---:|---:|---:|
| Horsetail reported occurrences in trackbed | 7 of 16 | 17 of 54 | **24 of 70** |
| Horsetail reported occurrences on slopes | 8 of 247 | 14 of 395 | **22 of 642** |
| Horsetail share of trackbed species reports | **43.75%** | **31.48%** | **34.29%** |
| Horsetail share of slope species reports | **3.24%** | **3.54%** | **3.43%** |
| Ratio of trackbed vs slope species-report shares | **13.51×** | **8.88×** | **10.01×** |

Among all microhabitat reports pooled, horsetail accounts for **98/1334 = 7.35%**; among trackbed reports it is **24/70 = 34.29%**, an **observational compositional overrepresentation ratio of 4.67×**. Both regions have the same direction. These ratios were calculated **independently by the new code**, not stated as causal effects in the original paper.

**Biological hypothesis generated:** Railway microhabitat occupancy is not a single "all grasses benefit" process; trackbed conditions/management may favor certain vigorous resprouting or persistent taxa (horsetail) while other vegetation is concentrated in embankments. However, source categories conflate cultivation/control treatment, region, accessibility and survey frequency; horsetail's physical tolerance, wind resistance or selection mechanism **was not measured**.

**Important contrast:** *Solidago altissima* / セイタカアワダチソウ has **0** reported trackbed occurrences and **46** slope occurrences in the published Table 4, despite **90** pooled occurrences in all rail microhabitats. This is a reported microhabitat composition, **not evidence that the species was absent from all railway rights-of-way**. *Setaria* spp. / エノコログサ類 have **3** reported trackbed occurrences across the two regions, not the predominant trackbed group. This cautions against oversimplifying railway habitats as selecting only short annual grasses.

There are **no p-values or confidence intervals** because the original source aggregates counts across not-necessarily-independent railway segments, management categories and species. Pseudoreplicated pairwise p-values would make the study misleading.

## Independent treatment literature supports the filtering mechanism, not our estimates

- Ito, Ueki & Sakamoto (1982), **[doi:10.3719/weed.27.41](https://doi.org/10.3719/weed.27.41)**: 37 national railway depots plus approximately 367–370 km observational Osaka train-window survey, documented *Setaria* spp. more under historical herbicide regimes than cutting, whereas *Solidago* and *Imperata* tended lower under herbicide. Observational, not a randomized herbicide/mowing experiment.
- Hayakawa (2007), **[doi:10.2219/rtriqr.48.240](https://doi.org/10.2219/rtriqr.48.240)**: 2002–2004 Okayama railway embankment field comparison: 5 treated/control field plots (untreated plot subdivided into A/B), selective repeated **bispyribac-sodium growth retardant**, not an ordinary "mowing-vs-herbicide" trial. Author reports approximately **30% lower above-ground biomass by year 3** under treatment and a transition from *Solidago altissima* and *Artemisia princeps* dominance to smaller *Setaria* and *Digitaria* grasses. One experimental field with limited spatial replication: cannot generalize as a modern national effect.
- Abe et al. (2012) surveys are **not a temporal repeat** of Ito et al. (1982): different surveyed geography, methods, taxon grouping and denominator. Do not call either a 30-year trend estimate.

## Proposed original field/data research

The sharpest new independent questions are now:

1. **Microhabitat filtering:** at verified railway *trackbed* vs *slope* vs matched *road verge* and unmanaged sites, do horsetail and other resprouting perennial taxa dominate on strongly managed ground? Measure repeated permanent quadrat taxon cover. **Do not equate distant georeferenced iNaturalist photos with trackbed plants.**
2. **Management type and regeneration:** separately record `cutting`, `herbicide`, `growth_retardant`, `sheeting`, `unknown`. A growth retardant is NOT equivalent to a herbicide spray or physical cutting. Compare same plant species' **vegetation projected ROI recovery and flowering fraction**, using fixed-view photos and dated real treatment logs.
3. **Slipstream:** with operator-authorized anemometers, compare actual peak gusts along safe rail locations to ambient wind, independently of management/ballast microclimate. No wind adaptation inference from these source counts.
4. **Propagule dispersal (independent main hypothesis):** rail network/freight/ports and citizen/GBIF species geography; a transported seed may arrive but be filtered by trackbed management before flowering. Separate arrival and establishment in ecological interpretation.

## Data/analysis safety

Curated numerical rows are a **factual, attributed reconstruction** of published results; the original PDF/images are not republished in this repository. The code writes:
- `results/literature_rail_microhabitat/summary.json` — headline overrepresentation and all source consistency flags.
- `results/literature_rail_microhabitat/species_habitat_composition.csv` — per taxon/region and habitat descriptive shares, with repeated-measure caveat.

The output files are CI artifacts, not ecological effect-size or selection tests. Use provenance and the source DOI for publication, and obtain raw independent sampling units before any significance testing. **No newly observed railway treatment response has yet been established by tetsu.**
