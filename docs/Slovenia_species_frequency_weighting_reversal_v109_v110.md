# tetsu v10.9–v11.0 REAL Slovenian railway source: same species, opposing wind trait contrasts when abundance opportunity changes weights

**Original study source authors:** Šabić, Grašič & Rovanšek (2026), *Vegetation patterns, functional traits and invasion dynamics across railway microhabitats in Slovenia*, DOI [10.1007/s44473-026-00155-z](https://doi.org/10.1007/s44473-026-00155-z), CC BY 4.0 supplement, original SHA256 `26be7ca9513dd54136d4ec499d8413c1cb288045f87fe3cf8814731f2edac99f`. This is tetsu's **secondary reanalysis** of field data from the original authors, **NOT** original field collection by tetsu.

**Actually successful source GitHub Actions:** [v10.9 reproducibility-corrected species-equal occupancy/permutation #43566361231](https://github.com/zuizui0223/tetsu/actions/runs/43566361231) (fully source-code-sorted deterministic species-label draws); [v11.0 exactly same source taxon set under different weights and rarity interaction #43566562331](https://github.com/zuizui0223/tetsu/actions/runs/43566562331). All source workbook files are re-downloaded from locked earlier successful original Dryad/publisher artifacts and SHA checked. Synthetic math tests and full generic CI have passed; no invented observed seeds, flower times or train gusts.

## Source and estimands, not a literal Simpson causality claim

From **89 authentic Slovenian railway sites / 173 matched Ballast Shoulder (BS) versus Ballast Toe (BT) sides across 15 original rail sections**, first take **union across eligible left and right sides WITHIN the same location** (no independent side pseudoreplication). For each recorded species with both source-known strict wind-diaspore strategy and lifespan coding (n=**261** taxa observed at least once; strict wind Epilobium/Phragmites **53**, other strategy **208**), define its railway microhabitat preference

**q_s = (n locations species present only in BS − n locations species present only in BT) / n locations present in BS and/or BT.**

This species preference is bounded [-1,+1], pertains to original field observed presence, and is not a plant's genetic, physiological or wind trait itself. If present in both BS and BT, this adds an opportunity but zero difference. Observed study BS and BT zone plot areas are NOT independently confirmed equal; differences therefore remain an area/detection- and habitat-confounded observational pattern. Source species codes and all original author strategy codes are read exactly, missing traits remain UNKNOWN.

**Two explicitly DIFFERENT estimands on the IDENTICAL source taxa and lifespan groups:**

1. **Equal species weight**: Within each strict annual/other-known-lifespan stratum, compare arithmetic average of q across wind-classified taxa against other taxa. Combine stratum differences proportional to **numbers of source taxa**.
2. **Occupancy-opportunity weight**: For each group, average the EXACT SAME species q values, weighted by **their original number of rail locations where observed in either BS or BT**, then combine lifespan strata proportional to total summed source site occupancy opportunities. This places more weight on frequently observed species but **does not measure individual abundance, biomass or seed flux**. It does NOT algebraically equal v10.6's site mean of composition fractions, and no quantities are interchangeable across those estimands.

### Genuine v11.0 sign reversal on ALL identical 261 taxa

| Minimum source sites with BS/BT presence per species | Source species retained | Equal-species source wind − other BS preference | SAME species weighted by occupancy opportunities |
|---|---:|---:|---:|
| **1+** | **261** | **+0.189388** | **−0.098125** |
| 3+ | 170 | +0.073379 | −0.113407 |
| 5+ | 120 | +0.066597 | −0.127400 |
| 10+ | 75 | −0.154450 | −0.204091 |

The **1+ row is a pure weight-scheme sign reversal on EXACTLY the same 261 observed taxa**. Removing each of **15** original rail sections in turn keeps **equal-species contrast positive in 15/15** and **site-opportunity contrast negative in 15/15**. This is not a positive universal wind filter or a proof of actual seed dispersal. The robust result is that **railway habitat composition and wind-related traits interact with source species occurrence frequency, and that estimated direction depends on the biological unit/weighting**.

Original v10.6 plot-composition result from jointly trait-known plants remains **strict wind BS−BT −0.037128** [−0.072362,−0.004974], with **within-lifespan strict wind component −0.040424** [−0.073793,−0.008518]. These are valid descriptively but they do **not** say that *each equally weighted wind species* prefers BT. v11.0 explains this non-equivalence; it does not replace an abundance-weighted source result with a contradictory false negative.

### Species-level permutation and rarity interaction, restricted validity

- **v10.9** with one original site per species minimum: n=261, equal-species **+0.189388**; lifetime-stratified wind-strategy label permutation 9,999 shuffles gives two-sided p **0.0357** (deterministic Monte Carlo null). However for **>=3 sites** n=170, +0.073379 with p **0.3530**; for >=5 n=120 +0.066597 p **0.3764**; >=10 n=75 −0.15445 p **0.0463**. Four cutoffs were assessed and *not* multiplicity-controlled; **DO NOT** claim independent/replicated 5%-level hypothesis tests.
- **v11.0 continuous species-frequency interaction** after adjusting original species strict-annual-vs-other-known lifespan label: **wind × log(source site occurrence count) = −0.193820** in a descriptive species-unit OLS. Wind category assignment permuted only within lifespan × four pre-created frequency bands (1–2, 3–5, 6–9, 10+ source observed sites), **3,999** label randomizations, centered two-sided **p=0.0165**. This is **post hoc** because bands and interaction were motivated by v10.9. The original same-species source p/coef is a sensitivity **not** a causal or phylogenetic association test; small per-band guild sample sizes (e.g., strict-annual wind n=1 at rarest band, only n=4 at 10+).
- **Sign robustness**: leaving out any single source original section (15/15) preserves the sign of BOTH weighting estimands. Separately, leave-one-original-taxon-out in the >=3 taxon-equal contrast retains the initial positive sign in **170/170** original taxa; one rare wind taxon does NOT cause the sign inversion. This corrects a reporting-code-only sign-change-count bug from the initial v10.9 exploratory run; **the corrected v10.9 #43566361231 has zero taxon-sign flips**, not 170.

**Do not promote a mechanistic claim:** Rare original source detected taxa may reflect sampling, zone widths, source landscape composition, plant competition and taxonomic/polyclade structure, not arrival on train carriages. More frequent species' stronger negative guild preference is a source occupancy/detection association, not direct evidence of wind-mediated mortality, mowing adaptation, heritable phenotype, or freight seed export.

## Next macro-ecological replication: ACTUAL independently sampled Kashmir railway transects source, currently access-gated

We found independent published railway surveys by Rashid et al. (2021) *Railways redistribute plant species in mountain landscapes*, DOI [10.1111/1365-2664.13961](https://doi.org/10.1111/1365-2664.13961): **31 Kashmir Himalayan railway T-transects surveyed in 2014 and 2017**, including track, embankment and interior locations (the biological design is DIFFERENT from Slovenian BS−BT). Their discussion explicitly describes rare pioneers being replaced by more abundant competitors and an evolving community composition. Independently sourced real Dryad source DOI is **10.5061/dryad.0p2ngf20h** (NLM dataset catalog and live Dryad API **v11.1 successful metadata #43566940782**). The article bibliography erroneously truncates it to **10.5061/dryad.0pngf20h** (missing digit 2).

- **Actual Dryad live API:** dataset id **60988**, version URL id **99735**; independently verified original **Matrix_plot.xlsx**, **122,187 bytes**, SHA256 `529f3b97521933018e12623b61f531f012c67241cae2602a6675dbd534192474`, file endpoint `https://datadryad.org/api/v2/files/546003/download`. Dataset title, DOI and file inventory live-verified; raw XLSX itself has **NOT** been successfully downloaded.
- In **v11.2 Actions #43567065358**, the genuine direct Dryad raw file HTTP download returned **401 unauthorized**, so file parsing/species-year check was skipped; **NO replication estimate was made**. This is a source access blocker, not a biological null. Do not report original Kashmir source raw data as acquired or imply 2014/2017 paired station counts exist in tetsu until access is resolved legitimately.
- French RER independently published Rakotoson et al. (2025) DOI [10.1111/wre.70025](https://doi.org/10.1111/wre.70025) includes **394 authentic railway/embankment 2×1m surveys and 255 taxa**, but authors' source says **data available on request**. It is a strong future original-data request target, not a completed tetsu cross-country species+wind replication.

**Priority genuine next tests** when at least one original independent country matrix is legitimately available and the source life-history/wind taxonomy can be harmonized: (i) repeat equally sampled locality-level BS-vs-BT/track-vs-embankment species preference with source-area adjustment, (ii) compare both species-equal and occurrence-opportunity-weighted lifespan-conditioned wind contrasts, (iii) test country-by-frequency interaction with country-level uncertainty, (iv) measure true train traffic, physical propagule catch and maintenance independently if causal transport/adaptation is desired. Keep existing PR as DRAFT; main branch unchanged.
