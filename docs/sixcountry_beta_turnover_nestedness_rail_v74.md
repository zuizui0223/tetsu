# v7.4 — REAL six-country railway macro: species turnover versus nestedness

**Actual source success on 2026-11-11:** [workflow #43035013575](https://github.com/zuizui0223/tetsu/actions/runs/43035013575) · [real JSON artifact #13861753546](https://github.com/zuizui0223/tetsu/actions/runs/43035013575/artifacts/13861753546) · [analysis](../scripts/analyze_rail_beta_turnover_nestedness_v74.py) · [tested equations](../tests/test_rail_beta_turnover_nestedness_v74.py).

## Biological question
Does apparent similarity of country-listed introduced flora along relatively short physical RAIL routes (versus ROAD) correspond to taxonomic **replacement** or **nestedness-resultant dissimilarity**, and is either relationship robust to source dataset provider removal and indexed-record count rarefaction?

The frozen source framework is **six prechosen countries**, one 200×200 km area per country, fixed 4×4 equal-area 50km cells, 2010–17 versus 2018–25. GRIIS-listed introduced species are matched to exact public GIFT 1.2.1 herbaceous/woody groups, then 2026 Overture physical RAIL/ROAD path lengths. Country blocks, not pairs that reuse cells, are the independent geographic summaries. Source providers are split FULL vs WITHOUT largest all-Plantae GBIF dataset chosen independently from the ecological response. Both provider variants and BOTH growth groups use the *same exact physical matched cell-pair cohort* within period×route-anchor; partial ranks hold geodesic distance, road route, all-plant GBIF record ratio/scale fixed.

Based on Baselga (2012), for shared species a, and unique species counts b and c:
- Jaccard total beta_jac=(b+c)/(a+b+c).
- Turnover beta_jtu=2min(b,c)/(a+2min(b,c)).
- Nestedness-resultant beta_jne=beta_jac−beta_jtu.

The **RAW** version uses recorded species incidence; the K=20 rarefied version uses expected inclusion probabilities (hypergeometric source-record drawing) and partitions a *ratio of expectations*, **not** the expectation of biological beta diversity. Neither corrects true observer/dataset bias or establishes species absences.

## Main actual results: 2018–2025, nearest-in-cell mapped rail anchor

After insisting on ≥8 unique cell sources, ≥20 physical source pair routes and ≥3 taxa/≥20 indexed records in EACH growth form and EACH source-provider set, only **4/6** country blocks have estimable coefficients (US Midwest, FR Burgundy, GB Yorkshire, DE Saxony). Canada retains 8 cells but only 16 valid matched rail/road pairs; Australia has 7 matched source cells. Do not convert unsupported countries to biological nulls.

**Woody GRIIS introduced flora:**

| Observational beta response | FULL provider: compatible positive / 4 blocks | Without largest provider: compatible positive / 4 blocks | Country median rho FULL → WITHOUT |
|---|---:|---:|---:|
| RAW total Jaccard beta_jac | 4/4 | **4/4** | +0.091107 → **+0.141035** |
| RAW species turnover beta_jtu | 2/4 | **3/4** | −0.003779 → +0.132464 |
| RAW nestedness-resultant beta_jne | 3/4 | **1/4** | +0.071566 → −0.001979 |
| Rarefied20 total beta_jac | 4/4 | **2/4** | +0.091771 → **+0.024878** |
| Rarefied20 turnover | 3/4 | 2/4 | +0.021850 → +0.010308 |
| Rarefied20 nestedness | 3/4 | 2/4 | +0.045506 → −0.015756 |

Positive beta-dissimilarity versus log(RAIL/ROAD path) is railway-hypothesis-compatible (more direct RAIL paths correlate with more similar recorded species). This reverses the signed convention for prior similarity correlations.

**Individual woody FULL-EXCLUDED source partial-rhos** (2018–25, nearest physical rail anchor):

| Country | RAW beta_jac | RAW beta_jtu | RAW beta_jne | rarefied20 beta_jac |
|---|---:|---:|---:|---:|
| US Midwest | +0.065273 | −0.182020 | +0.345308 | **−0.154987** |
| France Burgundy | +0.055597 | +0.103113 | −0.002880 | **−0.089708** |
| UK Yorkshire | +0.225203 | +0.328263 | −0.114606 | +0.445059 |
| Germany Saxony | +0.216798 | +0.161816 | −0.001077 | +0.139463 |

For herbaceous flora, FULL→WITHOUT total raw beta_jac has compatible directions 4/4→3/4, versus rarefied20 3/4→2/4, with latter country median +0.070427→−0.038920.

**Alternate source anchor:** Woody WITHOUT provider RAW total beta 4/5 positive (median +0.055597), rarefied20 only 2/5 positive (median −0.089708). **Earlier 2010–17 period:** only 2 qualified country blocks. Neither supports temporal or historical freight generalization.

## Verified integrity and inference boundaries

The workflow finished **SUCCESS**, all 8 source-math/provenance unit checks passed, and a local independently coded Jaccard decomposition matched the source results. The v7.4 rarefied total dissimilarity partial-rhos are exactly the negative of **52 comparable existing v7.2 similarity-rhos**, with pair support identical and maximum coefficient discrepancy **0.000000**. No pseudoreplicated p-values were estimated.

**Supported:** the tempting across-country RAW botanical association is metric-dependent, and species replacement and nestedness components have distinct signed behavior. RAW total +4/4 does not prove a robust biological effect: even after largest-provider exclusion, the *same cells and network paths* give rarefied20 only 2/4 compatible regions.

**NOT established:** actual propagules carried on trains; freight operating frequency, direction or seed release; local railway-verge plant establishment; selection by mowing, train wind, heat or herbicide; plasticity or evolved adaptation. The data are positive GBIF indexed records, not standardized site surveys, and the 2026 physical infrastructure postdates 2010–25 source records.

**Next macro identification gate:** obtain verified active historical freight route/volume and true railway seed-vector data, independent original GBIF event/provider/observer samples in multiple country×ecoregion blocks, and independently measured rail-versus-road verge vegetation + mowing/wind exposure. Then test *transport delivery × disturbance × real species dispersal and regrowth traits*. Growth form is only a broad functional class.

Reference: Baselga (2012), Global Ecology and Biogeography 21:1223–1232, https://doi.org/10.1111/j.1466-8238.2011.00756.x. CRAN betapart Jaccard beta.pair is the canonical observed-incidence decomposition.
