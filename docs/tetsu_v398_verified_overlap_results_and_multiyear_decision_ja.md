# v398 比較可能性監査の実測結果と多年度解析の判定方針
2026-11-25. Original source run [45549118059](https://github.com/zuizui0223/tetsu/actions/runs/45549118059) SUCCESS, source commit 4cbb1e6b640317b347080720244b91a7a649e8c3. Original v396 maximum-cardinality pair counts reproduced in every panel. Outcome-blind 9-feature 0.75 SD caliper, exact GloNAF region × continent × 1000km block. No raw records downloaded or reclassified for this audit.

## 100km all-effort support
| Boundary guard | Land fraction | Mapped-zero cells | Zero with any partner | Positive with any partner | Admissible edges | Optimal pairs | Blocks | Zero ESS | Positive ESS | Maximum weighted imbalance |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.3° | 0.5 | 52 | 15 | 23 | 28 | 15 | 8 | 12.3871 | 19.5254 | .220351 |
| 0.3° | 0.8 | 51 | 14 | 21 | 26 | 14 | 8 | 11.2941 | 16.3986 | .231610 |
| 0° | 0.5 | 84 | 34 | 56 | 76 | 33 | 16 | 24.5760 | 31.7139 | .159574 |
| 0° | 0.8 | 83 | 32 | 52 | 69 | 31 | 16 | 23.0977 | 30.7979 | .179884 |

Only 0/0/1/1 additional distinct mapped-zero cells have a partner beyond optimal 1:1 matching. The combinatorial pairing algorithm is not the principal bottleneck under the fixed restrictions. In high-effort conditions eligible zero cells with a partner are 2/2/11/11, block counts 2/2/7/7, zero ESS 2/2/9.1875/9.1875. No admissible edges in any 200km all-effort panel; high effort no mixed strata. Every v398 panel fails at least one prespecified support criterion, including zero ESS>=30. Therefore no v398 biological estimate, not a zero estimate or a negative biological finding. ESS is a weight concentration diagnostic and does not guarantee spatial independence. Strict conditions also fail the mean imbalance <=.20 gate.

## Scientific disposition
- v396 relaxed one-to-one comparisons remain descriptive, +1.233256/+1.450385 pp for observed GloNAF-documented fraction in 33/31 pairs, intervals positive; they do not become stronger through all-edge weighting.
- v394 simultaneous population and local control substantially attenuates count associations. v395 positive-only railway length slopes all cross zero; binary rail presence association under relaxed region support only.
- Thus no established rail-specific dispersal mechanism or globally representative invasion effect. The observed association is a geographically restricted landscape pattern. Nonlisted species are not native; mapped zero is not verified absence; no time-of-establishment information.

## Multiyear decision
v397b run [45548846372](https://github.com/zuizui0223/tetsu/actions/runs/45548846372) is at the exact original snapshot scan/model step as inspected, no source output artifact yet. Retain frozen 2023/2024/2025 and pooled effort/geometry/matching/support gates. Only annual or pooled results passing the original reproduction gates and specified support should be interpreted. Separately report annual newly matched and frozen 2025-pair retests (not independent ecosystems), and 3-year cumulative unique species/observer sets (not sum of species counts). If strict/high-effort support remains inadequate, record nonidentifiability rather than relax conditions in response to outcomes. Historical observations use one frozen rail geography and GloNAF classification; cannot identify opening-year effects or temporal invasion rates. No additional parallel outcome analyses before inspecting the multiyear source. Main unchanged.
