# v351–v355：無作為全球陸地格子における鉄道網と植物観測努力の真正データ統合（2026年11月24日）

## 研究上の問題と観察的な結論

従来のチェコ/イタリア/スロベニア/クロアチア・ボスニア/スロバキア/カシミール/西豪州の鉄道植物研究は、生データが取得できた研究の事後的な集積である。これだけから「世界中の鉄道で暖地性植物が多い」と一般化する根拠はない。

これに代えて、鉄道・植物の出現データを見る前に、全球の陸地・気候帯だけから自前で有限標本を抽出した。現段階で真正データによる結論は次の３つ：

1. **自前抽出は実行可能。** Natural Earth 110m地理情報とWorldClim2.1 1970–2000年BIO1の真正原本を結合し、EPSG:6933の100km等面積格子49,640から、陸地中心・気候値有効13,437格子を母集団として確定した。6大陸×3気候帯の18層のうち16層が非空。各層2格子、計**32格子**を気候層内の同確率無作為抽出で固定。
2. **真の地図鉄道線分は標本32格子中13格子に存在。** 元の100kmセル内へOverture 2026-11-18.0の線分を正確にクリップし、元sourceの国や植物観測に依存せず全32格子を完走。13格子が本研究のrail subtype/class/flagsに一致し、19格子はソース取得完全でも該当路線無し、**ソース欠損0**。抽出セルに含まれる有効鉄道線分長の合計は**1,245.517km**。この合計を「全球鉄道総延長」と誤読してはいけない。
3. **植物観測の地域偏在は極端。** iNaturalist公式APIで2025年の公開・野生・research-grade植物検出件数を**全32格子の外接緯度経度矩形**で取得。元source全32API成功、矩形の元報告件数合計**43,025**、0件の矩形**11/32**、100件以上の矩形**15/32**。大陸別の機会的観測合計は欧州23,084、北米10,568、オセアニア7,144、南米2,060、アフリカ164、アジア5。矩形は元EPSG:6933セルと同じ形・面積ではなく、矩形の重なりも未除外なので**43,025は元のユニーク個体数や全球の植物豊度ではない**。

今回結合した32 sourceセルでは、鉄道あり13格子中10格子の源植物API矩形で100件以上が報告されていた。しかし観測者のアクセス/写真のGPS精度・植栽・法律上の鉄道敷地は完全未検証であり、**100件は鉄道近傍の植物観測数ではない**。

## Source-authenticated original analysis runs

| Source | Original analysis Actions | Important verified state |
|---|---|---|
| v351 Natural Earth map-units＋WorldClim BIO1 | [#45330702677](https://github.com/zuizui0223/tetsu/actions/runs/45330702677) SUCCESS | 13,437 source population, 32 original sampled centers, sample hash b65a13e32d65e8fa818bbefe2d7fd1f0e15c6c5cb97ff6f6212aea23cf138bf2 |
| v352 six continents warm vs nearby? NO, rail presence only | [#45331051469](https://github.com/zuizui0223/tetsu/actions/runs/45331051469) SUCCESS | first six source landcells (temperate) 3 mapped rail positive, 3 absent |
| v354 32-cell Overture fixed release | [#45331600841](https://github.com/zuizui0223/tetsu/actions/runs/45331600841) SUCCESS | all 32 actual cells source complete; 13 present, 19 absent, 1245.517km; first six rail status/km identical to v352. Aggregate original artifact ID 14891363991 |
| v353 iNaturalist official source count | [#45331986000](https://github.com/zuizui0223/tetsu/actions/runs/45331986000) SUCCESS | all 32 original source APIs, 43,025 cumulative bbox reports, 11 zeros, 15≥100; aggregate artifact ID 14892187868 |
| v355 original source runs and sample-ID join | [#45332475164](https://github.com/zuizui0223/tetsu/actions/runs/45332475164) SUCCESS | exact 32 original grid keys, source-shas, π and BIO1 values cross-verified; six continent and 18 climate strata totals rechecked, no source GPS/species names re-exported |

Source Natural Earth map-unit GeoJSON pinned SHA256 718ac51d68ed2bf5de4f563fff2e2e045db706c1d3222805f3e317140a64213b, 2 polygon make_valid repairs area drift <0.000001%; continent correction records French Guiana geographically South America and a 60E Russia split to distinguish Europe and Asia. WorldClim2.1 10arcmin BIO ZIP SHA256 00513224583665ec0f2f955a4ec252730c4deb2004cce9e793492a3f26df4dcf, annual mean temperature source units °C checked Tokyo15.578 and Prague9.345; not species thermal niches.

## Source design-based inference only

Unweighted sample mapped-rail occurrence = 13/32 = **40.625%**. But source climate strata have very unequal population N and were drawn n=2 within each eligible stratum. Appropriate **sample-only finite-population stratified estimator** for mapped rail presence in an original cell-center-on-land 100km cell is

P-hat = [sum over h: N_h * (n_source_rail_positive,h / n_sampled,h)] / 13,437 = **0.2106496986** (21.06%).

This describes mapped rail presence in the **defined finite land grid cell population**, not active train usage, true trackbed, rail-network length, all land area or the global rail-biome vegetation. The 2 source cells per stratum make uncertainty unstable, especially if both happen to have the same source rail status, so the estimate is **exploratory**. No narrow confidence interval or broad railway universality claim is justified. The Russia geographical band correction occurred before the v351 source rail/plant results and was documented; v351 original fully geography-corrected sample SHA is immutable.

### Six continents original unweighted opportunity (not continent-representative)

| Continent | Original sampled 100km cells | Actual mapped rail-present cells | Mapped rail km summed in sampled cells | 2025 original iNat bbox report sum | Original bbox report zero cells |
|---|---:|---:|---:|---:|---:|
| Africa | 4 | 1 | 0.098 | 164 | 2 |
| Asia | 6 | 1 | 105.844 | 5 | 4 |
| Europe | 4 | 4 | 785.645 | 23,084 | 0 |
| North America | 6 | 3 | 305.542 | 10,568 | 2 |
| South America | 6 | 1 | 11.405 | 2,060 | 2 |
| Oceania | 6 | 3 | 36.983 | 7,144 | 1 |

Individual source rail segments and source original plant observations were deliberately not re-exported. Source code reads full original data privately within Github CI; aggregates-only artifacts, raw author plant × plot matrices and users' geotags are excluded.

## Source method warnings and unresolved biological inference

- v351 cell presence selection is a **global LAND-CENTRE grid probability sample**, NOT a worldwide railway centerline length-proportional sample. Better research of the railway-plant flora specifically needs a full source rail-line length sampling frame and probability sampling within actual physical railway right-of-way. Non-rail land cells are legitimate source zero-rail presence when original Overture scan complete.
- The original 2025 iNaturalist API bbox count is a public observation opportunity indicator; without exact original record GPS/error, distinct observers, plot effort and rail-local background pairs it is not even a temperature-bias outcome. It is also not standardized survey species richness; zero source records is not zero vegetation.
- Overture source snapshot is Nov 2026 and is 1yr newer than 2025 public plant reports. Rail centerline geography does not certify a physically accessible or train-operated ballast surface.
- Site WorldClim BIO1 is **macroclimate covariate**; warm-affinity biological trait is species-level T across global range, not the physical temperature of selected cells. The Harris et al. 2022 EIDC actual 200k+ species global tmean table source full ZIP was not fully downloaded due nonresumable EIDC transfer; **no real global plant warm affinity effect is reported or estimated**.
- The strong observational effort heterogeneity is a genuine positive *data-quality* result. It argues AGAINST opportunistically choosing high-observation cities and treating them as a global representative ecological test.
- Next test before any biology: within the 13 source-confirmed rail-positive cells, select original rail centerline segments with reproducible length-based sampling, then collect nearby vs 300–1000m local nonrail comparable plant data with original observation-effort and GPS accuracy, preferably fixed-area independent botany; obtain truly global temperature niche per plant species and correct spatial sampling bias. A rail-negative control site should stay in the frame rather than be replaced post hoc.

## Reproducible files
- Original global random source frame: data/derived/v351_global_land_climate_100km_probabilistic_source_verified.json
- Original genuine source first six rail: data/derived/v352_original_global_temperate_landcell_Overture_rail_first6_aggregate.json
- Original 32-source reconstructed aggregate from successful source-run logs: data/derived/v355_source_matched_global32_rail_plant_opportunity_verified.json (crosschecks source final real-run artifact totals)
- Regression validator: scripts/validate_v355_global_rail_plant_source_join.py
- Guard tests: tests/test_v355_source_verified_global32_join.py
- Scientific study designs: docs/tetsu_v351_global_equal_area_worldclim_land_frame_frozen_ja.md, docs/tetsu_v352_stratified_global_original_Overture_rail_availability_frozen_ja.md, docs/tetsu_v353_global_random_landcell_iNat_sampling_opportunity_frozen_ja.md, docs/tetsu_v354_full32_random_world_landcell_rail_source_frozen_ja.md.

DRAFT PR #57 stacked on PR #56 and #55. Main untouched, original author/plant rows and source map geometries not rehosted.
