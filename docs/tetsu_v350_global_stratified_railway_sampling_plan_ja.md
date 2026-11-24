# v350 — 鉄道植物の全球標本抽出設計：原著選別から確率抽出へ

## 既存研究7つを主母集団にしない理由
これまでのCZ/IT/SI/HRBA/SK/Kashmir/WA研究は生データが公開された順の convenience sample である。地域数・環境数・鉄道の運行状態・調査面積・参照の自然生息地を網羅する抽出枠ではない。同研究の駅数も独立大陸数ではなく、多くは欧州の研究。**その7研究の平均や符号数をglobal railway thermal effectへ変換してはいけない。**

## 自前一次データで作る二つの異なるデータベース

**A. 全球 railway-proximity species detection**（個体/位置観測、調査区植生ではない）
- Overture Maps 2026-11-18.0 official global transportation `type=segment`、`subtype=rail`、`class=standard_gauge/narrow_gauge/unknown (unknown is a sensitivity only)`。RailFlags `is_tunnel`, `is_bridge`, `is_abandoned`, `is_disused`, `is_under_construction`, `is_covered` は保守的に除外。Source provenance: [official transport docs](https://docs.overturemaps.org/guides/transportation/) / [global S3 source](https://docs.overturemaps.org/getting-data/cloud-sources/)。
- 全球の個々のmapped rail featuresを国・気候帯・都市化度の標本抽出枠に入れる。空間独立単位は**研究者が選ぶ駅**ではなく、全国レール線分を10km equal-area gridへ入れた線分集落/路線block。長い路線は長さに比例した標本選択確率、非常に短いマップ分割segmentを独立生物学的反復にしない。Overture dataset全体を取得可能/選択支持を監査。
- サンプル候補は continent(6) × Köppen major climate class or pre-frozen macro-temperature bins (e.g. cool/temperate/warm) × urbanization strata(urban/rural) で事前層化し、全カテゴリを含む地域を**温度効果を見ずに決定論的SHA seedで抽出**。各層のrailway length support、採否確率、未観測地数を記録し、世界の鉄道網全体へのデザインベース重みを使う。海洋島・線路の存在しない気候層は構造的NA。最低国数/独立路線数の閾値は事前凍結。
- 各元source Overtureの選定segment midpoint近傍と、同じlocale/year/season/landcover/道路距離/人口密度で比較可能な300–1000m離れたreference windowを定義。Source Overture rail lines 2026は古い観測との時間一致を保証しない。rail近傍記録は鉄道権利区域内植物とは限らない。
- Source iNaturalist/GBIF public 2018–2025 wild plant species exact georeferenced record with original positional accuracy and eventDate, exclude nonpublic/obscured, cultivated, coordinate centroid/institution, camera duplicates and incompatible dates. 1 data provider at a time; iNat→GBIF重複に注意。基準は `distance_to_rail + source accuracy ≤ 80m`、reference `distance-accuracy ≥300m` and `distance+accuracy ≤1000m`。元日付や観察者ごとのsource sampling effortをフィールドに保持するが元ID/GPSをGitHubへ転載しない。種の不在は観測できないので、見かけの頻度をpresence/absence occupancyへすり替えない。
- **標準化推定対象**：同じ original observer×season×10km source geographic block内で記録された植物のT_scoreとrail-proximity labelの関連。source observation‐effort model / target group background を含む multi-level logit or retrospective point-process model。Clustering whole source railway route / grid / original observer, not individual photos counted as independent ecosystems. 
- **主global ecological hypothesis** (analysis following data capture, not yet tested): `P(near rail | plant detected, matched sampling effort, local macroclimate)` vs independent species global T trait, with `species thermal preference × local macrotemperature` interaction. Potential warm-affinity positive in cold rail regions, but spatial selection vs preferential plant observation cannot be fully separated. Negative-control pseudo-rail lines within similar urban landcover and matching other linear corridors are mandatory. Use original global plant climatic affinity Harris et al 2022 EIDC species `tmean` if fully obtained, never European EIVE in tropics.
- **Claim level**: globally stratified *proximity observation association* ONLY. Even significant results cannot prove train seed transport, actual ballast hot microclimates, within-species adaptation, causal filtering.

**B. 真正の global railway/adjacent-nonrail vegetation study database**（論文で主張する機構）
- Independently sampled fixed-area vegetation relevés within operating railway ballast or documented railway land AND geographically paired nonrail natural/analogous substrate; original measured site year, area, method, raw plant×plot list, rail ownership, disturbance/mowing, macroclimate. Separate actually abandoned/active rail. Only these support habitat-level community mean difference (thermal position/width, taxonomic+functional composition). Collect from papers and open vegetation archives beyond seven original convenience studies; report excluded studies and bias.
- Global representativity requires >1 climatic region/biome independently, not country-level cross-reference to unrelated city-tree-pits or MIREN roads. WorldClim air temp not actual measured trackbed soil microclimate. Any global railway-specific causal filtering claim remains conditional on measured rail microclimate/experimental controls.

## v348 first empirical API opportunity; why it is not yet a representative global study
One **urban** window selected independently of species identity in each of Tokyo, Prague, Cape Town, Sydney, San Francisco and São Paulo (year 2025, +/-0.012°). True public iNaturalist source API returned all 6 counts [170,536,96,89,133,17] total **1,041** source observations BEFORE geoprivacy/accuracy/rail-buffer exclusion; [original CI #45328376103 SUCCESS](https://github.com/zuizui0223/tetsu/actions/runs/45328376103). These six hub windows are fixed and not a global random sample, nor do 1,041 occurrence records form 1,041 independent plant communities.
- Source rail geometry from OSM Overpass/official OSM map and independent Overture 2026 transport source audited by workflows v348 and v349. If an actual source request fails, report unavailable rather than 0 plants.
- v348b geographic source availability establishes technical feasibility only, not a temperature effect, plant colonization, species diversity or comparison of actual railway vegetation.

## Prioritised go/no-go tests
1. Need complete genuine near and background *spatial source cohort* in ≥4 continental anchor windows with ≥10 near individual plant observations and ≥30 background observations/5 and 10 different species, respectively. If fail, do not fit thermal response from opportunistic citizen observations; broaden predeclared representative rail-grid sample independent of favourable signs.
2. Need independently globally calibrated species T and verified source-match coverage >=50% for each within-window group; otherwise global warm-affinity remains **not calculated**.
3. Need spatial sampling intensity and local habitat confound diagnostics + equivalent roadside/pseudo-line negative controls. Any pilot result cannot carry an intercontinental causal claim.

## Sources
- [Overture global rail segment class and rail flags](https://docs.overturemaps.org/guides/transportation/segments-and-connectors/) (official latest release Nov 18 2026)
- [Geofabrik OSM regional extracts](https://download.geofabrik.de/) for locally complete PBF alternative.
- [iNaturalist official API usage limits](https://www.inaturalist.org/pages/api+recommended+practices) (200 records/page, ~1 call/s, 10,000 result page cap; bulk must use download workflow).
- [GBIF API paging and downloads](https://techdocs.gbif.org/en/openapi/).
- [Global EIDC plant climate envelope](https://doi.org/10.5285/ca339c86-3674-4030-b891-35326e71141e); its actual full source CSV fetch had technical HTTP timeout, so no global species T CWM has yet been observed.

This is a frozen **future sampling design**, not a completed six-continent macroecological result. Main unchanged; new Draft PR #56 stacked on research PR #55. No individual raw species observation IDs, geotags or rail geometry republished.
