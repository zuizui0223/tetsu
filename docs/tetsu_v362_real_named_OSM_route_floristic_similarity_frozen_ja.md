# v362 — 本物の名前付きOSM鉄道路線関係を共有する植物調査地点は、線路網距離と直線距離を揃えても似ているか

2026-11-24。v360で元OSM route membership coverageを把握後、v362の生物学的結果を見ずに固定した**探索的事後設計**。本当に「路線単位」の群集類似性を検証するため、著者の15 `section`ラベルではなく、OSM `type=route`かつ`route=railway/tracks`で名前またはrefのある**source relation member way**で同一路線を定義する。旅客運行の`route=train`を混同せず、同路線関係が列車の種子運搬を意味すると主張しない。

## Source authenticity locked

1. 2026-10-01 static monthly Geofabrik Slovenia `slovenia-261001.osm.pbf`（313417178 bytes, published MD5 `703113982dc57e356f08dfccb1ddbf1a`; original source SHA256 `cce34d7bdfd310238dd302f50b9ad61c3ba9e66bc376856bf2cf56e002dbd9f7`）。本物v360 source [#45345176407](https://github.com/zuizui0223/tetsu/actions/runs/45345176407) で**73名前付きインフラroute relations**、元OSM way membership source 7981 distinct ways。86/89 botanist sites mapped physical Overture; OSM named relation **exactly one43、multiple42、none1、physical ambiguous3**。通常class source rail 6613 GERS segments中named relation match722（10.9%）。この不均一なsource被覆率をそのまま結果へ残し、個別の路線名で植物種を恣意的に選ばない。
2. 植生：Šabić et al 2026原著source original XLSX 401279 bytes、SHA256 `26be7ca9513dd54136d4ec499d8413c1cb288045f87fe3cf8814731f2edac99f`、89真实 site/15作者 rail section、植物BS/BT元種観測集合と2022–23 source年/単線複線。
3. 物理鉄道網：2026-11-18.0 Overture source `subtype=rail` standard/narrow/unknown、disused/abandoned/construction除外、bridge/tunnelは実接続を維持するため残す。元86 siteをGERS physical connectorにマッチ。25km等の距離条件は新たにチューニングせず、v359 mainの**40km以内、直線距離差≤2.5km、true rail connector graph経路長差≤5km、ネットワーク経路長≤300km**で固定。

## Two strictly separated source line identity definitions

- **Primary: unambiguous named-infrastructure-line**。原著候補anchor, same reference j, different ref kは3者ともOSM名前付きsource route relationが**ちょうど1つ**。same referenceは同一の一つのsource relation ID、other referenceは別の一つ。Source original line IDを別に自動推定せず、元source route provenanceを使用。Minimum 12 actual anchor botanical sites and >=4 distinct original anchor named OSM route relation IDs。未達は**NOT_IDENTIFIABLE**で、数値Jaccard差を公表しない。これは複数routeをもつ約42/86 sitesを除く保守的なサブセットで、大量のsource data利用不可能が起こりうる。
- **Sensitivity: genuine multi-membership any-overlap**。source named relations>=1の3者で、same refは元route ID集合のintersection非空、differentはdisjoint。元route-setごとのalias/route hierarchyを勝手に縮約せず、multi-membership route setを記録。>=12 anchorsかつ>=4 distinct original named **anchor route-set signatures**。元route間の論理的関係は未監査であるため、この感度分析は**同一営業路線と同義ではない**。元source particular nameの選択で肯定結果を得ることを禁止。
- いずれも候補j,kは植物の種名・類似度を使わず、辞書式の`(|Δgeo|,|Δnetwork|,max_geo,source site index)`で1トリプル/anchorに固定。計算source部位主BS+BT、感度BS/BT別。
- 両路線の元original年、線路layoutまで一致するstrict sensitivity。同一次元の標本支持不足時は未識別。
- 類似度primary taxonomic Jaccard。各独立アンカーの繰り返しを独立nとせず、元route relationまたはanchor route-set signature等重み。地理・ネットワーク距離の平均とsource match support、source source rankを報告。距離差のsource caliperを満たさない比較を0と数えない。
- Source null control: 元89 botanical sitesの**種出現頻度とsite種数を完全保存**した2×99 correlated checkerboard MCMC、上記のsource当該same/different named-route comparison trioを固定して評価。主推論を支えるのは更に調査年×SINGLE/DOUBLE layout内の元種頻度を保つnull。Match support不足でnullを走らせない。Technical rankはfield p値でない。検定は事後探索。
- Overture sourceでcoverageが低いのは地図source availability問題であり、この1国サンプルを全球全鉄道ネットワークとして結論しない。元OSM PBFと原著サイト位置、学名、source路線名・ID・wayID・完全rail geometryをGitHubへ転載しない。

## Source success criteria

Original v360 73 named-relations / 86 mapped sites / 43 unique+42 multi+1 none /3 physical ambiguousをexact source pinで再現。元v359 27-anchor 7section source comparison and original v28 59anchor baseline stay separate, because route membership defines a genuinely different **estimand**. If source match sparse, output supports and scope only, no biological difference. Report independently against authored-section result v359b which shrank from +0.0270 to residual +0.00848 under source year×layout species-frequency preserving null.

Files `scripts/analyze_Slovenia_true_OSM_named_rail_route_floristic_similarity_v362.py`, synthetic tests, GitHub source flow. Research Draft PR #58 main unchanged.
