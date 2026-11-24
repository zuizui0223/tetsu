# 3ヶ国の本物の鉄道網距離と植物相類似性：v364–v367真正データ統合（2026-11-24）

## 結論：物理鉄道網の存在・接続そのものは植物相の世界的均質化をまだ支持しない

**全球GISの全件取得は達成済み**。Overture 2026-11-18原データ128/128 Parquet source、2,054,855元source鉄道線分、2,892,198.162 km physical line centerlineをSHA検証し、通常鉄道候補1,536,226線分・2,237,043.481km・44,710真のGERS connector-connected componentsを構築。物理接続が巨大成分へ集中するが、これは運行実績／本当の営業路線名／列車の種子運搬ODではない。[PR #58 source-verified full global analysis](https://github.com/zuizui0223/tetsu/pull/58)。今回、現地原著の植物群集データを付けた**3独立国の研究別検証**まで進んだ。

## 独立した本物の植物ソースと物理鉄道接続情報

| 国 | 真正原著植物調査の単位 | 本物の元植物データ | 実GERS接続元調査地点 | source connected n dependent pairs | 実路線数の制約 |
|---|---|---|---:|---:|---|
| Slovenia, Šabić et al. 2026 | 100m field site, BS ballast shoulder / BT ballast toe | 89地点/15著者定義線区。2022–23 | 86/89 uniquely mapped、3曖昧 | 2,167/3916 connected undirected pairs <=300km | 原著15 source sectionsは本物のOSM named-line route関係と一致とは限らない |
| Czechia, Kutlvašr et al. 2024 | **whole railway station flora**, different surveyed areas | 39 author stations/761 contemporaneous taxa/6941 positive incidence、**35駅** source SR70代表点あり、4駅不明 | **35/35**／5 original author stretches | **568/595** true connected railway paths | 2026 third-party official-station **centroid**, NOT original botany GPS |
| Slovakia, Májeková & Limánek 2016 | **whole railway station flora** | **19 stations/318 original taxa/1763 positives**, author original DMS GPS, 2012–14 original field period | **19/19** originally printed DMS station point | **171/171** connected railway paths | **ONE original true Poprad–Plaveč line**, 19駅≠19独立路線 |

Czech original station centreline offsets中央値0.152m／max0.933mはSR70 derived mapsに由来し、原著野外植物調査地点の測量精度ではない。SK元駅DMS coordinate-to-GERS offsets中央値6.378m/max21.486m。Global GERS graph 2026と植物調査2012〜2023の非同期にも注意。元学名×GPS・OD等はrepoへ再配布せず、workflow本物元source内だけで集計。

## Slovenia：線路網距離＋直線距離を揃えた小さい正の記述的類似性、独立network効果は未支持

[Original SI WGS84 source metric v364 SUCCESS #45358570723](https://github.com/zuizui0223/tetsu/actions/runs/45358570723)：source path4,334 directed元original pairsについて、UTM投影距離→真のGERS WGS84 geodesic-at修正で元経路長中央値39.446m／max100.929m変化したが、元距離二重マッチ27 original anchors/7 author rail sectionsすべて同じ。BS+BT source Jaccard SAME author section−other +0.027033。元station site richnessとyear×SINGLE/DOUBLE内species source occupancy完全保持technical nullで残差 **+0.008484**、technical upper-rank **0.312**。厳密に同年／線路構成をsource matchingにも要求すると6元anchors/3 sections、元最低12/4不足で未識別。[SI source-original strict species-null #45345069976](https://github.com/zuizui0223/tetsu/actions/runs/45345069976)。

[Original SI 15 author-rail-section HELDOUT source prediction v365 #45359255096 SUCCESS](https://github.com/zuizui0223/tetsu/actions/runs/45359255096)。真の物理network linked2,167 unordered source field-site pairs、1,749 GERS経路無しはmissing。M0=geo log+square+same survey year/layout+source BS+BT plant richness、M1=M0+true rail network shortest path log、M2=M0+author same-section、M3=M2+network。λ=.05 fixed ridge / source whole-line 15-fold holdout。

| 元original railway habitat | M0 only geography/source covariates MSE | M1 + true network distance MSE | absolute M0−M1 MSE | direction |
|---|---:|---:|---:|---|
| Combined BS+BT | 0.003851991 | 0.003857378 | **−0.000005387** | **0.14% WORSE** |
| BS ballast shoulder | 0.005679457 | 0.005703874 | **−0.000024417** | **0.43% WORSE**、0/15 source folds improved |
| BT ballast toe | 0.004983032 | 0.004941408 | **+0.000041623** | 0.84% better, 10/15 source folds, but source omitted-richness and line-label sensitivities unstable |

したがって、**同じ線区内にやや種の共通性が見えても、実線路距離を加えるだけでは独立の留保線区の植物組成を広く改善できない**。線路の本物の物理接続はsource GISの事実だが植物種子運搬とは異なる。

## Czechia：元35駅・5線区、真のGERS距離まで揃えると元原著区間類似性ほぼゼロ

[**v367 source SUCCESS #45362326266 (Czech job)**](https://github.com/zuizui0223/tetsu/actions/runs/45362326266)。原著2024 source S1/S2とSR70 derivative source SHA照合し、source元地理距離のみの旧v282 matching 33 anchor/5 original stretches raw Jaccard SAME−OTHER**−0.0017311806**を完全再現した後、新規source matcher:
- source station proxy35/35 network geolocated、595 possible pairs→568 true network connected within300km。
- 元直線距離最大160km、same/other original source geographical gap≤10km、実GERS true route path差≤20km、1 source triad/anchor（植物結果を未使用）。
- **31 station anchors /5 original stretches**がsource両距離matching eligible；same vs other geographic mean44.675 vs45.697km、physical GERS network mean58.877 vs65.091km。
- **Jaccard SAME−OTHER +0.003225**、5 source stretches中2正方向、one-source-stretch-LOO gap範囲[−0.013329,+0.010495]。
- 全35 source駅に対してoriginal station richness+individual species exact occupancy保存2-chain×99 technical null mean**−0.002293**、observed−null**+0.005518**、upper technical rank**0.1608**。原5線区ごとの種頻度まで固定するOVERCONTROL null mean+0.026550、observed−null**−0.023325**；このnullは被験の線区種プールを先取りするのでprimaryではない。
- 本物のmap rail network距離を含めても**チェコの原著各線区間には頑健な種類似性の優位が観測されない**。元 station-area, surrounding urban habitat, Czech SR70 GIS centroid-to-field GPS offset等は未調整。

原original source auth：S1 SHA df6bc0b949e7f6a5cb57755dc8f9c86187223f66ddfbef959cf0e03ef7c46b06、S2 SHA 00c7cfdf5047854ba99d47dc8ec8c55b1b7df72731ce647c22b286ddef24a5a2、SR70 SHA fc43a1969af64f17285c7c6e3a30b00c3994d0ba1cfda6db09d8e4baaa166019。結果 [v367 CZ genuine original study aggregated JSON](../data/derived/v367_source_verified_Czech35_real_railnetwork_same_other_Jaccard.json)。

## Slovakia：本物の線路沿いの距離減衰はあるが地理距離を超える新規予測力はない

[**v367 source SUCCESS #45362326266 (Slovak job)**](https://github.com/zuizui0223/tetsu/actions/runs/45362326266)。元2016原著PDF249,373B SHA99465bdd7ed56985e65b1163f16eca4f702fc821619d84c89bf986d4d3736ae1。19原著作者DMS駅座標，318本物分類群、1763元駅×種positive marksと著者公表richness 5点原値をsourceチェック。真のGERS 2026 graph **全171/171 original station pair paths**。Network km range **1.1008–60.9542 km**、直線距離range **1.0721–46.7372 km**、network/geo length ratio中央値**1.1342**。

| 2016 SK original paired station Jaccard taxonomic distance decay | 原著グラフ上の線路距離(km) | WGS84直線距離(km) |
|---|---:|---:|
| original 19 station slope | **−0.001651745 /km** | **−0.001999674 /km** |
| exact original site richness and each species occupancy-conserved checkerboard source null expected slope | **−0.000251879 /km** | **−0.000313089 /km** |
| observed−technical null | **−0.001399866 /km** | **−0.001686586 /km** |
| 2-chain×199 technical null lower rank (NOT field p) | **0.002506** | **0.002506** |
| leave-ONE whole original station out slope direction | **19/19 negative**, source range [−0.001952,−0.001461] | existing old source shows 19/19 negative |

Source geodesic paper v310 previous original slope **−0.0020021839** exactly re-produced using published DMS+original Haversine, independent original WGS84 geodesic sensitivity −0.001999674. **True source line has clear botanical spatial turnover**; but **not network-specific**: ridge M0 `geography log+square` vs M1 `geography+GERS network log distance`, leave-one-entire-source-station-out 19 folds (153 training pairs/18 test pairs each), mean test MSE **0.003796289 → 0.003798021** (M0−M1=**−0.000001732**, 0.046% WORSE). Although 12/19 folds individually improved, mean prediction gets no independent lift. One source Poprad–Plaveč line; **no between-line category** so cannot claim cross-line railway dispersal.

### CRITICAL ORIGINAL DATA-AUDIT ERRATUM
Early v367 SK produced an apparently stronger network distance slope using erroneous Python `network.get(i,j)` rather than tuple-keyed actual `network.get((i,j))`; Python implicitly returned numeric station index j as the lookup's default. This was caught because the spurious source network-distance range (1–18km) implied physical paths shorter than geodesic endpoints. **All biological figures from that incorrect run invalidated and NOT SAVED**. Source was corrected with synthetic regression that station index and true GERS WGS84 path km must differ, bidirectional symmetry gate, original source 19/19,171/171 QA. Corrected true path 1.10–60.95km, median network/geo1.134, no invalid shortest-path below original true site chord once checked against source original GPS-to-track offsets. This report cites **only corrected genuine original SOURCE SUCCESS #45362326266**, and the before/after source incident is documented in frozen source design erratum [v367](tetsu_v367_CZ_SK_real_rail_network_floristic_similarity_frozen_ja.md). No artificial confirmation of network dispersal.

Original corrected SK study aggregate [v367 SK genuine original study JSON](../data/derived/v367_source_verified_Slovakia19_real_railnetwork_distance_decay_Jaccard.json)。

## 論文としての落とし所

自前全球全原本GERSネットワークは存在し、局所の植物群集に有意な地理的パターンはある。しかし、**地図の線路網自体の接続だけでは、独立した（地域外へ移せる）分類学的同化を予測できていない**という国横断（Slovenia/Czechia/Slovakia）に整合的な限定的結論が得られた。3国は完全に独立した標準化全球試験ではなく、**source field grain differences**が大きい。特にスロベニアはtrack ballast region、チェコ・スロバキアはstation-wide floraである。世界一般のnull network effectを証明する「非支持ではなく反証」言説は避ける。

**次の意外性のある仮説**: 線路網がpotentialな生物地理的移動路になっても、元沿線地形・生育適地・複数の接続方式・駅ごとの人為攪乱・種子繁殖機会というlocal establishment filterが効かなければ、flora identityは物理的接続だけでは均質化しない（Connectivity–establishment decoupling）。これはデータと整合する**新しい仮説**であり、確立した機構ではない。検証には複数国で同じ面積・同じ部位（BS/BTなど）・路線ID・year・local nonrail controls・種子散布形質・列車ODと実種子移動を組み合わせる。

第2の機能的な可能性: species identities are locally idiosyncratic while functional profiles (lifeforms, dispersal strategy) converge along rail management. Prior CzechVeg original within-publication guild heldout AUC lifeform+strategy .6035 vs source full taxon AUC .5051, but original CzechVeg rail/road labels text not confirmed active ballast and 1country observation; separate independently standardized functional convergence test needed.

No worldwide plant thermal-affinity result or train-carried seed deposition / adaptation mechanism proved.

Research Draft [PR #59](https://github.com/zuizui0223/tetsu/pull/59) stacked upon [PR #58](https://github.com/zuizui0223/tetsu/pull/58), no main merge. Genuine source original Excel/PDF botanical full data, station GPS, OSM source rail GIS geometries and per-source rare original taxonomic data not redistributed in repository.
