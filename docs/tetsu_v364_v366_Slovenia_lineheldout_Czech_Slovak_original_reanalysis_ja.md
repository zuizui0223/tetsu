# v364–v366 — 実線路網距離と植物相の独立地域検証（2026-11-24）

## v364：Overture 線路長の測地線単位修正

[真正原著 v364 CI成功 #45358570723](https://github.com/zuizui0223/tetsu/actions/runs/45358570723)。スロベニア Šabić et al.2026、原著89地点・15著者線路区間、実物Overture GERS 6,612 normal candidate rail segmentsの12,346 nodes /14,181 edges、86実線路に対応＋3曖昧。Overture source connectors[].at は元WGS84測地線長の割合。以前のUTM投影線路長とそのsource atをそろえ、各segmentと元地点の位置重みを測地線で修正した。真の対応graph 4,334 directed pathsの距離は絶対中央値**39.446m**、最大**100.929m**変わったが、元27地点・7区間の地理距離＋network距離二重matchは**27組全て不変**。原著BS+BT Jaccard差**+0.027033**と、元2022/23年×単線/複線の植物種頻度を保持した帰無残差**+0.008484**（技術的順位.312）は完全再現。生物学的に独立な鉄道種子散布の証拠は得られない。[原ソース出力](../data/derived/v364_source_verified_Slovenia_true_geodesic_network_reanalysis.json)。

## v365：真の線路網最短距離の独立予測力

[真正原著 v365 CI成功 #45359255096](https://github.com/zuizui0223/tetsu/actions/runs/45359255096)。同じスロベニア89地点のうち、本物の鉄道GERS graphで300km以内の経路が得られた**2,167独立でない地点ペア**を標本、残る1,749ペアはnetwork path不明であって生物学的ゼロではない。15個の著者定義原著線路区間を**丸ごと１区間ずつholdout**。訓練に留保した原著地点が入るペアは一切使わない。

モデルM0：元WGS84直線距離のlog＋二乗、同じ調査年／同じ単線複線構成、元各地点の植物種数mean/差。M1：M0＋本物のGERS鉄道網経路最短距離log。M2：M0＋原著同一区間ラベル。M3：M2＋鉄道網距離。固定 ridge λ=.05、同条件同test site-pairを15原著区間で比較。BS＋BT源 union種数をeffort共変量とし、元その種数を使わない感度も実施。

| 実source植生パネル | M0 test MSE | M1 test MSE | M0−M1 source delta | 正の改善が出た原著区間 |
|---|---:|---:|---:|---:|
| BS＋BT | 0.003851991 | 0.003857378 | **−0.000005387**（約0.14%悪化） | 9/15 |
| BS ballast shoulder | 0.005679457 | 0.005703874 | **−0.000024417**（約0.43%悪化） | **0/15** |
| BT ballast toe | 0.004983032 | 0.004941408 | **+0.000041623**（約0.84%改善） | 10/15 |

さらにsource同一区間ラベルを足したM2→M3の鉄道網追加gainはBS＋BT −0.000004932、BS −0.000031900、BT +0.000024428。元植物種数共変量を外すとM0→M1 gainはBS＋BT −0.000027293、BS −0.000016408、BT **+0.000006958**（小さくなる）；M2→M3ではBTは**−0.000012586**と悪化。したがって、**鉄道網経路が植物群集を距離や既知の管理proxy以上によく予測するという広い仮説は未支持**。BTに現れる微小な好転は微生息地依存性の探索的仮説であり、本物の列車種子散布・群集同化と断定してはいけない。原著89地点は一国一研究、2,167ペアは独立サンプルではない。

[元source集計JSON](../data/derived/v365_source_verified_Slovenia_true_GERS_network_distance_lineholdout_Jaccard_summary.json)、完全fold別元source artifact ID14904395436、[固定解析設計](tetsu_v365_true_rail_network_distance_leave_line_out_frozen_ja.md)。

## 既存の独立原著研究を発見：チェコ・スロバキア

- チェコ原著 Kutlvašr et al 2024 Preslia：**39本物の鉄道駅、761現生分類群、原著5線路区間**。先行PR #14では、公式station registryに基づく2026年の代表駅点（**植物調査plot原GPSではない**）35駅＋4駅不明を厳密に分離。33駅/5原著線区の地理距離対応で原著同一線区Jaccard優位は**−0.001731**、頻度帰無残差+0.003165、technical rank.357。頑健な同一路線の植物均質化なし。原著S1 79,261B SHA df6bc0b949e7f6a5cb57755dc8f9c86187223f66ddfbef959cf0e03ef7c46b06、修正S2 13,533B SHA 00c7cfdf5047854ba99d47dc8ec8c55b1b7df72731ce647c22b286ddef24a5a2、SR70代表点のoriginal SHA fc43a1969af64f17285c7c6e3a30b00c3994d0ba1cfda6db09d8e4baaa166019。
- スロバキア原著 Májeková & Limánek 2016：**19駅・318種、原著著者が発表した19駅DMS座標**、元source 1,763 station×taxa。本物の**1つ**のPoprad–Plaveč線路沿いの19駅の地理距離減衰−0.002002 Jaccard/km、種頻度/駅種数固定のsource帰無残差−0.001667 Jaccard/km。元原著PDF SHA 99465bdd7ed56985e65b1163f16eca4f702fc821619d84c89bf986d4d3736ae1。本物の独立営業路線は19本ではなく1本。
- 3ヶ国の研究は地点grain・source survey時期・GPS精度・駅面積が異なるので、地点ペアをまとめて全球群集メタ解析とは呼べない。

[v366 source設計](tetsu_v366_CZ_SK_independent_rail_station_true_graph_source_frozen_ja.md)では、原著CZ35駅代表点とSK19駅DMSを植物結果を見ずに同じ実Overture physical GERS graphへsource別対応し、鉄道からの距離・Ambiguity・本当のnetwork shortest-pathの有効地点ペア数を確認中。**まだ本物のCZ/SK network花相効果は算出せず、現研究は全世界の種子輸送因果証拠ではない**。original raw taxon×plot、site GPS、train rail geometryをGitHubの公開repoへ転載しない。

Research Draft [PR #59](https://github.com/zuizui0223/tetsu/pull/59) に保持、main未マージ。
