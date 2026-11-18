# tetsu v27.7–v27.9 — 鉄道植物群集の局所性は、数kmの地理的近接だけで説明できるか

更新：2026-11-18 JST。原著者真正補足データに基づく探索的な再解析。mainは未変更。成果は研究ブランチ research/tetsu-route-distance-local-pool-v277-20261118 に限定する。元植物種×個別地点、元GPS、元路線距離標と学名リストは派生JSONへ再配布しない。

## 1. 一次資料と本来の独立単位

Šabić et al. (2026) Slovenia RAIL original published plant source: https://doi.org/10.1007/s44473-026-00155-z 。元author supplement XLSX 401,279 bytes SHA256 26be7ca9513dd54136d4ec499d8413c1cb288045f87fe3cf8814731f2edac99f 。2022–23年、真正89個体群現地地点・173同側BS道床肩/BT法尻ペア、15実国有線路線区間。元地点の左右のtaxaを同じ現地地点内でunionしてから比較する。14路線区間・88地点で同国有線内に別地点があり、残り1路線区間1地点では別地点比較を作らない。

作者元の区間内100m transectのstart_km / end_kmの中点を鉄道距離標として使用、start_lat/lon・end_lat/lonの中点を実緯度経度として使用。地球表面の2点間距離はhaversineで計算。これらは実列車の出発→到着のサービスグラフや貨物ODを含まない。元SOURCE89地点を14路線に整理すると異地点組合せは340 unordered / 680 directed site pairs。うち10 directed source pairsでは同線区内距離標の差が直線距離より0.5km以上小さく、実際の軌道走行距離とroute km差が完全に一致しないことをSource geometry caveatとして明示する。作者の空間地点配置をゼロ距離の実輸送と解釈しない。

実データGitHub Actions:
- v27.7 original source, mathematical 5/5 SUCCESS https://github.com/zuizui0223/tetsu/actions/runs/44256989504
- v27.8 author original 5–50km source frequency fixed species exchange, mathematical 5/5 SUCCESS https://github.com/zuizui0223/tetsu/actions/runs/44257548836
- v27.9 original same-source route vs geographic exact matched target habitats, mathematical 6/6 SUCCESS https://github.com/zuizui0223/tetsu/actions/runs/44258051883

全Actionsで元source XLSX bytesとSHAを照合。独立野外研究3件ではなく、同じSlovenia鉄道原著1件からの3段階仕様の再解析である。

## 2. v27.7 本物の最寄り別地点を参照した初回の空間的検証

各BS_i同地点のBT_iとのJaccard類似度（同site）と、元同じ国有線内で距離標の絶対差が最小の別地点BT_jとの類似度を比べた。jの選択は植物種共有を見ず、距離標だけで決定。元サイトを1票、区間を1票で平均し、他の線区や複数年のsourceサイトを偽造しない。

| 最寄りBT別地点の許容距離 | 元線区・本物のBS元地点 | 同地点BS-BT類似度 | 最寄り別地点のBT類似度 | 差 | 正の元路線区間 |
|---|---|---:|---:|---:|---:|
| 5km以内 | 11区間/60site | 0.353266 | 0.205078 | **+0.148188** | 11/11 |
| 10km以内 | 11区間/77site | 0.352708 | 0.217333 | **+0.135375** | 11/11 |
| 20km以内 | 13区間/83site | 0.342639 | 0.218101 | **+0.124538** | 13/13 |
| 距離上限なし | 14区間/88site | 0.363236 | 0.217983 | **+0.145253** | 14/14 |

5km内で選ばれた最寄り別地点のoriginal source line-equal中位距離は約2.30km、10kmでは2.70km、全88地点では3.80km。同年・同SINGLE/DOUBLE track typeの候補だけに限定した感度でも、同方向。v25.7 source exact original same-site Jaccard類似度0.36323635は元SOURCEから再現。距離0と距離2.3kmに差があるので、地理的近接・土壌/地形・検出差・地域種供給の区別はまだできない。

## 3. v27.8：5km内だけでBT植物種を交換する厳密な元source incidence null

BS本物のtaxon集合を完全固定し、BTの種出現の2×2checkerboard swapsをsourceの元同一路線区間の距離標差5/10/20/50km以内に限定。BT各元siteの種数およびBT各種の同source線区内occupancy frequencyは厳密に保存。遠隔路線区間へ元学名や植物種群集を飛ばさない。99 MCMC source drawsは連続した相関サンプルであり野外ランダム化試験p値でない。

5kmでは11元線区・81本物のsource地点が対象線区内にあるが、**これら81site全部に5kmの交換相手がいるわけではない**。元source11線区内に交換可能な原siteペアが52個、v27.7で5kmの最寄り別地点があるのは60地点。交差可能なsite以外は固定のままsource matrix marginsには残る。このためv27.8の5km効果は保守的に薄められている。全距離範囲で同じ元11線区・81siteを使った厳密比較：

| 交換相手の路線距離上限 | 実測βjac | 元種数+種頻度完全固定のnull βjac | null−observed局所共有超過 | 正の源線区 |
|---|---:|---:|---:|---|
| 5 km | 0.645700 | 0.718592 | **+0.072892** | 10/11 |
| 10 km | 同じ0.645700 | 0.746158 | **+0.100458** | 11/11 |
| 20 km | 同じ | 0.757540 | **+0.111840** | 11/11 |
| 50 km | 同じ | 0.770075 | **+0.124375** | 11/11 |
| 同源路線全体 | 同じ | 0.768581 | **+0.122880** | 11/11 |

同source初回調査年のBT交換にさらに制限しても5kmのexcess +0.070634（10/11）。ここでは、遠くのsource地点へBT種を入れ替えたから見かけ上の局所共有が出た、という説明だけでは足りない。なお元11線区の1つで5kmの差が帰無と逆符号なので、全面的な普遍則にしない。

## 4. v27.9：鉄道路線距離と直線距離のどちらが局所共有を説明するか

v27.8の「交換可能な5km内の対象でないsiteも集約に残る」問題を解消。対象を**source距離標5/10/20kmの近隣が存在し、同時に元GPS直線距離5/10/20kmの近隣も存在する同一元site**に先に限定。元線区内のBT全体の種数×taxon頻度は保存したまま、路線距離ベース・直線距離ベースの二つの別MCMC交換グラフを使い、**完全に同じ元siteのobserved βjac**へ対してそれぞれ帰無を評価する。これは因果的な「線路を通った種子」と「道路/風で直進した種子」の実験ではなく、同一source field location supportによる幾何学的対照である。

| 近接候補範囲 | 共通源線区／site | 実測site source βjac | 路線距離条件null−実測 | 直線距離条件null−実測 | 路線−直線の超過差 |
|---|---|---:|---:|---:|---:|
| 5km | 11線区／60site | 0.646734 | **+0.093658** | **+0.099306** | **−0.005648** |
| 10km | 11線区／77site | 0.647292 | +0.102945 | +0.108623 | −0.005677 |
| 20km | 13線区／83site | 0.657361 | +0.100736 | +0.105550 | −0.004813 |

5kmに元路線沿いの候補ペア52組、直線距離の候補ペア59組。5km源Siteの路線・直線の双方で正のlocal excessを示すのは10/11源線区であり、10kmでは11/11、20kmでは13/13。元同survey yearに制限しても5kmのroute +0.090050、geo +0.099442、両者差−0.009392。v27.7の同一5km源site集合60地点での同地点Jaccard類似度0.353266と今回の1−βjac=0.353266は完全に一致。

**科学的判定**：
- たとえ同じ元国有線上の5km内に参照を制限しても、その現地のBS/BTは元種数・source種頻度の期待以上に共有する分類群を持つ。
- 地図上の直線距離を使っても同方向・同程度のlocal effectとなり、距離標の方が**特有に強い**結果にはならなかった。現段階で「鉄道サービスネットワークの接続性が局所共有を作った」とは支持できない。
- どちらの幾何学でも同site vs近隣siteは環境・管理・採取条件・地元の種子供給を分けていない。原BS/BT双方が鉄道内の微生息地で、非鉄道の同面積/同観測努力の強い環境対照ではない。元sourceのBS/BT survey areasの正確な等価性も未保証。鉄道路線距離標の絶対差は複線分岐/実貨物ODのグラフ距離そのものではない。

## 5. 投稿可能な研究主張と次段の必要条件

SOURCEで反復するのは「共通する人工基質に生育しながらも、近隣数kmの土地固有の植生・プール/環境と結合した群集の形成」である。分類上の種の置換（先行v26.9のSloveniaのβjac）と同地点における強い共有は両立している。これは単純な「鉄道が全球的に全種を均質化する」より、複数のスケールを分けた群集形成モデルと整合する。

一方、実列車本数・貨物OD・有効種子の実運搬・管理/除草・風力・送粉・個体の繁殖/生存・遺伝選択は元資料にない。前v27.5/27.6の本物の2015/2020貨物強度は、源9物理区間で局所共有の負の相関を見せるが48路線内貨物ラベル交換に耐える固有の予測差ではなかった。今回の徒歩/車両/風の輸送起点推定には使わない。

次の確証的な問いを短く定義すると「局所種プールと気候・現地微環境を固定した後、*鉄道網の実接続および実交通量*は植物種共有に独立の説明力を持つか」。そのためには同地域・同時期・同面積の真正鉄道vs道路/荒地/自然対照、独立複数路線地域、実列車OD/線区別通過量、原source植物種出現と撮影/検出条件、必要なら現場種子採取が要る。現在のv27.9だけでは鉄道網固有の種子運搬効果の検出に成功していない。

再現ファイル：
- scripts/analyze_Slovenia_rail_route_nearest_plant_microhabitat_v277.py；tests/test_Slovenia_nearest_route_km_local_flora_v277.py
- scripts/audit_Slovenia_rail_nearby_route_distance_fixed_frequency_null_v278.py；tests/test_Slovenia_route_local_species_frequency_v278.py
- scripts/audit_Slovenia_rail_route_vs_geo_distance_source_null_v279.py；tests/test_Slovenia_route_geo_distance_reference_v279.py
- Corresponding original author SHA-gated GitHub Actions workflows tetsu-v277, v278 and v279.
