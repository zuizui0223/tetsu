# tetsu v27.4–v27.6：スロベニア鉄道の局所種共有、幹線／地域線と実列車本数の検証

更新2026-11-18 JST。真正一次資料の検証を通過した探索的な二次解析。mainは未変更。元個別現地地点×種のリスト、元GPS、元公式官庁資料の全文をリポジトリに再配布していない。

## 1. データ源・観測の独立性

- 植物：Šabić et al. 2026（DOI https://doi.org/10.1007/s44473-026-00155-z）、2022/23年の89真正現地地点・173同側BS道床肩/BT法尻ペア・15原国有線区間。元作者XLSX 401,279 bytes、SHA256 26be7ca9513dd54136d4ec499d8413c1cb288045f87fe3cf8814731f2edac99f。左右は元同地点内で先に統合し、元173左右ペアを独立地点としない。
- 幹線G／地域線R：スロベニア政府2024年改訂の国家鉄道計画7/36頁、https://www.gov.si/assets/organi-v-sestavi/URSZR/Datoteke/Drzavni-nacrti/dn_zelezniska_azuriran_2024.pdf 。原著国有線15区間（G:10,20,41,44,50,64、R:21,31,32,34,42,70,72,80,81）。G6線区34現地点、R9線区55地点。これは法定カテゴリであり交通量でない。
- 2023年インフラ公式駅距離標（Network Statement Annex2B、真正16頁PDF）https://infrastruktura.sz.si/wp-content/uploads/2023/02/NS_2023_1.2_Annex_2B_Stations.pdf、SHA256 83f12dfdb91784a653b9a39a8d927a027e7eeb5d8af18ce17b850b5318796425。
- 実列車本数：Eurostat公式方向別2015 https://ec.europa.eu/eurostat/api/dissemination/statistics/1.0/data/rail_tf_ns15_si?lang=EN （source SHA256 44edbd7a2fc857fe363612ce44c25de6c77d4f83c9d915b3455846c548efe2c3）、2020 https://ec.europa.eu/eurostat/api/dissemination/statistics/1.0/data/rail_tf_ns20_si?lang=EN（SHA256 05bc92c3555a47947f4acb8cf5aa1ebfe72607487e81713177570cc0b9d5114d）。元国有線番号、両端実駅距離標、方向別区間コードと両端駅ラベルを全て合わせて、逆方向の実本数を二方向ともある場合だけ加算。

元89植物地点のうち、この実交通区間に厳格に対応するのは27地点／9独立物理区間／5国有線区間。62未対応地点の列車本数を0としない。2022/23植物と2015/2020列車本数には年差があり、これらは同時の曝露を測っていない。

真正資料再取得・SHAと数学テストに成功：
- v27.4 source SUCCESS https://github.com/zuizui0223/tetsu/actions/runs/44254444702 （4件のQA試験成功）
- v27.5 source SUCCESS https://github.com/zuizui0223/tetsu/actions/runs/44255016472 （5件成功）
- v27.6 source SUCCESS https://github.com/zuizui0223/tetsu/actions/runs/44255375169 （4件成功）

## 2. v27.4：法定G/Rを調べるが、単複線交絡を分離する

道床肩BSの真の種リストを固定し、同じ元国有線内の法尻BT種リストだけcheckerboard交換する。各BT実地点の種数と元国有線内の各BT種の出現地点数を完全固定。元地点のsource Jaccard βの帰無平均との差を、局所的な種共有の超過量 local_excess = null β−observed β として定義。正ならBS/BTは同元路線の種頻度から期待以上に種を共有。各地点をまず元国有線の中で等重みとし、線区等重みの平均でG/Rを対比した。これは実地点ラベルを無作為操作した効果ではない。

| genuine source cohort | G(main) | R(regional) | G−R |
|---|---:|---:|---:|
| 元14線区・88地点、G6線区とR8線区 | +0.110093 | +0.125253 | **−0.015160** |
| 元の単線のみでBT species-frequency nullを再計算、12線区・64地点、G4線区とR8線区 | +0.105823 | +0.125378 | **−0.019554** |

元の線区を1つずつ除いても全線区差は −0.02788〜−0.00314、単線だけでも −0.04167〜−0.00698と負方向を保持。ただし路線選択・地形・地域植物相・車両交通と管理が交絡する。Gには複線が元24地点/2線区集中、Rに複線はない。複線のG/R比較は元R複線支持域0で識別不能。線区数の少なさゆえ独立国研究CIや貨物本数効果とは呼ばない。

## 3. v27.5：2020実列車本数と種頻度を固定した局所種共有の超過

最初に全89元植物地点のもとの国有線内で species-incidence reference nullを作成してsite別 local_excess を出し、その後27site/9真交通区間のみを抽出した。源site植物βは国有線内のBT種数・各種のBT出現頻度を固定して調整。列車量との相関は真9物理区間に集約した。G/R交通区間の支持はG7/R2、元5国有線のうち線内に異なる交通区間をもつのは2線だけ。

| source 2020交通量と局所種共有超過 | 原値 |
|---|---:|
| 実貨物本数との9区間Spearman rho | **−0.316667** |
| 実旅客本数とのrho | **−0.583333** |
| 国有線内でdemeanしたlog(1+貨物列車数)の局所超過への傾き | **−0.033879** |
| 国有線固定＋旅客列車数調整の貨物傾き | **−0.060724** |

元5国有線内の貨物本数ラベルを完全48通り入れ替えたSource診断では、無調整貨物傾きの絶対値が観測以上となる割当は16/48、旅客調整では10/48。この比率は実地無作為化のp値ではなく、路線固有の交通量差が観測固有の明瞭な信号かを問い直す技術対照。強い生態機構を肯定しない。

## 4. v27.6：2015/2020の真9区間は完全一致、植物応答は同一

2015元Eurostat directionalコード＋元ラベルが2020の同一区間と完全一致する場合だけ2015交通量を使用し、真9物理区間・27植物地点・5国有線を固定できた。元植物local_excessは同じ source であり2015と2020を独立生態学試験として数えない。2015/2020元交通9区間の順位Spearmanは貨物0.866667、旅客0.933333。

| source歴史的交通量の年 | 貨物rho | 旅客rho | within line貨物傾き | 48割当中同等以上の絶対傾き | 旅客調整傾き | 同等以上 |
|---|---:|---:|---:|---:|---:|---:|
| **2015** | −0.433333 | −0.716667 | −0.022429 | **24/48** | −0.037515 | **21/48** |
| **2020** | −0.316667 | −0.583333 | −0.033879 | **16/48** | −0.060724 | **10/48** |

元2015貨物列車の年間両方向本数514〜27,157、2020年262〜28,639。両年で負方向を再現するが、路線内の貨物ラベル再配置で同程度の傾きが頻繁に生じる。貨物列車の実輸送が局所種共有を減らした証拠と断定するには不十分。2020年のCOVID外乱、地理、旅客列車、単複線、管理が交絡。元2015・2020とも植生の現地観察2022/23と同時期ではない。

## 5. 論文としての主張・停止条件

新たに支持された記述的パターンは「幹線Gは地域線Rより局所種共有超過が小さい傾向があり、元単線だけでも同方向」「2015/2020の実貨物・旅客本数が多い9交通区間では、局所超過が小さい負の相関」。しかし、**同路線内に貨物量を割当て直すと観測と同等の効果量が容易に生じる**ため、路線内で列車数が特異的な局所植物相の共通性を弱めるという因果仮説は現状未支持。

元Eurostat区間別実貨物本数は物理区間を通過した本数であり、積載した植物種子の採取・実際の出発→到着ODや生物の流入経路を表さない。線路接続性のグラフ指標・同時年列車本数・同地域同方法の非鉄道植物群集・実草刈り・風ストレス・遺伝変異・個体生存を今回直接得ていない。強い因果同定を避ける。次に必要なのは2022/23年と整合する真正区間別交通、真正origin–destinationと外来種側のtransportable pool、別国・別研究の同粒度の鉄道生息地対照。

実装と数学QA:
- scripts/analyze_Slovenia_official_main_regional_local_coupling_v274.py、tests/test_Slovenia_official_GR_local_species_coupling_v274.py
- scripts/analyze_Slovenia_historical_train_frequency_local_species_coupling_v275.py、tests/test_Slovenia_real_train_freight_local_species_coupling_v275.py
- scripts/audit_Slovenia_2015_2020_real_train_local_coupling_temporal_v276.py、tests/test_Slovenia_historical_2015_2020_local_source_coupling_v276.py

対応するGitHub Actions source XLSX/PDF/Eurostat JSONを検証し、元taxon×地点・現地GPS個別資料は新たに同梱していない。mainは未更新、Draft PRのみ。
