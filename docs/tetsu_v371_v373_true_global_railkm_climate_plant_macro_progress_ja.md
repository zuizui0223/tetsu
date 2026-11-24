# v371–v373：実際の全球鉄道網を母集団とする植物マクロ解析と鉄道の気候帯偏り（2026-11-24）

## 研究の主解析対象を全球へ戻す

ローカルなスロベニア・チェコ・カシミールの花相比較を追加し続けるのではなく、**世界中の真正鉄道網と、別途確立した全球植物記録の確率・観測モデル**を本筋にする。鉄道全セグメントの地理的輸送sourceと植物観測sourceを結合し、区域内の記録植物種数および隣接区域の分類群Jaccardを、観測努力・地理・気候の下で説明できるかを大陸を丸ごと伏せて確認する。**全球鉄道敷そのものの定面積植生群集を測ったものではない**。

### 検証済み全球鉄道母集団

- Overture Maps **2026-11-18.0** の全128/128元Parquet（元全rail subtypeを含む2,054,855 GERS線分・2,892,198km）から、元鉄道class standard/narrow/unknownかつsourceのabandoned/disused/constructionを除いた、**1,536,226通常鉄道候補セグメント、WGS84測地線2,237,043.481km**を確立。地図中心線であって営業線、列車運行本数、実バラスト生育面積ではない。
- [v368b元全球rail-km長さ比例サンプル 2,048 点 SUCCESS #45373590287](https://github.com/zuizui0223/tetsu/actions/runs/45373590287)：世界の鉄道延長に沿い、植物記録を見ずに地点確率抽出、同じ区間の複数抽出もsource設計どおり保持。
- [v370元全球100km/200km鉄道km全セル走査 SUCCESS #45373985228](https://github.com/zuizui0223/tetsu/actions/runs/45373985228)：全128元source・source WGS84線分をEPSG6933へ地理クリップし、全球の真の線路中心線kmをsource格子に配分した。原人口計と残差を照合。全市町村・独立の個別列車routeではない。

## 新たな真正データの肯定的「地理的事実」：全球鉄道の分布は温帯に集中

[**v373 原鉄道2,237,043km×WorldClim2.1×Natural Earth 原source SUCCESS #45378280660**](https://github.com/zuizui0223/tetsu/actions/runs/45378280660)。**8/8合成QAテスト成功、元EIVEや元植物観測ではなく真の物理鉄道位置の環境基線を測定**。真正WorldClim 2.1 BIO1（1970–2000の平均気温）ZIP SHA256 `00513224583665ec0f2f955a4ec252730c4deb2004cce9e793492a3f26df4dcf`、Natural Earth 110m原map-unit source SHA256 `718ac51d68ed2bf5de4f563fff2e2e045db706c1d3222805f3e317140a64213b`。

100km等面積格子の**気候が与えられる陸地中心**は13,465セル。そのうち鉄道kmのあるセルは5,650。全全球鉄道源2,237,043.488kmのうち**2,146,220.175km（95.940%）**を陸地中心にWorldClim気温のあるcellへ分類。残り90,823km（4.060%）は海・島・沿岸のcell-center問題などでmacroclimate分類外。**海岸/島の鉄道が無いという意味ではない**。

| BIO1気候帯 | 100km陸地セルの構成比 | 100km気候同定された鉄道kmの構成比 | 鉄道km（元WGS84 geodesic） |
|---|---:|---:|---:|
| cool: T<8°C | **34.304%** | **29.871%** | 641,091.688km |
| temperate: 8°C<=T<20°C | **22.020%** | **54.200%** | 1,163,248.470km |
| warm: T>=20°C | **43.676%** | **15.929%** | 341,880.018km |

100km源鉄道延長による加重気候平均は**11.9595°C**、全陸地中心セルの等重み気候平均は**13.0617°C**。**鉄道のマクロ気候は、世界中の陸地の気候中心と同分布ではない**。200kmでも鉄道延長加重11.6408°C、陸地中心13.0989°Cで同方向、ただし200kmの気候分類支持は鉄道長の93.57%に減る。局所線路敷の気温や暖地性植物を観測したのではない。

**注意**：v373のEPSG6933気候分布の全球原グリッド外接範囲は、従来のv351確率抽出時の元land-center frameと完全同一ではない（v373 100km中心source 51,504に対しv351固定frame 49,640）。したがって13,465 vs 13,437を同じ母集団として統合/差分解釈してはいけない。全球格子を整合させたbio-rail 2025 source macro解析だけを比較する。

**この方向が何を意味するか**：もし鉄道沿いで暖地性植物が多いのかを世界で検証するなら、鉄道建設の温帯偏重および大陸・気候・土地利用・地理的植物供給源の偏りを、まずSource sampling frame段階で調整しなければならない。**気候による植物の選別、微気候的な生育可能性の拡大、列車による種子輸送の証拠ではない**。

## 2026-11-24 追記：v369 真正全球植物データ全件抽出SUCCESS

[**v369 original iNaturalist full world source SUCCESS #45375453375**](https://github.com/zuizui0223/tetsu/actions/runs/45375453375)、源観測gzip 13,697,505,636 bytes（source ETag `4539360677c8190d9aa12dbf6982ee72-27`、完全SHA256 `1c48a9646f3bf6d56ffdfdc91a1a176c173c98014830ad48fa162b2803c3b279`）、源taxa 40,184,871 bytes（ETag `0a6a162c57bb2f24ea4f00c76c20080d`）を元S3で全件検証済み。

| source真正集計 | 完全原sourceの実結果 |
|---|---:|
| 原世界2026-10-27 photo-linked観測レコード全行 | **288,721,938** |
| 原source分類体系の全taxa | **1,675,070** |
| original active Plantae species taxa | **284,070** |
| 2025年・source研究用／座標条件後の公開植物観測 | **7,827,936** |
| 100km等面積格子に2025年植物が1件以上 | **9,029** |
| 100km格子：≥50観測/≥5観測者/≥20植物種 | **4,219** |
| 100km格子：同支持条件の真正隣接2格子Jaccard | **5,944** |
| 200km格子：source植物観測あり | **3,355** |
| 200km格子：同支持条件の植物cell | **1,923** |
| 200km格子：source真正隣接pair | **2,722** |

原観測者・元species×GPSの個別地点はCI上でのみ使用し、集計成果物に転載しない。**source研究用写真観測の豊富さと真の植生plotの被度・未検出種は別物**。全球のコース100km Jaccardは個別の鉄道路線に沿った種子散布・同一生息地の局所群集類似性ではない。

Source artifacts： [v369 全球100/200km source cell植物種数・努力・Jaccard集計 #14912667253](https://github.com/zuizui0223/tetsu/actions/runs/45375453375/artifacts/14912667253)、[元S3 source QA #14912582489](https://github.com/zuizui0223/tetsu/actions/runs/45375453375/artifacts/14912582489)。

この真正成功を確認して、[v371 全球鉄道km×全球植物記録の実データ統合分析](https://github.com/zuizui0223/tetsu/actions/workflows/tetsu-v371-world-rail-plant-Jaccard-source-macro.yml)を源原データ完了条件つきで再実行依頼した。**予測改善の実データ結果は、このv371実ジョブの成功と結果を確認する前に肯定／否定とも主張しない**。

## 全球植物相の真の解析条件（未完了と結果を区別）

- [v369 2025年元iNaturalist全世界13.70GB研究用source](https://github.com/zuizui0223/tetsu/actions/runs/45375453375)：元S3 iNat observations.csv.gzとtaxa.csv.gzのETagを固定し、全世界の2025年research-grade野生植物の100km/200kmセル内記録件数、unique observers、種/属/科 richness、隣接セル対のJaccardを直接抽出する。**公式source全件のstreamがSUCCESSとなり実際のchecksum・aggregateが揃うまで全球植物の効果量は未計算**。iNat sourceは種位置をobscuredするgeoprivacyフラグを持たないため、近傍80mの真の鉄道植生には使わない。全年/季節/分類群を完全に観測したfloraでもない。
- [v371 全球railkm×2025植物richness／regional Jaccardモデル基礎CI #45377420094](https://github.com/zuizui0223/tetsu/actions/runs/45377420094)：**9/9合成テスト成功**、sourceデータ実計算は当初実行時はv369未完了のためスキップ。v369は後続のsource実行でSUCCESS。baselineは種数を元観測努力・元distinct observers・WorldClim BIO1・緯度でモデル化。M1は元格子内鉄道kmを追加。隣接セルの源Jaccardは種記録数の幾何平均と差、元気温、緯度、方向を調整してrailkm平均＋差を追加。1000km元地理block等重み、source大陸を1つ丸ごと伏せて予測loss改善を判定。研究の主材料は**鉄道網 vs 周辺地域全体の植物記録**であり、実線路内植物群集ではない。
- 既存別実装 `scripts/model_global_railkm_plant_landscape_macro_v371.py` とworkflow `tetsu-v371-world-railkm-global-botany-inference.yml` も同じ実全球入力を対象に存在する。片方の実装で良い結果が出た方を都合よく採用せず、結果を得た後はsourceキー、モデル説明変数、foldと推定数値を**相互監査**する。v369結果が使えなければ両実装とも肯定/否定の全球生態学的推定値は報告不可。
- WorldClimは**現地気候**であって種レベルのtemperature preference Tではない。植物の暖地選好を直接検定するには独立全球種別温度ニッチ（EIDC等）と元学名の真正クロスウォークが必要。

## 全球機能形質のデータゲート

NHM London 2026 ADEPT/WFO元データ DOI `10.5519/p3dm31kc` は**59,369被子植物レコード・84 trait columns**（source official CKAN DataStore `combined.csv`、life form/habit/dispersal mode/height等）を持つことが[original v371a source QA SUCCESS #45376217612](https://github.com/zuizui0223/tetsu/actions/runs/45376217612)で確認済み。**これらのレコード全てで散布様式と生活形が両方測定されているとは限らない**。元CSV直アクセスがHTTP403、公式DataStore小分け読取でもページング上の制約を修正中。原source59,369全件のjoint trait completenessが確証されるまでは全球機能的収斂や未取得lifeform=otherといったカテゴリー埋め込みを主張しない。

## 科学的な区別

**確立済**：全球鉄道中心線の母集団、物理ネットワーク接続、距離、気候と地図上の延長、地理/気候層化の標本フレーム、同一原source全世界グリッド曝露。

**未識別**：全球の同一調査デザインを持つ真の鉄道敷内植物群集と隣接自然群集（現存種/不在を揃えたもの）、種別生理的温度ニッチ/全球生活形・散布ギルドの十分な対応、駅/主線/支線実運行頻度と実種子輸送。全球macroregional associationで植物記録の予測が改善しても列車種子散布、バラスト環境適応の機構にはならない。

### 再現
- 世界元rail kmの気候帯分布：`scripts/analyze_world_railkm_climate_distribution_v373.py`, `tests/test_world_railkm_climate_distribution_v373.py`, `.github/workflows/tetsu-v373-fullworld-rail-climate-exposure.yml`
- 真正全球植物×鉄道km区域マクロ：`scripts/analyze_original_world_railkm_plant_richness_Jaccard_v371.py`, `tests/test_original_world_railkm_plant_richness_Jaccard_v371.py`, `.github/workflows/tetsu-v371-world-rail-plant-Jaccard-source-macro.yml`
- 全球functional trait source gate：`scripts/audit_original_NHM_global_functional_trait_coverage_v372.py`, `tests/test_original_NHM_global_functional_trait_coverage_v372.py`, `.github/workflows/tetsu-v372-NHM-global-59k-functional-trait-source-coverage.yml`
- 元著者の個別植物観測者ID・稀少植物GPS・位置、現地floraや鉄道網の個別原元線分の生ファイルは公表リポジトリには転載しない。
