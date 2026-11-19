# tetsu v29.7 — 75%系統樹単系統群で残る草丈の空間構造は道路・内陸でも成立するか

2026-11-19 JST、source-native新しい効果量を見る前の分析設計。元植物site×taxon出現・LEDA高さ＋種子重のpaired2D形質値・野外調査単位を変更しない。実測値はGitHub source workflowが成功してから追記する。

## 科学的な判定の軸

v29.6：元の独立した欧州鉄道4原著（チェコ39駅、イタリア28駅、クロアチア＋ボスニア15駅＝1出版、スロベニア89個100m道床地点）では、真正V.PhyloMaker2（LCVP/WP/TPL）の既存系統樹末端種を使い、根から先端の枝長中央値75%で切った単系統群の中で草丈・種子重の2D形質を組み替えても、観測された草丈CWMの地理的分化が4研究全部で帰無より大きい。元原著等重み平均LCVP+0.056695、WP+0.055429、TPL+0.058998。これは「鉄道への適応」そのものではない。

新しい独立な参照はMIREN山岳道路調査。道路際と同じSite×Yearで同じ100m2を調べたInterior植生を持つ。RoadsideとInteriorの両方に同じ草丈系統構造があれば、鉄道特異的選択でなく地域種プールと地理・系統の共変動による一般植物群集パターンでも説明できる。Roadside正でInteriorなしなら線状の人工生息地に関連する追加の群集構造候補。RoadsideとInterior両方正、ただしRoadside−Interior追加超過距離が負なら、実群集の非ランダムな草丈の地域配置と相対的な回廊均質化の両立を示唆。

## 元sourceとfield独立性を事前固定

MIREN公表原CSV DOI https://doi.org/10.5281/zenodo.15749536、元17,899,241 bytes MD5 b640195cc44f252ccb6ddb28eb6fff13、真の植物出現167,804行、2007–2023年。元MIRENのRoadside Plot 1と同位置Interior Plot 3をペアにし、元Region×Road×Siteの最初の有効年だけ使う。同一siteの後年の反復記録を新しい独立サンプルとしない。v29.1以来の真source2パネルを変更しない：

- Primaryは元各PlotでLEDA草丈＋種子重両方が既知の植物5種以上、元生物全リストの25%以上。8元Regionコード、22実Road、272実Site。
- Strictは元各Plotで10種以上、50%以上。5元Regionコード、13実Road、192実Site。

これは同じMIREN field programの入れ子感度比較であり、独立した2実験や鉄道5本目のfield研究ではない。元全5,739植物種中LEDA 2形質既知1,096種（19.1%）に限られ、元SourceのTrait coverage偏りを常に記載する。MIREN元データ自体の2D source z軸で先行v24.9道路際−Interior差−0.10381093795211047を再現できなければ新たな帰無計算を停止。

今回の高さと種子重の座標はMIREN datasetから再fitせず、欧州鉄道原著で凍結済みの同じz(log height, log seed mass)平均（−0.580236575,+0.293611319）とSD（1.380100055,2.165766056）を使用。

## 本物のpublished系統樹と元sourceの支持域

Jin & Qian 2022 V.PhyloMaker2（DOI https://doi.org/10.1016/j.pld.2022.05.005）、LCVP/WP/TPL元著者の公開R-data original phylo3樹を固定GitHub commit 7af3fb5152f691af2e4ec9d5e2e467d1b50505e9 と3 SHAで照合。これらは三命名法版の同じMegaTreeに基づく感度解析であり、独立野外論文3本にはしない。元作者のtreeに存在しない種は植生データから除外するが、系統樹へ新しいtipを推測で接ぎ木しない。treeのtip掲載はその種のDNA直接配列がある証拠とは限らない。

各source Region×tree×基準panelを元sourceの生態効果を見ずに審査：
1. 元LEDA両形質既知相異taxaの60%以上がtreeの元exact original tip。
2. 元Region内Roadside＋Interiorの元site×species出現数の60%以上が元direct tip。
3. 元sourceと同じ1:1 Site×Year pairのRoadsideもInteriorも、元source primaryは5種、strictは10種の元Tree直接tip既知taxaがあること。Roadで少なくとも4地点、Regionで元Road2本以上。
4. 元source選定されたfield siteの80%以上が条件に残ること。

これらが不成立ならNOT_IDENTIFIABLE_ORIGINAL_DIRECT_TIP_FIELD_COHORTとし、残存植物数・site数・road数のQAを報告、草丈0という値は作らない。元植物の種リスト・GPSは出さない。

## 真枝長の条件付き形質ランダム化

同じtreeで原original tipのみの種subset・同じ元site×taxon incidence・同じEurope固定zlog形質座標で4帰無を比較する。

- 元Source Region内の実種出現頻度五分位でheight+seed paired trait vectorを組替える基準。
- 根から先端までの元Tree枝長中央値に対し0.25で切った実単系統群内のみ組替え。
- 同じ0.50 depth単系統群内でのみ組替え。
- 同じ0.75 depth単系統群内でのみ組替え。鉄道v29.6で支持した真の枝長cutと同じ定義。

source region内で道路site間の高さ、種子重、2D CWM距離をRoadごとに平均して各Road等重み。RoadsideとInteriorには同じ一回のtraitラベルの交換を適用、元現地ペアの相関を壊さない。各Region内、2種以上の交換可能true-cladeが5以上、元Plot1/Plot3の両方で交換可能taxaのサイト出現25%以上、元Road2本以上で両habitat各25%以上あるときに限り、2seed×199=398 technical randomizations。条件未達のモードはNOT_IDENTIFIABLE_ROAD_INTERIOR_PHYLOGENETIC_CLADE_NULLとする。技術的順位は野外field p値ではない。

各tree×Source primary/strict panelで3元Region以上の適格siteがあれば等Region記述平均と1Region除外範囲を出す。元regionコードを独立論文や国として扱わない。v29.6の鉄道4原著は元plot survey grainが異なるため、MIRENとの値を一つのランダム効果メタ平均へ直接投入せず、まず方向・源tree tip欠測支持域・群集系統構造の再現性を比較する。

## 因果的限界とコード

ここで検証するのは地理に沿う植物種と種中央値の草丈という群集のパターン。線路バラストを使った稼働中railと同等面積・同調査年の開放砂礫非rail植生の真source matched対照は未取得。列車の種子運搬、風、草刈り、植物個体の形質可塑性・遺伝適応の機構は示せない。

実装：
- scripts/validate_MIREN_original_roads_true_phylogeny_clade_negative_control_v297.py
- tests/test_MIREN_road_interior_true_tree_source_v297.py
- .github/workflows/tetsu-v297-source-original-MIREN-ROADS-three-megatree-clade-placebo.yml

元植物の個別source site×taxon、元LEDA生データ、phylo tree raw種名と位置をGitHubに再公開しない。元source原著MD5/SHAと元Region単位のQA、効果量のみartifactに残す。
