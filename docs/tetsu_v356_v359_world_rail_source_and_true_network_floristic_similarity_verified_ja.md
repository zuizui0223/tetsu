# 全球鉄道網全件抽出と「同一路線」植物群集類似性の検証 v356–v359（2026-11-24）

## 主題
世界の鉄道そのものを生物地理のネットワークとして扱い、鉄道の任意の切片や一国の「駅数」ではなく、**実際の線路接続と、名前のある原著区間内外で植物相の類似性を比較する**。同じ地理距離で鉄道に沿う移動距離が異なるなら、植物群集のtaxonomic Jaccardが鉄道網距離と関連するかを明示的に検証する。

**元の全球鉄道路線のsource**: [Overture 2026-11-18.0](https://docs.overturemaps.org/blog/2026/11/18/release-notes/) official transportation theme schema v2.0.0. 発行元による全rail-subtype（地下鉄・トラムその他の全クラス混在、元source線路モデル上のcenterline総延長）**2,892,198km**。世界の独立した「名前付き営業線」の総距離、列車の運行頻度、権利地の土地面積ではない。

## v356 元全球輸送sourceインベントリ真正確認 — SUCCESS
[GitHub Actions #45342025751](https://github.com/zuizui0223/tetsu/actions/runs/45342025751)。公式 `s3://overturemaps-us-west-2/release/2026-11-18.0/theme=transportation/type=segment/` 匿名 S3 全ページ一覧：
- 元全sourceオブジェクト128個/128 Parquet（元source全量 72,770,647,616 compressed bytes）
- full original inventory SHA256 `6237da737bb783f23f352638cd5aef287bbb3d7aa8f0fe0bed8c1cc4558f6b34`、listing全ページ取得確認
- プラハの初回線路元source preview 947 rail-subtype segments：`standard_gauge`441、`unknown`157、`tram`317、`subway`31、`narrow_gauge`1。接続IDあり947/947、観測ユニークconnector 3323、複数元segment共有connector 776。source `names.primary` 付与63/947（約6.65%）、source `routes[]`付与0/947。Source OpenStreetMap `w<ID>@<version>` way型参照も確認（1,953 ref entries、元segmentは複数source refs可）。
- 生物学的出力を作る前に確認：**1段階で物理GERS connector network、2段階で名前付き路線を独立に復元しないと「本当の路線」と混同する**。単に地図座標が交差するだけではconnected componentに入れず、`connectors[].connector_id`の共有を条件にする。巨大な物理的連結成分は多くの営業線を包含しうるため、`component` ≠ `route=railway`。

## v357 128ファイルの全球rail-only Parquet抽出 — 実行状況
[World source scan #45342573085](https://github.com/zuizui0223/tetsu/actions/runs/45342573085) を16分割（各source8 Parquet）で開始。該当subtype全てを抽出し `id,class,connectors,names,sources,rail_flags,geometry` を固定、元sourceは保存・再配布せず新しいrail-only ParquetをGitHub Actions artifactとして出力する。**全16完了＋全128元ファイルの同一性＋release公式rail length検証が通らない限り「全球全件抽出完了」と表記しない**。

- Streaming pilot [#45342353780](https://github.com/zuizui0223/tetsu/actions/runs/45342353780) SUCCESS、元S3の実GeoParquet/connector field schemaを実確認。Pilot sourceの最初65,536 transport rowsにはrail1つで、**source geographic partition selection**を意図したものではない。
- 最初の完成source shard #1は **8/8 source files**・元transport rows **22,144,979**件・`rail` **94,946件**・元中心線測地長 **125,163.125km**、rail-only parquet47,651,690 bytes SHA256 `729da5544db3ac0dac238d0566360123381cb51a7526abf8689633b72b225977`、source full checksumとID維持、正式な「全球全件」ではない。

## v358 真の線路connector network × スロベニア原著植物相 — SUCCESS
[Authenticated source CI #45343086076](https://github.com/zuizui0223/tetsu/actions/runs/45343086076) は真正のŠabić et al.2026、原著サプリ401,279B、SHA256 `26be7ca9513dd54136d4ec499d8413c1cb288045f87fe3cf8814731f2edac99f`、原著89地点/著者15線路区間と固定先行v28.0のoriginal59 geo-matched anchor/8 source sections、原著BS+BT Jaccard差+0.019670941578を実データで再現。そのうえでOverture物理ソース：
- Overture source regio transport 861,914件、`rail` 8,886件、通常線路クラスと放棄・廃止・工事中除外後、**6,612真のsource rail segment、12,346接続点、14,181graph edges、286 regional physical connected components**（地方ウィンドウによる切断が含まれうる）。
- 原著89植物地点のうち**86/89**を実source線路位置へ割当；残り3は異なる物理成分を同時に近くに持つため割当保留。原著地点→地図線路の最近接距離中央値**1.865m**、最大報告距離18.97m。地図上の交差点をconnector共有無しに偽接続しない。
- 元59 source geodesic-matched triples中、真の物理経路がboth sidesとも完結した**48 anchor /8 publisher sections**。Geo same vs other平均 **21.331 vs 21.413km**、真のconnector-network geodistance同一section vs other **27.691 vs 55.812km**。類似度original Jaccard同一区間−別区間**+0.020759**、5/8 sectionsが正。
- このpatternは単純な同じ水平距離であっても、別の著者区間へ行くのに線路網を迂回することを実物グラフで示す。しかし「線路上の短い経路だけが原因」という因果解析ではない。両者の元群集は2022/23調査、線路マップは2026/11で時代差もある。

## v359 **地理距離と真の線路網距離の両方を対応**した源floraの強い反証テスト — SUCCESS
[Authentic original-source Actions #45343455821](https://github.com/zuizui0223/tetsu/actions/runs/45343455821)。原著89地点/15 source section、元site source geocoordinates/真のOverture rail network mappingをv358どおり再現。
- マッチングsource-original事前条件：アンカー→同じ原著sectionの別地点j、別sectionの別地点k。両方の水平距離 ≤40km、水平距離差**≤2.5km**、実railway connector graph上の距離差**≤5km**。植物種名やJaccard値を元対象の抽出に使わず、各アンカー1 source triple。完全original physical paths無しは「不明」扱い。
- **27 original anchors /7 source sections** がjoint matching可。section等重みのoriginal geodesic距離は同一source20.188km vs 別19.624km；実線路網距離同一 **25.059km** vs 別**27.329km**。地理・network距離の条件を両方揃えられた。
- BS+BT同一section Jaccard **0.264197** vs 別**0.237164**、差 **+0.027033**、7 source sections中5正/2負。Leave-one-author-section-out差範囲 **[+0.008934,+0.051121]**。方向はこの条件付き小標本ではpositive。単にnetwork距離が短いだけではsource line-species identity clusteringを消せない可能性はあるが、**実際のsource section-specific seed dispersalであることを示す訳ではない**。
- BS shoulder-only差+0.035867（7 section中3正/4負、LOO[−0.004831,+0.060375]）。BT toe-only差+0.009538（4正/3負、LOO[−0.007486,+0.031165]）。部位を分けると符号は安定しない。
- さらに同じ調査年＋単線/複線レイアウトまで源referenceを揃えると**6 anchors/3 sections**で事前最低12/4に達せず、「NOT_IDENTIFIED」。効果が負やゼロになったわけではなく、推論できるsource支持がない。
- 以前の本物v28.1年/線路構成ごとの種出現頻度固定nullでは同一sectionの距離マッチ後の余剰効果は+0.007693、technical upper-rank 0.226だった。v359の27地点とは異なる元集合なのでその数値を流用した`p`を作ってはいけない。今回の強い帰無（原著植物の種頻度、観測努力、周辺土壌、非鉄道群集、交通OD）も未評価。

## 何が分かったか
**直接観察**：固定路線内の物理的距離と植物の種組成を、実際のGERS connector graphに接合可能。スロベニア原著で、地理距離と実線路網距離の二重マッチ後も一部では同一路線の分類学的類似性が上乗せするという小さい記述的傾向が確認できた。一方で厳しいsource年/線路様式条件では標本不足。環境・管理・地域プール・線路沿い交通流と市民科学観測機構は未分離。既存自然植生matched controlsがないので線路固有の効果を証明していない。

**まだ言えない**：世界中の路線で同じ結果、列車が同じ植物種を運搬している、線路長があれば地理的隔離を消せる、鉄道の微気候選別、C4や暖地植物の適応。全世界の同一名前付きrailway relationを照合するにはOSM route=railway/tracks relationをOpenRailwayMap等から別取得する必要。実物線路のconnected componentを実際の独立営業路線と同義に扱わない。

## 再現
- Source global inventory script: `scripts/audit_global_Overture_all_rail_parquet_inventory_v356.py`.
- ALL 128 global original Parquet streaming scan: `scripts/extract_entire_world_Overture_rail_parquet_v357.py`、源source 16 shard aggregator `scripts/validate_entire_world_Overture_rail_shards_v357.py`、`.github/workflows/tetsu-v357b-entire-world-rail-source.yml`.
- v358 network+botanical original `scripts/analyze_Slovenia_Overture_connector_network_botany_v358.py`, sources and QA `.github/workflows/tetsu-v358-original-Slovenia-real-connector-network-botany.yml`.
- v359 original joint matched `scripts/analyze_Slovenia_geod_and_rail_distance_matched_v359.py`, QA `.github/workflows/tetsu-v359-SI-original-botany-joint-rail-network-distance.yml`.
- [Frozen global original railway source design](tetsu_v356_global_rail_all_segments_line_network_source_frozen_ja.md), [v358 fixed graph botanical design](tetsu_v358_Slovenia_original_rail_graph_similarity_frozen_ja.md), [v359 fixed dual match design](tetsu_v359_source_Slovenia_real_network_distance_matched_flora_frozen_ja.md).

Original source publishers/OSM ODbL/copyright respected, Overture full geometric original data not committed directly to GitHub repo. Draft PR #58 stacked on PR #57, main unchanged.

## v359b：本当の線路距離と地理距離を揃えた27地点の**同じ比較集合**に種頻度保存帰無を適用 — SOURCE VERIFIED（2026-11-24）

[**Original-source and Overture-GERS network run #45345069976 SUCCESS**](https://github.com/zuizui0223/tetsu/actions/runs/45345069976)。元データŠabić 2026原著401,279Bを元SHAと照合し、v359真正27 anchors/7 author source sections のBS+BT source Jaccard差**+0.027033312**を1e-10以内で完全再現。生態的なネットワーク・地理距離二重マッチのアンカーと植物記録は**一切入れ替えていない**。

同じ89原著調査区の種×地点行列を元に、各地点の植物種数と各種の出現頻度を完全に固定したcheckerboardシャッフル。2 independent seeded chains×99 retained correlated MCMC draws。全89地点での種頻度のみの帰無H0-Aと、2022/2023調査年×SINGLE/DOUBLE線路構成内での種頻度を保持したより厳しい帰無H0-Bを、**新しい元27地点セットに**適用した。両familyは元source結果を調べる前に固定（ただしv359観測値を見たあとの事後頑健性分析）。

| v359b source original same 27 matched anchors / 7 publisher sections | H0-A 全site種頻度保存 | **H0-B 年×線路構成別site種頻度保存** |
|---|---:|---:|
| 元観測の同一区間−別区間植物Jaccard | +0.027033 | +0.027033 |
| 技術帰無の平均 | +0.004555 | **+0.018549** |
| 観測−技術帰無 | +0.022479 | **+0.008484** |
| 技術帰無2.5–97.5% | [−0.03238,+0.03874] | **[−0.01895,+0.05987]** |
| 上側技術順位（**フィールドp値ではない**） | 0.0955 | **0.3116** |
| 原著7 sectionで観測＞各section帰無平均 | 5/7 | 5/7 |

**最も強い結果解釈：** source network＆地理距離を両方揃えても、小さい同一区間類似性+0.027は観測される。しかし、元の調査年・単線／複線別の種頻度を保持する技術帰無の平均+0.01855を差し引くと残差**+0.00848**となり、技術帰無でも観測と同程度の効果が普通に現れる（rank0.312）。したがって、**線路網の独自の輸送や分散が植物相を同化する独立した証拠は得られていない。** これは「鉄道植物の同一路線類似性が存在しない」の意味ではなく、観測された記述的差の説明機構を厳格に限定するもの。原著89地点/15sectionは1か国であり世界普遍性/列車ODは不明。帰無シャッフルは地理的な原点種プール、同じ管理、列車ODを完全制御していないので、フィールド因果仮説のp値を生成しない。

[元の技術検定凍結](tetsu_v359b_Slovenia_network_geo_matched_fixed_species_frequency_null_frozen_ja.md)、[原著集計JSON](../data/derived/v359b_original_27site_7section_real_network_matched_species_frequency_null.json)、[解析コード](../scripts/analyze_Slovenia_networkmatched_species_frequency_null_v359b.py)、[真正元original source workflow](../.github/workflows/tetsu-v359b-SI-rail-network-species-null.yml)。これにより先のv359の「可能性」は**機構としては未支持**へ格下げした。
