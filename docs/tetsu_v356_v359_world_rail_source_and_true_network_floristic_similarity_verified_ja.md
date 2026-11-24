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

## v357 世界原本抽出の追認済み途中経過とsource route識別（2026-11-24、完了版ではない）
[Official all-source workflow #45342573085](https://github.com/zuizui0223/tetsu/actions/runs/45342573085) のoriginal completed source shards 0,1,2,3,4の**40/128** official Parquet fileを、各アクションログとoriginal source file countのQAで照合。元original交通segment 110,628,831件のうち真正 `subtype=rail` segment **516,926**件、source元線路ジオデシック長合計 **760,285.070 km** をrail-only Parquet **245,500,601 bytes**に圧縮して保存済。source full worldwide railway 128/128は未実行完了であり、760,285kmを単純に (128/40)倍して世界総延長を推定することは禁止。別のsource 11 shardについてはAction実行中／未着手。SOURCE-COMPLETE判定は128/128を行った後だけ。

- original OSM `route=railway/tracks`と実物GERS接続成分を混同しない独立のsource coverage設計 [v360](tetsu_v360_named_real_OSM_railway_relation_membership_frozen_ja.md)。最初の2026-10-16臨時PBF原本はGeofabrik HTTP404で解析前停止。公式source-only [#45344719906](https://github.com/zuizui0223/tetsu/actions/runs/45344719906) により入手可能な月次固定版 **2026-10-01** PBF 313,417,178 B、公式MD5 `703113982dc57e356f08dfccb1ddbf1a`をSOURCE数値を見る前に固定してroute relations解析を実行中。route名を元GERS地理近接から勝手に補完しない。
- 真の全球物理鉄道connector成分構築コード [v361](../scripts/build_global_Overture_physical_connector_components_v361.py) と [元source完全性を要する手動実行workflow](../.github/workflows/tetsu-v361-full-world-real-rail-connector-graph.yml) を追加。**元128Parquet全件完走前には実行・結果宣言しない。** コンポーネントの数は名前付き営業路線数ではない。

## v357 original global source full 128/128 SUCCESS (2026-11-24)
[Source-authenticated all-128-file Actions #45342573085](https://github.com/zuizui0223/tetsu/actions/runs/45342573085) completed all **128/128** original transportation Parquet files as 16 source-verified rail-only GeoParquet shards. Exactly **354,631,439** original transportation segments, and **2,054,855** subtype=rail segments. Geodesic source global rail centerline sum **2,892,198.162 km**, within **0.162 km** of publisher's 2,892,198 km. Full original global source inventory SHA256 6237da737bb783f23f352638cd5aef287bbb3d7aa8f0fe0bed8c1cc4558f6b34. Source classes count: standard_gauge 1,483,865; narrow_gauge 30,874; unknown 401,716; subway 57,042; tram 58,962; light_rail 17,863; monorail 3,143; funicular 1,390 (sum exactly 2,054,855). All original rail segments have connector info, 7,630,212 source connector incidences; 310,851 original GERS segments have names.primary. Source flag incidences abandoned231,134, disused142,612, under construction10,731, bridge180,576, tunnel80,497; overlaps possible. These are ORIGINAL SEGMENT features and centerline lengths, **not operated named-railway route counts or active track length**. [Verified original aggregate](../data/derived/v357_original_global_full128_rail_source_aggregate_verified.json). True global connector-component topology is a separate v361 physical graph step and not a route-name result.

## v360 original true named OSM railway route membership source verified
[Original Geofabrik monthly Oct 1 2026 PBF×Overture rail source×89 authentic Slovenian botanical sites SUCCESS #45345176407](https://github.com/zuizui0223/tetsu/actions/runs/45345176407): Original PBF 313,417,178 bytes, source MD5 703113982dc57e356f08dfccb1ddbf1a, source SHA256 cce34d7bdfd310238dd302f50b9ad61c3ba9e66bc376856bf2cf56e002dbd9f7. The genuine OSM scan read 115,044 original relations including **73 named infrastructure railway route=railway/tracks** relations, 42 different passenger route=train relations excluded. Original route-to-way member incidence 9,587, unique original OSM ways covered 7,981, and 1,414 OSM ways have MULTIPLE named railway route relations. Source physical-rail subtype mapped region 8,886 rail segments; normal class eligible 6,613 of which 722 (10.92%) matched at least one named OSM infrastructure railway route, including 272 with multiple named route membership. Of 89 author-original field botanical locations, 86 source physical tracks uniquely matched: 43 with exactly one named route, 42 with multiple, 1 with no named relation; 3 ambiguous physical railway components. Thus a physically connected GERS component does NOT define a single original named railway route, and forced one-to-one route IDs would misrepresent 42/86 genuine field locations. Route coverage of original botanist sites is highly nonrandom relative to global/regional rail segments. Source Oct1 vs Overture Nov18 temporal mismatch also remains. [Source-verified OSM infrastructure route aggregate](../data/derived/v360_original_Slovenia_OSM_route_railway_true_line_membership_verified.json).

## v362 actual named infrastructure routes and original species-Jaccard test pending
[Design locked](tetsu_v362_real_named_OSM_route_floristic_similarity_frozen_ja.md): real named source OSM route relations, one-to-many true memberships, authentic original 89-site botanical incidence and Overture physical graph distances. Primary route uniquely identified in all three contrasted field locations, sensitivity any overlapping named route vs disjoint; geographic and railway-network distances jointly matched, genuine survey year/layout and exact source species-frequency null. No named-route botanical effect asserted before authentic source run succeeds. Previous v359b physical-network source section-effect was nonexceptional after stricter source frequency null.

## v361 全球physical rail connector network 全件構築 SUCCESS（2026-11-24）
[**v361 TRUE original full worldwide physical GERS connector graph Actions #45354737204 SUCCESS**](https://github.com/zuizui0223/tetsu/actions/runs/45354737204)。元v357全128/128 Parquetとrail-only16 blobsを原SHAで照合し、元GERS physical connector UUIDを同じなら連結するunion-find物理グラフを2枚構築。地図線路の単純な視覚的交差を真の線路接続と取り違えない。Original worldwide source 2,054,855 all rail-segment identity and 2,892,198.162km lengthを再現。実線路名・実際の乗車可否は別source。

| 原本 physical network 指標 | All Overture rail subtype | 通常鉄道の地図上active candidate（standard+narrow+unknown class, abandoned/disused/under-construction除外。橋/トンネルは維持） |
|---|---:|---:|
| actual source GERS segments | 2,054,855 | 1,536,226 |
| source rail centerline km | 2,892,198.162 | 2,237,043.481 |
| **true source connector nodes** | **4,970,641** | **3,547,570** |
| **physical connected components** | **142,792** | **44,710** |
| physical component with only 1 original segment | 114,552 | 35,704 |
| LARGEST component original segments | 1,208,759 | 1,009,974 |
| LARGEST component centerline km | **1,548,769.636** | **1,369,687.050** |
| LARGEST component share of original subject total rail km | **53.55%** | **61.23%** |
| Top10 source components centerline km | 2,421,827.420 | 2,030,263.105 |
| Top10 source components share of subject km | **83.74%** | **90.76%** |

**科学的な意味**：数十万のsource small isolated componentsと非常に少数の巨大なsource connected networksが同時に存在する。**巨大なsource physical graph component＝一つの営業路線ではない**：路線の分岐/乗り換え/連結管理・複軌間・軌道ゲージ・通行禁止などと、複数営業路線が接続する可能性を原GIS connectorのみでは区別できない。Graphのsource conn UUIDは図形/地図上の物理的連続の表現であって、貨物OD、旅客運行、列車走行可能性、種子散布実験を表さない。地図source境界欠損や誤りはcomponent countを変える可能性がある。この結果は**ネットワークの地理学的/インフラデータ構造の新事実**であり、**種組成による全球均質化が観測されたという意味ではない**。

検証可能な集計source [v361 derived global component metrics JSON](../data/derived/v361_source_verified_world_all_rail_physical_connector_components.json)、実装 [v361 connector build](../scripts/build_global_Overture_physical_connector_components_v361.py)、原GitHub artifact ID14902715264。

次の真の生物地理的問いは、（a）植物調査地点が各source named infrastructure routeに割当可能か、（b）同一路線・物理network距離・地理距離・地域種プール・生息地の観察頻度を同時に含むモデルが距離減衰や独立留保した路線を予測できるか、（c）route=railway関係被覆の地理的不均一を補正できるか。Sloveniaでは元89地点の一部がnamed relation複数所属し、名前付き路線を無理やり単一ラベルに統一してはいけない。追加独立地域の同法による植生定面積プロットが必要。
