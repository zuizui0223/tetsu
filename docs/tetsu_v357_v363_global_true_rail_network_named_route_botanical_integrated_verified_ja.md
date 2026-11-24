# 鉄道網は植物相を似せるか？ 全球128/128線路抽出・巨大物理接続網・真の路線別帰無検証（2026年11月24日）

## 要約と最終判定

**全球鉄道ネットワークの実線路原本が本当に全件抽出され、物理接続・路線名の階層が判別できた。** しかし、名前付き鉄道路線を共有した植物相の類似性が強い独立の種子運搬信号であるとはまだ言えない。唯一真正の植物×OSM路線対応があるスロベニアでは、差は正でも元調査年・線路構成ごとの種出現頻度を保つ帰無モデルで例外的ではない。全球のrail-specific floraの類似性そのものは未推定。

### 全球の真の鉄道元source

[Overture original all-128 source Parquet #45342573085 SUCCESS](https://github.com/zuizui0223/tetsu/actions/runs/45342573085)
- 全世界 transport segments source 354,631,439、うちSubtype railのGERS線分 **2,054,855**、original geodesic centerline **2,892,198.162 km**（公式2,892,198kmとの整合は0.162km）。元source128/128検証、16 rail-only source artifacts、元source checksumゲート成功。
- GERS class counts: standard_gauge1,483,865、narrow_gauge30,874、unknown401,716、subway57,042、tram58,962、light_rail17,863、monorail3,143、funicular1,390。**地図上の物理線分数**であって営業路線数や運行ネットワーク長ではない。
- 公式元sourceのrail_flagsにはis_abandonedやis_disusedを含む。運行線路を地図から確証できない。元source raw transport72.77GBをGitHub repoへコピーせず、source rail-onlyをGitHub artifactに置いた。

[Original global physical connectors graph #45354737204 SUCCESS](https://github.com/zuizui0223/tetsu/actions/runs/45354737204)
| 真のOverture物理接続グラフ | 全railクラス | 通常鉄道の地図上active candidate（標準・狭軌・unknown、abandoned/disused/construction除外） |
|---|---:|---:|
| 使用segment | 2,054,855 | 1,536,226 |
| 使用線路中心線km | 2,892,198.162 | 2,237,043.481 |
| original connector IDs（node） | 4,970,641 | 3,547,570 |
| original connected components | **142,792** | **44,710** |
| 1 segmentから成るconnected components | 114,552 | 35,704 |
| 最大1 connected component中心線km | **1,548,769.636** | **1,369,687.050** |
| 最大1 componentの対象線路長比率 | 53.55% | **61.23%** |
| 最大10 componentの対象線路長比率 | 83.74% | **90.76%** |

**科学的に重要**なのは、sourceの物理的な連結性が極端に大きなcomponentに集中していること。世界の植物類似性を「同じconnected componentなら1、違えば0」という粗いbinary変数で予測すると、世界の線路の多くを同じカテゴリーに置き、ネットワーク距離や経路の役割を見失う。代わりに、実線路距離、結節点・分岐、路線relation/運行、周辺生息地を説明変数へ分解すべき。しかしGERS physical connectorは**移動可能な旅客／貨物の運行路線ではなく、同じ列車が全大componentを走行可能という主張ではない**。名前付きOSM路線名もGERSのconnected componentとは異なる。

[Original geodesic connector spatial verification #45356279797 SUCCESS](https://github.com/zuizui0223/tetsu/actions/runs/45356279797)：元世界segment2,054,855と元connector incidences7,630,212を再走査、UUID SHA固定サンプリングで19,375 original connectors抽出。source 2段以上の線路でUUIDを共有する**8,807 connectors**において、元source WGS84 LineStringを**楕円体geodesic distance fraction**で補間して同source UUIDの座標を照合すると、**5m超の不一致は0**。Source位置中央値と最大は出力の小数第3位で0m。元コードの平面経緯度補間では734/8807>100mの**偽の不一致**が出たが、[Overture公式線形参照仕様](https://docs.overturemaps.org/guides/transportation/linear-referencing/)の測地線補間に直すと消えた。これは解析方法ミスの校正、元source地図の破損ではない。サンプル検証であって原全connector4,970,641の絶対完全保証ではない。元GERS UUIDを公開しない。ソース検査資料 [v363 design/erratum](tetsu_v363_world_GERS_connector_spatial_consistency_frozen_ja.md)。

### 「本当の路線」≠ GERS connected component

[OSM genuine route relation × Overture source #45345176407 SUCCESS](https://github.com/zuizui0223/tetsu/actions/runs/45345176407)：Geofabrik 2026-10-01（original MD5とSHA256認証）のスロベニアPBFから、OSM type=route & route=railway/tracks & name/ref非空の元関係が**73件**。元route=train旅客サービス関係42件は除外。元OSM wayが複数本の元実路線relationに所属するのは**1,414 ways**。v358真の89植物野外地点のうち86が物理線路に確実に対応、そのうち名前付きOSM route relationが**一つ43・複数42・無し1**、他3は物理接続が曖昧。対象source通常線路全6,613GERS segment中本当のOSM路線source relationが覆うのは**722 segmentsのみ約10.9%**。植物調査地点の高被覆率から全路線の元route名前対応が充実していると誤ってはならない。

### 真のOSM named-route membershipによる元植物相の検証

[**Authentic named railway route original botanical outcome #45355198378 SUCCESS**](https://github.com/zuizui0223/tetsu/actions/runs/45355198378)。スロベニアの著者本物89植物地点・原著元標本（2022–23、ballast shoulder+toe）、OSM本当の73 named route infrastructure relations（2026-10-01）、物理鉄道ネットワーク（Overture 2026-11-18）で対応。植物種組成を知らずに、元GPS直線距離差≤2.5km、実線路網距離差≤5km、≤40km内の1 anchor-2 referencesを対応した。author元section IDと本当のnamed OSM relation IDを混同しない。

| 真のsource路線relationの扱い | 対応できた元anchor地点 | 元named route signature群 | Jaccard |
|---|---:|---:|---|
| primary: 3者とも元OSM路線IDを**ちょうど1個**持つ | **4** | **2** | **標本不足で未識別** |
| より厳しいsurvey year×track layout（single route） | **0** | 0 | 未識別 |
| sensitivity: 元relation setsを1つ以上共有/別は互いにdisjoint、複数路線所属を許可 | **26** | **9 route-set signatures** | same0.277402 / different0.245569、差**+0.031833** |
| 同じyear×track layout（any-overlap） | **7** | 5 | 最低支持12/4に達せず未識別 |

any-overlapの元26地点・9 source route-set signature単位（独立営業路線9本とは限らない）で、地理距離same22.96 vs different22.47km、actual network shortest path same29.65 vs different31.28km。素の植物相Jaccard差+0.031833は元named route set signature9群中7正方向。しかし元source89地点のsite richness/species frequencyを保存したtechnical checkerboard 2chains×99 draws：
- 全original sites species prevalence null mean+0.006613、観測−期待+0.025220、technical upper rank0.0503（**field pではない**）。
- **survey year×SINGLE/DOUBLE**による厳密なspecies prevalence保存null mean**+0.017226**、観測−期待**+0.014607**、technical upper rank**0.1809**、null technical 2.5–97.5%範囲[−0.013739,+0.050331]。**真のnamed route relation共有と植物相類似性が独立して増加するという強い証拠はない。** Source元固定路線の構成や管理・地域種プールの違いで説明できる可能性が高い。元source路線ID署名が相互に重複していることも独立性の制約。

### 先行v359 / v359bとの整合
元**著者定義15 section**で、真のrail network距離と地理距離を両方揃えた27地点・7sectionのJaccard差+0.027033。元survey year×track layout別の固定original species frequency nullを差し引くと+0.008484（上側technical rank0.3116）、独立の交通運搬特有の信号ではない（[source #45345069976](https://github.com/zuizui0223/tetsu/actions/runs/45345069976)）。本v362の実路線membershipで記述的差は類似して正、帰無モデルで弱いという傾向が再現したが、元source sampleと路線membershipが違うためこれを2国/独立2studyの再現とは数えない。厳密条件では標本不足。

## 科学的な次の一手

新たな中心問：「鉄道網が局所の植物群集へ及ぼす影響は、同じ沿線にあるかより、**実際の移動経路の短さ、途中の分岐・管理の連続性、近接した周辺種プール**により予測できるか」。
- 原source rail 128/128全球network/真のtopologyは取得完了。ただしこれだけでは元flora種組成がない。iNaturalistの点観測は区画完全種リスト・不在を証明しない。
- 真正な全球路線別植物相は、独立した各国line IDとraw source section plot-by-species、同一調査年・定面積・近傍非鉄道対照、source原地気候を追加してleave-one-line/leave-one-region-heldoutの同じ効果を測る必要。植物のglobal thermal niche（EIDC200k原データ全件）は未取得のままなので「全球暖地化」推定も未完成。
- Connector全世界の物理連結が成立することは「列車の種子散布を証明」しない。貨物ODや列車頻度、ハネの曝露、草刈りストレス、原著の温度偏りも未観測。

## 検証sourceの位置
- [v357 original full global rail source census](../data/derived/v357_original_global_full128_rail_source_aggregate_verified.json)
- [v361 original true worldwide connected components](../data/derived/v361_source_verified_world_all_rail_physical_connector_components.json)
- [v363b true WGS84 geodetic connector spatial QA](../data/derived/v363b_original_global_GERS_geodesic_connector_spatial_QA_verified.json)
- [v360 real OSM named infrastructure route membership](../data/derived/v360_original_Slovenia_OSM_route_railway_true_line_membership_verified.json)
- [v362 actual original true named OSM route botanical Jaccard and strict source null](../data/derived/v362_original_Slovenia_true_named_railway_route_floristic_similarity_verified.json)
- [v359b previous original strict source species-frequency null](../data/derived/v359b_original_27site_7section_real_network_matched_species_frequency_null.json)

すべて [Draft PR #58](https://github.com/zuizui0223/tetsu/pull/58) 内、main未変更。元sourceロウ研究者の植物位置、元OSM ways/route IDsと原生植物種名・線路GISをrepoに再配布していない。元データの出典とODbLを尊重。
