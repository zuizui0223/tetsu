# v356: 全球の線路そのものを母集団にした物理鉄道網と植物類似性（固定設計 2026-11-24）

## 科学的な問い
「鉄道の同じ路線で離れた区間の植物相は、単なる直線距離・周辺気候・観測努力を揃えた別路線より似ているか？」および「鉄道**網上距離**（network shortest-path）が直線距離より植物相の変動をよく説明するか？」。既存7原著の機会的抽出やv351の32ランダム**陸地格子**とは異なり、鉄道物理ネットワーク全体を研究対象の抽出枠にする。

## 使う唯一の真正全球交通sourceと完全性の確認
- Overture Maps 2026-11-18.0 schema v2.0.0（11月6日OSM cutoff）。公式release notes https://docs.overturemaps.org/blog/2026/11/18/release-notes/ は`subtype=rail`の報告全長**2,892,198 km**（全クラス・物理モデルの重複を含むsource線分長で、営業路線長と同義ではない）。Original all transport segment source `s3://overturemaps-us-west-2/release/2026-11-18.0/theme=transportation/type=segment/*`、connector別source `.../type=connector/*`。ODbL attributionと元公開地図のsource条件を遵守。
- *最優先*: AWS S3 original global Parquet file indexを再帰的に列挙し、`file count`、全sourceバイト、原本キーlistのSHA256、ページ全件取得を検証。匿名ListBucketが403などで拒否なら、DuckDB read_parquet globで全source enumerationが可能か独立確認する。ファイルインベントリの完了なしに「全球全線路抽出完了」と呼ばない。
- 全sourceファイルを列挙して、その全ファイルについて`subtype='rail'`だけを抽出。対象を植物種・国・暖地効果で絞らない。source固定元全件count、class別、flag別の全国数・元line lengthを計算する。少なくとも公式released 2,892,198 kmと計算値が合理的に一致するか照合。単純なLineString WGS84長（経緯度 degree）で距離を計算しない。官方sourceの全長を引用するときは正本published value、pipelineの幾何原source実測は別に表記。
- 元OSMの全`railway=rail`を混同しない。Overture `subtype='rail'` はtram/subway/light rail等も**全て**含む。主比較は `class in {standard_gauge,narrow_gauge,unknown}` として意図する通常線路、全railはsource totalの監査パネル。unknown種別が不明。公式v349bの「橋/トンネルを落とす」は連結を分断するため**ここでは橋梁とトンネルを残す**（列車が走るなら接続）。 `is_abandoned,is_disused,is_under_construction` のactive-candidate除外は別の感度であり元全rail censusには含める。使用中の保証・運行頻度は地図から導けない。
- Overtureの`connectors[].connector_id`を参照して無向グラフを作る。地図上交差・線路近接だけでedgeを繋がない。connector共有で複数segmentが連結。異なるsegmentの生態学的比較の独立単位をconnector辺の細切れ一個にしない。source `id`が公開GERSの一意な線分ID。source名 `names.primary`は任意で、**connector連結成分は「路線名」ではない**。交換駅で複数路線が一つの巨大な連結成分に吸収される。
- 真の「名前付きインフラ路線」にはOSM`type=route,route=railway`および`route=tracks` relationsをsource OSMapPBFから結合する。これは旅客列車サービス`route=train`と別物。OSM relation被覆が国によって不均等なため、被覆率を明示し、欠損を自動的に「単独路線」とみなさない。Overture `sources[].record_id` は元OSM way IDを含む可能性があり、源の実例で必ず確認。元のrelationへの多対多結合は路線名と区間の関係、直線的な走路名がない場合は`corridor_id`を欠測として保持する。

## 植物類似度の凍結エスティマンド（元独立植生リレベが必要）
1. 元source botanical plotを公開された位置精度・線路距離に投影し、確定したroute relation/connector-corridorに割当。列車駅全域のフロラを鉄道区間の定面積群集と混同せず、鉄道肩・バラスト・路肩など原著の部位を保持。観測由来ならphotos/個体観測は不在を含む完整種リストとせず、検出過程/effort補正。
2. 同一沿線 vs 別沿線の植物同定Jaccard・Sørensen、functional guildまたはGower distanceの距離。ペア群集は地理近接・気候類似・季節・原著/調査努力をマッチし、`same_named_line` と `network shortest path` の増分予測（距離、気候、土地利用、空間ランダム場を先に入れる）。接続した線路を通る観測距離だけで鉄道運搬を証明しない。
3. 帰無: 距離ビン/同一年・同生息地・同源論文/同気候帯に制限した群集の路線ラベル置換（route cluster）、道路のネットワークを使った偽corridorなどのnegative controls、地理距離だけで学習したモデルに対するleave-one-line-out/leave-one-region-outの予測利得。非独立pair 36/3916などを反復数として利用しない。
4. 既存真正原著での独立なpilot: Slovenia原著89 trackside BS+BT植生地点を15**原著source section**に分類済で、同一section内外のJaccardが算出済（v141）。これは真のroute line IDではなく著者定義section、実際のconnector graph shortest-pathではない。Czech原著39 stations/5 railway stretchesも別の駅全域のフロラである。**すでにsource section 類似性を計算したことを全球network因果検証済と言ってはいけない。**

## Global extraction practicality / abort criteria
全球original transportationは巨大（roads ~97% of segments）；全世界全source ZIP/GeoParquetをGitHubへ投稿しない。Source inventory→partitioned streaming `subtype=rail`→aggregate n/length, plus source **rail-only** Parquet shards as GitHub Action artifacts if licensing/size permits. Original full raw official network data are linked to authoritative publication instead of duplicated in repository. All original file source objects have provenance and per-shard source completeness; failure in ANY shard prevents describing the subset as complete global census. S3 anonymous read, geography partition, source latest release must be pinned, otherwise no silent fallback.
The source v351–355 land-grid source random sample and the six urban-window observation pilot remain separate materials; neither contains enough plant plot-by-source-rail-line data for universal similarity estimate yet.

## Phased real outputs
v356 SOURCE INVENTORY: all S3 Parquet object count, total bytes, names SHA and per-file source schema sample with real connector/names/source_osm id. Source original S3 list and DuckDB cloud-glob fallback.
v357 FULL RAIL: every original parquet object scanned, `subtype=rail` and graph connector IDs materialized into shards; source complete status and all rail class/flags totals; original all-source rail n+coverage official check.
v358 TRUE LINE: OSM route=railway/tracks member relations on source OSM ways, coverage and multiparent associations; no invented named line in uncovered nations.
v359 BOTANY: Czech 39 stations/5 stretches and Slovenia 89 sites/15 sections as separately defined pilot plus real map-matched Overture corridor for source coordinates only if original mapping fit is verified. No global ecological claims until botanical data coverage/range is real.

Research Draft PR #58 stacked on PR #57 (main unchanged); fixed design before source results.
