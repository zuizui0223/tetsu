# v371 全球：実鉄道網2,237,043kmと2025年全球植物記録によるマクロ予測（2026-11-24）

## 位置づけ
**チェコ／スロベニア／カシミールなどの便宜的な野外事例の統合ではなく、真正の全球source 2本を同じ100km・200km等面積格子に結合した大規模探索的マクロ生態解析。** これは地域全体で公開・写真連動・research-gradeとして**記録された植物種数と隣接格子間の種構成**を対象とし、鉄道バラスト上の実植物種数、列車による種子輸送、個体群・現場適応を測った実験ではない。

### 2本の原source
1. **2025年全球植物記録**：iNaturalist original open-data S3 2026-10-27 monthly snapshot、圧縮観測源13,697,505,636バイト、元S3 ETag `4539360677c8190d9aa12dbf6982ee72-27`、全件SHA256 `1c48a9646f3bf6d56ffdfdc91a1a176c173c98014830ad48fa162b2803c3b279`。元写真連動の世界観測288,721,938行を全件処理、2025年の植物に絞ってsource位置精度<=1000m・quality=research・真の植物speciesへのcanonical化等の基準を満たす**7,827,936レコード**。[完全S3実計算SUCCESS #45375453375](https://github.com/zuizui0223/tetsu/actions/runs/45375453375)、source全球粗格子成果物[14912667253](https://github.com/zuizui0223/tetsu/actions/runs/45375453375/artifacts/14912667253)。源元taxonomy1,675,070分類群、active Plantae species284,070で、後者が2025年7.8m観測で全て出現した意味ではない。
2. **2026年全球鉄道GIS**：Overture Maps 2026-11-18.0 source全128/128元Parquet、候補usual rail1,536,226 GERS線分、WGS84測地線**2,237,043.481km**。[全球全km格子確定SUCCESS #45373985228](https://github.com/zuizui0223/tetsu/actions/runs/45373985228)。鉄道全延長の**実営業状況、輸送列車の本数・草刈り、種子運搬は推定できない**。鉄道の格子内source曝露は物理線分km。

両sourceは同じEPSG:6933 100km・200kmの粗格子に集計した。WorldClim2.1元BIO1気温（原ZIP SHA256`00513224583665ec0f2f955a4ec252730c4deb2004cce9e793492a3f26df4dcf`）、Natural Earth110m大陸境界（元SHA256`718ac51d68ed2bf5de4f563fff2e2e045db706c1d3222805f3e317140a64213b`）をsource照合。iNatの公開種位置にはgeoprivacy/obscuredフラグが無いため**100/200kmより細粒の線路近傍分布は推定しない**。真の標準化した植物調査区のabsence情報やcoverもない。

## 解析単位と主検定
元植物源100kmには実観測のある9,029格子、記録≥50・distinct observer≥5・植物種≥20の条件を満たす4,219格子、隣接Jaccard 5,944組。元land・WorldClim・continent境界を通して最終のモデルは**100km 3,455有効格子と4,736隣接組**、200km **1,442有効格子と1,953隣接組**。欠測・offshoreや低努力cellを「植物0種」や鉄道無影響と読み替えない。

- **M0種数**：`log(1+観測植物種数)`を`log(1+観測記録数)`、`log(1+distinct観測者数)`、WorldClim現地BIO1と2次項、絶対緯度と2次項で予測。
- **M1種数**：M0＋実際のsource鉄道格子内`log(1+railkm)`。
- **M0隣接格子のJaccard**：観測件数・観測者数・元植物種数の隣接2格子平均と差、元BIO1平均・差、緯度、隣接方向で予測。
- **M1隣接Jaccard**：M0＋2格子のrailkm平均、railkm対数値の差。
- 地理1000km source blockは各blockで合計の重みを等しくし、**全blockの合計重みを1へ正規化したWeighted Ridge λ=0.1**で学習。予測ごとの特徴標準化は**学習大陸のみ**で行う。大陸を丸ごと1つずつ伏せて、6外部foldのblockweighted MSE/MAEを評価。**種数/Jaccardの記録やsource誤差を独立した野外studyとは扱わない**。

## 重要な監査修正と二実装の比較
初版[CI #45379245943](https://github.com/zuizui0223/tetsu/actions/runs/45379245943)では、ridgeペナルティλ=0.1に対して元1000km blockの学習重みを正規化せず使い、実質的には格子block数に反比例する**はるかに弱いridgeペナルティ**だった。初版の「6大陸中5大陸で種数予測改善」は**修正後の結果ではないので採用しない**。

[source blocksを1へ正規化して実データを再実行・既存の独立実装で照合した CI #45380248289](https://github.com/zuizui0223/tetsu/actions/runs/45380248289) **SUCCESS、9/9合成QAテスト**。別の既存コード`scripts/model_global_railkm_plant_landscape_macro_v371.py`は、真source植物cellのblockフィールド`original_equal_area_1000km_block`を扱わず`block`だけ読むバグを修正し、源データ実計算が通ることを確認した。正規化を一致させても２実装は学習blockなど一部の詳細が異なり、**同一モデルの完全な数値複製ではなく独立実装による方向・平均の頑健性検証**として扱う。

## 修正版の真正SOURCE結果（各6大陸完全ホールドアウト等重み）

以下のΔMSEは**大陸を丸ごと伏せて評価したM0−M1**。正なら鉄道kmの情報を加えたモデルの予測が改善。括弧内は6foldのうちΔMSE>0のfold数。

| 全球モデル | 新実装・正規化済みridge ΔMSE | 同独立源上の既存v371実装 ΔMSE | 主要判断 |
|---|---:|---:|---|
| **100km 植物記録種数** | **+0.004463745 (4/6大陸)** | **+0.004313222 (4/6)** | 平均の追加予測改善は２実装で再現 |
| **200km 植物記録種数** | **+0.006614244 (3/6大陸)** | **+0.006349818 (3/6)** | 平均は改善するが大陸別に方向が揺れる |
| **100km 地域全植物Jaccard** | +0.000147667 (6/6) | +0.000087581 (4/6) | 小さい予測上積み、モデル実装・空間scaleの感度あり |
| **200km 地域全植物Jaccard** | **−0.000061200 (2/6)** | **+0.000001871 (3/6)** | ほぼ無し／マイナス、安定した均質化の証拠なし |

参考：既存独立実装の基準モデルの大陸等重みMSE=100km種数**0.107848**／200km種数**0.120466**、鉄道モデルMSE=0.103534／0.114116で、予測誤差の相対減少**4.00%／5.27%**。ただし大陸単位では100kmは4/6、200kmは3/6のみ改善。**「全球鉄道が植物種を増やす」とは言えない**。実地の野外absence/observabilityの差も未統制。

100km新実装の種数ΔMSE大陸別は、アフリカ+0.00352、アジア**−0.00659**、欧州+0.02114、北米+0.01496、南米+0.00498、オセアニア**−0.01122**。200kmではアフリカ−0.00137、アジア−0.00350、欧州+0.02202、北米+0.01953、南米+0.00445、オセアニア−0.00144。**欧州・北米の改善が大きく、全世界に同じ向きのモデル性能が一般化しているわけではない**。

### 生態学的解釈の位置
**何が肯定的か**：物理鉄道kmは、公開植物観測努力・気候・緯度の基準モデルに無い**広域的な人為地理情報を追加**する可能性がある。とくに欧州と北米で見え、2実装でrichnessの正の平均予測上積みは維持された。

**何が未支持か**：鉄道による線路敷の植物種数の増加、鉄道による全球群集の分類学的均質化、実際の鉄道路線ネットワークに沿った種子散布、暖地植物の気候制約解除。200km隣接群集Jaccardは改善しないため、単純な全球群集収斂の命題は今回の広域データから支持できない。

**最大の交絡**：世界の鉄道は歴史的な人間活動・人口密度・鉄道以外の道路延長・土地利用・観測者アクセスと一緒に発達している。公衆植物観測の数/観測者を調整しても、それらの社会経済・土地利用共変量はまだモデルに含まれていない。**鉄道kmは人間居住のproxyかもしれない**。またiNatの観測記録は現地植生の豊富さ・完全な種相や全土の長期的多様性と異なる。

## 追加反証の設計
本結果を見てから、同じ大陸×BIO1気候帯×絶対緯度10°内で元鉄道kmだけをラベル入替する、**鉄道地図の偽配置negative control**を別ジョブで実施する設計（v375）。真の鉄道地図の予測上積みがその条件付き無作為配置を上回らなければ、単なる人為インフラの共存地理proxyに留まる。真地図が上回っても、都市・道路・管理・鉄道敷植物との因果を示さない。**このv375が原source全件SUCCESSするまでは結果を引用しない**。遠隔の元地域をそのまま訓練に混ぜるランダムセルCVや植物ゼロの補完をしない。

## 研究の次の方向（全球範囲を維持する）
1. iNatデータを別の元全球植物標本・分布DBの種選好/形質と結合できる率を測り、形質欠測や種の異名を監査する。iNatの現在のsource成果物は**cell別species IDsを保存していない**ため、今回のrichness/Jaccard集計だけでは実植物形質を直接復元できない。生活形や種子散布ギルドの全球解析には、直接source taxon ×位置を安全に取り直す別のデータ経路が要る。
2. 都市・道路・土地利用・人口などの真正全球説明変数を同じ100/200km元source格子で照合し、偽鉄道地図と競争するproxy説明を識別。既存のWorldClimだけでは都合良い結果を鉄道固有の予測力へ帰属できない。
3. 鉄道敷内の全球種組成・同一GERS線路内のJaccard・機能的収斂・実列車頻度・直接種子運搬は**地域全体の撮影植物数や鉄道延長だけでは推定できない**。全球マクロ予測を第一階層、真正元鉄道敷内の定面積floraの世界各地での反復を第二階層として区別。

## Source reproducibility links
- [Original full global 13.7GB photo-linked 2025 iNat records SUCCESS](https://github.com/zuizui0223/tetsu/actions/runs/45375453375)
- [Original global railway source 2.237 million km EPSG6933 100/200km full grid SUCCESS](https://github.com/zuizui0223/tetsu/actions/runs/45373985228)
- [Corrected ridge plus second pre-existing implementation on same source SUCCESS](https://github.com/zuizui0223/tetsu/actions/runs/45380248289)
- Source model code: `scripts/analyze_original_world_railkm_plant_richness_Jaccard_v371.py`, `scripts/model_global_railkm_plant_landscape_macro_v371.py`, tests and workflows in matching paths.
- No original source individual plants/observers/rare original GPS or original rail UUIDs mirrored to public GitHub. 
