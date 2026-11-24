# v376c：全球鉄道網の植物相予測情報は都市化によって消えるか？（2026-11-24）

## 研究質問と科学的境界

**世界全線路の実際の地理位置は、周辺地域の公開植物観測種数を、現地マクロ気候・緯度・観測者数・記録数・建築面積割合を統制したあとも予測するか。**

原研究では「全球の鉄道が植物を増やす／種子を運ぶ」ことを問うていない。全世界の植物記録は、野外定面積plotではなく **2025年公開写真付き研究用植物観測7,827,936件**から100km／200km地域の**記録種数と隣接地域Jaccard**を集約したもの。線路敷内の植物生育・種子運搬・適応・鉄道熱の実測はない。

## 独立真正ソース、検証済み
- **植物**：iNaturalist 2026-10-27 original worldwide photo-linked S3 gzip 13.70GB、288,721,938原観測行の完全source処理から2025年7,827,936件の元研究用植物レコードを抽出。[v369原データSUCCESS #45375453375](https://github.com/zuizui0223/tetsu/actions/runs/45375453375)、元OBS SHA256 `1c48a9646f3bf6d56ffdfdc91a1a176c173c98014830ad48fa162b2803c3b279`。
- **鉄道**：Overture2026-11-18.0 original 128/128世界鉄道GIS源。通常鉄道候補1,536,226 GERSセグメント、真のgeodesic 2,237,043.481km、100/200km全地理格子に分割。[v370世界全源SUCCESS #45373985228](https://github.com/zuizui0223/tetsu/actions/runs/45373985228)。
- **都市化**：European Commission JRC GHSL R2023A Built-up Surface 2020 [DOI:10.2905/9F06F36F-4B11-47EC-ABB0-4F8B7B1D72EA](https://doi.org/10.2905/9F06F36F-4B11-47EC-ABB0-4F8B7B1D72EA)。原publisher ZIP 152,527,564bytes、全源SHA256 `d61ce2e88cf7441b062c0fbd9624c98618c682f4f6bee5522663669d4cd11fe0`、元Mollweide ESRI:54009 1km各pixelのm² built-up（source真NoData=4,294,967,295）をEPSG6933等面積100kmへ`rasterio Resampling.sum`でaggregateし、200kmを元100kmの2×2source面積合算で作成。**21,895有効100km・6,143有効200km**、100/200km全源built-up m²保存確認。[v376b元GHSL全源SUCCESS #45382671629](https://github.com/zuizui0223/tetsu/actions/runs/45382671629)、[source artifact #14914622910](https://github.com/zuizui0223/tetsu/actions/runs/45382671629/artifacts/14914622910)、digest `sha256:4d5ce1c05d0e426f4b0e6142848a9cc9dcb1c91a870b5bc7341b2a2d4df4c443`。
- **現地気候・地理**：WorldClim2.1 1970–2000年BIO1 source原ZIP SHA256 `00513224583665ec0f2f955a4ec252730c4deb2004cce9e793492a3f26df4dcf` と Natural Earth原行政map SHA256 `718ac51d68ed2bf5de4f563fff2e2e045db706c1d3222805f3e317140a64213b`。

**sourceモデル実行・チェック成功**：[v376c ORIGINAL GLOBAL RESULT SUCCESS #45384121307](https://github.com/zuizui0223/tetsu/actions/runs/45384121307)、合成元schema/加重モデル9/9テスト。source派生aggregate [artifact #14915810785](https://github.com/zuizui0223/tetsu/actions/runs/45384121307/artifacts/14915810785)、digest `sha256:a64fd1757d01d416f083cd2e1f080fae80aa8b2853e6909e7557aa53f4ebdac6`。出力に元個別植物・観測者・GPS・個別線路ジオメトリは再配布しない。

## 厳密に同じoriginal global field support

元GHSLデータ結合ではsource v371の元地理/観測努力を通過した **100km 3,455植物cell／隣接4,736ペア、200km 1,442植物cell／隣接1,953ペア**の全てが都市化データありとして残った。すなわち**都市化が不明なため標本が減ることによる擬似的効果はない**。source城市建築率＝0は100km 3 cell、200km 1 cell、元NODATAをゼロへ置換しない。

比較モデルの同一原-source・6大陸全折・1,000km原block等重み（総和1）・Ridge λ=.1（各訓練大陸のみ標準化）：

- **M0**：`log1p(植物記録種数)`については元記録数、元distinct観測者数、WorldClim BIO1と二乗、絶対緯度と二乗。隣接Jaccardは両セルの観測努力・記録種数・気温・緯度・隣接方向も使用。
- **MR**：M0＋実GERS鉄道km（元v371再現用）。
- **MU**：M0＋元2020 GHSL built-up割合（`log1p(1000×fraction)`）。隣接Jaccardは2セル平均・差。
- **MUR**：M0＋GHSL同一説明変数＋GERS同一説明変数。

**解析上の重大バグと修正**：v371の`continent_holdout(base, extended)`は`extended`に**追加列だけでなく全説明変数**を渡すAPIである。途中のv376c試験的コードでGHSLが二重に加わると誤って推定し、追加列だけを送る誤りが生じて元QAで停止した。修正版は**MU/MURそれぞれの全説明変数リスト**を渡し、各列一度のみ・MU列集合がMURへ包含される検査と9/9回帰テストを追加。修正版の元source CI #45384121307のみ採用し、失敗runは一切結論に使わない。

## 全世界の実データ結果

ΔMSEは**大陸を丸ごとholdoutした平均`M0−MR`、または`MU−MUR`**。正なら鉄道kmを追加したモデルのheldout予測誤差が小さい。

| 元地域全景観floraの指標 | 鉄道単独の追加改善 M0−MR | **都市化調整後の鉄道追加改善 MU−MUR** | 都市化調整後、改善した大陸 |
|---|---:|---:|---:|
| **100km 元植物記録種数log1p** | +0.0044637455 | **+0.0025802022** | **4/6** |
| **200km 元植物記録種数log1p** | +0.0066142438 | **+0.0042858862** | **3/6** |
| 100km 隣接景観Jaccard | +0.0001476674 | +0.0001268978 | 5/6 |
| 200km 隣接景観Jaccard | −0.0000612002 | **−0.0000242899** | 3/6 |

全6大陸等重み平均のM0/MR/MU/MUR元MSE:
- 100km種数：**0.108133 / 0.103669 / 0.107607 / 0.105027**。MUからGERS追加は相対**2.398%**のMSE減少。
- 200km種数：**0.121283 / 0.114668 / 0.120343 / 0.116057**。MUからGERS追加は相対**3.561%**のMSE減少。
- 100km Jaccard：0.0048162 / 0.0046685 / 0.0048467 / 0.0047198。
- 200km Jaccard：0.0049042 / 0.0049654 / 0.0048000 / 0.0048243。

**2014〜現在の鉄道植物研究との解釈の差を維持**：このglobalな元flora richnessは各鉄道線路敷に存在する生物種数ではなく、元2025年source写真連動の地域植物レコード種数。「鉄道網が植生の多様性を増やした」「実列車種子散布が起きた」「鉄道のバラストの暖地化」は支持されない。

### 大陸ごとの不均質性（都市化を加えた後の種数追加改善ΔMSE）

| 大陸 | 100km MU−MUR | 200km MU−MUR |
|---|---:|---:|
| Africa | **+0.007983** | **+0.005807** |
| Asia | **−0.003931** | **−0.004581** |
| Europe | **+0.011152** | **+0.012886** |
| North America | **+0.008181** | **+0.013265** |
| South America | +0.000571 | −0.000340 |
| Oceania | **−0.008475** | **−0.001322** |

**最も重要な肯定的結果**は、全球のGERS鉄道延長が、元iNat植物種数の大陸外予測に**GHSLだけでは代用できない追加的な人為地理情報**を持つこと。これは100km/200km平均で同方向だが、**Asia/Oceaniaでは改善しない**。「全球すべての大陸で再現した」という主張は不可。

## 交絡と次の全球優先検定
1. **全球道路網、人口・土地利用・標本採集アクセス**はGHSL built-up割合と異なる。本結果は都市化proxyのみで調整しており、鉄道が残ることは**鉄道特有の因果機構を意味しない**。次の独立原source global road kmあるいはGHSL population densityを同じ元source格子に結合し、列車輸送の経路情報との因果推定と分ける。
2. 2020GHSL・2025植物観測・2026鉄道の**年代不整合**を明示する。元true railroad geometryは全候補通常線路の地図中心線で、現役運行・貨物・列車風の計測値ではない。
3. 2025iNaturalist open metadataに真正obscured/geoprivacyフラグは無く、元individual plant GPSの80m以内という近接検査をここから導出していない。source粗格子生物群集の種同一性・機能形質も出力しておらず、全球機能的収斂/H2はまだ未検定。
4. **v375**大陸×source気候帯×絶対緯度10°固定の39fake rail mapsでは、100/200km種数と100kmJaccardの真GERS地図の予測上積みに匹敵するものが0/39。これに今回の源GHSL調整を加えると「気候位置だけの地図擬似効果」と「都市開発proxy」の両方を一部除けるが、*道路網の地理配置*に対する反事実ではない。
5. 巨大観測者数のある格子のsource観測機会バイアスを検証する固定上位5%/10%除外の全球source解析 v377 は別の事後感度であり、まだ結果を確認するまではその効果を引用しない。

## 再現コード
`scripts/aggregate_original_GHSL2020_global_builtup_v376b.py`、`scripts/analyze_global_original_railkm_biodiversity_GHSL_control_v376c.py`、対応テスト・真正元Source Actions。解析用GPSや個別植物種名、元line UUID、元植物観測者は本研究GitHubのsource結果に含めず、公開元提供元から取得する。研究Draft [PR #60](https://github.com/zuizui0223/tetsu/pull/60)、main未変更。
