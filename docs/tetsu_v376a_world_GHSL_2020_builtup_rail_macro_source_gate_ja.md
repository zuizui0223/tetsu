# v376a：全球鉄道km×植物観測種数の独立都市化交絡データGHSL、真正ソース取得ゲート

2026-11-24。ユーザーの方針は**全球鉄道の全線路・植物の全球観測のマクロ解析を優先**すること。国別の鉄道plot研究はサブ解析に留める。

## 先行実データ結果
元iNaturalist公開観測2025年7,827,936、元GERS2026全鉄道候補2,237,043.481km。修正版v371全球100km景観記録種数の大陸留保6-fold平均MSE差は+0.00446374549（4/6大陸正）、200kmは+0.00661424383（3/6正）。v375元continent×気温帯×緯度帯で鉄道kmの場所だけ入替え39回は100/200km richness、100km Jaccard実地図の予測改善を39回とも下回った（run 45380848733）。しかし**人口・市街化・道路・写真アクセス**の同時代理効果は未統制。

## 新たな本物の独立ソースの定義
- [European Commission JRC Global Human Settlement Layer GHS-BUILT-S R2023A](https://human-settlement.emergency.copernicus.eu/ghs_buS2023.php)、2020年実績の全球 built-up surface（元1kmセルに含まれるbuilt-up面積m²）、標準世界Mollweide 54009（本研究EPSG6933の等面積100/200kmへ面積加重合計で再投影）。元 official DOI `10.2905/9F06F36F-4B11-47EC-ABB0-4F8B7B1D72EA`。**1km格子値は「当該セル内建築地の㎡面積」であり、人口や鉄道交通量を表さない**。
- 元大容量全球ZIPを推測ファイル名としてダウンロードする前に、実 publisher directory/HEAD/HTTP Rangeによる原ソース取得可能性を監査する。対象固定resource:
`https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL/GHS_BUILT_S_GLOBE_R2023A/GHS_BUILT_S_E2020_GLOBE_R2023A_54009_1000/V1-0/GHS_BUILT_S_E2020_GLOBE_R2023A_54009_1000_V1_0.zip`
- 最初のv376aは**HTTP source metadata/schema/zip signatureの取得ゲートであり生物学的予測改善の結果ではない**。成功したsource ContentLength、ETag、LastModified、ContentType、Range supportをsource監査に残す。source HTTP403/404/5xxのとき「source unavailable」と明示し捏造フィットをしない。

## v376bの分析規約（sourceが本当に入手可能になった時のみ）
1. 真正JRC ZIPファイル全量SHA256を確定、内側唯一のoriginal `*.tif` GeoTIFFのCRS、nodata、pixel size、最大値、東京都/プラハ周辺の数値をチェック。global target原100km/200kmはv369/v370源EPSG6933 `x0=-17400000,y0=-7400000`とする（v351別原点の土地格子と混ぜない）。
2. 元1kmセルの built-up m²を100km/200kmセルへ合計・面積規格化する。単に**100km格子中心の1km built-up値**を全球市街化率とみなさない。GeoTIFF nodataがsea/unknownであれば都市0と同一視しない。得られる値は2020建物面積の代理変数で、樹冠や舗装路面の総面積を含むとは限らない。
3. v371修正版と同じ元source全球植物種数・近隣Jaccard 100km/200kmの格子と努力しきい値を保持し、**M0a（既存観測数・観測者・BIO1・緯度）→M0b（＋GHSL建築地面積率）→M1b（＋実際のGERS鉄道km）**を元大陸留保・train weighted ridge λ0.10で比較。M0a→M0bとM0b→M1bを分けて報告し、都市化による予測改善と鉄道kmの条件付き残差を区別。
4. 引き続き全球**景観植物観測**をアウトカムと呼び、線路上の真植生plot、生理的温度選好、気候制約の解除や列車種子散布を識別したと主張しない。市街化proxyだけでは道路km・人口・農業・観測生息地アクセスは未統制。
