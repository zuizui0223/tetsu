# v305 — スロベニア鉄道調査地点の周辺人工被覆：取得可能性と計測監査

確認日：2026-11-20 JST。植物群集との関連を見ずに、89地点すべてへ同じ処理を適用した。

## 採用した一次資料

ESA WorldCover 2021 v200（DOI [10.5281/zenodo.7254221](https://doi.org/10.5281/zenodo.7254221)）の公式公開S3 COGを利用した。公式データは3度四方のEPSG:4326 GeoTIFFで、部分読出しが可能。[公式配布仕様](https://esa-worldcover.org/en/data-access)、[開発元の配布案内](https://github.com/ESA-WorldCover/esa-worldcover-datasets)。

分類50は建物だけでなく道路・鉄道を含む人工被覆である。都市の公園はこのクラスに含まれず、森林下の人工物も樹冠に隠れ得る。したがって「都市人口密度」「市街地の行政区分」ではなく、**周辺景観の人工被覆率**として扱う。[公式PUM v2.0、分類表p.15](https://esa-worldcover.s3.eu-central-1.amazonaws.com/v200/2021/docs/WorldCover_PUM_V2.0.pdf)

2021年地図は植物調査2022–23年に先行する1–2年前の景観である。地図上の変化を同時に観測してはいない。2020版と2021版の差をそのまま土地利用変化として扱わない。現地土壌、草刈り、風、種子供給の直接計測ではない。

## 実装と取得確認

- コード：[extract_Slovenia_worldcover_buffers_v305.py](../scripts/extract_Slovenia_worldcover_buffers_v305.py)。既存の `src/tetsu/worldcover_landcover.py` の公式URL・クラス辞書を再利用。
- 原植物位置は、著者公開Excelの各100m区間の両端GPS中点。既存v103の元ファイルSHA検証とv277の89地点／15路線検証を通した。
- 外半径500m、1000mの円、およびそれぞれ内側100mを除いた環状域を取得。中心の鉄道構造そのものが人工被覆率を押し上げる問題に対し、環状域を感度分析として残す。ただし円の外側へ延びる線路や駅施設は除去されない。
- 円への所属はWGS84楕円体の中心―ピクセル中心距離で判定。各行のピクセルの楕円体面積で重み付けする。カテゴリの補間はしない。タイル境界をまたぐ候補も複数タイルから読み、重複しない画素を合算する。
- 欠測値0を人工被覆なしとはしない。分類済み面積比を別に保存し、クラス率の分母は分類済み面積とする。
- local-onlyキャッシュ：`../v305/urbanization.json`。元地点IDで結合できるが、GPS、植物種リスト、原画像は出力しない。元ExcelのSHAと読出し窓のピクセルSHAを記録する。ピクセルSHAは全タイルファイルのハッシュではない。

公式タイルURL：

- [N45E012](https://esa-worldcover.s3.eu-central-1.amazonaws.com/v200/2021/map/ESA_WorldCover_10m_2021_v200_N45E012_Map.tif)
- [N45E015](https://esa-worldcover.s3.eu-central-1.amazonaws.com/v200/2021/map/ESA_WorldCover_10m_2021_v200_N45E015_Map.tif)

ローカルPython3.13のrasterio/pyprojでリモート部分取得を完了した。89/89地点、各4範囲すべてで分類済み面積比1.0。円・環状域の画素中心近似による面積と幾何学的面積との差の絶対最大は0.1952%。500m円の人工被覆率0–0.841619、1000m円0.001056–0.674204。この幅は都市化の測定に利用できるが、幹線・地域線の共通支持域を保証するものではない。

人工の小配列で、面積重み、欠測の分母除外、100m内側除外を確認した。現地での分類精度検証は行っていない。全地点の取得完了と、都市化交絡の解消は別問題である。

## キャッシュの結合仕様

`sites` は89行。各行は `source_id`, `section`, `buffers`, `raster_windows` を持つ。`buffers` のキーは `r500_inner0`, `r500_inner100`, `r1000_inner0`, `r1000_inner100`。各範囲には `built_up_fraction`, `tree_cover_fraction`, `cropland_fraction` 等の11クラス率、`valid_area_fraction`, `sampled_area_m2`, `n_pixels`, `n_valid_pixels` が入る。`complete: true` を要求してから生態解析へ渡す。

主共変量は1000m円の人工被覆率、500m円と100–1000m環状域を感度分析に使える。これらの選択は解析計画で固定し、結果が好ましい半径に切り替えない。森林・農地率は都市化以外の景観差を表すが、89地点・15路線、実交通量ではさらに少ない支持数で、多数の景観変数を同時投入して過適合させない。
