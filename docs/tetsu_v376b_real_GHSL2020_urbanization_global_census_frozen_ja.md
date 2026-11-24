# v376b : 元JRC 2020年全球Built-up 1km rasterを真正取得し100/200km鉄道・植物格子へ集計（実行前固定）

v376a真正元source [CI #45381799454](https://github.com/zuizui0223/tetsu/actions/runs/45381799454)で欧州委員会JRC GHSL R2023A全世界2020 built-up 1km ZIPのHEAD HTTP200、Content-Length **152,527,564 bytes**、ETag `91762cc-622dc4d319ac9`、Range HTTP206、original ZIP magic `PK03\x04`を確認。元ソース DOI `10.2905/9F06F36F-4B11-47EC-ABB0-4F8B7B1D72EA`。[公式仕様](https://human-settlement.emergency.copernicus.eu/ghs_buS2023.php)。

## 真正源GIS変換
1. 元publisher 2020 1km `GHS_BUILT_S_E2020_GLOBE_R2023A_54009_1000_V1_0.zip`、先行HEADの`Content-Length`・ETag（GET If-Match）を照合し、全圧縮ZIP SHA256を新たに凍結。唯一の`*.tif`を抽出して元WGS84でない**World Mollweide ESRI:54009**、ピクセル幅約1000m、`m2 of built-up area within original pixel`という源単位を確認。元NoDataは実source確認で **uint32最大値4,294,967,295（0xFFFFFFFF）** と判明し、0（実際のbuilt-up無し）とは明確に区別する。初期の「NoData負値」という想定は誤りのため、この確定した源GeoTIFFを優先する。
2. 元GERS（v370）・元iNaturalist全植物（v369）は独立のEPSG:6933 100km grid原点`x=-17,400,000, y=-7,400,000`で完全一致。**同元100km CELL=10,000km²**の中に重なる元1km source pixelsの built-up m² を`rasterio.warp.reproject(Resampling.sum)`で空間合計し、`builtup_m2/(1e10)`をbuilt-up fraction（0〜1）として保存。**単なる中心の1km市街化の値は面積10,000km²の人為被覆率を代表しないため採用しない**。
3. 元200km結果は、**同じ元100km 2×2格子のm²加算/4e10**により作る（元GHSLを別にWarpすると微細な数値差が出るのを避け、総量の保存を検証する）。元nodata/海域セルは負のまま残し、測定不能地域にbuilt-up 0を代入しない。粗格子への推定変換は目視だけでなく、統計（GHSL no data数・builtup fractionsの最大/平均・両スケール元built-up m²整合）を監査する。
4. このv376bは**独立した2020年人間活動の説明変数ソースのみ**で植物や鉄道とまだフィットしない。出力は全球粗格子100kmと200km built-up fractions、原文献・schema・源全ZIP SHA256とaggregateのみ（原GeoTIFF生画像や居住地原個票のGitHubへの再配布なし）。
5. source CRS/単位/Extentが予想と異なる場合は fail closed で架空の都市化係数を報告しない。独立ソース真正取得後v376cで元v371源同一植物全大陸保留・λ0.10 ridgeで`M0: geoclimate+observers`、`M0b: M0+GHSL built-up`、`M1b: M0b+GERS railkm`を同じsource対象に対し予測し、予測差を判定。
