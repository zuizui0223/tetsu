# tetsu v28.7 — チェコ鉄道14駅 vs 周辺CzechVeg植生：HABITAT明示・調査年代・面積・標高記載の厳格対照

更新：2026-11-18 JST。ユーザー研究 tetsu の真正原データ再解析。**本ファイルの数値結果は元GitHub Actionsの原資料再取得・元14/70比較の既報数値再現後に追加する。** Source SHA/MD5を通過しない限り肯定・否定を記録しない。独立調査原sourceの個別学名×site/個別GPS/plot IDを再配布しない。

## 研究の問いと既報の出発点

以前v28.5で4独立著者出版の欧州5か国鉄道源植物は、元site×taxon配置を完全に固定してLEDA源Canopy height＋Seed massを2Dペアでtaxonに付け替える帰無と比較すると、**地点間の高さCWM距離が4原著全てで正方向に分化**していた（CZ +0.05725、IT +0.03823、HR+BA +0.09533、SI +0.05653）。これは鉄道内部の**source traitsとtaxon位置の対応**の結果で、非鉄道原著対照との差ではない。

v28.6のCZ1出版社元Preslia14 station-wide植物群集 vs 真のCzechVeg他原著の70原参考plot（近傍10km、元source70別ID＋70別WGS84 100m rounded cells、元14駅×5元plot・21参考原著）では、駅間高さCWM距離差Rail−周辺=−0.28140、種子重=−0.29723、2D=−0.42673で全100源technical species equal-sampling drawsで負、元4source railway stretchesのleave-one-outでも負。つまり内部の分化と周辺群集に対する相対的な類似は両立できる。

しかし鉄道調査は元駅**全域**の複数生息地合算（各原著元station area 2–90ha）、CzechVegは元別原著の小植生調査区×5で、季節・管理・生息地・植物poolと現場面積がそろわない。このまま「鉄道が草丈を均質化した」と因果主張できない。

## v28.7 の事前固定した source raw metadata の候補選択

v28.6で原CzechVeg植物2,510,340出現記録と117,739元header、原Preslia source S1/S2真正xlsx、SR70第三者駅代表点、Oldenburg LEDA height/seed元種中央値をいずれもSHA/MD5で照合する。元16本物のトレイト合格駅・元candidate reference plot opportunityを固定する。

元HABITAT欄の記述文字列のみに開放草地・メドウ・ルーデラル・砂礫・岩石・採石場という植物**分類群に依存しないテキストの指標**がある場合に「開放地proxy」候補とする。森・湿地など明示的に両義性のある場所は曖昧として除外。元SYNTAXONという**種組成から事後的に決まる植生分類**は選定に用いない。これは現地で物理的に基質を確認したことの代用にはならない。

8源パネル：
1. original source broad nontransport context（元v28.6の14/70基準、完全数値再現で元原著比較固定）
2. original reference YEAR>=2020（鉄道元2020–22のeraへ絞る）
3. source \`HABITAT\` original text explicitly OPEN/BARREN/RUDERAL（元species outcome blind）
4. open+YEAR>=2020
5. open+元CzechVeg reference \`SURF_AREA\` 16–400m²
6. open+YEAR>=2020+SURF_AREA 16–400m²
7. stronger explicit sand/gravel/rock/quarry original \`HABITAT\` source
8. open＋reference \`ALTITUDE\`元記録あり（標高がsource元駅と一致したという意味では**ない**）

すべてのpanelに対し元16観測可能Preslia駅は保持し、元control候補は本物の別原著、元駅の10km以内、元reference2015–23/0＜PRECISION≤250m、種数5以上／既知2D trait coverage ≥55%、各駅>=5参考plot・>=2元source reference papers、参考70 plot全体のID・WGS84の0.001°セルを重複禁止し、整数MILPの最適性が確認できた場合だけsource元駅数を確定する。比較可能な元Preslia線区が3本未満、元駅10未満なら**結果0ではなく判定不能**とする。選定する際に元植物のtaxa identity、CWM、seed/heightを使わない。

採択supportがある場合だけ、元v28.6と全く同じ100種数均等化technical drawsを使い、同じ元観察駅の真Rail vs近傍5plot composite 2Dとheight/seed 1D βCWMを再計算。元Preslia採集線区内ペアの4/3元線区等重み比較、元線区leave-one-outを報告し、元繰返し抽出を100独立野外実験と数えない。

## 識別上の停止条件

- **元Preslia側に個々の鉄道駅源実測標高がない**ので、CzechVeg元参考plotのALTITUDE値があっても、駅との標高一致は検証できない。周辺5plot内部の標高ばらつきは記述できるが、駅との実差を推定してはいけない。
- 元Preslia駅全域の調査面積はha規模、CzechVegの単一relevé source SURF_AREAはm²規模。元plant taxaの種数を希薄化しても野外観察面積・生息地の混合を同じにしたとは言えない。
- CzechVeg source ROCK coverは元source missing→0補完が疑われ、元0を実測0%岩盤と呼ばない。
- 元HABITATにrock/sand/gravel/grasslandが書かれていても実鉄道と同じ稼働中のバラスト・剪定・乾燥ストレスを物理確認したわけではない。
- 現地源研究はCZ rail1著者原著、元CzechVeg source複数reference papersは別植物生息地の観測であり、これを21独立RAIL研究の再現と数えない。観察地選抜・種の元LEDA存在数・地理／標高／管理差は依然として交絡。
- 同一次元のfunctional β距離の差の変化を、元supportが変わるパネル間で純粋な「生息地マッチング効果」と解釈しない。source cohort/元site×plot grainの変化を常に併記する。

## 実装と出力

- \`scripts/audit_Czech_original_14station_reference_habitat_era_area_altitude_v287.py\`
- \`tests/test_Czech_reference_original_HABITAT_area_year_v287.py\`
- \`.github/workflows/tetsu-v287-CzechVeg-original-habitat-area-era-controlled-rail-reference.yml\`

すべてのsource元Preslia/CzechVeg/Oldenburg/SR70を厳格ハッシュ確認し、元v28.6の2Dsource baseline −0.42673326131926104が1e−9の精度で再現しなければそこで停止。生物の結果はsource field supportがあるパネルにのみ記載する。
