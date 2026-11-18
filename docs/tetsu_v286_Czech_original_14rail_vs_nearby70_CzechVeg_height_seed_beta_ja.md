# tetsu v28.6 — チェコ14鉄道駅と周辺70独立CzechVeg植生区：草丈と種子重を分けたβ機能多様性

2026-11-18 JST 新規事前固定設計。真正Source GitHub Actionsの原資料再取得・厳格元v24.4結果再現が通った後で、生態結果を追加する。元site×taxon・GPS・生学名や参考元調査区の座標を新規repoへ再配布しない。

## 問題設定

既存の新v28.5で4独立著者原著（CZ, IT, HR+BA, SI）では、元駅・道床サイトの**種配置を固定しLEDA草丈と種子重の二形質のペアをspeciesへ付け替える帰無**より、元著者地理区内のCWM草丈差・2D距離が4/4正だった。だが、同じ草丈構造は同地域の**非鉄道植生**により強く現れる可能性がある。各研究内の種の形質地理対応は、鉄道特有の選択の原因を示すわけではない。

既存v24.4 [真のsource GitHub Actions SUCCESS #44006279169](https://github.com/zuizui0223/tetsu/actions/runs/44006279169)では、CZ元Preslia 39駅のうち厳格な位置・参照plot重複禁止ができる**14本物の駅**と、各駅10km以内の独立CzechVeg植生plot**5つ/駅・計70本**を比較した。全参考plotはsource IDがユニークで、緯度経度を丸めた約100mセルが異なり、CzechVeg本来の元位置誤差条件は \`0 < PRECISION ≤250m\`。元reference publication数は21。元生態結果は同一駅14×6?違う、14真駅の91依存駅ペアでsource両形質CWM Euclideanの鉄道−局所reference距離**−0.426733（鉄道0.270484、参考0.697218）**。これは**「鉄道植物相が全て同じ草丈に収斂した」という因果効果ではなく、元鉄道の駅全域植生対CzechVegの小plot5つ和集合という不揃いな調査単位**。

今回のv28.6はこの元14駅・70元独立reference plot・100mセル制約・各真駅の標準化種数のmatchingを**完全に固定**したまま、既存2D距離を**草丈単独・種子重単独**に分解し、元source4著者rail stretches内のみの駅ペアにも同じ指標を算出する。

## 一次資料

- CZ railway: Kutlvašr et al. 2024 Preslia DOI https://doi.org/10.23855/preslia.2024.247。元author S1 SHA256 df6bc0b949e7f6a5cb57755dc8f9c86187223f66ddfbef959cf0e03ef7c46b06、元source corrected S2 SHA256 00c7cfdf5047854ba99d47dc8ec8c55b1b7df72731ce647c22b286ddef24a5a2。元source 2020–22 39 stations, 5 survey stretches, 761 current recording taxa。元2026 independent SR70 station centroid SHA256 fc43a1969af64f17285c7c6e3a30b00c3994d0ba1cfda6db09d8e4baaa166019（**元現地plot GPSではなく駅代表点**）。
- CzechVeg original true plots: Chytrý et al 2026, https://doi.org/10.5281/zenodo.19231152 version 1.1 (2026-03-26)。元117,739 real original header rows (MD5 0bc920462d1f9eac473e7b2a231e6782), 2,510,340 original species-plot rows (MD5 b296e5e19f03fb3bad0a8347c317f8a3)。元公開CSV約460MB、全データを取得しverify、部分記録やnative vegetation statusを生成しない。2015–23年、元sourceのHABITAT/LOCALITY等に交通語がない現場記録を調査対象とする。ただし非交通のテキスト表示は実際の鉄道・道路の物理的非存在の完全保証ではない。
- Oldenburg Univ. LEDA元canopy height（m）と元seed mass（mg）の元観測種代表分布、源本物のSHA256 cbf5e761390d78a843512dfc1139acfc1154d46d6c1995a827354840a9f78a92と34aa12c333fc1ba6de2e560a495ac5d693aa0722a09ecf68418ad94ac1c19cef。CZ14 source rail+local reference両方の元source種のunique matchのlog traits平均・SDを**一度だけ統一推定し同じ2D座標**にする。v28.5の四研究前3研究622種による共通欧州軸とはスケーリング母集団が異なるため、観測効果量を無条件に合算しない。

## v28.6 実行仕様・検定停止条件

1. 元source原植物調査駅39→真SR70代表座標35→2015–23 reference10km内plot数5/元出版研究数2以上19→source species+LEDA結合後16本物の元駅→MILP0/1最適性証明の14駅。元参考plotは計70異なったIDと小数3位丸めWGS84セル、真正4元rail stretchesを保持。生物の形質によるreference plotの手動再選択はしない。
2. 各元rail station植物taxaと5非交通文脈original relevéの和集合taxaの両traits既知種集合から、元各駅で厳密 \`n=min(n_rail,n_local_ref)\` 種を100technical drawsランダム抽出し、元駅・元referenceのsource種数を一致させる。元Source selection seed \`20261117+244+291*draw\` とv24.4と完全一致。全100drawの2D結果が **鉄道0.270484／周辺0.697218／差−0.426733** と絶対誤差0.001以内で一致しなければ、結果を保留してCI fail closed。
3. 同じ元14駅の2D CWMで\`|height_Ri−height_Rj|\`/ \`|height_Ni−height_Nj|\`、\`|seed_Ri−seed_Rj|\`/\`|seed_Ni−seed_Nj|\`を別々に計算。元14駅91依存ペア全体での差と、元4rail stretches**内だけ**の駅ペアの差、および元stretchを等重みした差を別々に計算する。元stretchを1つ除いても方向が一致するかを点検する。
4. α site trace（2D形質幅）は別応答として保持し、駅間のβ距離と混ぜない。元source連続種中央値は現地植物個体の測定や進化した草丈を意味しない。
5. 推定対象の重要な違い：v28.5は**完全に鉄道内のsite植物群集を固定**した種→形質帰無。今回v28.6は元駅と別原著・10km内局所CzechVeg 5plot poolの地理ペア比較で、**現地調査面積・群集タイプ・季節・観察努力、都市化・標高・土壌・地元種プール**が完全対応しているわけではない。元観測体数調整の同種数希薄化でも環境選択の因果効果は特定しない。
6. 著者の元CzechVeg位置誤差最大250mだから、100mセルの重複ゼロは実際の現場個体集団の完全独立を保証しない。民間列車貨物OD、実現場種子移送、草刈り日・風圧・遺伝・可塑性の元データはない。

## 持続的出力

\`scripts/audit_Czech_original_rail_vs_local_CzechVeg_height_seed_beta_v286.py\`、\`tests/test_Czech_14station_local_nonrail_height_beta_v286.py\`、\`.github/workflows/tetsu-v286-Czech-14rail-70nonrail-height-beta.yml\`。全元素材真正MD5/SHAチェックとMILP/既知2D結果厳格復元、100drawのサイト・元文献調査単位に応じた効果集計。元生plot×speciesと元reference座標・原元CSVは出力しない。
