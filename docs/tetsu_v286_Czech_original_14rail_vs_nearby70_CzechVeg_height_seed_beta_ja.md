# tetsu v28.6 — チェコ14鉄道駅と周辺70独立CzechVeg植生区：草丈と種子重を分けたβ機能多様性

2026-11-18 JST **真正Source GitHub Actions SUCCESS #44298978338**。著者元植物・独立CzechVeg・LEDA形質をすべて原MD5/SHA照合し、数式QA5件合格、元v24.4の既知2D数値を1e−10精度で再現。参照14駅・70源植生plotのMILPも最適性証明（status0・gap0）。元site×taxon・GPS・生学名や参考元調査区の座標を新規repoへ再配布しない。

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

## v28.6 真正source結果・比較対象の確定

[Original GitHub Actions #44298978338 — Success](https://github.com/zuizui0223/tetsu/actions/runs/44298978338)。原referenceの117,739 headers / 2,510,340本物の植物記録を元出典 2026 Zenodo v1.1 から再取得し、元sample本数と公開MD5を通過。14本物駅×5独立source reference relevésは元70異なったsource IDと70異なったWGS84小数3位（約100m）座標セル、**計21本物の元CzechVeg出版研究**から成る。元鉄道路線支持駅は I 4, III 2, IV 4, V 4。元Rail/参考併合1,224 source taxa名からLEDA2Dが両方真正一致するのは782元ソースtaxa名（rail側434、CzechVeg参考側644）、**train/noneでzを別々に決定しない**。元CzechVeg参考の学名947種中644元LEDA一致という欠損があり、特定の地理や生活形への偏りは未解消。

### 1．同じ実14駅の全91依存ペアにおけるbeta距離

100回の元reference同一species-rarefaction seedと元のsite matched 14駅を保持。結果は元v24.4の二次元β CWM差を **−0.42673326131926104** で厳密再現し、同じdrawでheight vs seedを分解。

| 指標（source zlog軸での絶対CWM距離） | Rail 14実駅 | 同一駅10km周辺source5plot結合・14組 | **Rail−周辺** | 100種数均等化technical draw の差2.5〜97.5%幅 |
|---|---:|---:|---:|---:|
| **Canopy height のみ** | **0.191842** | **0.473240** | **−0.281398** | −0.351113〜−0.204919 |
| **Seed mass のみ** | 0.164166 | 0.461395 | **−0.297230** | −0.348574〜−0.241076 |
| 2D Euclidean（元v24.4値） | 0.270484 | 0.697218 | **−0.426733** | −0.501741〜−0.343705 |

全3指標について100/100の技術的種数均等化drawで負。これは元野外観察研究の1つを百回独立に再調査した信頼区間ではなく、元のtaxa含有群集から任意に種を選ぶことによる不確実性幅である。

### 2．元の4国有線採集区間の内部だけで比較

14元rail駅のうち同じOriginal author source stretch内にある真の駅ペアだけを対象。全source元区間内駅ペアを合わせたとき：
- Canopy height: rail 0.204378 vs local reference 0.456065、差 **−0.251686**。
- Seed mass: rail 0.174328 vs local reference 0.445979、差 **−0.271651**。
- 2D: rail 0.286008 vs local reference 0.667498、差 **−0.381490**。

元**路線区間を1票**にして再重みづけした比較：
| 元4 railway source stretch平均 | Rail | Local source reference | Rail−Local |
|---|---:|---:|---:|
| height | 0.188121 | 0.445538 | **−0.257417** |
| seed mass | 0.169687 | 0.380810 | **−0.211123** |
| 2D | 0.270464 | 0.617361 | **−0.346897** |

元4source stretchesの各1区間を丸ごと除いても、残る11〜12駅の全pair height差（**−0.339714〜−0.199784**）、seed mass（**−0.371118〜−0.174881**）、2D（**−0.525804〜−0.294665**）のそれぞれで負方向維持。**これは原4/4の別原著を4つの独立野外研究として新たに確認したのではなく、1Czech original study内部の4原採集路線分割での感度検証**。

### 3．同駅内のα機能幅は頑健ではない

元source同駅species-equalの2D trait variance trace：Rail 2.092008、局所source参考2.233032、差 **−0.141024**。technical drawによる2.5〜97.5%幅 **−0.340861〜+0.108040** で0をまたぐ。つまり、**駅単独の機能幅が鉄道で必ず縮小する**という強い主張はここからできない。高さ・種子重の**地点間β重心距離**とは別の応答である。

### 4．v28.5 の新しい4独立原著結果と矛盾しない理由

内部の形質ラベル交換帰無（v28.5）はsite×plant speciesを完全固定し、元taxonに付けられた2D高さ・種子重をsource元taxa間で移す。鉄道内の実源元地理区分ペアではheight CWM分化の observed−source trait shuffle はCZ +0.057253／IT +0.038231／HRBA +0.095326／SI +0.056525 と4/4正。この結果は**周辺非鉄道植生との比較でない**。

一方、元v28.6はたったCZ1 original researchでも、本物の付近10km内CzechVeg参考5plot/駅と同種数（形質既知植物）のペアを別の参考グループとして置いて、**絶対地点間beta**を比較する。周辺プールに現れる高さ・種子重の地域差がより強ければ、元鉄道CWMは内部のランダム性よりも空間分化している一方、周辺元植生ほどは離れないことが同時に成立する。**siteから見た二つの分母が違うので、正と負を一つの都合の良い帰無として合算しない。**

### 5．次に必要な反証・原因同定

「鉄道の草丈・種子重は地理的に均質化する」あるいは「鉄道で均質化の方向と草刈り強度・列車風圧が関連する」という強い因果主張は保留。元CzechVeg参考は近傍別出版21件の本物の植物記録だが、**原観察面積/種類と生息地の違い**が残り、元鉄道駅は様々な駅内部zoneを含む、元参考は少数plotの和集合。希薄化は元種数を揃えるが、原plot面積・開放度・都市化・土壌・剪定・観察努力を揃えない。選択・可塑性・遺伝・列車実散布を同定していない。また referenceの WGS84 最大250m位置誤差は100m-cell間の現場重複可能性を完全には消せない。

次段優先：真正同じ100m² rail ballast × 非鉄道道路/礫地/荒地の原現地plotを、同一季節・元流域/土壌/地理/管理条件で少なくとも複数研究に拡張し、元local plant species pool とtraitが同じsource taxonomyで対応する母集団を定義する。その場合初めて rail specific environmental filtering を検証できる。元v28.6は「鉄道が地域種poolと同じように変化するか」の**直接の比較候補証拠1原著**である。

---

## 持続的出力

\`scripts/audit_Czech_original_rail_vs_local_CzechVeg_height_seed_beta_v286.py\`、\`tests/test_Czech_14station_local_nonrail_height_beta_v286.py\`、\`.github/workflows/tetsu-v286-Czech-14rail-70nonrail-height-beta.yml\`。全元素材真正MD5/SHAチェックとMILP/既知2D結果厳格復元、100drawのサイト・元文献調査単位に応じた効果集計。元生plot×speciesと元reference座標・原元CSVは出力しない。
