# tetsu v28.7 — チェコ鉄道14駅 vs 周辺CzechVeg植生：HABITAT明示・調査年代・面積・標高記載の厳格対照

更新：2026-11-18 JST。ユーザー研究 tetsu の真正原データ再解析。**原著真正SOURCE CI 成功（2026-11-19）**： [Actions #44433101484](https://github.com/zuizui0223/tetsu/actions/runs/44433101484)・[派生集計artifact](https://github.com/zuizui0223/tetsu/actions/runs/44433101484/artifacts/14487129986)。元CzechVeg 117,739 metadata + 2,510,340 source species-plot occurrences MD5、Preslia S1/S2・SR70駅代表座標・Oldenburg LEDA元trait SHAと5数学テストすべて成功。**非鉄道元現地情報を厳格化すると必要なsource駅/出版元plotのsupportがなく、対応する新しいheight/seedの生態効果は判定不能**。独立調査原sourceの個別学名×site/個別GPS/plot IDを再配布しない。

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

## v28.7 新しい真正原資料・生息地対照の支持域

[原著再取得・数式テスト・元基準値再現を含めたSOURCE SUCCESS #44433101484](https://github.com/zuizui0223/tetsu/actions/runs/44433101484)。5元合成テスト成功、Preslia S1/S2、Czech Veg元header 117739・出現行 2510340、LEDA height/seed、SR70原registry derivativeがチェックサム認証を通過。元v28.6の観測2D CWM betaは**−0.42673326131926104**を1e−9以内で厳密再現し、元MILP source optimal status \`EXACT_SOURCE_MILP_GEO_UNIQUE_OPTIMAL\`・14真駅/70真CzechVeg source plot/原4 source railroad stretches・reference 21原著を再構成した。初回v28.7のWorkflowは元−0.426733261を本文表示用の小数6桁−0.426733と比較したため意図せず失敗しており、全桁のsource v28.6真正値へ修正後に元生物同定・原著 source は一致して成功した。元情報を捏造して修正していない。

### 原CzechVeg候補（本物の元出版原植生）

元のtrait QA合格16rail source駅に対応する元参考植生は**339個の元distinct plot IDs**、これが複数駅の10km圏に重複して現れるため**560 station×plot候補機会**。元339本の実 \`HABITAT\`、\`YEAR\`、\`SURF_AREA\`、\`ALTITUDE\` を、種組成・CWMを一度も対照選定に使わずに検査。

- 元HABITAT欄が欠損ではない：**163/339**。元HABITATに開放地/草地/ルーデラル等の明示：**22/339**、砂・礫・岩石/採石場等の明示：**8/339**。両義的な森・湿地と開放地の語を同時に含む元記載：3件（開放地候補から除外）。**記載がない**ことは現地が開放地ではないことと同じではない。
- 元現地観察YEAR>=2020：153/339。元reference SURF_AREA 16–400m²：228/339。元reference ALTITUDE明示301/339、元傾斜数値270/339。**元rail駅側標高が著者元sourceにない**ため、ALTITUDEが記録されても駅側との差を計算できず、標高一致とは呼べない。
- 原70reference plotのsource \`HABITAT\` 記載は **16/70**、うち**開放地明示2/70、明示的砂礫岩石0/70**。元採取YEAR>=2020は20/70、元plot面積median 25m²（range4〜625m²）、標高記録49/70（median287m）、傾斜記録42/70、GEOLOGY文字記述1/70、SOIL文字記述0/70。つまり元14駅と周辺5plotを組んだ比較は**鉄道の基質と同等の砂礫・乾燥・開放対照を作ったことにならない**。

### 事前に固定したsource PANELを厳格化したときの本物の比較支持

Source前処理は16候補rail駅、5個/駅の元植生参照plot、最低2つの独立したCzechVeg原出版研究/駅、元reference座標precision0–250m、原10km圏、70元referenceplotが全て異なる元IDと0.001°rounded cellを保つMILP。以下「事前条件を満たす候補駅」の数は**元source metadataのみの上界**であり、各plotの100m-cell衝突が全て解けた実際のMILP採択駅数とは区別する。

| 原source参照候補を選ぶ非植物基準 | 元候補 station×plot機会 | 源5plot・2研究を満たす元駅（採択駅数の上界） | source3元rail stretches以上＋10駅以上の検証 |
|---|---:|---:|---|
| 従来nontransport文脈（reference全地域） | **560** | **16/16** | 原MILP **14駅・4路線・70plot実現** |
| YEAR>=2020 | 268 | **6/16** | **判定不能** |
| HABITAT開放地/ruderal明示 | 37 | **2/16** | **判定不能** |
| 開放地＋YEAR>=2020 | 33 | **2/16** | **判定不能** |
| 開放地＋source SURF_AREA 16〜400m² | 33 | **2/16** | **判定不能** |
| 開放地＋YEAR2020+＋面積16〜400m² | 30 | **2/16** | **判定不能** |
| HABITATに砂/礫/岩石/採石場明示 | **11** | **1/16** | **判定不能** |
| 開放地明示＋reference ALTITUDE 記録 | 36 | **2/16** | **判定不能** |

後の厳格7パネルは源実駅・元rail stretches支持閾値を満たさないため、**CWM height/seed mass βの符号も大きさも再計算して『効果ゼロ』にしていない**。この結果は「開放地同士で鉄道群集が均質化しない」を示す否定生態効果ではなく、**厳格な対照でその生態効果を識別するための元HABITAT×trait×geography同時supportが足りない**というsource情報の否定である。確認済みの先行v28.6「従来広域非交通文脈での鉄道−reference β差」自体は変わらない。

### 研究の結論・投稿で守るべき境界

- 元14/70で「元鉄道駅間の草丈・種子重のCWM地域差が元参考CzechVegで選ばれた異なる植物群集間より小さい」ことは再現性がある。しかし元referenceはsourceのHABITAT open補足がわずか2/70で、明示砂礫は0/70。**鉄道バラスト環境特有の機能均質化**という因果推論はまだ支持不能。
- 原\`SYNTAXON\`で乾燥草地や砂地を選べば見かけの参照supportは増えるが、その分類は元plant communityから定義されており、同じ植物種・形質の群集類似性を検証する参照選定に使うと循環する。既存v19.6の同じRAIL/非交通元plotにおける\`SYNTAXON\` reference vs abiotic\`ALTITUDE\`参照の不一致はこの問題を先行指摘した。今回は元HABITATに限定するため仮説の検定は正しく\`not identified\`となった。
- 元\`HABITAT\`の空欄は本当に閉鎖森林・湿地の測定ではない。今後は独立の現地基質・オープンハビタット分類（元々の植物群集から定義されない外部GISまたは現場確認）が必須。鉄道側現地駅標高の真正DEM抽出が可能なら、地形差の一部を独立に制御できる。ただしSR70駅代表点2026は元植物調査区GPSではない。
- 次の優先データは元鉄道100m²バラストplot×本物の非鉄道100m²開放・砂礫地plotが元同年度・同季節・近傍・同観察努力で独立した2以上のpublication source。路線内外の交通種子散布・現地風・草刈り・遺伝適応を元CWM数値から逆推定しない。

---

## 実装と出力

- \`scripts/audit_Czech_original_14station_reference_habitat_era_area_altitude_v287.py\`
- \`tests/test_Czech_reference_original_HABITAT_area_year_v287.py\`
- \`.github/workflows/tetsu-v287-CzechVeg-original-habitat-area-era-controlled-rail-reference.yml\`

すべてのsource元Preslia/CzechVeg/Oldenburg/SR70を厳格ハッシュ確認し、元v28.6の2Dsource baseline −0.42673326131926104が1e−9の精度で再現しなければそこで停止。生物の結果はsource field supportがあるパネルにのみ記載する。
