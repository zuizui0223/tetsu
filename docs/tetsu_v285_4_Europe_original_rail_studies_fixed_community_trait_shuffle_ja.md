# tetsu v28.5 — 四独立鉄道植物調査で「機能的分化」は植物種の空間的配置そのものを固定しても残るか

2026-11-18 JST original source-verified continuous-trait 4-field-study investigation. 研究branch \`research/tetsu-fourstudy-trait-label-null-v285-20261118\`、\`main\`は未変更。**ここまでのv24.5–v24.6の肯定結果と、新しいv28.5の未評価/評価後を分離する。**

## 背景：既にsource真正検証済みの4原著

- CZ Kutlvašr et al. 2024, DOI https://doi.org/10.23855/preslia.2024.247：39駅・5元採集路線区間、2020–22年、公開S1同時期の761taxa、LEDA高さ＋種子重の両方真正結合した種507、元出現対応80.81%。元S1/S2それぞれ source SHA \`df6bc0...\`/ \`00c7cf...\`。
- IT Adorni & Ghillani 2020, 28駅・438taxa、平野11駅vsその他17駅、2012–13年（追加調査一部）。原出版社PDF SHA \`c09bd96ff7307a936030195cc4510bf2591576a713de54019625c04aa9333787\`、LEDA両形質結合271種、出現対応71.09%。
- HR＋BA Jasprica et al. 2017 DOI https://doi.org/10.20302/NC.2017.26.21：15駅・359taxa、元TWS/TSS/MED気候分類、2013–14年。クロアチア＋ボスニアを**1著者研究**と数える。原PDF SHA \`90d489ca8cbc29464d3e6a2b1d8b5d9568c9940e953b1ce720964991366ab956\`、LEDA両形質190種、出現対応58.32%。
- SI Šabić et al. 2026 DOI https://doi.org/10.1007/s44473-026-00155-z：本物の89個**100m道床調査地点**・15本の元国有線区間・173同側BS/BT実植物対、2022–23年。89駅全域ではない。原XLSX 401,279bytes SHA256 \`26be7ca9513dd54136d4ec499d8413c1cb288045f87fe3cf8814731f2edac99f\`。LEDA両形質真正結合224種、元source植物数と各siteの実分布を保持。
- Oldenburg大学LEDA原canopy_height (m) SHA \`cbf5e761390d78a843512dfc1139acfc1154d46d6c1995a827354840a9f78a92\`、原seed_mass (mg) SHA \`34aa12c333fc1ba6de2e560a495ac5d693aa0722a09ecf68418ad94ac1c19cef\`。原学名を二名法で厳格結合、種・変種の重複カウントを防ぐ。3研究で本物の622種が両形質に一致したときの \`z(log(height))\` / \`z(log(seed_mass))\` スケールを**Sloveniaが未知の時点で事前固定**。本v28.5も同じ欧州共通z座標を使い、国ごとの再スケーリングで差を恣意的に消さない。

真正過去original source [v24.5 CI SUCCESS #44010794675](https://github.com/zuizui0223/tetsu/actions/runs/44010794675) / [v24.6 CI SUCCESS #44011602082](https://github.com/zuizui0223/tetsu/actions/runs/44011602082)。

過去source **「現地site×種」をcheckerboardで交換し、各元site種数と研究元地域内の種頻度を固定する**帰無との比較では、地点間LEDA2D CWM距離が元source帰無を超えるのは4原著すべて正（CZ +0.044116 / IT +0.038670 / HR+BA +0.074704 / SI +0.064977）。これは本物の鉄道内での**種の配置の異質性**に対する条件付きベースラインであり、地域自然植生より機能が異質という証明ではない。

## v28.5では異なる帰無を追加する

前の帰無は植物種自体を別のsiteへ動かす。今回は元著者のplant occurrence matrixを完全に固定し、**個々の元speciesに付属する草丈・種子重の2変量traitのペアだけを別の元speciesへ割当て直す**。研究全体の2D trait平均・共分散・値の分布、各source siteの種数と植物分類群同一性・種出現頻度・その群集の種ペア共有・元著者地域群は**完全固定**。元taxonomic Jaccardは観測sourceそのまま。

2つの別の事前帰無：
1. 元研究の真正LEDA両形質既知taxa全部の間でpaired 2D trait vectorを独立にシャッフル。植物種の地理的配置と形質ラベルの無関係を仮定。
2. 元speciesの元site frequencyから5分位を作り、分位内に限りpaired trait vectorをシャッフル。「多くの駅で普通に出る種ほどどの形質か」という大まかな関係を保存してsource頻度−形質の結合を制限する。ただし完全な厳格な頻度条件付きモデルではない。

別々の2つの固定seed×各299独立label draws。MCMCの長期記録ではなく、各ラベル付替えは明示的に独立したpseudo-random permutation。ただし**同じ原研究×元siteという植物生態学のデータ自体は独立に収集し直していない**。Permutation quantilesは「無作為現場実験の信頼区間」ではなくtechnical assignment null。

真の応答：
- 種のpresenceを等重みした各siteのCWM2Dベクトル間Euclidean平均距離。主解析は元著者の**同地域区分内だけのsiteペア**で、地理をまたぐ差ではない。
- 研究全体の全siteペア2D CWM距離（副解析）。
- height/seed massそれぞれのsource CWM差（2D方向の分解）。
- 各真siteの2D traitsのpresence-equal variance trace（**α機能幅**であり、地点間β距離とは別指標）。

### 科学的判別

- 新帰無でも観測CWM間距離−Trait-label nullが正：固定されたtaxonomic assemblagesの元**traitsと種の場所**に、ランダムなtrait-species対応を超える空間的構造がある可能性。ただし気候、近接原種プール、生物地理、系統、観測欠損を含み得るので「鉄道だけの機能的適応選択」ではない。
- 新帰無で負：実際のtaxonomic assemblagesの間で、species trait valueの分配はランダムより均質なCWMを作っている。**観測配置を残した条件付き機能収束**だが、鉄道vs非鉄道の外部差ではない。
- 4原著中の符合一致を報告しても、4出版・非無作為に選ばれたヨーロッパ研究である。HRとBAは国2つでも一元研究、SIは駅ではなく100m地点であり、\`171 independent stations\`としない。未知の熱帯・乾燥地域・別種プールへ一般化するには同条件の真正外部非鉄道plotが必要。
- LEDA北西欧由来のcoverage偏り（HRBA 58%）が無視できない。元の欠損種は\`0 height\`や\`1 mg\`を補完しない。二名法レベルの種中央値traitsは現地植物個体の表現型/遺伝子変化でない。

この試験は、別の公開研究が用いるtrait-label permutationと同じ「種の出現構造を固定しtraitとsite occupancyとの対応を切る」概念と整合するが、観測単位・source種プール・traits欠損を一致させない限り他研究の推定係数と混ぜない。たとえばtrait randomizationでコミュニティ構成を保持するnullは Swenson 2014などの機能多様性研究で使われており、\`https://doi.org/10.1111/ddi.70270\` の2026年公表研究にもtaxon-swap方式が記載されている。

## v28.5 真正SOURCE完走後の新しい肯定結果（2026-11-18）

**[真正4独立原著SOURCE GitHub Actions SUCCESS #44296227526](https://github.com/zuizui0223/tetsu/actions/runs/44296227526)**。新しい数学QA5/5・Czech original S1/S2 SHA・Italian original source museum PDF SHA・HR/BA原著PDF SHA・Slovenia original publisher XLSX SHA・Oldenburg LEDA height/seed原表 SHAがすべて成功した。過去のsource v24.5／v24.6に報告済みの4つの本物の全siteペアCWM値（CZ 0.16501788905535814、IT 0.2014200675821101、HRBA 0.2232461322998969、SI 0.35003008736045554）を、作者原sourceから **1e−10以内で全て再現**。国ごとにtrait z-scoreを再標準化せず、CZ/IT/HRBAの622種を用いた前回のzlog平均・SDも **1e−12以内**で元値と照合できた。初期のCIは過去Github Actions集計artifactをダウンロードする権限エラーで失敗したが、今回は**全元著者sourceからtraitの初期較正を再生成する方式**へ変更し、別の元集計artifactを取得できなくても独立source検証を成立させた。失敗は元生物解析の不一致ではなく、旧集計artifactの認証経路だった。

### 1．全species/site植物群集を固定し、種に付くLEDA二形質のペアだけをランダム化する

| 真正独立field study | 本物の駅・道床地点 | 両形質既知種 | 2D CWM距離source元実測（同著者地理内） | species trait-label null期待値 | **実測−trait-label null** | 草丈1Dの実測−null | 種子重1Dの実測−null |
|---|---:|---:|---:|---:|---:|---:|---:|
| CZ Czechia Kutlvašr2024 | 39駅 | 507 | 0.154620 | 0.093251 | **+0.061368** | **+0.057253** | +0.026440 |
| IT Parma Adorni2020 | 28駅 | 271 | 0.203252 | 0.154865 | **+0.048388** | **+0.038231** | +0.027471 |
| HR+BA Jasprica2017 | 15駅（2国・1原著） | 190 | 0.209951 | 0.146172 | **+0.063779** | **+0.095326** | −0.028017 |
| SI Šabić2026 | 89実100m道床source地点 | 224 | 0.336127 | 0.286862 | **+0.049265** | **+0.056525** | +0.004771 |

**全4独立原著で元地点内地理群の2D CWMβ距離が、元site/taxon assembly完全固定・paired height/mass trait-label shuffleより大きい**。注目すべきは、**草丈軸は4/4原著で強い正方向**であり、種子重はHRBAで負、SIでほぼゼロという点。これは地理に沿って種々の大きさの植物が入れ替わることが、元研究内でランダムな種→形質対応を超えて成立することと整合する。

技術的なsource permutation上側順位（各モデルの2固定seed×299draw=598本の独立**形質ラベル組替え**であり野外実験p値ではない）：2D距離はCZ 0.00167, IT 0.00668, HRBA 0.02170, SI 0.04007。草丈のみCZ 0.00167、IT 0.01836、HRBA 0.00167、SI 0.01503。種子重はCZ 0.0050、IT 0.0518、HRBA 0.8915、SI 0.3990と非一様。これは確証的なmultiple-testing補正済p値ではない。元地理group数はCZ5、IT2、HRBA3、SI15であり、それぞれの元siteペアは互いに依存する。

### 2．出現頻度の近いspeciesでのみ二形質を入れ替えるより保守的な帰無

各元speciesの本物のsource研究内出現site数で上位/下位5等数binを一意に定義し、同じbin内でだけ2D height/seed paired trait vectorsを交換する。元サイト種数・species identity・種ごとのexact site occupancy、元site×species co-occurrence、元2D trait相関はそのまま保つ。**出現頻度とtraitの関連を分位レベルでは保存する**が、5分位以内の厳格frequency同値は固定していない。これはsource群集の種名の局所構造がどのくらい形質との特異的対応をもつかの補助確認。

| 原著 | 源元同地域内2D CWM距離実測−frequent-bin trait-shuffle null | 草丈1D実測−null | 種子重1D実測−null | 2D technical上側順位 |
|---|---:|---:|---:|---:|
| CZ39駅 | **+0.058534** | **+0.054819** | +0.025337 | 0.00167 |
| IT28駅 | **+0.048722** | **+0.038810** | +0.026291 | 0.01169 |
| HRBA15駅 | **+0.070119** | **+0.099189** | −0.024478 | 0.00501 |
| SI89道床地点 | **+0.080191** | **+0.078781** | +0.023874 | 0.00167 |

この別定義でも**2D・草丈は4/4原著で正**。特にSIでは全源種を一緒にshuffleした初期null +0.0493からfrequency-bin +0.0802へ増加。これは多くのサイトに出る種と稀な種の間の元形質−占有率関係だけでは、元現場の形質空間的分化を説明しきれないことと整合する。元の頻度binの定義を固定しているだけなので、真の地域種pool・plant climate/soil/management・phylogenetic similarity・缺測biasを消した結果とはならない。

### 3．同一サイト内部の2D形質幅αはβと別の応答

源site内height+seed massの2D variance traceとtrait-shuffle null平均との差：

| 出版研究 | 全species trait-label nullへのsite α trace差 | 出現頻度5分位内shuffleへの差 |
|---|---:|---:|
| CZ | −0.02823 | −0.10306 |
| IT | −0.25224 | −0.03038 |
| HR+BA | −0.31447 | −0.03674 |
| SI | −0.42385 | −0.00716 |

Unrestricted trait shuffleではα traceの符号が4/4負で、一見機能幅縮小に見える。しかし**site出現頻度bin別にすると、特にIT、HRBA、SIは負の効果量が大幅に小さくなる**。種子/高さと通常のspecies occupancyの相関によるαの見かけの狭さの可能性がある。従って**鉄道で個々のplant communityが普遍的に機能幅を狭める**という強い一般則にはしない。各研究の多重分類群とsite集団構成の交絡・trait missingnessが残る。

### 4．新しい論文の中核と未解決

高い信頼度の観察：鉄道に共通した人工環境が反復されても、その場にある植物種の**草丈の地理的対応**は研究内でランダムな種→形質割当を超える。鉄道のどの場所にも同じ植物の背丈が再現する、という単純な「均質化」予測とは異なり、**regional plant pool × site filters / plant life-history sorting**による*局所特有の機能状態*という解釈と整合する。この解釈は主張可能だが、どのフィルターが草丈を変えたかという具体的因果機構はまだ未同定。

現地plant height実測や遺伝子、剪定日、車両風圧、種子実測到着は取得していない。LEDA元種中央値のheightを現地の可塑的草丈や進化的耐刈取り形質と誤認しない。また、**自然植生vs鉄道の同季節・同面積の対照**の分布を取り込めていない。したがって「鉄道特異的に草丈の地域差が生まれる」「鉄道で進化的適応する」という因果主張はできない。

*先行研究との結合：* 元4/4のβCWM分化は、v27.0–v27.9の「鉄道内での同実現地site BS/BT植物種の強い局所共有」、v28.0–v28.4の「別地点同一路線の追加taxonomic均質化は弱く、実散布能力や実貨物ODの効果は未認証」と同時に成立する。異なる空間/群集単位の機構を一つの効果量に混ぜない。

**暫定投稿の一文案**：
“Repeated railway habitats do not imply uniform functional states: spatially differentiated plant stature recurs across four independently surveyed European railway floras, even when plant community memberships and species occurrence frequencies are held constant.”

ただし表題/要旨でより大きな全球の因果フィルタと見せかけず、真の4原著/5国、各publicationを独立単位とした探索的scopeを明示する。

---

## 外部真正独立鉄道＋非鉄道data入手監査

Penone et al 2012 France71実site：元著者論文は2線、実71 site×5m²/186種で接続/跨線橋断絶を検証した。ただし公開MNHN博士論文 Appendix Bは種別のNumber of observationsであり、その意味が71真site出現頻度かは未認証。元71site×186taxa＋GPS＋駅/跨線橋の中断位置と植生環境のfull matrixは取得できず、tetsu既存v14.5/14.6では6種の周辺度数から同線差の符号すら一意に推定できないことをsource conditional mathematical witnessで検証済み。これを今回の「4独立研究」へ5個目として捏造しない。

Májeková et al 2025 Slovakia30実駅、https://doi.org/10.1016/j.flora.2025.152780：413実植物taxa（258 Native/155 Alien）、元論文がCargo trainsとAlien station α-richnessの肯定的な関連を報告。これを駅間種名Jaccardや元列車貨物の種子arrivalと混同しない。元全30site×413taxaと駅ごとの個別貨物列車本数の同一source・利用許諾は未確認。

新たに判明した2026年Czech rail–agriculture 4site×5habitat study (Winkler et al, Land DOI https://doi.org/10.3390/land15040523)は、83植物taxa・10m² plot・2021年3回調査で railway yard/embankment/rail land/field margin/arableの同時現地対照という比較として非常に有益。しかし出版元のData Availability Statementに**「data are not available to the public」**と明記、同時刊行Supplementary Table S1は植物taxon名と略語であり、真正plot×speciesの出現行列ではない。元4site×plant×habitatを取得・配布したかのように主解析へ加えない。

## 実装ファイル

- \`scripts/audit_Europe_4independent_railstudies_fixed_community_trait_label_null_v285.py\`
- \`tests/test_Europe_4study_fixed_taxa_paired_trait_label_null_v285.py\`
- \`.github/workflows/tetsu-v285-four-original-rail-studies-trait-label-null.yml\`

各publisherの原資料SHAとv24.5/v24.6過去original baseline CWM値を厳格に再現できなければ失敗を返す。元XLSX/PDF、EPPOコード・生の元site×species / 元GPSとLEDA元のtrait一覧は出力artifactへ含めず、各original publicationの現地単位ごとの集計JSONとnull QAだけ保持する。
