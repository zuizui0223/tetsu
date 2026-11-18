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

## 外部真正独立鉄道＋非鉄道data入手監査

Penone et al 2012 France71実site：元著者論文は2線、実71 site×5m²/186種で接続/跨線橋断絶を検証した。ただし公開MNHN博士論文 Appendix Bは種別のNumber of observationsであり、その意味が71真site出現頻度かは未認証。元71site×186taxa＋GPS＋駅/跨線橋の中断位置と植生環境のfull matrixは取得できず、tetsu既存v14.5/14.6では6種の周辺度数から同線差の符号すら一意に推定できないことをsource conditional mathematical witnessで検証済み。これを今回の「4独立研究」へ5個目として捏造しない。

Májeková et al 2025 Slovakia30実駅、https://doi.org/10.1016/j.flora.2025.152780：413実植物taxa（258 Native/155 Alien）、元論文がCargo trainsとAlien station α-richnessの肯定的な関連を報告。これを駅間種名Jaccardや元列車貨物の種子arrivalと混同しない。元全30site×413taxaと駅ごとの個別貨物列車本数の同一source・利用許諾は未確認。

新たに判明した2026年Czech rail–agriculture 4site×5habitat study (Winkler et al, Land DOI https://doi.org/10.3390/land15040523)は、83植物taxa・10m² plot・2021年3回調査で railway yard/embankment/rail land/field margin/arableの同時現地対照という比較として非常に有益。しかし出版元のData Availability Statementに**「data are not available to the public」**と明記、同時刊行Supplementary Table S1は植物taxon名と略語であり、真正plot×speciesの出現行列ではない。元4site×plant×habitatを取得・配布したかのように主解析へ加えない。

## 実装ファイル

- \`scripts/audit_Europe_4independent_railstudies_fixed_community_trait_label_null_v285.py\`
- \`tests/test_Europe_4study_fixed_taxa_paired_trait_label_null_v285.py\`
- \`.github/workflows/tetsu-v285-four-original-rail-studies-trait-label-null.yml\`

各publisherの原資料SHAとv24.5/v24.6過去original baseline CWM値を厳格に再現できなければ失敗を返す。元XLSX/PDF、EPPOコード・生の元site×species / 元GPSとLEDA元のtrait一覧は出力artifactへ含めず、各original publicationの現地単位ごとの集計JSONとnull QAだけ保持する。
