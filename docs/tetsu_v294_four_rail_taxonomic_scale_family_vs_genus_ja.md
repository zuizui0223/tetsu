# tetsu v29.4 — 鉄道群集の草丈分化が「どの分類階層」で失われるか

2026-11-19 JST。**新しい実測値はSOURCE CIとTaxonomy QAが成功してから追加**。4真正出版研究の元植物出現 site×species 行列と各研究元地理区分、元LEDA2形質を変更しない。

## 背景：比較対象は3つの異なる帰無

- **v28.5 種出現頻度5分位内でplant height/seed massの2D paired vectorsを交換**：4独立欧州RAIL研究すべてで元地理グループ内の草丈観測CWM β距離−nullが正。4publication等重み平均+0.0678997。ただしこの帰無は属をまたいでtraitを交換する。
- **v29.3 同じ植物属の中だけで2D traitを交換**：全4元source原著を再取得してv28.5を完全再現後に検証。4出版等重みの高さ超過が+0.0062803、Sloveniaの符号は−。元source元trait-known speciesのうち属内交換できない1種属もある。
- **本v29.4 FAMILY（科）の中だけで2D traitを交換**：属が違っても同じ科なら交換できる。科をまたぐtrait差を保持して、属の差は動かす。従来の「属をまたぐ広いランダム化」と「同属の中だけ」の間に、科を保持する中間帰無を入れる。

属や科は **植物の実分子系統樹ではない**。属内/科内randomizationで効果が変わっても「系統遺伝率」「選択圧」の因果推定にはならない。頻度分位を固定するv28.5とfamily/genus固定する帰無は全て同一の現実のcommunity/LEDA profileを対象とするが、ランダム化空間が異なるので、超過の差を「属が何%説明する」式の分散寄与率にしてはいけない。

## 原source publisher・観察単位・外部分類データ

既存の原著 source-authenticated スクリプト \`validate_Europe_original_4rail_studies_genus_trait_label_null_v293.py\` の処理を再使用する。

- CZ Kutlvašr2024 Preslia原著 S1+S2：39駅全域（2020–22）。
- IT Adorni2020：28駅全域。
- HR+BA Jasprica2017：15駅、**2国で1著者出版**。
- SI Šabić2026：**89個の100m鉄道道床区画**、駅全域ではない。
- 真正LEDA canopy_height[m]、seed_mass[mg]元2テーブル（原SHA）、元欧州3原著622speciesに基づくz(log height,seed mass)標準化を固定。元全4studiesのfield siteは171だが**171 equivalent stationsではない**。region/line/groupは元著者定義を保持。

新しい科分類は **GBIFの公表 Species API** https://techdocs.gbif.org/en/openapi/v1/species の \`GET /v1/species/match\` で元植物真正 **GENUS** のみを \`rank=GENUS&kingdom=Plantae&strict=true\` の条件で照合。返答が \`EXACT\`、rank GENUS、Plantae、confidence≥90、照合genusが元の同一名、familyが明記のときだけfamilyを採用。別の属へ書き換わるaccepted synonym、fuzzy、HIGHERRANK、非植物、family不明は**不明**として固定し、勝手に分類を創作しない。GBIF taxonomyは実行時点の参照でありバージョン固定された分子系統樹ではないため、照合日・元属数・exact分類属数・照合マップの正規化SHA256 digestを保存する。元学名一覧/元source地点×種はartifactへ出さない。

## 分析前source taxonomy support gate

4実研究ごとに以下を凍結：
1. 元site×trait-known plant incidenceの **70%以上** が厳格照合familyへ割当て可能。
2. その研究の**元2属以上が含まれる科が8科以上**。
3. 元各site×species出現の **30%以上** が2属以上の科にある（family conditionで交換可能）。
4. 元著者区分する地理グループで、(1)70%と(3)30%の双方を満たすgroupが **2以上**。

満たさなければ \`NOT_IDENTIFIABLE_ORIGINAL_FAMILY_TAXONOMY_SUPPORT\` で止める。未分類の元speciesはもとのsource観測に残し、元trait値固定の1種単位の帰無にして群集から消さない。これにより、同じsite×speciesの元community CWM observed値を維持する。分類できない数をsourceで公開しないまま推定を行うことは禁じる。

## 計算と事前解釈

原v29.3の真正集計ファイル \`data/derived/tetsu_v293_4European_original_rail_genus_null_summary.json\` を照合し、各研究で元source v29.3のwithin-genus source結果（完全な真元tableと同じsource 種の順序）を再現できなければ今回の新解析を**停止**。

科内だけ元height+seedの2Dペアを独立2seed×各299回=**598技術的シャッフル**。科内種でないものは元source値固定。各siteのheight/seed CWMと、同元著者地理内の全siteペア距離を再計算し、現実の元observed−family-nullを求める。2D / height1D / seed1Dを別々に集計。4独立publicationが少なくとも3研究source taxonomy supportを満たす時のみ、研究等重み記述平均・1研究除外平均を出す。technicalシャッフル順位は野外統計のp値ではない。

- **Family-nullでも4研究の大きな正値が残る**場合：属間（同科内）plant height差と元地理配置に情報があり得るが、family trait差自体は固定。バラストで選択が起きた証拠ではない。
- **Family-nullで弱くなるがgenus-nullより残る**場合：元trait構造は科の組成を固定すると縮小するが属の分布を移動させる余地はある、という分類群階層への依存と整合。
- **Family-nullとgenus-nullの両方がほぼ小さい**場合：元広いtrait-null正効果は科以上の種プール/環境/生活史の共変動に依存する可能性があり、「鉄道内で独立に進化した低草丈」という解釈は弱まる。
- **GBIF照合が十分でない**場合：元studyを削除/効果ゼロとせず、分類データ入手がボトルネックとして記録する。

## 区別して残す科学的問い

今回測るのは「同一の人工ハビタット内で、**どの尺度の分類学的な先験構造**が種と形質の地理的分化を担うか」である。列車運行頻度、種子輸送、風圧・草刈りへの進化的適応は現在の群集種中央値から測れない。鉄道の因果的な機能収斂には、実稼働バラストと真の同等面積・同調査年・開放非鉄道砂礫地という対照が引き続き必要。

## SOURCE実装

\`scripts/validate_Europe_original_4rail_family_genus_taxon_null_v294.py\`
\`tests/test_Europe_original_family_classification_null_v294.py\`
\`.github/workflows/tetsu-v294-original-EU-four-rail-family-trait-null.yml\`

元植物種×site、正確な地点と所有者に関わるGPS情報、LEDA原測定値、taxon listは元author原sourceからGitHubへ再配布しない。分類APIレスポンスの元全文も配布せず、classificationsのaggregate品質とSHAだけを記録する。
