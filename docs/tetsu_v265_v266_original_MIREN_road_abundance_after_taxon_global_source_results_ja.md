# tetsu v26.5–v26.6 — 公開MIREN山岳道路再調査による初期優占度の反証検定（道路は鉄道ではない）

原データ実行・検証日：2026-11-18 JST。対象は現地調査された MIREN 道路（ROADS）であって鉄道（RAIL）ではない。著者元データ DOI https://doi.org/10.5281/zenodo.15749536（2025年v3、2007–2023年）、方法とデータ論文 https://doi.org/10.1111/geb.70284、標準化調査プロトコル https://doi.org/10.1002/ece3.8590 。原調査観測はMIREN研究者のものであり、tetsuが野外で新たに採集したのではない。

真正一次資料：原CSV 17,899,241 bytes、MD5 b640195cc44f252ccb6ddb28eb6fff13、source実記録167,804行。元v9.8の署名付き成功run #43212464976 の原ファイルをGitHub Actionsで再取得して毎回MD5とbytesで照合。原CSV・元plot×taxon・元GPSはこのPRへ再公開しない。

- v26.5 元source [Actions SUCCESS #44235004684](https://github.com/zuizui0223/tetsu/actions/runs/44235004684), [派生結果artifact #14397041743](https://github.com/zuizui0223/tetsu/actions/runs/44235004684/artifacts/14397041743)。合成数学5/5。
- v26.6 元source [Actions SUCCESS #44236248506](https://github.com/zuizui0223/tetsu/actions/runs/44236248506), [派生結果artifact #14397880667](https://github.com/zuizui0223/tetsu/actions/runs/44236248506/artifacts/14397880667)。source4件合成検算・厳格な2012年以降コホートのチェック通過。
- 本PR研究ブランチ research/tetsu-global-road-abundance-control-v265-20261118、main未変更。

## 1. 鉄道の独立反復へ加算しない理由と生態学的問い

カシミール元鉄道2014/17年のv26.1–26.3では、初期種別量的順位は種同一性・種×habitatの記録傾向を加えても、その種が2017年に再記録されるかを多少改善した。一方、その結果が鉄道固有のストレス/適応に由来するとは未証明。そこで独立の人工交通回廊である**真正MIREN山岳道路**でも、道路際の初期量的多さが後年再記録の予測に寄与するか、同地点の内陸対照との「初期量的順位への依存の差」を検定した。

ReSurveyGermanyの公式公表データは2022年の論文 DOI https://doi.org/10.1038/s41597-022-01688-6 に記載された大規模反復種別被度データだが、2026年11月18日時点のiDiv原配布ページ https://idata.idiv.de/ddm/Data/ShowData/3508?version=0 は **"Currently not available"**。元科学論文は7,738反復plot、23,641 observation、RS_PROJECT/RS_PLOT/YEAR/ Cover_Percを記載。対応する公開GitHubのRコード https://github.com/idiv-biodiversity/Read_ReSurveyGermany は読込方法のみで原二つのCSVを含まない。元正規CSVが得られず、さらに真のrailway ballastと確認できる区画ラベルも未検証。ReSurveyGermanyを取得済み鉄道実証として追加しない。初期2022年preprintに記された外部cloud公開URLは現在の再配布元として確認できておらず、原本にアクセスできたとも書かない。

## 2. MIREN v26.5：原調査の「初期量的記録→後年再記録」の同一地点での比較

source Original Region×Road×Site×Plot×Yearを厳格に固定し、記録年ごとに Plot1道路際 vs Plot3内陸（ともに100m²）を使用。Plot1A/3Aの置換された調査区やPlot2は主解析に混ぜない。元種名 Accepted.Name.MIREN の完全同一表記で全種を比較。全167,804原source行から厳格Plot1/Plot3に対応する120,506観測行、完全重複したsite×plot×year×taxonを5,152行まとめた。**1,876同じSite×YearでPlot1+Plot3の両方がある組、484地理地点で複数回調査が存在、11元Region code**。

**決定的な検出制約**：484再訪地点すべてで初回→最終回の元Protocolバージョンが変わっている（v2007/2012/2017→v2012/2017/2022等）。初期v2007/08は外来種を優先して記録したため、主解析は**初回観測年が2012以降**の270地点だけを候補とする。2Plotそれぞれに初回の再記録種と未再記録種の両方があり、初回Abundance全種が数値で確認できるのは252地点。4地点以上が揃う本物の道路を各Regionで最低2本保持する支持基準で、**7 Regionコード・20道路・251地理地点**を主結果とする。

毎個体の生残ではなく、初回plotに登録され後年同じ元plotで再記録された種と未記録種の*初期数値 Abundance=1,2,3（個体数階級）*について、plot内順位AUCを評価：AUC= P(再記録種の初回Abundance>非再記録種)+0.5P(同値)。実測植物個体の生死や繁殖成功は追跡していない。Region内-road内-site内 AUC -> そのROAD内のpaired difference -> Region内roads等重み -> Regionコード等重み。ペアの各Siteを一票とする。原sourceの地域codeを7つの独立な国・無作為実験とは扱わない。

| 原SOURCE計算 / 响応 | 調査単位 | 平均AUC差（Plot1道路際 − Plot3内陸） | 正のRegion / 負のRegion |
|---|---:|---:|---:|
| **主解析：2012年以降初調査、Abundance 1–3階級** | 7Region・20道路・251地点 | **−0.023337** | 2 / 5 |
| 同一248地点：Abundance階級 | 7Region・20道路・248地点 | **−0.021181** | 2 / 5 |
| 同一248地点：Cover 1–8階級 | 同上 | **−0.040503** | 2 / 5 |
| 初期2007/08も混合した全初回–最終 | 11Region・31道路・406地点 | −0.008071 | 6 / 5 |
| **初期2007/08だけ**（当時は外来種中心） | 5Region・13道路・145地点 | **+0.010290** | 3 / 2 |

主解析のRegionを一つずつ除くと等Region平均のAUC差は −0.03267〜−0.01162。これはロバストネス診断であり、世界の交通回廊を母集団とした推定CIではない。初回プロトコル時代を混ぜると符号の集計が動くので、旧プロトコルを独立した全球道路結論に混ぜない。初回階級値は科学的に解釈できても時代間の検出・方法変化が残る。

**観測判定**：一般的な初期量的順位→後年再記録は道路と内陸の両方に存在しうるが、単純AUCで「道路際により強い」を支持しない。ただし道路と内陸の植物相・基礎種数・元測定階級は異なり、鉄道同条件との効果量差を直接求めていない。

## 3. MIREN v26.6：種同一性を固定したあとの量的追加予測力

厳密同じ上記**251地点・20真正道路・7Region・初回2012以降**のsource母集団。元初回に出現したtaxon×plot riskイベントは**14,682**、初回登録TaxonNameは**1,925種類**。元道路全区画を丸ごと1本ずつleave-one-original-road-out、同じ道路の全source地点・両plotの後年再記録ラベルを学習から遮断し、source元Road内のSite等重み→Region内Road等重み→Region等重みで予測log-lossとBrierを評価。連続数値の学習データ正則化は固定L2 logistic C=1、試験データでチューニングしていない。

事前固定した比較：
- M2 = Region+plot+source初回年+調査間隔+taxon identity。
- M3 = M2 + 初回plot内数値 Abundance の種間midrank [-1,+1]。
- M2Z/M3Z = 元taxon×Plotの特異性まで固定してから初回rank追加。
- M4Z = M2Z + そのtaxon×Plotの **train-ROADSだけ**から計算した平均rankを差し引いた、留保道路での「通常の多さからの偏差」。訓練に元留保道路の種別未来ラベルを使わず、未知taxon×Plot組み合わせのセンタリング残差を0とする。M4ZはM3Zと異なる推定対象であり「全種で普遍的に道路際ストレスを測定した」ではない。

| 元7Region等重みのleave-ROAD-out log-loss利得 | 道路際 Plot1 | 内陸 Plot3 | 道路際 − 内陸 | 右側差が正の地域 |
|---|---:|---:|---:|---:|
| M2→M3 初期raw順位追加 | +0.019004 | +0.022622 | **−0.003618** | 3/7 |
| M2Z→M3Z 種×zoneの上にraw順位追加 | +0.018119 | +0.021411 | **−0.003292** | 2/7 |
| M2Z→M4Z 種×zoneの上にtrain-centered順位追加 | +0.013227 | +0.010699 | **+0.002528** | 6/7 |

Brierでも同じ符号：
- M2→M3：道路際 +0.007638 / 内陸 +0.008672。
- M2Z→M3Z：道路際 +0.007237 / 内陸 +0.008413。
- M2Z→M4Z：道路際 +0.004857 / 内陸 +0.003778。

**解釈：初回量的順位は種同一性・生息地を越えて道路でも後年の再記録を予測する。** したがって Kashmir RAIL v26.1–v26.3 の「種同一性を越える初期量的順位の追加予測力」は、鉄道だけの固有の生態選別過程であるとは言えない。道路でraw rankの差は内陸優位で、train-centered rankの差は道路際優位。目的変数は**観測された再記録**であって個体死・繁殖・真の定着は識別していない。3モデルの差の範囲は小さく、事後的に正の設定だけ選んで世界一般則にしない。

## 4. 科学的にはどこが新しく、次に何が必要か

Kashmir RAILの同一地点区画4群での初期量的順位の再記録予測（v26.1–26.3）と、MIREN ROADの7地域20道路での独立source留保（本v26.5/26.6）を並べると、**「初期量的優占度依存の再記録は人工回廊一般でも見られる」**という非鉄道ベースラインが得られた。ただしKashmirは1原著に過ぎず、MIRENは25地域を記述する別研究、量的元階級、観測年月差、source protocolが異なるため原logloss差を同一計量のmeta効果にしてはならない。MIREN元研究での roads×regions は独立の国・鉄道試験ではない。

Kashmir元原著自体が稀少な先駆植物から優占種への遷移を既報している。本研究の追加貢献は「種名を調整した初期量的順位の地点外予測」という比較可能な検証枠組みと、その**非鉄道の標準化されたフィールド再調査での確認・仕様依存性の提示**にある。

次に必要なデータ/検証：
1. 鉄道で独立した複数地域の同一地点×複数年×種別量的 abundance source、同時期の隣接非鉄道対照を探す。ReSurveyGermanyは現行raw不明で物理railway labelingも不明。
2. Kash RAILとMIREN ROADSの比較は**同一分類群の検出確率、調査季節、年差、source変数の共通定義**を評価しない限り、鉄道の効果強度差としない。
3. 「列車による運搬」「刈取り・風への適応」「個体の真の生存率」を主張するなら、輸送起点終点・種子検出・実際の剪定/突風・個体識別が必要。

再現：scripts/audit_MIREN_original_road_resurvey_initial_abundance_auc_v265.py、tests/test_MIREN_initial_abundance_road_resurvey_auc_v265.py；scripts/analyze_MIREN_road_abundance_after_taxon_out_of_road_v266.py、tests/test_MIREN_road_abundance_taxon_roadholdout_v266.py、対応するGitHub Actions workflows。ソース再配布でなく正規著者原データ／MD5 source artifactsから同じ結果を実計算する。MIREN公開データを投稿論文に使用する前には MIREN 提供者との連絡・研究計画共有を推奨する元ガイダンスを尊重する。
