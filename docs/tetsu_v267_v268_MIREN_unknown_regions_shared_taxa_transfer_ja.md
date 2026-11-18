# 鉄 tetsu v26.7–v26.8：全球道路の「初期優占度による再記録予測」は、未知地域と同一種の両生息地対比に耐えるか

最終元データ監査日：2026-11-18 JST。独立した**MIREN山岳ROAD**研究の source-authentic 再解析。鉄道のもう一つのサンプルと数えない。Zenodo https://doi.org/10.5281/zenodo.15749536 の原CSV v3（17,899,241 bytes、167,804行、MD5 b640195cc44f252ccb6ddb28eb6fff13）を、source-verification済みGitHub Actions #43212464976 artifactから取得して原バイトを照合。原 source 地理・個別種・species×plotデータはGitHubへ再配布していない。元MIREN consortiumと2026年データ論文 DOI https://doi.org/10.1111/geb.70284 に帰属する。

**v26.7 SOURCE VERIFIED:** [Actions #44244622350](https://github.com/zuizui0223/tetsu/actions/runs/44244622350), [derived artifact #14401330553](https://github.com/zuizui0223/tetsu/actions/runs/44244622350/artifacts/14401330553)、synthetic数学・source-foldテスト5/5成功。

**v26.8 SOURCE VERIFIED:** [Actions #44245397243](https://github.com/zuizui0223/tetsu/actions/runs/44245397243), [derived artifact #14401023100](https://github.com/zuizui0223/tetsu/actions/runs/44245397243/artifacts/14401023100)、synthetic same-taxon / site-matching / nested-weightテスト5/5成功。

## 1. 一次資料とモデリングの真の支持域

前v26.5–v26.6と全く同じ一次母集団を固定する。MIREN 484真正反復道路地点のうち、初回2012年以降・source Plot1道路際とPlot3内陸の両方が調査・後年再訪され、初回Source Abundance階級1–3が利用可能で、地域内2道路以上・各道路4地点以上という事前supportを通過した **7 Regionコード／20本の実道路／251実地地点**。初回にsource種が記録された正の原risk eventsは **14,682 taxa×plot**、元文字列speciesは1,925学名。100m²両plotの年次初回/最終回のSource Protocolのversionは全484再訪地点で変化しており、検出変化・種同定基準・調査員変化を因果的に除外していない。

予測するもの：初回2012年以降にその元site×plotで実際に記録された植物種が、最終回の同じsite×plotで**再記録されるか**。個体を標識した生残・真の局所消滅・植物の遺伝適応ではない。

L2 logistic C=1を固定。v26.6では同一Region内の別Roadを学習に残す leave-entire-ROAD-outであったが、v26.7からは **leave-one-entire-Region-code-out**。留保地域の全元Road・Site・初回の翌年結果・種記録を学習から除外。未学習の新地域ではRegion one-hotが学習不能なので、他の地域を含めた共通の plot/initial_year/elapsed_years/species のメタデータだけでモデルを構築する（元履歴のtargetは不要）。全モデルの学習と照合は地域単位で分離。元site内source taxonイベントを平均→Road平均→Regionコード平均という3段階集約を維持し、14,682件を独立実験と偽らない。

- M2→M3：plot・初回年・再訪年差・species identityから、初回plot内の源Abundance中順位 [-1,+1] を追加。
- M2Z→M3Z：species×habitat Plot1/Plot3別の履歴を学習した上で順位を追加。
- M2Z→M4Z：種×habitatの通常の初回順位を、**他地域のTRAINデータのみ**から計算して差し引いた中心化順位を追加する。留保Regionの中心化平均や将来結果を学習に使用しない。元税on×Plotの訓練記録なしの場合の中心化順位は0として、未知種の効果を捏造しない。

各モデルの差は、値が正なら新変数を追加すると元未知Region source地点のlog-lossが改善することを意味する。Brierは独立の予測誤差指標。サイト・Road・Regionの選ばれ方と2012以降のprotocolはランダムでないので母集団CIにしない。

## 2. v26.7 — 真の地域外予測で残るもの

| 7 Regionコード等重みのsource留保log-loss利得 | Plot1道路際 | Plot3内陸 | 道路際−内陸 | この差が正のRegion数 |
|---|---:|---:|---:|---:|
| M2→M3 species + initial rank | +0.018182 | **+0.031805** | **−0.013623** | 2/7 |
| M2Z→M3Z species×habitat + rank | +0.018594 | **+0.030589** | **−0.011995** | 2/7 |
| M2Z→M4Z training-only species×habitat centered rank | **+0.011957** | +0.005119 | +0.006837 | 4/7 |

BrierでもM2→M3差は道路際+0.009193／内陸+0.015623、M2Z→M3Zは+0.009148／+0.014443で内陸側優位。中心化M4Zでは道路際+0.005673／内陸+0.002364。**どちらの生息地にも一般的な初期多寡による再記録予測力は転移する**が、どちらが強いかはもとのtaxon平均による中心化と種のpoolの範囲に依存する。

### 予測先地域の種共有率そのものが大きな制約

1つのRegionコードを留保したとき、14,682 taxa×plot risk observationsのうち、他Regionで同一種に少なくとも1つ2012+ source訓練例があったのは**4,311件（29.4%）**、同じ学名が他Regionのtrainに存在しないのは**10,371件（70.6%）**。元7Regionのそれぞれにおける既知種risk event割合は17%〜43%程度。Regionコードを国や論文と数えない。この違いによって、種同一性モデルの適用可能な部分が地域で大きく異なる。

種で層別したときの原順位 M2→M3の改善（Region等重み）：
- **他Regionのtrainにも存在する種**：Plot1 +0.030620、Plot3 +0.017977（道路際が強い方向）。
- **他Regionで未学習の種**：Plot1 +0.010578、Plot3 +0.035429（内陸が強い方向）。
- 未学習種の種固有係数は訓練データでは推定不能であり、モデルが取得できるのは生息地/年次メタデータと源量的順位の関係。「新種の鉄道輸送を予測した」とは全く別。
- M2Z→M4Zの中心化項は未知taxon×Plotで構造上0に設定されるので、未学習種で近傍0の利得は「種内優占度差が存在しない」の直接測定ではなく**定義上の限界**である。

## 3. v26.8 — 同一元Site・同一2012+初回年・同じtaxonがRoadsideとInteriorの両方にいたケース

v26.7のroadsideとinteriorの比較では、元の地域内でも出現種の種類が異なることが交絡。より対称な比較として、MIREN真正2012年以降初調査sourceの同一Region/Road/Siteにおいて、**同じ学名taxonが元Plot1と元Plot3の双方に同じ元年で出現したときだけ**テスト損失を比較する。学習モデルはv26.7と厳格に同じ14,682 risk events、7 Region留保、同じ5モデル。テスト対象だけを固定された共通植物種へ絞る。将来の再記録はPlot1とPlot3で別々のsource観測であり、同じ個体の処理前後ではない。

- 真の251/251地点に少なくとも一つ両zone共通taxonがあり、総**3,280 matched-original-site×taxon pairs（Plot1+Plot3で6,560 taxon×plot events）**。どの実Siteでも元存在共有を観測した。本研究の7Region/20Road支持を全て維持。min2 sites/roadとmin4 sites/roadの両方で原結果が不変だった。
- あるsiteの同じtaxonの両zoneで初回順位が違う場合でもtaxon ID/原region/source siteと初回年を固定して後年再記録のprediction gainを比較。真の元植物相の組成差を大幅に縮小した反証用サブセットだが、両zoneにいるspeciesだけに選択するため全種群集のランダム標本ではない。

| 同一taxon・同一Site両zone、元7Region等重みLOHO log-loss利得 | 道路際 Plot1 | 内陸 Plot3 | 道路際−内陸 | 差が正のRegion数 |
|---|---:|---:|---:|---:|
| M2→M3 initial raw ordinal rank | +0.006471 | **+0.035891** | **−0.029421** | 1/7 |
| M2Z→M3Z taxon×habitat + initial raw rank | +0.007036 | **+0.033858** | **−0.026821** | 1/7 |
| M2Z→M4Z taxon×habitat + TRAIN-centered rank | **+0.011021** | +0.006553 | +0.004468 | 3/7 |

上位2つのraw rank利得は6/7地域で内陸側が強く、1 Regionずつ除いても平均差の負が残る（M2Z→M3Z −0.03273〜−0.02069）。中心化の結果はRegionを1個除くと −0.00063〜+0.00637と符号が安定しない。**同じ植物種・同じ元Siteへ統一しても「roadsideほど初期多さが重要」とは言えず、raw rankの相対予測力はむしろ内陸側で強い傾向。** ただし中心化の絶対利得は正でも、地域でその相対強度が反転するため、都合の良い仕様だけで一様な攪乱勾配の一般則を主張しない。

## 4. 鉄道植物マクロ論文への結論と停止条件

Kashmir RAIL source-authentic単一の2014/2017原著（v26.1–v26.3）では Trackの初期量的記録が後年再記録を、種同一性・taxon×微環境を超えて少し改善した（地域間訓練ではない地点留保）。MIREN ROADでは、それに対応する種同一性を越える初期量的記録の追加予測力が**未知のRegion**にも移る。ただし:
1. 初期量的順位による後年再記録の予測は**鉄道以外の一般的な道路回廊でも成立する**。鉄道固有の列車風圧・草刈りへの適応としては解釈できない。
2. 未知地域への予測で最も重要な交絡は**学習済みの植物種がそもそも別地域に少ない**こと。species pool共有は予測の受容域（transportability）を変える。
3. Siteとtaxonの両方を固定した対比では、raw rankの道路際優位がなく内陸側の利得が強い。一方でtaxon-normalized rankは地域によって不安定。**平均形質による環境フィルターと種同一性・順位による群集編成を同じメカニズムに安易に還元できない**。
4. 元MIREN観測のfirst/lastの全再調査では調査Protocol版が変わっている。検出過程、季節、原種数、調査年差、地域種poolの混合、同じ種でも異なる個体と生育条件は未同定。真の個体生存・列車種子運搬・刈取りや風への進化的選択を証明しない。
5. これらは**1つのグローバルROADデータセット内の7Region holdouts**。独立した7原著を再現したと数えない。鉄道の独立再調査は今もKashmir 1研究だけであり、新たな原著鉄道×対照同調査法データが必要。

再現コード・synthetic QA・原ソース認証CI：
- scripts/analyze_MIREN_source_region_holdout_abundance_v267.py；tests/test_MIREN_region_holdout_abundance_v267.py；.github/workflows/tetsu-v267-MIREN-whole-region-holdout.yml
- scripts/analyze_MIREN_source_taxon_matched_road_interior_region_holdout_v268.py；tests/test_MIREN_same_taxon_both_plot_region_holdout_v268.py；.github/workflows/tetsu-v268-MIREN-taxonomic-matched-region-holdout.yml

探索・感度検定であり、事前登録の確証的な比較研究ではない。原著と配布者の利用・引用条件を順守し、MIRENデータを主要投稿材料に用いる場合には MIREN consortium の事前研究計画共有・提供者関与のガイダンスを尊重する。
