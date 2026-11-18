# tetsu v25.8–v26.0 — 鉄道植物の局所供給・初期優占度・時間的再記録：原データ検証結果

最終照合 2026-11-18 JST。研究ブランチ research/tetsu-local-source-2014-2017-v258-20261118 に保存し、main は変更しない。**同じ Kashmir 原著の独立した再解析3仕様**であり、3つの新しい野外研究・3つの独立国の検定ではない。元種×地点の生データ、実所在地・種別 abundance 数値はGitHubへ転載しない。

## 1. 厳格な元原著資料

Rashid et al. (2021), Journal of Applied Ecology, DOI https://doi.org/10.1111/1365-2664.13961、original Dryad https://doi.org/10.5061/dryad.0p2ngf20h。原著 version 99735の Matrix_plot.xlsx 原本122,187 bytes、SHA256 529f3b97521933018e12623b61f531f012c67241cae2602a6675dbd534192474 を毎回GitHub Actionsで再取得・照合。31真正地点（16駅・15駅間）、2014/2017両年、100m²の4ゾーン、123元分類群、真の欠損セル4つ（0に補完しない）。全4ゾーン×両年に完全な同じ現地地点は駅15・駅間12。

著者の本文では種別 abundance を各plot内の12個の1×1 m quadratから記録したとある。**元Excel数値を一律「被度%」「個体数」だと断定しない。** 今回は同一site×year×plotで元数値の順位に限って使用する。

- [v25.8 SOURCE SUCCESS #44217504159](https://github.com/zuizui0223/tetsu/actions/runs/44217504159)（合成QA 7/7）。[集約元結果ZIP artifact #14389535828](https://github.com/zuizui0223/tetsu/actions/runs/44217504159/artifacts/14389535828)。
- [v25.9 SOURCE SUCCESS #44218411342](https://github.com/zuizui0223/tetsu/actions/runs/44218411342)（QA 5/5）。[同artifact #14389855269](https://github.com/zuizui0223/tetsu/actions/runs/44218411342/artifacts/14389855269)。
- [v26.0 SOURCE SUCCESS #44218855802](https://github.com/zuizui0223/tetsu/actions/runs/44218855802)（QA 4/4）。[同artifact #14389237770](https://github.com/zuizui0223/tetsu/actions/runs/44218855802/artifacts/14389237770)。

## 2. v25.8 — 2017年新規記録種は、2014年の同地点の内陸植物相に予め記録されていたか？

各source地点で2017 Trackには存在し2014 Trackには記録されていなかった種を「新規記録」（真の到着・移入とは区別）とする。2014 Interior1+Interior2に存在した割合を、同じ駅/駅間群の他の地点の2014植物相による期待割合と比較。真の地点の数と種リストを完全固定して999回site identity derangement。真のサイト別期待は厳密に全他地点平均で算出。原zone自身の2014分類群を潜在供給元には入れない。

| 元source層 | 元地点 | 2017 Track新規taxon×site | 2014局所内陸に記録済み（地点等重み） | 他地点内陸の期待 | 局所追加予測 | 技術上尾順位 |
|---|---:|---:|---:|---:|---:|---:|
| Station | 15 | 24 | 0.1410 | 0.1355 | **+0.0055** | 0.558 |
| Between | 12 | 17 | 0.2417 | 0.1568 | **+0.0848** | 0.177 |

2014 Embankmentだけをsourceとすると局所追加予測は駅 −0.0311、駅間 −0.0318。Embankment+Interior1+Interior2なら駅 −0.0110、駅間 +0.0045。両内陸を固定したEmbankment新規記録の局所予測差は駅 −0.0245、駅間 +0.0052。**局所source優位は新規記録に頑健ではない。**

これに対し「2014 Trackに記録され2017 Trackには記録されなくなったtaxa」に対する2014局所内陸の情報増分は駅 +0.1086、駅間 +0.0737。Station 15地点では2014 Trackに116 taxon×site記録、2017まで継続22、消失記録94、新規記録24。Between 12地点では2014 Trackに54件、継続19、消失記録35、新規記録17。「新規記録」は検出漏れも含み、列車からの真の種子加入ではない。

## 3. v25.9 — 初期の量的出現順位と2017年再記録

2014の各元plotで記録されたtaxaのみ母集団とする。同じsource plot内の2014量的数値（順位）で、2017にも記録された種が、2017未記録の種より大きな数値をもった確率（同順位は0.5）をsite別AUCとして計算。2014 site内のtaxonセット、再記録/未記録の件数を厳密固定した999 label permutation帰無を確認。0.5はその期待値。各siteを1票として等重み。

| 元sourceゾーン | Station（AUC; eligible sites） | Between（AUC; eligible sites） |
|---|---:|---:|
| Track | **0.7659; 14** | **0.8535; 12** |
| Embankment | 0.6989; 15 | 0.7094; 12 |
| Interior1 | 0.5636; 15 | 0.7036; 12 |
| Interior2 | 0.5553; 12 | 0.6412; 12 |

Track sourceのlabel permutation上尾順位は両層とも0.001。Trackの元地点1つ除外AUC範囲は Station 0.7479–0.8008、Between 0.8402–0.8856。これは**量的な初期優占度と後日の種記録の条件付き予測**であり、野外で同じ個体が生き残ったと追跡した結果・遺伝的選択の確認ではない。

重要な先行研究上の留保：Rashid et al. 2021 は元論文中ですでに「初期に希少だった先駆植物が失われ、優占種に置換されてゆく遷移の可能性」を示した。**新規性は初期種別数値の再記録AUCと、同一地点・異環境の比較による追加定量化であり、遷移そのものを本研究が初発見したわけではない。**

## 4. v26.0 — 完全に同じ現地地点で Track の「初期優占度による再記録予測」は内陸より強いか

両方のplotで2014数値のある種の「後年再記録あり／なし」が各最低1分類群ある真の同一地点だけを保持し、site内 AUC_Track − AUC_reference を計算。1999回原siteをクラスタ再標本化した技術分位幅で、field treatment 95%信頼区間ではない。比較ゾーンのtaxon母集団、2014種数、優占度分布は違うため、実列車運行由来の追加選択と直接解釈しない。

| 元source地点群 | 実plot比較 | 両ゾーンでAUC評価可能な同一地点 | AUC差 Track−reference | original-site resampling 2.5–97.5% 技術範囲 |
|---|---|---:|---:|---|
| Station | Embankment | 14 | **+0.0684** | −0.0664〜+0.2024 |
| Station | Interior1 | 14 | **+0.1977** | **+0.0270〜+0.3448** |
| Station | Interior2 | 11 | **+0.1872** | −0.0009〜+0.3593 |
| Between | Embankment | 12 | **+0.1440** | **+0.0335〜+0.2451** |
| Between | Interior1 | 12 | **+0.1498** | **+0.0339〜+0.2670** |
| Between | Interior2 | 12 | **+0.2123** | **+0.0743〜+0.3351** |

Trackの再記録予測は現地内陸に比べ高い傾向がある。駅のTrack−EmbankmentとTrack−Interior2は技術分位幅がゼロをまたぐ（後者はごくわずか）。Station vs Betweenは独立国の野外反復ではなく同じ調査の空間層なので、全球鉄道差と混同しない。

## 5. 本来の統合科学的位置付け

これまでtetsuでは、スロベニア道床肩の一年草率上昇、駅別植物の欧州間出現頻度順位の部分的な再現、真の鉄道微生息地間の大きな分類群置換、強い同一地点の地理的結合を観測している。本Kashmirの3回の原データ監査で新たに示したのは、**2017の新規種の局所植物相による追加予測は弱いのに、2014種別量的記録は2017の再記録を強く予測する**という区別である。

現段階の候補機構は「地域種プールによる加入可能性 × site別フィルタ × source初期優占度/記録可能性」の相互作用。元Kashmir原著で既に議論された succession とも整合的。ただし、列車本数・貨物OD、除草・強風・実際の個体寿命、遺伝変化の直接測定はない。検出漏れ、植生季節、前回記録数の違い、taxon identityによる効果を分離しない限り因果的な鉄道適応を主張しない。

次の段階は、①新しい独立鉄道研究の同年・同面積の zone×taxon×abundance再調査を取得、②**taxon identity固定後に初期量的 abundanceが leave-site-out で追加の予測力を持つか**、③他の攪乱立地と同方法比較、④実運行数や管理履歴を導入、という手順。MIRENの道路世界データは別機構の対照にはなるが、鉄道の独立実地サンプルとして扱わない。

## コードと検証

- v25.8: scripts/analyze_Kashmir_2014_local_flora_predict_2017_rail_gains_v258.py、tests/test_Kashmir_local_flora_gain_v258.py と対応 Actions workflow。
- v25.9: scripts/audit_Kashmir_2014_species_cover_persistence_2017_v259.py、tests/test_Kashmir_early_cover_persistence_v259.py と対応 Actions workflow。
- v26.0: scripts/audit_Kashmir_same_sites_early_cover_AUC_track_vs_adjacent_v260.py、tests/test_Kashmir_matched_zone_auc_v260.py と対応 Actions workflow。

原著 Dryad 元データ利用規約を順守し、出典・バイトdigest・各元単位を明記して、原種×観測の全文やGPSは再配布しない。
