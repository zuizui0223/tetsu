# tetsu v30.0 — 実測標高は未使用地点の草丈CWMを予測できるか

2026-11-20 JST。初回実データCV推定前の固定設計。実生物試験・標高を予測できるという事実は source Actions 成功まで主張しない。

## 既存の真正実証と今回の推定対象

[原著source検証v29.9 #44659650978](https://github.com/zuizui0223/tetsu/actions/runs/44659650978)までSUCCESS。MIREN Mountain ROADS（2007–2023、Zenodo v3元CSV 167,804植物観測行、Plot Elevation原著6,446地点区画）、元の同Site×YearでRoadside Plot1とInland Plot3（100m²）を観測。真V.PhyloMaker2 LCVP/WP/TPL系統樹 root-to-tip元枝長中央値75%で区切った単系統群と各種のsource-occupancy頻度五分位を同時に固定した形質ラベル帰無を使い、同じ元Road内の高さCWM（zlog canopy height）のElevation(km)への傾きが予想以上に負であることを確認（primary Roadside LCVP −0.267 / WP −0.296 / TPL −0.213、Inland −0.328 / −0.351 / −0.365 標高1,000m当たり）。ただしこれは同じsource地点に当てはめた関係である。

v30.0新規の主要推定量：同一実original Roadの未使用Original Siteの群集高さCWMを標高だけから予測し、訓練Siteの平均値しか使わない予測より改善するか。独立の別Region/鉄道への転移予測とは区別する。

## 真の Original-Site Leave-One-Out予測

元original sourceで同じRoadの調査に、source phylogeny tip×LEDA両形質既知の植物を残し、最初の適格調査年、元種出現行列、元ElevationとRoadside・Interior 1:1対を固定。各Roadで一つの真のSiteを除外し、残る全Siteを用いて (A) CWM height ~ intercept + Elevation(km) のOLS直線、(B) CWM height ~ 訓練Site平均のみの切片モデルを作り、除外SiteのCWMを予測。全Siteを一度ずつ除外し、予測残差RMSEを計算。

**primary held-out predictor metric**：ΔRMSE = RMSE_B − RMSE_A（z(log canopy height) の単位）。正のとき実標高線形モデルが未観測Siteの誤差を減少させたと判断。負なら標高を使う予測の方が悪い。「observed−ラベル帰無」も別に計算するが、実観測ΔRMSEが正でなければ予測の肯定としない。

一次物理標本単位はRoad。同じRoadのsite対を独立研究としない。Region内は元Road等重み、Region間は元Region等重みの記述平均。元published study programはMIREN一件。

## 原sourceの解析前support gate

- 元v29.9に同じ元author botany/source年/野外調査を再結合し、元traitとメガツリーの源SHAを認証。v29.9 original tree×primary/strict 6組のRegion数・真Road数・真site数・高さ×Elevation観測−phylo×occupancy label-null斜率を **1e−10以内**で再計算する。元sourceが一致しなければ新しいCVを発表しない。
- 元Trait known source cohortは primary 8 Region/22 Roads/272 Site (min5 taxa & 25% flora coverage)、strict 5/13/192 (min10 taxa & 50%)。新CVでは元v29.9のaltitude-eligible Roadから **6 Site以上**の原Road、元Elev幅150 m以上、各留置訓練モデルが特異にならない本物のRoadのみ使用。元Regionには少なくとも2真Road、元v29.9標高解析eligible siteの70%以上を保持。効果符号でSiteを選ばない。
- 年次感度：同じRoadに実採用Yearが1つしかないRoadのみ残し、元Region内2roads以上、元v29.9標高eligible sitesの60%以上を保持。対象が少なければ NOT IDENTIFIABLE、効果量ゼロとはしない。
- 源データMD5: MIREN 17,899,241B b640195cc44f252ccb6ddb28eb6fff13。LEDA height/seed: 既存のSHA pin維持。V.PhyloMaker2元LCVP/WP/TPL: 固定作者commit 7af3fb5152f691af2e4ec9d5e2e467d1b50505e9 と元3 Git blob pin。**新tip graftなし**。元種中央値を現地個体の表現型として扱わない。

## 元系統×植物頻度を同時に固定した帰無

同じoriginal Region内のSource直接tip植物を、元真系統樹の75%切断でできたmonophyletic cladeと、元site出現頻度五分位の交差群に割付け、その交差群内のみで **元草丈＋種子重のpaired2D species trait vectorを置換**。元全地点のplant×site membership、元Site Elevation、Roadside/Inlandペア・survey Yearを一切変えない。各描出は2固定seed×199=398技術ラベル配置で、毎描出同じRoadで新たに除外Siteへの予測利得ΔRMSEを計算し、Observed-NullMeanを記録。Technical permutation rankは野外field p値ではない。Region等重みの集計は3Region以上支持があれば提示し、LOO（一Region除外）と地域ごとの符号を示す。3本のtreeは同じ生物データのtaxonomic感度であり独立3論文ではない。

## 解釈の分岐

実LOOCV利得も元sourceに対するconditional null超過利得も正で、元3樹・厳格source trait cohort・同一年Road感度に維持：標高は同Roadの未観測SiteのCWMに予測情報があり、元系統と種頻度を条件にしても真元植物配置と環境の結び付きが残る。ただし同Roadに限定した内挿的予測で、別道路・別国・別鉄道への転移能力は未知。

元source v29.9の標高傾きは負でも、CV利得が0か負：観察回帰はあっても未使用Siteへの予測機構としては不十分。元clade×frequencyの置換性が不足する場合：源データで未識別、効果ゼロと呼ばない。

RoadsideとInlandのCV利得の差は、同じOriginal Road/Siteに対する観察上の違いであり、非ランダムな管理や土壌/森林開放性・気温・降水は残る。**鉄道固有の草丈適応、列車種子輸送、風・剪定の因果機構を証明しない**。独立したRail4原著に今回のMIREN環境係数を流用しない。

## 実装

- scripts/evaluate_MIREN_native_elevation_heldout_stature_predictability_v300.py
- tests/test_MIREN_true_site_elevation_LOOCV_v300.py
- .github/workflows/tetsu-v300-original-MIREN-LOOsite-Elevation-phylo-occ-source.yml

元site elevation/GPS/種ごとの高さ/植物の現地観測行はActionsの原資料元でのみ使用し、GitHub artifactにはRegion/無名Road数と予測利得の集計のみ公開する。
