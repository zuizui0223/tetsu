# tetsu v30.1 事前登録：山岳道路を丸ごと伏せても、標高―群集草丈の関係は移せるか

2026-11-20 JST。**現段階：解析コード・事前登録。原著ソースの新規 v30.1 結果は CI 成功まで未確定。** この計画は MIREN の山岳**道路**データ1研究プログラムであり、鉄道ではない。

## 背景と新規性

真正原著 v29.9 は同一 Road 内で実測標高が高いほど、群集構成種の LEDA 種中央値に基づく CWM 草丈が低いことを示した（実系統樹75%単系統クレード×植物出現頻度5分位の同時制約下）。真正 v30.0 は同じ Road 内の未使用 Site を標高から予測すると、訓練 Site の単純平均に比べる RMSE 利得が Roadside +0.035156 (LCVP; 系統×頻度帰無超過+0.028757) と示した。しかし**実観測利得が正なのは元5地域のうち3地域**。Roadを丸ごと未使用にした予測は未検証。

次の反証可能な問い：**別の本物の Road から習得した群集草丈―標高の傾きは、1 Site だけ草丈を測った新規 Road の残りの Site へ移せるか。**

これは**1地点較正付きのRoad外推**であり、「まったくデータのない未知Road」の cold start ではない。

## 事前固定するデザイン

1. v30.0 と同じ真正 MIREN 2007–2023 原著 CSV (167,804行)、元Site×Year の Roadside Plot1/Inland Plot3 対、実測 Elevation、LEDA 原著 canopy height + seed mass、V.PhyloMaker2 LCVP/WP/TPL 三原樹を使用。元原資料の MD5/SHA と v29.9 の真正6集計を再検証し、さらに**公開済み v30.0 の系統樹×primary/strict計6集計を1e−10以内で再現**できない場合は新結果を報告しない。
2. 生息地、種×Site、各種占有頻度、元標高、現地調査年、原著道路を一切変更しない。採用道路は元v30.0の先決条件（各Road ≥6個の実地点、実標高幅≥150mなど）を**そのまま再使用**。各Regionで少なくとも2本の実Road、記述的地域等重みの統合には最低3 Region。満たさない組は効果ゼロではなく NOT_IDENTIFIABLE。
3. **Leave-one-original-Road-out**：対象の1 Road の全Siteを訓練から外し、ほかの各 Road 内で CWM z(log canopy height) に対する標高(km)のOLS傾きを推定、訓練元Roadを等重み平均する。傾きの学習には伏せたRoadのCWMを一切使用しない。
4. 新Roadには現実的な**1 Site較正**を許す。Siteの草丈を見ずに標高最低のSiteを決め、そのSiteの実CWMだけを基準点にする。残りの実Siteで「較正点の高さ＋他Road学習傾き×標高差」のRMSEを、「較正点高さの一定値予測」のRMSEと比較。主要指標は **ΔRMSE＝RMSE(較正点一定値)−RMSE(他Road傾き外挿)**。正なら未使用Siteで実予測改善。各Road→各Regionの順に等重み集計。
5. **事前固定の感度解析**：較正Siteを中央値標高に最も近いSiteへ変更。元v30.0と同じ単一調査年だけのRoad集合。Roadside・Inlandで別々に検定し、同一採用Road集合でRoadside−Inland利得差も報告。
6. **条件付き帰無**：本物の同一系統樹75% depth単系統クレード×species occupancy quintile交差群内で、LEDAの二形質ベクトルを対のまま置換。事前固定の2乱数列×199＝398技術置換で**同じRoad外学習・較正・RMSE**を再計算。実ΔRMSEと、実−帰無期待の超過を併記する。398置換は独立生態研究や野外p値ではない。
7. 地域ごとの符号、原Region一つ除外感度、元道路数、採用実Site数を必ず開示。LCVP/WP/TPLは同じ植物群集の系統名寄せ感度であり、3研究と数えない。

## 判定と失敗時の解釈（結果を見る前に固定）

- 観測ΔRMSEが複数の原Regionで正、最低標高・中央標高どちらの較正でも維持：**調査済み同一Regionの別Road間で、標高に沿う群集草丈の傾きが移せる可能性**。
- 元v30.0同Road CVのみ正、今回Road外の実ΔRMSEが0以下：**局所的な環境―形質対応はあるがRoad間転移できない**。
- 帰無より良くても実ΔRMSEは0以下：植物の系統と出現頻度だけでは説明できない配置はあるが、**実用上の未使用Road予測は改善しない**。
- Roadside−Inlandの地域符号が混在：**一般的な回廊の緩衝効果**を主張しない。
- いずれの結果でも、鉄道固有のフィルター・列車種子散布・草刈り・列車風圧・遺伝的適応を識別できない。元MIRENは独立研究**1プログラム**であり、地域削除範囲は全球予測区間ではない。

## 実装と元データの権利

スクリプト：[v30.1 Road外検証](../scripts/analyze_MIREN_original_road_out_transfer_v301.py)、[数理・漏洩テスト](../tests/test_MIREN_original_road_out_transfer_v301.py)、[真正source CI](../.github/workflows/tetsu-v301-original-MIREN-LOORoad-calibrated-transfer-source.yml)、[元v30.0凍結集計](../data/derived/tetsu_v300_original_MIREN_heldout_Elevation_predictability_summary.json)。元の種×Site表、個別実Elevation/GPS、種別LEDA中央値、系統樹tip名は公開成果物に再配布しない。

## 鉄道側へ戻るための外部データ候補（未解析）

フランスの Rakotoson et al. (2025), Weed Research, DOI https://doi.org/10.1111/wre.70025 はRERの真の線路内128区画・法面266区画、合計394区画(各2 m²)、255種を測定し、線路内の一年草・攪乱性種と法面の多年草・競争性種という差を報告する。**元個別調査行は公開ダウンロードでなく著者へ依頼制**と明記。原著データを得るまでは第5の再解析研究と数えず、線路内 vs 法面の原著ペア情報も未確認とする。

スロバキア Májeková et al. (2025), Flora DOI https://doi.org/10.1016/j.flora.2025.152780 は貨物列車頻度と外来種**α種数**の正の関連を報告。これは種組成の駅間β均質化や列車による種子の直接運搬と同義ではなく、原著個票の検証なしに meta 効果量へ混ぜない。
