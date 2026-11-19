# tetsu v29.0 — 離れた鉄道群集で生活史選別の方向が再現するか：区画参照と種同一性の反証

更新: 2026-11-19 JST。`research/tetsu-kashmir-gradient-taxon-influence-v290-20261119` 上の追加 SOURCE-based 解析。原著SHA・形質SHAを検証した GitHub Actions 実行に先立ち、数値を事後的に選んで普遍性の支持に転用しない。main は変更しない。

## 既存出発点
- Slovenia: Šabić et al. (2026) 89実地点・173同側BS/BTペア、strict annual fraction BS−BT +0.054703、15路線区間 whole-cluster technical bootstrap範囲 [+0.017619,+0.086547]。BS=バラスト肩、BT=法尻。Source DOI https://doi.org/10.1007/s44473-026-00155-z 。
- India Kashmir: Rashid et al. (2021) Dryad DOI https://doi.org/10.5061/dryad.0p2ngf20h、元31本物地点、4区画×2年完全共通27地点。Nature 2023（https://doi.org/10.1038/s41586-023-06644-x）元23.6万種life cycle dataとKew厳格2種の補完後、元 `Track−Interior1` の一年草率2年平均 sharp bounds [−0.202726,−0.053184]、元27地点cluster bootstrap最悪上限のtechnical q025/q975 [−0.124482,+0.022116]。Source run https://github.com/zuizui0223/tetsu/actions/runs/44298354325。標本内の形質欠測への頑健性であってインド全体のCIではない。
- Kew POWO 原著species：`Cynodon dactylon` perennial/rhizomatous（https://powo.science.kew.org/taxon/urn%3Alsid%3Aipni.org%3Anames%3A397065-1/general-information）と `Breea arvensis` が `Cirsium arvense` の homotypic synonym（https://powo.science.kew.org/taxon/urn%3Alsid%3Aipni.org%3Anames%3A186181-1）；受容名はperennial（https://powo.science.kew.org/taxon/urn%3Alsid%3Aipni.org%3Anames%3A195034-1/general-information）。

## 今回の5つの元source microhabitat比較と反証
1. Track−Interior Plot 1：原著同一サイト近側の本来の内陸植生参照。
2. Track−Interior Plot 2：より遠い本来の内陸植生参照。
3. Track−Embankment：どちらも元鉄道影響区内。rail-vs-nonrail対照と呼ばない。
4. Embankment−Interior 1。
5. Embankment−Interior 2。

すべて同じ27の元 site を2014/2017の2年にわたり保持。Station n=15 と Between n=12を別々に記述し、train-frequencyを測ったかのように改称しない。各比較は「全元観測植物のうち厳密な一年草の存在種割合」のsite-equal mean。LifeSpan未分類17 source種を同一original speciesなら全区画・地点・年で一貫して annual=0/1 と扱い、線形 partial ID の sharp lower/upper を算出。未知の二年草はstrict annualでは0。元species incidenceを再配布しない。

### Primary falsification gates
- **参照の頑健性**：Track−Interior1だけが負なら、一般的stress gradientではなく参照配置依存として解釈。Track−Interior2/Embankmentと符号が違った場合は明示。
- **場所の頑健性**：元27 physical siteをひとつずつ丸ごとLOSO。2014/17の同じsiteを両年留保し、分母の元植物種数を保持。各groupは別にtrue nで集計。全27LOSOが負でも未調査国の普遍的予測力は示さない。
- **種への依存**：既定のoriginal species Cynodon dactylon / Breea arvensis / Taraxacum officinale / Lactuca serriola / Senecio chrysanthemoides / Lotus corniculatus を元source全plotから1種ずつ除去して fractionの分母も真正に変更。この sensitivity は「確認済みperennial labelを unknown に戻す」という別の実験と混同しない。除去によって元zoneが空になったsiteがあれば同地点比較を減少させ、reportで記す。
- **不確実性と効果の区別**：Cynodonのsource perennial labelが無ければnegative sharp upperを確定できない可能性がある。その場合は「生物学的結果が一種の確認済み生活史ラベルに条件づけられている」と誠実に報告。原種の出現・器官測定が因果的な除草・耐風適応を識別することはない。
- **外的検証**：これを Euro4著者原著の trait-label shuffle null と同じ平均量へpoolしない。Slovenia BS/BT、France track/embankment、Kashmir Track/Interiorは真の物理勾配が異なる。大陸で学習した種の life form 方向予測は未評価。

## 他の大陸への募集候補
- France Rakotoson et al. (2025): 128 track /266 embankment、著者原著はannual/ruderal track enrichment、元完全行列は request only。
- Perth/Beijing Liu et al. (2026): DOI https://doi.org/10.3390/land15050714、3 sites per city ×5 quadrats 1m²、年と季節の違い、active/abandoned/redesigned typology、published annual count but plot raw unavailable（Data Availability: not publicly available）。1研究、2都市を別publicationとして数えない。
- China–Laos Li et al. (2026): DOI https://doi.org/10.1016/j.jenvman.2026.130380、14駅84plot、植物・土壌・菌類対照があり機構探索に有望だが元source個別taxon×plot未取得。
- Brazil Mendes et al. (2026): DOI https://doi.org/10.1590/2175-7860202677029、railway４区間400 herb quadrats、rail外真正対照なし。

追加原著の選抜台帳: `data/literature/tetsu_v290_independent_rail_annual_gradient_evidence_registry.csv`。同じ publication の複数都市/季節は1原著。真の meta effect のPIを出すには、同じ推定対象で効果・field SEが複数独立著者原著から必要。予測区間ゼロ非交差だけで「全地点で普遍」とは書かない。

## Reproducibility
- 原source実行スクリプト `scripts/analyze_Kashmir_gradient_taxon_influence_v290.py`
- 合成数学/元source spatial groupingテスト `tests/test_Kashmir_gradient_taxon_influence_v290.py`
- CI `.github/workflows/tetsu-v290-Kashmir-original-gradient-taxon-influence.yml`
- GitHub workflowは Dryad `Matrix_plot.xlsx` SHA256 529f3b97521933018e12623b61f531f012c67241cae2602a6675dbd534192474 を確認し、Figshare Nature global life-cycle CSV SHA256 83b23e7874d73387a0c24aeae3e2931c50ec7c78cad52ff26cdfc99d1b6623b1 を読み、Kew証拠URLと合成テストを使う。元source×site×species raw table・GPSを artifactに含めない。
