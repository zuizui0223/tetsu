# v29.6 — 4独立欧州鉄道原著と枝長付き植物系統樹の条件付き帰無

2026-11-19 JST、元植物の本物のsite×species観察を固定したまま、真の系統樹の内部枝に対応する単系統群を使う。新しい実測値はsource GitHub Actions完走後にだけ記載する。

## これまでの結果と独立性

元4原著：CZ39駅全域、IT28駅全域、HR+BA15駅（2国で1原著）、SI89個100m道床地点。野外物理単位合計171は独立原著171件ではない。v28.5の広い種頻度5分位帰無では草丈CWM β超過が4/4研究で正、等重み+0.06790。v29.4科内traitシャッフルは4/4正、+0.03434。v29.3属内traitシャッフルは3/4正、+0.00628。3つは異なる帰無空間であり差は系統が説明した分散率ではない。

v29.5 original-source CI SUCCESS #44517543086 では真の3元 mega-tree LCVP/WP/TPL の既存tipへ完全一致できたLEDA両形質既知taxaが元原著種の82.1–91.5%であり、12/12 original study×treeが事前に定めた対応率support gateを通過した。ただし元mega-tree自身が著者により推定・増補された系統樹であり、tip掲載は当該種のDNA直接配列による検証ではない。今回さらに欠落した種を接ぎ木しない。

### 重要なSOURCE QA修正

初回v29.5派生表の「地点ごとのtip割合の中央値」は、Rのlengths(data.frame)で元サイトの種数でなく列数6を分母とする不備があった。一部に20など1を超える値が現れたため無効。対応種数・5種以上の元地点の数・全12組のsupport判定・元CWM observedには影響しない。v29.5 Rコードは元サイトの観察種数をnrow()で数える方式へ修正し、出力率が0–1に収まらなければ停止する。本v29.6ではこの修正済み source audit を同一 CI の中で再計算し、先行QAの分母が不正なら解析を行わない。

## 今回の枝長付きsource original

Jin & Qian (2022) V.PhyloMaker2 DOI https://doi.org/10.1016/j.pld.2022.05.005。著者GitHub公開のoriginal mega-tree LCVP、WP、TPLの原R-data（固定コミット7af3fb5152f691af2e4ec9d5e2e467d1b50505e9）を真正Git blob SHAで照合。元published phyloのtip.label、edge、edge.lengthを読み込み、rootから各tipまでの実枝長を計算する。長さのmin/median/max/CVを公表して、元枝長の根からの深さが揃っていない可能性も記録する。

系統樹内部枝を root-to-tip枝長中央値の **0.25 / 0.50 / 0.75** の深さで切り、その位置にある既存の単系統cladeの中だけで、本物のLEDA草丈と種子重の **paired 2D trait vector** を交換。これは真のbranch lengthを用いた内部node群への割付けで、科名・属名だけのgroup nullではない。相対深度を勝手に100万年前などの年代へ換算しない。Brownian trait conditional simulationやPGLSに相当すると誤称しない。

### 同じ標本で比較する4つの帰無モデル

1. source各treeで直接一致した元taxaのみを使い、元site×taxon出現率5分位内でtrait vectorラベルを交換する元source基準。
2. 実系統樹の相対深度0.25に対応するmonophyletic clade内のみ交換。
3. 相対深度0.50のmonophyletic clade内のみ交換。
4. 相対深度0.75のmonophyletic clade内のみ交換。

各研究・treeについて対象の元site・taxa・trait標準化軸は完全同一、先行v29.5の元study×tree直接tip種数と元CWM height observed 距離に一致しなければ新帰無を走らせない。今回の基準系列は元tip全種ではなく直接tipが一致したsubsetなので、v28.5のfull source betaと混ぜてメタ平均しない。

## 帰無の識別可能性ゲート

対象の真tree×元独立原著×相対深度の各比較について、実sourceの観測に基づいて事前に以下を固定する。

- 元観測植物の2種以上を含む交換可能な真系統cladeが5群以上。
- その交換可能なcladeに属する元site×taxon incidenceが全体の25%以上。
- 元地理区分の2群以上で交換可能trait-known source incidenceが各25%以上。

条件未達は NOT_IDENTIFIABLE_PHYLOGENETIC_CLADE_EXCHANGEABILITY、負効果や真の0として記録しない。十分なoriginal独立field出版（3本以上）がsupportを満たす時だけ等元出版重みの記述平均と1研究除外範囲を計算する。2固定seed×199=398種ラベル技術配置の順位は野外独立検定のp値ではない。3本の系統樹は同じ野外sourceに対する命名法感度解析で、12独立研究ではない。

アウトカムは元研究地理区分内の異なるsource siteの草丈CWM絶対距離の観測とtrait-label nullとの差。種子重・2D jointも併記し、元site×taxon組成、species occupancy、元field geographyを固定する。使用した系統樹自身に人工的なバックボーン補完が含まれうるため、種ごとの元DNA配列の直接根拠と区別する。

## 実装と原資料の秘匿

先行 v29.5 original-source Python helperで原author site×plant taxon×LEDA sourceをGitHub Actionsの private RUNNER_TEMP のみに復元し、元本当のplant list、現地GPS、各元site×taxon原行列、LEDA元の植物種別値、R tree原バイナリを公開 artifact・Gitに追加しない。source QAとして出版別×tree別×conditional nullの集計48行だけを公開する。

対象：
- scripts/analyze_original_4rail_phylogeny_clade_depth_trait_label_null_v296.R
- scripts/validate_and_aggregate_four_rail_true_tree_clade_source_v296.py
- tests/test_Europe_four_rail_true_tree_clade_source_v296.py
- .github/workflows/tetsu-v296-real-four-rail-phylo-clade-depth-null.yml

結論の限界：これは4欧州独立元出版のtraitと植物群集地理配置の関係であり、鉄道適応、列車による種子散布、草刈りや車両風圧という因果的メカニズムの実証ではない。稼働中rail ballastと同一面積/年/開放度の非鉄道砂礫地点原sourceは別途必要。
