# v29.5 — 鉄道4原著の植物に「実際の系統樹の枝長」を対応づけられるか

2026-11-19 JST、元原著×3元mega-tree の真正Source取得と照合率を **結果を見る前に固定**する。新しい系統進化の肯定/否定を source original action 完走前に書かない。

## ここまで明らかになったこと

- v28.5の欧州RAIL独立4原著（元 CZ39駅、IT28駅、HR+BA15駅＝**1出版**、SI真89個100m道床区画）で、元サイト内 plant assemblage 固定の占有率5分位 trait-label null から見た草丈CWMβの非ランダムな地理構造は4/4正。
- v29.3の属内限定shuffleでは平均超過が小さく、3/4正（Slovenia負）。
- v29.4の **科内限定** shuffle は4/4正、ただし属をまたげるので、科内の属の地理的分布と現地環境を切り分けられない。科/属のラベルは**枝長を含む分子系統樹ではない**。3種類の帰無モデルの平均差を直接「科/属が何%説明した」という分散寄与としない。

## 今回の新しい公開真正系統リソース

Jin & Qian (2022) *Plant Diversity*, [V.PhyloMaker2](https://doi.org/10.1016/j.pld.2022.05.005)。著者公開 [jinyizju/V.PhyloMaker2](https://github.com/jinyizju/V.PhyloMaker2) のmasterではなく **固定commit 7af3fb5152f691af2e4ec9d5e2e467d1b50505e9（2022-11-17）** に元から添付の phylo R-data original backbone を取得する。3標準化体系は：
- **LCVP** — GBOTB.extended.LCVP、約73,420 tip
- **WP** — GBOTB.extended.WP、約72,570 tip
- **TPL** — GBOTB.extended.TPL、約74,529 tip

原Git blob SHAの固定：LCVP \`0d761cb0bf36392cc13759a7ac9026db564217eb\`、WP \`8d3777762cb32a80e8c6ea097a0456e2b63a94c6\`、TPL \`d7bc9b69bf3cb11c3c3d090282723b6e73b1de55\`。GitHub Actions は **元binary \`.rda\` そのまま**取得し \`git hash-object\` でこれらの元publisher hashを照合。元R \`phylo$tip.label\`, \`edge\`, \`edge.length\` を基本Rで読み込み、**負/欠損枝長・重複tip/構造変化を許さない**。

**区別が絶対必要**：巨大系統樹の末端に種名が存在することは、その種自身のDNA配列を用いて枝が推定されたことを意味しない。末端が元treeに「ある」こと（source original direct backbone tip）／著者アルゴリズムが分類情報で種を「追加する」こと（graft）／分子配列の真正カバレッジを絶対に同一視しない。今回の監査はdirect tipの有無と真元treeの枝長・どの元author siteにそのtipを対応づけられるかだけを測る。**graftはしない**。

## 著者field communityと形質対象を固定

原author真正葉群植物元資料、Oldenburg LEDA元 \`canopy_height.txt\` と \`seed_mass.txt\` の先行SHA、v28.5欧州3-study 622 trait-known植物で固定した z(log canopy height, log seed mass)基準、実pub元geo/元site共同出現行列を一切変更しない。

本研究branchのPython helper は **Actions runnerの一時ディレクトリにだけ** \`(原出版, site番号, 元地理group, 元学名, zlog height, zlog seed)\` を書き、そのファイルを **GitHub artifact と gitに絶対含めない**。Rはこの元source植物を全3mega-treeに照合し、元研究×treeの12行の集計**だけ**出力。追加したsite×speciesは一切捏造しない。

出力する各出版×3treeメタデータ：
- Source元両LEDA対応speciesから真正backboneに**完全二名法一致でtipがある種数/比率**、現場元plant×site記録数に対する一致比率
- 元著者実際のsiteそれぞれに、直接tipが5種以上存在するsite数、元site正確母数、元geo区分に元source直接tip十分siteが2つ以上ある区分数
- 一致したtipの terminal phylogenetic edge長が正かゼロか、tree全体のsource edge品質
- 元site全体の**source草丈CWMβ**を既存v29.4真値で再照合。一致tip5種以上のsite cohortだけに制限したときの、元真全taxaとtip taxa双方のCWM height β（selection biasの大小を検査）。ただし計算上の候補値はまだ**新しい系統的効果量ではない**。
- 種名とsource genusが元backboneにいるかどうか。元genusがtreeにあるだけで、元speciesがtreeにあるとみなさない。
- 各treeとstudyのsource statusを、効果の方向を見ずにtaxonomy/coverageで判定。

## 事前に固定した「直接tipのみ」真phylogenetic inference gate

各独立publicationで、（1）元LEDA両形質既知の**60%以上の相異種が真正backboneに完全一致**、（2）元site×species出現の**60%以上が真正backbone tipである**、（3）実siteの**80%以上がtip対応5種以上**、（4）本来の地理groupの**2群以上に直接tip対応siteが各2以上**。満たしたoriginal publisherだけを \`DIRECT_BACKBONE_TIP_SOURCE_SUPPORT_PASS\` として、今後の真枝長patristic距離の解析に投入可能と判定する。満たさない元sourceは \`DIRECT_BACKBONE_TIP_SOURCE_SUPPORT_FAIL\` だが、「環境形質の効果が0」という意味は一切ない。

条件未達の場合、V.PhyloMaker2のS1/S2/S3の人工tip補完を許すかは、厳密に異なる次の**感度解析として別に事前登録**する。欠落種をすべて補完して「実際の配列がある種の系統樹」とすることは禁止。これらの閾値も最終の効果符号を見て変更しない。

## 本研究の到達点

本v29.5は**新しい生態係数を作るためではなく、初めて本物の枝長付き植物系統を4原著に結合できるかを判定する**。生物学的独立は「4 source field publications」であり、3系統樹を3独立研究/171地点を171独立メタ研究として水増ししない。taxonomic nomenclatureの違いは感度分析であり別の進化的実験でない。

元mega-treeや元source植物種×site/GPS、LEDA個別原データ、原著taxonomy listを本repoへ再配布しない。source author originalsとtreeはCI一時フォルダのみ、派生結果はsource QA数字だけ。

実装：
\`scripts/prepare_Europe_4original_rail_phylogeny_audit_source_v295.py\`
\`scripts/audit_Europe_original_4rail_VPhyloMaker2_three_megatrees_v295.R\`
\`tests/test_original_4rail_phylogeny_audit_source_v295.py\`
\`.github/workflows/tetsu-v295-four-original-rail-VPhyloMaker2-backbone-source-audit.yml\`
