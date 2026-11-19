# tetsu v29.6 — 実系統樹の枝長条件付き植物草丈空間分化：真正source検証済み

検証日：2026-11-19 JST。

**[元4原著資料・LEDA形質2表・V.PhyloMaker2元LCVP/WP/TPLの公開sha・数学15テスト・純Python原著解析 SUCCESS #44527225368](https://github.com/zuizui0223/tetsu/actions/runs/44527225368)**。実source派生artifact **[14532331216](https://github.com/zuizui0223/tetsu/actions/runs/44527225368/artifacts/14532331216)** が全study×tree×modeの完全なsource aggregateの正本。原著元個別site×taxon群集/植物種名/GPS、元tree生TIP一覧、元LEDA個票は公開artifactへ含めていない。

## 1．今回の本物の独立性と元source QA

4出版原著はCZ39駅全域、IT28駅全域、HR+BA15駅（2国で1独立出版）、SI89地点100m道床（駅全域でなく100m物理地点）。総数171はfield sitesであって独立研究171件ではない。3本の元樹は同じfield生物調査を別命名法で照合したもので、12独立研究でもない。元著者植物×site incidence、元author geography、種の占有頻度、元LEDA草丈・種子重量のpaired2D trait、元の欧州で凍結したz(log height,log seed mass)軸を保持した。

真正系統樹の既存tipに完全一致するtaxonの元field siteに5種以上が全171地点で存在。元tree×field studyの12組合せの元SOURCE tip種数・高さCWM observed値を先行v29.5真正sourceと再照合した。旧v29.5の地点一致率分母にあったRの data.frame列数誤りをPythonでは元site真正植物種数を使って訂正し、原source各siteの率が0〜1を外れたら計算拒否する。

元作者V.PhyloMaker2固定commit `7af3fb5152f691af2e4ec9d5e2e467d1b50505e9` のRData(.rda)本体を変更せずPython rdataで読んだ。3樹の枝長は全て正。root-to-tipの中央値は **400.787656**（source branch units）、最小〜最大差は約4×10^-6、CV約2×10^-10のほぼultrametricな木である。ただしoriginal backbone tipは各植物種のDNA配列が直接解析された証明ではない。追加のgraftはゼロ。

## 2．事前に固定した枝長条件付きtrait帰無、解析結果

同じtree内で同じ元研究の一致tip種集合を固定。元speciesがどのsiteにいるか、元plant height/seed massの対そのものは保持。帰無は元study内の元植物種出現頻度五分位でtrait vectorを交換する**基準null**と、実系統樹を原根→先端の枝長中央値の25%・50%・75%で切った単系統群の中のみでtrait vectorを交換する**3段階のtrue-clade null**。事前条件は2種以上の元LEDA両形質既知種を含む交換可能cladeが5以上、元site×taxon incidenceの25%以上が交換可能、original地理グループ2群以上で25%以上交換可能。2 fixed seed×199=398 technical shuffles。真のfield p値や独立地域実験ではない。

| 元original tree | 同じ一致tip種に対する出現頻度5分位null：草丈観測−帰無距離 | 25% clade内null | 50% clade内null | **75% clade内null** |
|---|---:|---|---|---:|
| LCVP | +0.067062（4/4正） | **不識別 0/4元出版** | **不識別 0/4** | **+0.056695（4/4正）** |
| WP | +0.062424（4/4正） | **不識別 0/4** | **不識別 0/4** | **+0.055429（4/4正）** |
| TPL | +0.066168（4/4正） | **不識別 0/4** | **不識別 0/4** | **+0.058998（4/4正）** |

表は同一研究の実地理内にある異なるsource site間の草丈CWM絶対距離についての元4出版等重み平均。25%では元植物がほぼ1巨大cladeに収まる、50%でも1〜2大cladeしかなく、「交換可能で異なるcladeが最低5群」というpredeclared識別性基準を満たさなかった。**0/4は検定不能であり、効果0ではない**。結果に合わせて25%/50%を弱い制約で再解析しない。

75%ではsource真正元独立field studyごとに14〜24の植物を含むclade群があり、元植物出現の97.6〜99.7%が2種以上のtrait-known植物を含む交換可能cladeに属する。以下は75%の新たな実source草丈観測−系統clade限定帰無の値（表示用丸め、完全精度と技術順位はartifactにある）。

| 独立真元鉄道field出版 | LCVP 75% | WP 75% | TPL 75% |
|---|---:|---:|---:|
| CZ39駅 | +0.030778 | +0.022670 | +0.034796 |
| IT28駅 | +0.035413 | +0.031801 | +0.035237 |
| HR+BA15駅（1出版） | +0.091946 | +0.105198 | +0.101976 |
| SI89地点100m道床 | +0.068644 | +0.062047 | +0.063980 |

各元作者研究×3代替樹=12表示のすべてで符号正。技術ランダム化の上側順位は0.0025〜0.0451。ただし原著研究は**4件だけ**、樹3本は非独立な感度解析で、技術順位は生態学的地域間のp値、無作為野外実験、全球普遍性の統計的検定ではない。

## 3．どこまで生態学的に解釈できるか

**成立する知見**：元植物群集の場所と本来のtaxonを全部保ったまま、実published treeの枝長に基づき高さ・種子重の2D traitを同じ深さ75%のclade群内で交換しても、元RAIL群集の草丈の地点間分化は帰無より4/4独立欧州出版で大きい。これは「植物群集の草丈の地理的構造は、かなり深い系統的グループのtrait値を保存しても完全にはランダムにならない」という観察的支持である。

**同時に制限される強い解釈**：v29.3のより強い「属内のみ」制約では同方向の一貫性が3/4に弱まった。したがって、系統樹に条件づけたpositiveが残ったことを「鉄道への適応的草丈変化が生じた」と呼ばない。75% cladeでも元属間の形質を交換できる余地が大きく、群集の属/系統の地域的構成と地理との共変動で生じ得る。源樹はメーカーのbackboneを採用し、個々の源種に新しい分子配列証拠を追加したわけではない。

非鉄道の同面積・同時期・開放/砂礫地対照が必要な「鉄道固有の機能的収斂」と、列車による種子散布・草刈り/車両風圧・集団内遺伝適応の因果検証は今回行っていない。source研究を増やすまで正式なrandom-effects予測区間も出さない。

## 4．再現可能性

Python-onlyの再実行コード：
- `scripts/analyze_four_rail_actual_3phylo_clade_null_python_v296.py`
- `tests/test_4rail_true_phylogeny_clade_null_python_v296.py`
- `.github/workflows/tetsu-v296-Python-noR-four-original-rail-source-phylogenetic-clade-null.yml`

R版はGitHub-hosted apt-getの503で停止したため、モデルの別実装として残した。新しい実証結果の正本は純Pythonの真正source成功 #44527225368。4 original field source×3 alternative treesの原資料条件を変えていない。元植物全サイト×種の個票は公開しない。