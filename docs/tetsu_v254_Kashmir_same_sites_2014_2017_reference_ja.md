# tetsu v25.4 — カシミール31地点の同一年・同一地点を固定した鉄道β多様性と時間的置換

更新：2026-11-17。研究ブランチ `research/tetsu-kashmir-temporal-reference-v253-20261117` の追加検証。`main` 未変更。

## v25.3 の明確な訂正

前回の「当該集約JSONでは2014と2017が同一IDかは確認できない」は **JSONファイルだけを見た場合** は正しいが、元の計算に用いた `research/validation_20261117/kashmir_paired.py` まで読むと、`complete` は `all(year, zone for both years)` を要求しており、`all_species` の各比較内では **同じ空間地点を2014・2017年の両方で使用している**。したがって、前回「2014と2017の採用地点が異なるかもしれない」と推測した部分は撤回する。一方、`Track–Interior Plot 1` と `Track–Interior Plot 2` **どうし** は使用地点集合が異なり得るため、両対照の結果差は地点構成と参照区選択の効果を混ぜていた。

**根拠**：当該コードの `complete = [s for s in population if all(np.isfinite(plots[s,y,z]).all() and (plots[s,y,z]>0).any() for y,z in product(years,[rail,ref]))]` と、その `complete` を両年の taxonomic Jaccard 計算に共通使用するループ。

## 新しい主解析（結果は真のsource Actions成功まで未確定）

真正一次資料：Rashid et al. (2021), JAE DOI https://doi.org/10.1111/1365-2664.13961 ; first-party Dryad DOI https://doi.org/10.5061/dryad.0p2ngf20h。公式 version 99735 の ZIP にある `Matrix_plot.xlsx` を取得し、真正 SHA256 `529f3b97521933018e12623b61f531f012c67241cae2602a6675dbd534192474` で **fail closed**。元123種×248 site-year-zone 行、31空間地点（16駅・15駅間）、4欠損セルを厳格に照合する。予めsource認証済みの `research/validation_20261117/kashmir_results.json` は Git blob SHA `75dd78dd2b66a67543c248ff686c3a9dc016ea73` で検証し、その12群の **raw Jaccard 観測差を 1e−10以内で再現しなければ結果を出さない**。本研究では原生物行列・個別GPS・地点ラベルを再配布しない。

### A：2種類の内陸対照を「完全に同じ地点」に制限

`station`, `between` それぞれで、同一空間ID × {2014, 2017} × {Track, Interior Plot 1, Interior Plot 2} が全て有効な地点の **積集合** を事前に抽出。内陸1・2で変わるのは参照ゾーンだけ、観測Track・原空間地点・年・種の全体リストは変わらない。

各群・年・ゾーンの異地点ペア間の平均 Jaccard 非類似度を `B(group,year,zone)` とする。`ref=I1/I2` について、

`D_ref = {[(B_station,2017,Track − B_station,2017,ref) − (B_station,2014,Track − B_station,2014,ref)]} − {[(B_between,2017,Track − B_between,2017,ref) − (B_between,2014,Track − B_between,2014,ref)]}`。

**基準変更の差 `D_I1−D_I2` では Track が完全に相殺され、内陸1と内陸2の地点間βの年次差だけが残る。** したがって基準変更後に符号が逆転しても、それを列車効果の違いと呼ぶことは数学的にできない。この不変性をコードで 1e−11以内で検証する。

### B：各地点自体の時間的入れ替わり

同じ地点 `i` の2014対2017 の種組成の Jaccard 非類似度 `T(i,z)` を、Track/I1/I2 で別々に計算。各地点の `T(i,Track)−T(i,ref)` を群内平均し、`station−between` を取る。これは **地点間収斂β** とは全く別の **同一地点の時間的組成変化** を測る。

### C：種数と疑似独立性への感度

- 地点ペアの生データでの分類群Jaccardを主解析にする。
- Aについて、同一地点・同一年にある Track/I1/I2 を3区共通最小種数へ199回希薄化し、元ペアのJaccardを平均。異なる参照で群集母集団・Track希薄化drawを変更しない。
- Bについて、同一地点における2年×3区の計6種集合を共通最小種数へ199回希薄化。
- 元空間地点IDを丸ごと再標本化する1999回の **技術的クラスターbootstrap** を実施する。異なる地点の組合せだけを積算する多重度重み付き U 統計を使用し、同一地点を二重抽選してできた自己比較 Jaccard=0 を標本のように追加しない。さらに元地点1つ除外の影響度を点検。
- Bootstrap分位幅は原標本内の安定性診断であり、ランダムに割り当てた鉄道交通の95%因果信頼区間ではない。

## 分析の停止条件

原XLSX SHA不一致、公開先ダウンロード不可、原source行・列・欠損セルの支持域差、既存12 Jaccard値の非再現、対象群の共通サイトが各5未満なら **新しい生物学的結果は公表せず** ワークフローを失敗させる。`Station/Between` の立地差は列車頻度・経路・貨物OD・管理履歴と同義ではない。外部対照のないこの再解析から世界的鉄道均質化、列車種子輸送、草丈の進化、刈取り・風への適応を主張しない。

実装：`scripts/analyze_Kashmir_same_sites_reference_year_v254.py`、合成テスト `tests/test_Kashmir_same_site_reference_v254.py`、原データActions `.github/workflows/tetsu-Kashmir-original-site-reference-year-v254.yml`。CIには元原著/LEDA元行列を保存せず、aggregate summary JSONだけを保存する。
