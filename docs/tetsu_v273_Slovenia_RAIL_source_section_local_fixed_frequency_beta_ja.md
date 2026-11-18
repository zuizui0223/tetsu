# tetsu v27.3 — スロベニア真正鉄道原著89地点・15路線区間：同一地点の道床植物相は元頻度以上に結びつくか

2026-11-18 JST source-verified。真正一次資料は Šabić, Grašič & Rovanšek (2026), DOI https://doi.org/10.1007/s44473-026-00155-z。作者が公開した source original XLSX 401,279 bytes / SHA256 26be7ca9513dd54136d4ec499d8413c1cb288045f87fe3cf8814731f2edac99f。元92現地地点・482 zone行、実際に同じ線路側でBS道床肩とBT法尻の植物記録が一致する89現地点、173 source同側ペア、15元線路区間を原データから再構成。元EPPOコード312taxaと作者コードブックの一つの既知不一致LOTHe対DOYPHは原まま保持し、同種に勝手に修正しない。

[真正資料のGitHub Actions #44250304215 成功](https://github.com/zuizui0223/tetsu/actions/runs/44250304215)。数式・源側別・現地地点・元section階層の合成試験4/4成功、原XLSX SHA認証も合格。原元資料・種×現地plot・元座標はrepoに再配布していない。結果派生のartifactはGitHub Actionsに保存。

## 問いと元作者の現地設計

「共通する鉄道微生息地は植物群集を均質化するか」には、実際には元植物の種同一性と同じ地点の現地植物相の両方が重要。v25.7はSlovenia原本を元siteごとに左右のBS/BT植物リストへまとめ、元source同一元section内で「同siteの共有種」と「別の同section siteへ参照を付け替えた場合の種共有」を比較した。元14の2site以上のsection（88現地地点）で同site Jaccard similarity **0.3632363548846349**、異site再割当平均 **0.2003803322207602**。しかし、元参照の割当先を変更すると各siteのreference種数も変わるので、この比較だけでは現地固有の豊富さを固定できない。

今回v27.3は **corridor BS実植物種集合を固定**したまま、BTのoriginal source site×species二値行列をsource同一線路区間内の2×2 checkerboard swapで組み替える。BT各現地地点の実種数は完全保存し、BT各元種の同一source line section内出現頻度も完全保存。各実original siteの右/左のペアはまずsite内でunionしてから分析し、sourceの173左右ペアを173独立現地点と誤らない。元section中1siteしかない区間はsource randomization不能なので除外、元14区間88siteが主比較。source4site以上の9区間76siteを厳格感度条件とする。

同siteでBS/BTの共通種a、BS-only b、BT-only cとすると、beta_jac=1−a/(a+b+c)であり、site平均を等重み→source sectionを等重み。各sectionの元species rich/frequencyを完全固定したまま99 Markov draw（burn-in3000、draw間提案350）で期待値を計算。MCMCは相関する技術サンプルであり、source sectionや植物種を無作為実地実験したp値ではない。

## 真のsource結果

| Slovenian source同section fixed-margin original plant beta | 実測 | 帰無平均 | 実測−帰無 |
|---|---:|---:|---:|
| 14 original sections・88 actual sites：Jaccard非類似度 | **0.636764** | 0.755724 | **−0.118961** |
| 同じ88site：共有taxa/site | **11.6210** | 8.5827 | **+3.0382** |
| 厳格9 source sections・76 sites：Jaccard非類似度 | **0.654320** | 0.789981 | **−0.135661** |
| 同じ76site：共有taxa/site | **10.9104** | 7.4921 | **+3.4183** |

14/14 source sectionsで実際のBS/BT Jaccard非類似度は元頻度固定帰無より小さく、実際の共有種数は期待より多かった。厳格化しても9/9 sectionsで同方向。元source 1site sectionはnullを作らず、ゼロ種数や欠測を偽って入力しない。原14 sectionsはいずれもaccepted checkerboard swaps>0（nullの機械的自由度あり）。特にv25.7で得た原source Jaccard similarity 0.3632363548846349 は、今回元本から **1−βjac=0.3632363548846349 として1e−10以内で完全再現**され、独立に計算された過去結果との一致も確認した。

## 他のsourceとの比較の上限

v27.0のMIREN道路原本（2012年以降first survey、16Region/48Road/669現地site）でも、同siteのPlot1道路際とPlot3内陸が、Road内のsource site×taxon frequency fixed-margin randomizationより **βjac−null=−0.11727、共有種+3.9114**であった。Slovenia元Rail v27.3は **βjac−null=−0.11896、共有種+3.0382**。したがって「鉄道でも道路でも、現地微生息地同士の種共有には、元site species richnessと元区間内species frequencyで説明しきれない空間的な情報がある」という**独立sourceでの方向性再現**が新しく支持された。

ただし、**Sloveniaの比較は道床BSとBTという双方RAIL内部のzone**で、MIRENの比較はRoadsideと52–102m内陸であり、元採取区の面積等価性もSloveniaでは保証できない。さらにRAIL源は15 source sections / 一国1論文、MIREN源は同一 consortiumの複数 Region / 一研究、空間ランダム化範囲も元sectionと元Roadで異なる。−0.119≒−0.117という近さは偶然あり得て、鉄道と道路の「同じ大きさの普遍効果量」と解釈しない。

以前のKashmir RAIL v27.0では駅 vs 駅間 stratum単位で元内陸種頻度を固定しており、実験条件が空間的に広く、2014/2017の相対的local signalが不安定だった。Sloveniaは元現地物理rail section内という、Kashmirより局所化したrandomization。両原著の真の現地調査法・対象年次・地理・zone配置が異なるため、差を「Sloveniaの鉄道輸送の方が強い」等とはしない。

## 生態学的な中心仮説と次段の判別

観測が整合する仮説：**交通回廊で種構成が変わっても、元地点・元周辺の植物相と共有する植物に基づき群集が組み立てられる。共通基質は植物種同一性の全球的均質化を保証しない。**

二つの独立原著で方向性が再現したのは「同site zone coupling」という記述統計であり、**種子源が同じこと、受粉、種子搬送、鉄道の貨物流動、維持管理、風選択、表現型可塑性、遺伝的適応を直接同定していない。** 同地点の環境・土壌・気候・観測努力・生育履歴もsite-level taxon-sharingを生む。次の科学的に必要な資料は、RAIL vs ROADを同地域・同面積・同季節で配置した真正非鉄道対照、現地列車運行量・管理実測、別の独立RAIL原著source種×zone行列。特に元線路通過数27site/9交通区間だけのSlovenia交通量sourceを大きな独立社会実験と誤らない。

再現：
- scripts/audit_Slovenia_original_rail_15sections_fixed_local_frequency_null_v273.py
- tests/test_Slovenia_original_fixed_section_frequency_null_v273.py
- .github/workflows/tetsu-v273-Slovenia-source-section-site-fixed-frequency.yml

これらは実原著データを取得/照合して導出した二次解析であり、原著者とMIREN consortiumそれぞれのデータ利用・論文引用ルールを遵守する。
