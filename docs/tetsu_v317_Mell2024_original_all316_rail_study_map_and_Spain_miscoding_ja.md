# v317追加：316研究の真正原表で発見された鉄道研究の全体像と新しい独立データ候補

2026-11-22 JST。Mell et al. (2024) *Environmental Evidence* 13:4 DOI [10.1186/s13750-024-00328-3](https://doi.org/10.1186/s13750-024-00328-3) の[Additional file 5](https://media.springernature.com/original/springer-static/esm/art%3A10.1186%2Fs13750-024-00328-3/MediaObjects/13750_2024_328_MOESM5_ESM.xlsx) を真正新規取得。**713,116 bytes、元XLSX SHA256 `85042f37984a6f8f316f1f578e7f3093334ded057ebdc0c72437ac8b34579c90`**。元出版review元5シート：README67行、A全文除外3,674行、B全文合格379行、Cバイアス審査404行、**Dレビュー採用研究317行（先頭ヘッダ1＋実316研究レコード）**。この元review採択316**研究レコード**は、205メタ効果量行ではない。

[出典原表の再取得・原316行の照合 Actions #45009957309 SUCCESS](https://github.com/zuizui0223/tetsu/actions/runs/45009957309)・[原著元レコードからの集計 artifact #14748943077](https://github.com/zuizui0223/tetsu/actions/runs/45009957309/artifacts/14748943077)。

## 独立の「研究」「出版論文」「メタ効果量」を混同しない

元Sheet D `LTI`欄で**Railwayと明記された採択レコードは10行**、その出版資料の`ID_doc`は**8件**。France Agence Gaïadomo (2012) 1文献がQ1/Q2/Q5の3研究質問として3行に現れるため、この10行を10独立出版論文や10独立野外調査とみなしてはいけない。原著レビュー本文は「railways 10 studies」と記述するが、定義の元デザイン行/文献IDは区別する。Railwayの文字が`LTI verge`にあっても`LTI`欄では他カテゴリーとなる**追加1行**があり、全体は「11対象行、9文献ID」になる（下記原表エラー）。この拡張行は鉄道が真の処置であるか原著で確かめたうえで限定的に候補へ載せる。

一方**補足8の205効果量行**は、Railway 3行、全て**2文献**：Agence Gaïadomo (2012)が1行、Everingham et al. (2019)が**2つの異なる応答**。したがって文献を追加発見したこと自体は、独立の新しい効果量を得たこととは違う。

## Railway LTIに明示された8つの出版資料（10レビュー行）

| 原著元source文献 | 原著採択Q・元バイアス | 原source対照と応答 | 元個票解析可能性 |
|---|---|---|---|
| **Everingham et al. 2019**, Australia, DOI 10.1111/aec.12719 | Q2・Low | 鉄道盛土・道路芝生等 vs intact bushland、在来種数・最大草丈・種子重量 | 元Appendix S4のspecies field raw data公開記載は確認。**Wiley直リンク取得403で未取得** |
| **Delbaere et al. 2015**, France/Belgium border, report `En marge...` | Q2・Medium | railway / road shoulder vs **adjacent woodland/grassland-like matrix**。植物種数・組成 | [公開研究報告書の全文PDF](https://likoto.hypotheses.org/files/2021/09/27a.pdf)の所在を確認。調査区の種組成が掲載されるか要source監査 |
| **Červenková et al. 2016**, Czechia | Q2・Medium | **営業中の鉄道盛土 vs 廃止線の盛土**、植物被覆・種組成 | 鉄道内の管理・利用状態比較として有望。ただし非鉄道自然植生の対照ではない |
| **Agence Gaïadomo 2012**, France | Q1/Q2/Q5・Medium | 2012 RFFの鉄道盛土、草刈り/除草剤管理、周辺の類似植生、景観文脈。植物種数と多様性 | 元全文は灰色文献。Q2に1元効果dがあるが`Pseudoreplicates` |
| **Ogawa & Mototani 1991**, Tokyo, DOI 10.1007/BF02347125 | Q2・Medium | **鉄道・道路・空地・庭園など複数土地利用のタンポポ出現** | 同一都市の土地利用比較。単一属で、全球全植物群集の類似性評価には向かない |
| **Prober 1996**, NSW Australia, DOI 10.1071/BT9960057 | Q5・Low | White box woodlandの**rail easement / 道路 / 墓地 / ストック移動回廊 / 森林保護区**、在来・外来の種数と被覆/植物組成 | **同国同生息地タイプ内の複数土地利用比較**で最優先候補の一つ。原著全文・地点別行列は別途確認 |
| **Reisch 2007**, Germany, DOI 10.1007/s10592-006-9244-4 | Q2・Medium | 鉄道・舗装目地 vs 自然岩場の *Saxifraga tridactylites* 遺伝変異 | 種集団の遺伝機構比較であり、群集の植物相・温度指標には直接使えない |
| **Tiébré et al. 2008**, Belgium, DOI 10.1007/s10531-008-9386-4 | Q2・Medium | 鉄道・道路 vs 森林・草地などで*Fallopia*の出現 | 特定外来群の生息地選別と景観比較。群集全体のmeta効果ではない |

## 元データの実質的な分類エラー：Rivera et al. 2014

元Sheet D **`WOS_25893` (Rivera et al. 2014, *PLOS ONE*, DOI 10.1371/journal.pone.0101413)** は原研究の`LTI`が **`Waterway`**、その同じ行の`LTI verge`が **`Herbaceous railway embankment`**。これは鉄道を機械的にLTIだけで検索すると落とす研究。

[原著PLOS ONE本文](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0101413)を一次ソースで確認すると、**スペインMadrid西部で施工中の鉄道路線の盛土**で、**5ブロック×3処置＝15区画**、各区画6×9 mの**無作為化区画試験**。処置は上層土なし（control）、深さ10cm、30cm。元研究は土壌呼吸、微生物酵素活性、植物被覆、種数、群集組成を同年2010/2011で測定。原著は上層土追加が植物被覆・種数の増加を後年助け、10cmと30cm間の被覆・種数差は顕著でないと報告する。**土壌管理が植物群集を選別する機構の直接操作**として、列車の種子運搬と切り分けて検討できる。原2014論文の単純なRailwayコード誤表記を元review行のまま放置した場合、機構別証拠マップから1件抜ける。元研究は元国別の研究でも交通量・貨物ODを操作しておらず、列車散布と呼ばない。

本研究の原論文は**Spain**だが、review原Sheet Dは **`Country=USA`** と記載。したがって`LTI=Waterway`だけでなく**国フィールドにも不整合**があり、元Source誤りを盲目的な気候/国コードとして使えない。原著[Materials and Methods](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0101413)とレビュー行を対照して確認済み。

## 新しい実行可能性と現段階の限界

- 同年同地点の**鉄道vs自然植生**を優先：Delbaere et al 2015報告書、Prober 1996 woodlandの元site別データ、China–Laos 2026 14駅×84区画。資料にデザイン情報があるだけでは元種×区画・地点気候が入手できたことにはならない。
- **鉄道盛土内の管理の直接操作**：Rivera et al. 2014 の元ランダム割付5ブロック×3上層土処置。原著の図表と補足Table S1は公開、ただし図表平均×SEだけで15区画の完全原観測値を復元できない場合は、メタ元生データへの再構成を偽装しない。元source報告の**species richness 2011 P<0.001**を独立再計算したp値と書かない。
- **鉄道内の利用状態比較**：Červenková 2016 active vs abandoned rail. 交通量や列車付着を直接測定した実験ではない。
- **種レベルの自然対比**：Ogawa東京タンポポ、TiébréベルギーFallopia、ReischドイツSaxifraga集団遺伝。これらを群集全体の同一Jaccard効果量へ混ぜない。
- Wiley原著Appendix S4取得試行[Actions #45009541422](https://github.com/zuizui0223/tetsu/actions/runs/45009541422)は公式の3URLともHTTP403。**Source unavailable**であり、解析データを取得した・元24対24地点で新しい効果をfitしたとは言えない。

**結果の位置付け**：外部の広域公開レビューを丸ごと鉄道メタとして使う構想は不成立。一方、元source設計を解剖すると、広域温暖種偏向よりも、同国内の近接対照（France/NSW/Spain）と、土壌・管理という操作的な生態機構へ前進する具体的な公開研究が見える。

## 再現と原著

- [原著316 study 元Excelの取得と区画/論文設計の照合（SUCCESS）](https://github.com/zuizui0223/tetsu/actions/runs/45009957309)
- [元316レビューstudy IDを識別するソースコード](../scripts/identify_Mell2024_railway_studies_316_source_v317.py)
- [元メタ205効果量3鉄道比較の別の原表検査](tetsu_v317_Mell2024_rail_only_original_meta_eligibility_results_ja.md)
- [真正原study Qレコードのみのaggregate artifact #14748943077](https://github.com/zuizui0223/tetsu/actions/runs/45009957309/artifacts/14748943077)
