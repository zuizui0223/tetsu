# v317 真正メタ原表：47研究205比較から鉄道は本当に何研究か（Mell et al. 2024）

2026-11-22 JST。出版社の原Supplementary Additional file 8を再取得し、XLSX本体の**SHA256 `0fd8ff9922c1138d5abd5289e743058b10be3b52c73d4d338cc9b9d2505def17`、64,630 bytes**を新規確定。元XLSXはREADME44行と`File`206行（ヘッダ1行＋真正205 case）の2シート、元`File`は**42列**。全研究の効果量d、分散、処置/対照と原研究IDが入っている。[第一原資料の照合成功 GitHub Actions #45009162011](https://github.com/zuizui0223/tetsu/actions/runs/45009162011)・[元XLSX原行由来集計JSONのartifact #14748957100](https://github.com/zuizui0223/tetsu/actions/runs/45009162011/artifacts/14748957100)。元研究を補助情報以上に転載せず、原著のXLSXはGitHubリポジトリへ再配布していない。

## 最重要監査：「205比較の鉄道メタ解析」ではない

Mell et al. (2024), *Environmental Evidence* 13:4 DOI [10.1186/s13750-024-00328-3](https://doi.org/10.1186/s13750-024-00328-3) は全交通インフラの植物調査を対象に系統的検索し、205効果量／47独立原著のQ2ベースを公開。真正Excelの`LTI`列の全行件数は：

| infrastructure（原著source LTIラベル） | 効果量行数 | 205行の比率 |
|---|---:|---:|
| Road | **127** | 62.0% |
| Powerline | **40** | 19.5% |
| Pipeline | **35** | 17.1% |
| **Railway** | **3** | **1.5%** |
| **合計** | **205** | **100%** |

元Excelの`Author_date`独立値は**47件**（元論文記述と一致）、**Railwayの3行は2原著だけ**：`Agence_Gaiadomo_2012`と`Everingham_2019`。後者2行は同じ24+24個体/plotの資料から種数と被覆という異なる応答を抽出。**3件を3研究扱いしない**。道路や電線の効果量に推論の重みを肩代わりさせない。原著でrailway占有割合を約1%と呼んだのは概数、実数は3/205≈1.46%。

## 鉄道3行の真正元効果と測定条件

処置がrailway verge、対照が別生息地。メタ原表`Hedges_d`と`Hedges_var1982`をそのまま転記し、**同じ植物種群・応答・対照ではないので統合しない**。

| 元論文ID／case | 国・原比較（鉄道vs参照） | 植物群・応答 | `Hedges_d` | `Hedges_var1982` | 原著バイアス評価 |
|---|---|---|---:|---:|---|
| `Agence_Gaiadomo_2012_1` | France 鉄道盛土20 vs adjacent similar vegetation20 | Plant_status Unknown（**全植物/在来外来の区別なし**）、種数 | **+0.643905** | 0.101403 | **Medium**、灰色文献、原表`Replicates=Pseudoreplicates` |
| `Everingham_2019_3` | Australia 鉄道盛土24 vs sclerophyll woodland fire trail24 | **Native**、種数 | **−1.633269** | 0.084288 | **Low**、元図5からの復元 |
| `Everingham_2019_6` | Australia 同じ24+24で別の原応答 | **Native**、被覆量（元`Cover_cm2`） | **−0.516972** | 0.084288 | **Low**、元Fig.S2Cからの復元 |

※元原著値を変更せず、独立研究から再算出した有意差と呼ばない。元メタ正規近似の個別95%区間を参照として計算すると、France種数[+0.020,+1.268]、Australia在来種数[−2.202,−1.064]、Australia在来被覆[−1.086,+0.052]。Franceは原著`Pseudoreplicates`と元バイアスMediumなのでこれらの形式的区間を独立野外CIと呼ばない。2本の母研究しかないので**研究間分散、予測区間、出版間メタ効果、緯度とのメタ回帰は不適格**。

**元研究の重大なソースラベル不整合：Franceの`Continent`が`North_America`と記録**されている。CountryはFranceなので地理結合時に原表のcontinentを盲目的に使わず、source typoとして監査に残す。別の研究に勝手に国を変更するものではなく、元著者の行ミスを追跡したもの。原分析に実際に使用したバイナリカテゴリの符号に関与したかは未確認。

## Everingham原研究の別の大きな進展：真正 raw species × site 行列が公開

*Everingham, Hemmings & Moles (2019)*, *Austral Ecology* 44:702–712, DOI [10.1111/aec.12719](https://doi.org/10.1111/aec.12719)の出版社[Supporting Appendix S4](https://onlinelibrary.wiley.com/doi/abs/10.1111/aec.12719)は**`aec12719-sup-0002-AppendixS4.xlsx`（約191 kB）、現場調査の原種・プロット個票**として公開されている。ニューサウスウェールズの森林/草地の生息地攪乱強度に沿うデータで、**鉄道盛土vs防火道ほか、道路芝生・道路肩・森林の比較**が原著同一地域で実施された。今回の3元効果のうちAustralia在来種数と被覆の元生物学的メカニズムを独立原表から復元可能か次に検証する。これはこれまでのWA(2008,120railway-only quadrats)とは別地域・別発表・別調査であり、まさに鉄道とnonrailの同地域比較の新原データ候補。

一方、Franceの`Agence Gaiadomo (2012)`は既存メタ原表からは処置側と対照側に種名の組合せがなく、近接対照の生物学的mechanismは要原原著。France 2025 RATP調査 [10.1111/wre.70025](https://doi.org/10.1111/wre.70025) は255種×394調査、同一地点で線路128・盛土266の間を部分的に対応できるが、元区画×種表は「Data will be made available on request」。**RATP盛土は線路敷の一部であり、鉄道外の自然植生ではない。**

## 研究判断

(1) **研究間の鉄道の普遍的種数増加は原公開meta datasetから支持不能**。観測符号が混在してもendpointとplantsのstatusが違うため「亜熱帯vs温帯の逆転」とは示せない。(2) 一方で、*Everingham et al. (2019)* は温帯オーストラリアの**元サイト×species＋生息地同地域対照**が公開されており、新しい実解析候補として優先する。(3) China–Laos Li et al. 2026（14駅×84区画、2022年のrail vs adjacent undisturbed habitat）も候補だが公開植物行列・土壌値は未取得。Figshare share link `https://figshare.com/s/0a7d857ba7c0a746ba4b` の中身を見たふりはしない。(4) 元meta表の原ケース数を増やしても独立対象研究は2のままなので、**原植物rawフィールドデータの新しい生息地対応比較**に方向を変える。

## 再現性

- [v317先に固定した抽出・メタ停止プロトコル](tetsu_v317_original_Mell_2024_rail_only_meta_source_plan_ja.md)
- [真正元Excel byte/hash+XML表構造+47原著×205cases+3Railway原行抽出](../scripts/audit_Mell_2024_railway_meta_original_xlsx_v317.py)
- [GitHub Actions #45009162011](https://github.com/zuizui0223/tetsu/actions/runs/45009162011)
- [原Excel3Railway元ケースのみ集約したaudit JSON](https://github.com/zuizui0223/tetsu/actions/runs/45009162011/artifacts/14748957100)
