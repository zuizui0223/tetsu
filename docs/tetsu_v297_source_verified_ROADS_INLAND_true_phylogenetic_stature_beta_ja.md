# tetsu v29.7 原著SOURCE VERIFIED — 道路際と内陸植生の草丈の系統的な空間構造

更新：2026-11-19 JST。元MIREN原著CSV、LEDA高さ＋種子重量、V.PhyloMaker2の公開R-data 3ツリーすべての固定ハッシュを再照合し、実元site×taxonの植物群集を固定した本物のデータ解析が GitHub Actions #44543008509 で成功した。実測値の正本は元Source artifact #14538279120（ https://github.com/zuizui0223/tetsu/actions/runs/44543008509/artifacts/14538279120 ）。元sourceの個別サイトの種名/座標/原author matrixやメガツリーraw tipは出力していない。

## 科学的に新たに分かったこと

欧州の元4独立RAIL出版のv29.6では、実published系統樹の枝長で切った75%相対深度の単系統群内に限って2形質を交換しても、草丈の地点間CWMβ距離は帰無より全4原著で大きかった（LCVP +0.056695 / WP +0.055429 / TPL +0.058998）。しかし鉄道固有かは不明だった。

今回は **鉄道とは独立した MIREN mountain ROADS 2007–2023 の元field program** で、同位置・同年・同100m²調査の道路際（Plot1）と内陸植生（Plot3）に同じ枝長75% clade-only trait-label nullを適用。**道路際・内陸の両方が、元属/種の地理的分布に伴う非ランダムな草丈の空間構造を持つ**ことが示された。したがって「真のphylogenyを固定しても鉄道4/4で正」というパターンは、**鉄道固有の因果選択の証拠とはならない**。

## 1．元source field/source taxonomyの成立性

Original MIREN 167,804植物出現行（Zenodo v3 DOI 10.5281/zenodo.15749536、MD5 b640195cc44f252ccb6ddb28eb6fff13）。5,739 unique binary species中1096 taxaのみOriginal Oldenburg LEDAの草丈＋種子重2形質が対応（19.1%、全球代表性には偏りあり）。独立した元野外研究は**1 MIRENプログラム**であり、8地域や5地域を8出版論文とは数えない。

元v24.9/v25.0の現地source選抜を変更せず再現：
- Primary source：5以上両形質種/元plot、25%以上の元taxa coverage、8元Region/22元Road/272元Site（Plot1/Plot3の同年1:1対応、同一Site初回適格年だけ）。
- Strict source：10以上両形質種/元plot、50%以上の元coverage、5元Region/13元Road/192元Site。同じ研究の入れ子感度サブセットで独立効果ではない。

系統樹の既存tipに植物種を完全一致させた後、元実Siteの80%以上が残る＋元同年両環境で各最低5種（strict10種）＋元2road以上・4site/road以上＋種源tip60%以上の事前source gateをかけた結果、解析対象の地域が縮小した。これは効果の向きに基づいて選別したものではなく、**系統への対応と実元plot richnessに基づくsource support**である。

| 本物の系統樹 | Primary source 8元regionから残ったコード / paired site | Strict source 5元regionから残ったコード / paired site |
|---|---|---|
| LCVP | AUT/CHE/IND/ITA/NOR **5地域・210サイト** | 同5地域・185サイト |
| WP | AUT/CHE/IND/ITA/NOR **5地域・210サイト** | AUT/CHE/IND/ITA **4地域・129サイト**（NORはsite保持80%未満） |
| TPL | AUT/CHE/IND/ITA/NOR/ORE **6地域・224サイト** | AUT/CHE/IND/ITA/NOR **5地域・186サイト** |

PrimaryでCLA/CLSは3元treeいずれもsite retaining不可能、OREはTPLのみ対応した。元data 8元regionを一律8地域とみなしたまま効果量を計算しない。元のMIREN2D道路際−内陸CWMβ −0.10381093795211047 も元完全SOURCEから再現した後にこの系統解析を実行。

## 2．正規化座標は鉄道と同じ

元LEDA log canopy height / log seed mass は前の鉄道3原著で凍結したmean（−0.580236575,+0.293611319）、SD（1.380100055,2.165766056）による z に変換。MIRENの種プールの平均を0に作り直していない。元植物種が道路と内陸に同時出現している場合、1つの技術ランダム化で両方に同じpaired height＋seed mass vectorを割り当てる。元source同じroad内の異なるsite同士の高さCWM絶対距離をRoadごとに平均し、同地域のRoadを等重み、最後に元Regionを等重みする。

## 3．SOURCE VERIFIED: 75%真枝長単系統群内の高さCWMβ 観測−帰無

系統樹の相対root-to-tip深度75%を切断して生じる元モノフィレティックclade内だけでLEDA草丈＋種子重paired vectorを交換。2固定seed×199=398技術的ランダムラベル。実系統樹への種対応を済ませた**同じSite×speciesで基準の頻度5分位null**も再計算（SOURCEに入る地域が違うので以前の鉄道原著とは単純加算しない）。

### Primary （元tree一致 5–6地域）

| tree | source直接tip集合の頻度5分位null：道路際 / 内陸 / その差 | 75%真系統clade内null：道路際 | 内陸 | 道路際−内陸 |
|---|---|---:|---:|---:|
| LCVP（元5地域） | +0.115676 / +0.174714 / −0.059037 | **+0.091204** (5/5正) | **+0.115357** (4/5正) | **−0.024153** |
| WP（元5地域） | +0.115416 / +0.174839 / −0.059422 | **+0.093536** (5/5正) | **+0.118308** (4/5正) | **−0.024772** |
| TPL（元6地域） | +0.089964 / +0.153331 / −0.063368 | **+0.069087** (5/6正) | **+0.104647** (5/6正) | **−0.035560** |

「道路際−内陸」75%残差の地域符号は、LCVP 2/5負、WP 3/5負、TPL 3/6負であり、**すべての地域で道路の分化が内陸より小さいとは言えない**。地域を1つずつ除外した等元地域平均の範囲は LCVP −0.039004〜−0.003555、WP −0.040603〜−0.007406、TPL −0.049986〜−0.021390。全て平均では負だが、元一研究プログラム中の地域の依存性・系統欠測・元種プールに依存しているため、global 95%CIではない。

### Strict （元tree一致 4–5地域）

| tree | 75%：道路際 | 75%：内陸 | 差 | Leave-one-region-outの差の平均範囲 |
|---|---:|---:|---:|---|
| LCVP（5地域） | +0.087338 (5/5正) | +0.102821 (4/5正) | −0.015483 | **−0.031714〜+0.004594** |
| WP（4地域） | +0.085510 (4/4正) | +0.107040 (3/4正) | −0.021530 | −0.046520〜−0.001389 |
| TPL（5地域） | +0.089739 (5/5正) | +0.102741 (4/5正) | −0.013002 | **−0.028921〜+0.006260** |

**最も頑丈な肯定結果は道路際の草丈空間構造の正方向**。内陸でも多数の元Regionで正だが、INDは負側になるなど異質性が残る。道路際−内陸の比較は平均では負であるものの、**厳格条件のLCVPとTPLのLOO感度ではゼロをまたぐため、全般的な道路上の機能収斂の普遍法則とは言えない**。

元tree切断深度25%・50%では交換可能な実matched taxonを含む単系統群が最低5つに届かず、**各tree・panelのすべてのregionでNOT IDENTIFIABLE**。効果がゼロ／有意差なしではなく元treeレベルで識別できない。事前閾値を後から緩めず、このsource supportがある75%のみ解釈する。

## 4．MIRENの新規結果をどう鉄道のマクロ解析とつなぐか

**生態学的な新しい反証**：鉄道原著で4/4正だった枝長75%clade帰無下の草丈の地理的構造は、鉄道のないMIREN山岳道路と自然側の内陸植生にも現れた。したがって共通する鉄道バラストの選択・遺伝適応だけを根拠とする**排他的な解釈は不支持**。このパターンは、各地域で供給される植物相の分類・系統組成と環境による空間的な選別のより一般的な結果と整合する。

**ただし正式な鉄道／道路／内陸の統合meta-effectではない**：MIRENは単一global道路研究、鉄道は4個の欧州出版。路線全体の「駅全域」と「100m道床地点」「100m²道路植生」は異なるfield grain、同地気候・植生位置・系統プールを外生的に一致させていない。標準化zの大きさをそろえても、同等の調査単位になったわけではない。

**機構候補を直接結果にするには**：同一研究内で道路際と内陸の地形・管理・基質・種プールを記録した元studyメタデータを使い、代表地で種の地理的配置と系統構造を同時に説明するメタ回帰や、real phylogenetic nullに環境と元taxon occupancyを両方条件付ける診断が必要。さらに鉄道特有性には稼働中バラストと同等の砂礫open nonrailの比較可能な原site×taxon field観察が必要。

## 5．SOURCE再現可能性・公開範囲

真正Source workflow https://github.com/zuizui0223/tetsu/actions/runs/44543008509
Full aggregate artifact https://github.com/zuizui0223/tetsu/actions/runs/44543008509/artifacts/14538279120

実装：
- scripts/validate_MIREN_original_roads_true_phylogeny_clade_negative_control_v297.py
- tests/test_MIREN_road_interior_true_tree_source_v297.py
- .github/workflows/tetsu-v297-source-original-MIREN-ROADS-three-megatree-clade-placebo.yml
- data/derived/tetsu_v297_MIREN_original_roads_inland_phylogenetic_height_summary.json

元のspecies×site/元GPS/元LEDA個別taxon値/元系統樹TIP名はGitHub artifactに出していない。source真正archiveに全地域等重みの解析出力はあるが、独立出版数は元鉄道4／MIREN道路1であり、感度版のtreesを独立研究として水増ししない。


## 6．追加の共通source地域による3系統樹間感度（事後の検証・主解析ではない）

元SOURCE artifactの元Region別符号付き効果から、**3系統樹すべての直接tip-support gateを通過した同じ真のRegionコードの積集合**を取った。この選抜は元tree source-supportという数理的条件から計算し、効果の符号で地域を選んだわけではない。ただし初期事前登録には含めていなかったため、**事後の補足感度**として主表から明確に区別する。地域が元MIREN一研究プログラムに属することも変わらない。

Primary：共通AUT、CHE、IND、ITA、NORの5元Regionのみ。

| 元tree | 道路際草丈 observed−真clade75% null | 内陸草丈 observed−null | 道路際−内陸の差 | Region除外平均差の範囲 |
|---|---:|---:|---:|---|
| LCVP | +0.091204 | +0.115357 | −0.024153 | −0.039004〜−0.003555 |
| WP | +0.093536 | +0.118308 | −0.024772 | −0.040603〜−0.007406 |
| TPL | +0.090792 | +0.113549 | −0.022757 | −0.037588〜−0.001844 |

TPL元Primary6元RegionのうちOREが5元Regionに含まれない。ORE単独はTPL道路際で負のsource効果−0.039438、内陸で正+0.060138。**原source6地域のTPL平均が他2樹と異なる主因は、元treeの生物学的枝が違うことだけでなく、元のRegion集合が違ったこと**である。同一5元地域では3樹の数値は互いに近い。

Strict：共通AUT、CHE、IND、ITAの4元Regionのみ。

| 元tree | 道路際 | 内陸 | 道路際−内陸 | Region除外平均差 |
|---|---:|---:|---:|---|
| LCVP | +0.081734 | +0.110339 | −0.028605 | −0.054620〜−0.006209 |
| WP | +0.085510 | +0.107040 | −0.021530 | −0.046520〜−0.001389 |
| TPL | +0.085258 | +0.110608 | −0.025349 | −0.050690〜−0.003783 |

ただしStrict共通4元Regionでも個別の道路際−内陸の差は負2/4と正2/4で**一致した向きではない**。regionの1つを除いた等重み平均が負でも、元4つの地域が全球の独立・代表的な無作為標本であることは保証されない。この感度の最も重要な意味は**選別されたSource region集合の違いを混ぜて系統樹間の効果を比較しない**ことであり、全球における「道路の草丈収斂の普遍性が統計的に確定した」とする根拠にはならない。
