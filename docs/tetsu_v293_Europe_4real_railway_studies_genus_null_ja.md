# v29.3 — 鉄道植物の草丈空間分化は、属を固定しても残るか

2026-11-19 source-author provenance gating. **この文章の新しい結果欄は真正原source CI完走後にのみ追加**。本コードは先行研究の実種×siteを乱数で発生させていない。

## 反証する強い主張

v28.5原著4研究（チェコ39駅、イタリア28駅、クロアチア/ボスニア15駅を**1論文として**、スロベニア89地点100m道床）が実真正site×speciesを固定したまま、種に付くLEDA草丈+種子重2Dベクトルを研究内全種間か元頻度5分位内でshuffleしても、草丈の観測CWM地点間β距離は帰無よりすべて正だった。しかしv29.1/v29.2の*独立したMIREN山岳道路とその周辺内陸植生*では、広いtrait種間shuffleで正方向だった効果の多くが、**属内のみのtrait swapを認める厳しい帰無**で縮小することが分かった。これは種群集の形質と地理との関係に、元種の属入れ替わり・trait保守性が寄与し得ることを示す。

そのため直接問う：
> 4独立鉄道野外原著の「地理に沿う植物の高さの非ランダムな配置」は、元植物**属**構成と元属のtrait分布を保存しても残るか？

元属のtaxon identityは真正元author species binomialのみ。LEDA1種中央値であって、現地個体の進化的適応形質ではない。属は系統の粗い代用であり、分子系統・遺伝子フローやエピジェネティクスは未検証。

## 分析前固定したsourceと母集団

source原著をすべて元publisher SHA/MD5確認（過去v28.5と同一）：
- **CZ Kutlvašr2024** Preslia original S1+corrected S2、39駅元2020–22 floras。39は駅全域、100m²のplotではない。
- **IT Adorni2020** original Parma pdf SHA、28元駅。
- **HR+BA Jasprica2017** original 15駅、3元気候群、2国でも**1元独立出版**。
- **SI Šabić2026** 元同一source XLSX 401279B SHA、真の89個100m道床植物地点、15元国有路線区分、駅全域ではない。
- Oldenburg LEDA原source \`canopy_height.txt\` SHA \`cbf5e761390d78a843512dfc1139acfc1154d46d6c1995a827354840a9f78a92\`、\`seed_mass.txt\` SHA \`34aa12c333fc1ba6de2e560a495ac5d693aa0722a09ecf68418ad94ac1c19cef\`。元3欧州鉄道source学名622種で固定したz(log height[m],log seed mass[mg])基準をそのまま保持、スロベニア元植物で標準化し直さない。

元source \`site×species matrix\` すべて保持。v28.5の元2×299 source occupancy-quintile paired-trait label swapを**新しいコードで完全再現**し、元公表6桁表示値との差が5e-6以下でなければ停止（source四原著原資料 driftによる偽positive防止）。著者の地理グループ内siteペアの2D、草丈1D、種子重1D距離のみmain endpoint。元siteペアを171独立研究とも新地域とも数えない。

## 今回のさらに厳しい「属内だけ」帰無

元研究ごとの元trait-known taxonを学名属でpartition。2種以上の属でのみ、元height/seed ペアをその属内の異なる種へ再割当て。**原source plant species×site incidence / 地域 / 種ごとの本来の占有頻度 / 共存 / 種数 / 元元属ごとのtraitの値とjoint covariance**をすべて固定する。元1種属は元trait値を固定し、変化可能な種数と元site×taxon abundance割合を必ず別QAとする。元source 2 chain×299 independent technical label draws。技術draw順位を野外p値と見せない。

事前の独立研究単位source support（効果の符号を見る前にtaxon/plotから固定）：
- 研究で**2種類以上の両形質既知種がある属が10以上**
- 元source研究の合計site×taxon植物元出現における、属内で交換可能な種の割合 **25%以上**
- 元著者source地理グループのうち、同割合25%以上の地理グループ **2以上**
- いずれか不足なら \`NOT_IDENTIFIABLE_GENUS_STRATIFIED_TRAIT_SHUFFLE\`、生態効果ゼロとせずsupportのみ出力。少なくとも3独立原著が通過した場合のみ、原著等重みの効果方向とleave-one-publication-outを記述する。

## 期待する識別

1. **両4研究で属内でも高さ観測−帰無が正**：属の転換だけでなく、各属内の種とheightの空間配置との関係が残る可能性。ただしfield気候、管理、系統内繁殖、trait欠損を分離した因果ではない。
2. **4研究で属内のみnullがほぼ0**：前の種ラベルを全属に跨ってshuffleして得た大きなpositiveは、属の生物地理や属に共有する特性から生じ得る。**ただし動かせないsingletonが多いので、何%属が説明したと主張しない**。
3. **研究で符号不一致/源taxonomy support不足**：全球普遍性を否定した結果ではなく、元原著の種層・元field grainや元LEDAの北欧偏りが主たる制約。原著間ランダム効果メタ解析の確証へ進めない。
4. **MIREN内陸も正・鉄道も正**の場合：鉄道plant height spatial sortingの正の向きは決定論的な一般群集形成と整合する可能性が高い。**鉄道特異的**なストレスや適応の証拠とは見なさない。

必要な次のフィールド原著は実際の *稼働中のrail ballast vs 同じ面積・時期・基質条件のopen nonrail/road/gravel*。現時点で元比較対象が厳格に満たせないことはv28.7源Czech Veg reference supportによって明示されている。

## 検証ファイル

\`scripts/validate_Europe_original_4rail_studies_genus_trait_label_null_v293.py\`
\`tests/test_Europe_4original_rail_study_genus_null_v293.py\`
\`.github/workflows/tetsu-v293-original-EU-four-rail-genus-trait-null.yml\`

source SHAが真正でも元著者site×species・taxon元表やGPSを再配布しない。source各出版別macro aggregates・元source taxonomy supportのみ保存。
