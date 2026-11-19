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


## v29.3 原著再取得・Source verified 実測結果（2026-11-19 JST）

**[元4原著再取得・SHA照合・元v28.5 source 2×299 occupancy-null再現・新genus null2×299・数学9テスト全成功 #44476119363](https://github.com/zuizui0223/tetsu/actions/runs/44476119363)**。原4出版author sourceには改変なし。真正集計 [artifact #14506416561](https://github.com/zuizui0223/tetsu/actions/runs/44476119363/artifacts/14506416561)、再利用用元出版等重み要約 \`data/derived/tetsu_v293_4European_original_rail_genus_null_summary.json\`。

元独立原著の属内交換可能性QA：CZ **507両形質既知種・100多種属・全元site出現の65.61%交換可能**、IT **271種・54属・56.04%**、HR+BA **190種・35属・48.42%**、SI **224種・42属・47.48%**。どの独立研究も事前の種数／25%incidence／原著geography 2群以上source gateを通過した。残りの出現は属内交換できないのでsource帰無の作用範囲が元v28.5より縮む。

### 1．元4原著すべて同じfield sourceで、従来の広い形質nullと属内trait-label nullを直接比較

| 実原著（HR+BAは1原著） | 本物のtrue field units | **草丈の観測−頻度5分位trait null**（v28.5 sourceを新runで完全再現） | **草丈の観測−属内trait null**（v29.3 新source結果） | genus限定の技術的上側順位（野外pではない） |
|---|---:|---:|---:|---:|
| CZ Kutlvašr2024 | 39駅全域 | **+0.0548185** | **+0.0078997** | 0.00668 |
| IT Adorni2020 | 28駅全域 | **+0.0388098** | **+0.0044320** | 0.22037 |
| HR+BA Jasprica2017 | 15駅全域（2国1出版） | **+0.0991891** | **+0.0172544** | 0.05843 |
| SI Šabić2026 | 89真100m鉄道道床地点 | **+0.0787812** | **−0.0044649** | 0.79132 |
| **元4出版等重み記述平均** | **4独立原著** | **+0.0678997** | **+0.0062803** | **属内は3/4正、厳格に全4独立研究同方向ではない** |

2形質同時CWMβも元広いtrait-label null +0.0643913（4/4正）から属内only +0.0050230（3/4正）へ小さくなり、Sloveniaでは属内 **−0.0081570**。種子重のみは属内平均 **−0.0009142**（1/4正）で、草丈と同等の普遍的形質選別と呼べない。

**元source4出版等重みgenus height残差** +0.0062803 は1出版ずつ除いても **+0.0026223〜+0.0098620**。しかし原著の数はわずか4で、離れた地域の無作為標本ではない。shuffle順位を野外実験p値や多重検定補正済p値に変換しない。特にItaly属内では+0.00443でもnull tail rank0.220、HR+BA +0.01725でも技術順位0.0584、Sloveniaマイナス。CZのみnull尾順位0.00668と比較的明瞭。**『属内分化が鉄道全4研究で再現する』とは言わない**。

### 2．科学的な構造：どこが再現したのか、どこが消えたのか

- 元v28.5の正結果 **4/4** は「実際の種組成を固定したときでも、属をまたぐ全体の高さの割当てはランダムではない」という結果として維持される。しかしその値は、**元属ごとのspecies traitsのまとまりを保持した条件**では同じ大きさを示さず、1/4研究では符号さえ異なる。前の『鉄道という反復した環境で草丈の空間構造が普遍的』を、**系統/属構成に依存しない独自選択の普遍性**と解釈するのは支持不能。
- 新たな外部の山岳道路・同所内陸MIREN v29.2では、元広いtrait-label nullでの正方向が、同じ属内のみの帰無でやはり大きく小さくなり、道路際−内陸の属内残差も基本/厳格で符号が変わった。これは、**植物が地域に供給されるときの分類群（属）構成と共有される生活史形質の重要性**という一般群集形成の有力仮説に整合する。
- **重要な限界**：属内だけの帰無では多数のsingle-species generaの元plant×plot traitが固定され、帰無変動の自由度・検出力は下がる。現地種の元LEDA matchingが系統偏りを持ちうる。表のnull差を『属が91%のCWM差を説明した』とすることは**厳禁**である。属の時間尺度が揃わず分子系統樹を固定した結果でもない。属と気候が共変動するため『遺伝的起源が分かった』とは言えない。

### 3．投稿できる主張の一文と新しい次段

**暫定一文**：
> *The apparent recurrence of stature differentiation across railway plant communities is largely contingent on the taxonomic scope of the null model: broad trait-label nulls show positive excess in four European railway studies, but genus-constrained nulls strongly attenuate the excess, including reversal in one study; analogous patterns occur in independently surveyed road and interior vegetation.*

ここまでの貢献は「鉄道の進化的適応」ではなく、**人工交通回廊で共通する機能的パターンを、地域種プール・属/系統構造と、現地フィルターに分ける**ための強い反証診断。次の外部原著で本当に必要なのは、物理的に同等の元調査面積・調査時期・開放度・基質で揃えた active rail ballast と open/gravel nonrail。初期の種数が少ないKashmir原区に無理やりtrait-missingnessを補完したり、非交通参考地を植物のSYNTAXONで後付け選択するのは避ける。別途全世界の植物系統樹を真のsource種にマッチして phylogeny-constrained null を作る場合でも、元source原著独立性は変わらない。

**解釈の判定**：分類学的共変動の可能性の支持を強めたが、鉄道固有の環境選択や進化を証明していない。mainへのマージはDraftの連鎖レビュー後のみ。
