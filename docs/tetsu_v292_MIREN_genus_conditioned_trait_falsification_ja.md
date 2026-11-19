# tetsu v29.2 — 山岳道路と内陸の草丈空間構造は属の置換だけで説明できるか

開始日 2026-11-19 JST。v29.1で新たに分かったsource-authenticatedパターンについて、**植物分類群の系統的な非独立性を考慮した強い帰無モデル**を事前定義する。v29.1は本物のMIREN道路のroadside/interior両方で高さ \`observed CWM β − regional trait-label null\` が正、Roadside−Interiorの超過差は基本8地域 **−0.06322014565930367**、厳格5地域 **−0.05122038893958619**（2 fixed seed ×299 technical assignments）。元欧州4鉄道原著v28.5のgrass statureの正方向も、鉄道固有の適応ではなく地域種プールとの対応かもしれない。

## v29.2で何を消せば判別できるのか

前のv29.1はRegion内の全植物種または頻度の似た種間で草丈と種子重を入れ替えた。その際、背丈の近い属同士が地理的にまとまっていると「実種の草丈が非ランダムに場所を占める」ように見える。**属内だけで2D高さ・種子重のtraitベクトルを交換**し、同一属に共通の傾向をそのまま保持する。

帰無の保持する要素：原著同じsource \`Region×Road×Site×Year×Plot\` 植物群集構成、種分布、元の地理/road/habitat、種数、すべてのspecies occupancy、真の属内のtrait組分布とheight–seed相関、**属間のtrait差**。交換する要素：同じ実属の複数の植物種の間のheight–seed 2変量の割当てのみ。

重要：属は真の時間較正済み植物系統樹の代用ではない。属内形質交換が未可能な1種属は**観測値をそのまま固定**するため、普通のtrait-label帰無よりはるかに変化できる状態が少ない。観測と帰無が似ても即「属進化のみが説明した」としない。

## SOURCE上のpredefined SUPPORT gate

元source自体から判定、事後に新たなCWM効果の符号を見てから選ばない：

- v29.1と**同じMIREN原著2007–23元CSV MD5 b640195cc44f252ccb6ddb28eb6fff13**、真正167,804元出現行、Oldenburg LEDA source 2元SHA、欧州3研究の622taxa由来の固定zlog高さ・種子重基準。
- v29.1と同じ各Region内 **original Plot1 roadside vs original Plot3 inland**、同元site、同年、同面積100m²、最初の採用年だけ、基本>=5/25%の8Region/22Road/272site、厳格>=10/50%の5Region/13Road/192site。元v24.9 CWM2D Road−Interior **−0.10381093795211047**、v29.1 Europe-z height raw Road−Interior CWM β **−0.0931493207367427**をsourceから完全再現できなければ停止。
- 源region内の**両形質既知植物が2種以上いる属が10以上**。
- 道路際と内陸の各方で、元trait-known植物 \`site×zone×taxon\` 出現の25%以上が**2種以上trait-known植物の属**に属する。
- これを満たさない地域は \`NOT_IDENTIFIABLE_GENUS_CONDITIONED_TRAIT_SHUFFLE\` として、原因（multispecies genera、固定singleton taxa、habitat別交換可能頻度）を公開し**効果ゼロとはしない**。少なくとも3つの元regionが閾値を通過した場合だけそのsubsetの等地域記述平均を出す。地域が変わればv29.1全8地域の値と直接比較しない。

帰無設定：一地域内の属ごと独立に元height/seedペアtraitラベルを置換、2fixed seed ×199 = **398技術label assignments**。各真road内の全組siteペアごとにCWM height、seed、joint2D距離を求め、road内平均→region内road等重み。placebo\`observed-null\`を道路際、内陸、両者の差で別々に算出。地域除外平均の範囲は95%CIでなく感度診断である。

## この検定の反証力と限界

| 新しいsource outcome（属内のみshuffle） | 説明 |
|---|---|
| genushuffle後も内陸・道路両者の草丈β超過が強く正 | 属の植物種構成を固定しても、**属内の種−形質割当て**に空間構造がある可能性。地理の因果効果ではない |
| genus shuffle後の両者の超過が小さい | 旧v29.1の広域traitシャッフルでの正効果は**属間の形質保守性を含む群集組成**と整合。ただしgenus singleton固定により帰無検出力は低い |
| 道路際−内陸の負の差が属内shuffle後も一貫して負 | 形質と属内の種配置の違いが道路・内陸の相対的分化へ寄与する可能性 |
| 属内交換可能な植物種が少ない | 生態効果は**識別不能**、属間・系統・形質被覆のデータ増量が必要 |

この追加のモデルも、MIRENという**同一山岳道路原調査**に対するmethodological sensitivityであり、独立出版研究数を増やすものではない。新しい関数の成功・GitHub Actions原資料再計算の成功を確認するまで実測の結果値を本文に書かない。

## 屋外鉄道との接続と次の独立データ

- 鉄道4研究の異なる種群集の草丈CWMβ分化(v28.5)は真のLEDA種中央値であり、現地個体の草丈可塑性や遺伝適応ではない。
- v29.2のroad−interiorは同実地点/同面積であるが、**鉄道バラスト−非鉄道砂礫地**の同条件比較ではない。
- 鉄道特異性を検証する主解析には、非鉄道referenceの物理基質・plot size・年季・標高・草刈りが植物組成と独立にわかる真正原調査を追加する。
- 既存候補はフランス2025 Rakotoson (10.1111/wre.70025, track/embankment 2m² data by request)、ポーランドDylewski2025 (10.1002/ldr.5383, 25 abandoned railway + 25 adjacent grassy meadow, data by request)、中国ラオスLi2026(10.1016/j.jenvman.2026.130380, 14 stations/84 plots; raw matrix availability not established)。後者2つも運行中バラストの同一管理条件を保証するわけではない。無断で「原データ取得済み」に変えない。

## ファイル

- \`scripts/validate_MIREN_genus_conditioned_trait_label_null_v292.py\`
- \`tests/test_MIREN_genus_conditioned_trait_label_null_v292.py\`
- \`.github/workflows/tetsu-v292-MIREN-within-genus-trait-source-null.yml\`

元著者植物種×site、正確な場所、LEDA元tableはGitHub artifactに再配布しない。regionコード別の支持域・平均効果などのみ集計。


## v29.2 実原著 SOURCE 完走結果（2026-11-19）

**[真正source Actions SUCCESS #44475032251](https://github.com/zuizui0223/tetsu/actions/runs/44475032251)**、原資料ハッシュ／9数学テスト／v24.9とv29.1の元field cohort baseline再現をすべて通過した。原元MIREN 167,804件、LEDA height/seed源2表、元source固定欧州z軸で **2 seed×299=598** 技術的な属内traitラベル再配置。新しい帰無の結果は原site×taxonを移動していない。

### 属内交換に必要な植物分類群の支持

- 基本5種・25%かつ8地域のうち、属内で交換できる2種以上のtaxaを持つ属が10以上、かつ両方の環境の源形質既知taxon×site出現の25%以上が属内交換可能だったのは **5/8地域：AUT, CHE, ITA, NOR, ORE**。その多種属の数は13〜89、対象サイト出現の交換可能割合は道路際55〜68%、内陸55〜66%。**CLA、CLS、INDはstrict source-support門で不識別**。これは『効果0』ではない。
- 10種・50%の厳格パネルでも **5/5地域：AUT, CHE, IND, ITA, NOR** が属内交換支持を満たした。基本5地域と厳格5地域は同一母集団ではなく、**OREとINDが入れ替わる**。厳格化でsource ecologyが変わったと解釈しない。
- 属は同じ種内系統を保証せず、元LEDA trait matchingはMIREN原5,739二名法のうち1,096 (19.1%) に限られる。原source multispecies属に同じ2D trait値があってもvectorの再配置は地理と関係なく事前条件から変えない。

### 1. 新しい genus-conditioned trait shuffle（ソース元regions等重み）

| 元site固定、属内traitシャッフル後の**観測CWMβ−帰無期待** | 基本（5元Region） | 厳格（別の5元Region） |
|---|---:|---:|
| 道路際草丈 height | **+0.012901**（4/5正） | **+0.023034**（5/5正） |
| 内陸草丈 height | **+0.016239**（4/5正） | **+0.014990**（4/5正） |
| **道路際−内陸の高さの追加超過** | **−0.003338**（2/5負） | **+0.008045**（1/5負） |
| 道路際2D(height+seed) | +0.013011 | +0.015563 |
| 内陸2D | +0.010697 | +0.009227 |
| **道路際−内陸2Dの追加超過** | +0.002313 | +0.006336 |

基本のheight追加超過差のRegion leave-one-out平均範囲 **−0.011145〜+0.007516** とゼロをまたぎ、厳格では **+0.003082〜+0.010838** の正。しかし基準のRegion構成が変わり、**Roadside−Interior の一貫した負のgenus内残差は識別されない**。technical traitラベルdraw598回は独立の野外field実験ではなく、Region leave-one-outは95%信頼区間でない。

### 2. V29.1 broad occupation-quintile trait-null と厳密に**同じ5地域**を再集計した比較

原著v29.1 SOURCE artifact #14500352634 の各region数値から元subsetを一致させた。*v29.1の公表全8地域平均*と新5地域平均を直接引き算しない。

| 元まったく同じMIREN field sites / 両LEDA trait、草丈 observed−null | v29.1 Region内全属横断、種占有率5分位内null | **v29.2 属内だけnull** |
|---|---:|---:|
| 基本条件の元5地域：道路際 | **+0.094161** | **+0.012901** |
| 同5地域：内陸 | **+0.179327** | **+0.016239** |
| 同5地域：道路際−内陸 | **−0.085166** | **−0.003338** |
| 厳格条件の別5地域：道路際 | **+0.116046** | **+0.023034** |
| 同別5地域：内陸 | **+0.167266** | **+0.014990** |
| 同別5地域：道路際−内陸 | **−0.051220** | **+0.008045** |

「元の空間配置を保持した高さ形質の非ランダムな関係」を表す広いtraitラベル帰無で大きな正方向は、属間のtrait差を固定すると**だいぶ小さくなる**。これは**属の種構成・生物地理との共変動が主要な候補**であることを示す。しかし属内nullはsingleton属や属に同じ値の種を動かせず、permutationsで変化できるtrait種が限定されるので、上表の差を「系統が何%説明した」といった寄与率に変換しない。

**反証された強い解釈**：「道路際の形質収斂が、属を一定にした上でも普遍的に発生する」という主張は不支持。元観測のRoadside−Interior CWMβ差自体は既存v29.1から消えたわけでないが、属内traitの再配置が示す追加的な差は小さく、条件間で符号が安定しない。

### 3. 主要な生態学的解釈と次の試験

鉄道v28.5は4原著で草丈 observed−**広域trait-label** nullが正だった。しかしMIREN道路・内陸でも同様の正は成り立ち、その多くが属間のtrait保守性・地理的な属配置に由来し得る。**従って鉄道特異的な草丈適応・剪定・車両風圧・列車散布を、現在のspecies-CWMデータだけで証明していない。**

**次の最も反証力のある実装**：MIRENのこのgenus-conditioned nullを、真正欧州鉄道4原著に全く同じ固定 source field \`site×taxa\` と厳格species→genus割付で再適用する。鉄道原著ごとの多種属数/交換可能source incidenceを先に検査し、条件未達なら判定不能にする。これで鉄道内部の草丈分化が、単なる同属・異属植物の地域置換を超えているかを比べられる。まだ非鉄道実測open/gravel対照は必要で、同じMIREN原sourceを5つの独立鉄道野外研究に水増ししない。
