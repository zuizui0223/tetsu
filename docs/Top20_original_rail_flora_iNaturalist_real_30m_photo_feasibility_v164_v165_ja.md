# 鉄 v16.4–v16.5：実原著から候補20種を固定し、iNaturalist写真×線路30mの可用性を監査

更新：2026-11-15。**v16.4元植物原表からの候補20種選定は本物の原著データで成功。v16.5の写真×30m実件数は真正API実行後に追記する（原本取得前に件数を捏造しない）。**

## 目標：種間の環境フィルターと種内の表現型シフトは同じ向きか

大きな問い：鉄道が強く選別すると予測される植物形質軸は、鉄道上に定着した**同じ植物種の集団内の形質シフト**にも表れるか。元データでスロベニア鉄道枕木の道床肩BSは、近接道床法尻BTより植物被覆が低く、原著の厳密annual比が+5.47ポイント高いことが確認されている。一方、wind diaspore原著形質の地理相関と列車風自体を混同しない。群集がannualへ収束することと、同じ種が鉄道で小型化/早期開花することは**異なる推定対象**。

種間 A / 種内 D を共通の reference trait を用いた exact additive CWM 分解にする：
- 対照Cと鉄道Rで、種iの真の元個体比・被度比をp_i,C / p_i,R、地理ごとの真の平均形質をz_i,C / z_i,Rとし、各種に共通の外部reference値z_i,*を定める。
- すると差ΔCWM は厳密に
  **Σ_i(p_i,R−p_i,C)z_i,* + Σ_i[p_i,R(z_i,R−z_i,*)−p_i,C(z_i,C−z_i,*)]**
  と書ける。前者は元種構成の入れ替わり（組成／種間）、後者は種内形質のずれ（ただし各生育地の種内実測が必要）。写真だけでは種内の可塑性と遺伝的適応を分けられない。各生育地に実出現していない種のzを任意に補完しない。元A/Dの分類ラベルは既存の「種子A」「管理D」のsource台帳とは対応が違うため、実装上は turnover / ITV と呼び分ける。
- 種間の鉄道環境におけるreference最適値θ_railを推定するには、**鉄道と非鉄道の同じ地域・努力・面積の比較フロラ**と独立の形質値が必要。現存するチェコの39鉄道駅フロラだけから「鉄道への選好順位」やθを独立に推定することはできない。

## v16.4 原著20種の候補定義（SOURCE VERIFIED）

[**成功元source #43692256719**](https://github.com/zuizui0223/tetsu/actions/runs/43692256719)。スロベニアの元著者2022–23年の89地点実地鉄道100m調査×312 taxon traits/codebook、チェコの元著者2020–22年39全駅REC調査（過去の非均等historyは使わない）を各原著の元SHAと外部原著付表で検証。元地点の**存在有無**を種ごとに数え、各原著で出現地点が多い順にrank、国ごと交互（SI奇数rank・CZ偶数rank）に候補を追加し、全20種は学名で重複させない。これは植物の**鉄道選好率**ではなく、鉄道出現記録が原著で比較的多い種のセット。非鉄道側の母集団出現数を事前に持たないため、source preferenceランキングとは呼ばない。

**重要な原著taxonomy修正**：スロベニアの元種出現は学名ではなく**EPPOコード**。元のspecies_codebookの明示的な学名欄からだけ照合し、原著のLOT​​HE / DOYPH（Lotus herbaceus）コード不整合を勝手な推測同定で結ばない。さらに「Taraxacum sect. Ruderalia」は分類学上の節であり、属名＋第2語があるからといって**1種（Taraxacum sect）として扱わない**。正規の二名法の種に限り候補を選び、同一種の別名による重複を可能な範囲で監査する。

**最終固定の原著20種**（奇数rankはSI、偶数rankはCZ；括弧内はその原著で観測された地点数）：

| 固定順位 | 元選定国 | 学名 | 元出現地点 |
|---:|:---:|---|---:|
| 1 | SI | *Convolvulus arvensis* | 71/89 |
| 2 | CZ | *Arrhenatherum elatius* | 39/39 |
| 3 | SI | *Erigeron annuus* | 69/89 |
| 4 | CZ | *Artemisia vulgaris* | 39/39 |
| 5 | SI | *Vicia cracca* | 68/89 |
| 6 | CZ | *Conyza canadensis* | 39/39 |
| 7 | SI | *Equisetum arvense* | 67/89 |
| 8 | CZ | *Daucus carota* | 39/39 |
| 9 | SI | *Lactuca serriola* | 65/89 |
| 10 | CZ | *Echium vulgare* | 39/39 |
| 11 | SI | *Galium mollugo* | 63/89 |
| 12 | CZ | *Hypericum perforatum* | 39/39 |
| 13 | SI | *Setaria pumila* | 62/89 |
| 14 | CZ | *Linaria vulgaris* | 39/39 |
| 15 | SI | *Rubus caesius* | 51/89 |
| 16 | CZ | *Lolium perenne* | 39/39 |
| 17 | SI | *Euonymus europaea* | 50/89 |
| 18 | CZ | *Medicago lupulina* | 39/39 |
| 19 | SI | *Geranium robertianum* | 50/89 |
| 20 | CZ | *Plantago lanceolata* | 39/39 |

*Conyza canadensis*はiNaturalistで*Erigeron canadensis*として認識されうることが、先行の別source監査で確認されている。その狭い別名だけをsource名と独立に許し、他の属の別種を混ぜない。*Equisetum arvense*は胞子を出すシダ植物の仲間で**花を作らない**ため、花期シフトの候補種として流用せず、肥沃茎・胞子嚢穂などの別形質軸のみ対象とする。

## v16.5–v16.6 真正iNaturalist公開写真×元実線路30mのSOURCE VERIFIED結果（2026-11-15）

[**20原著植物種を2国の真正2023/2024 OSM線路Polylineに30m距離で結合した実行成功 #43692490652**](https://github.com/zuizui0223/tetsu/actions/runs/43692490652)。2018–2025年iNaturalist research-grade・写真あり、公開座標の元位置誤差≤30m、真の線路Polylineまで点距離≤30m。元APIの20種×SI/CZ各1矩形の**40/40組を全件取得**でき、API件数超過による推定不能0組。**SI矩形40件、CZ矩形969件、2国合計1,009件**。うち元CC利用可ライセンス写真が少なくとも1枚ある観察は**732件**、さらに**点線路距離＋元GPS誤差半径≤30m**のより厳格な記録は**858件**。どれもその撮影個体が真の鉄道敷地内にいたことまでは保証しない。

続く[**v16.6 独立ソース再検算 SUCCESS #43693580409**](https://github.com/zuizui0223/tetsu/actions/runs/43693580409)で、2つの成功source artifactsを再取得して、元原著候補20種・元OSM PBF SHA・観察者・100m独立セル・100–1000m対照写真の一致を確認した。写真件数を見**る前に固定した**「同種×元国矩形で30m以内30件以上、観察者10人以上、100m独立空間セル10以上、100–1000m対照30件以上」の**最低写真機会ゲートに通ったのはCZ矩形の11種、SI矩形0種**。ゲートは写真が存在する必要条件で、全草姿・草丈・実際の管理ストレス/遺伝適応の検証ではない。

### 実20種すべての30m厳格写真数と全球写真総数の区別

**SI/CZ欄はあくまで各研究国を覆う固定矩形内で線路30m以内を示した公開位置記録であり、最右列の世界写真総数は全球の線路30m以内ではない**。元チェコでの*Echium vulgare*は30m以内178件、独立観察者126人、独立100mセル170地点、100–1000m写真725件と次段の実現性が高い。

| 原著20候補順位 | 種名 | SI 30m | CZ 30m | 2地域合計 | CZ写真ゲート | 全球写真付きresearch-grade総数（線路30mではない） |
|---:|---|---:|---:|---:|---|---:|
| 1 | *Convolvulus arvensis* | 2 | 23 | 25 | 未達 | 92,081 |
| 2 | *Arrhenatherum elatius* | 2 | 19 | 21 | 未達 | 8,837 |
| 3 | *Erigeron annuus* | 6 | 77 | 83 | CZ可 | 41,496 |
| 4 | *Artemisia vulgaris* | 1 | 64 | 65 | CZ可 | 109,519 |
| 5 | *Vicia cracca* | 2 | 41 | 43 | CZ可 | 81,062 |
| 6 | *Conyza canadensis* | 1 | 73 | 74 | CZ可 | 67,805 |
| 7 | *Equisetum arvense* | 6 | 90 | 96 | CZ可 | 75,149 |
| 8 | *Daucus carota* | 0 | 28 | 28 | 未達 | 144,306 |
| 9 | *Lactuca serriola* | 0 | 50 | 50 | CZ可 | 49,588 |
| 10 | *Echium vulgare* | 3 | 178 | 181 | CZ可 | 85,996 |
| 11 | *Galium mollugo* | 1 | 1 | 2 | 未達 | 25,042 |
| 12 | *Hypericum perforatum* | 2 | 27 | 29 | 未達 | 66,571 |
| 13 | *Setaria pumila* | 1 | 11 | 12 | 未達 | 13,692 |
| 14 | *Linaria vulgaris* | 6 | 81 | 87 | CZ可 | 90,133 |
| 15 | *Rubus caesius* | 0 | 10 | 10 | 未達 | 19,863 |
| 16 | *Lolium perenne* | 0 | 7 | 7 | 未達 | 18,879 |
| 17 | *Euonymus europaea* | 0 | 0 | 0 | 未達 | 33,191 |
| 18 | *Medicago lupulina* | 1 | 40 | 41 | CZ可 | 74,078 |
| 19 | *Geranium robertianum* | 3 | 88 | 91 | CZ可 | 126,734 |
| 20 | *Plantago lanceolata* | 3 | 61 | 64 | CZ可 | 151,658 |

同じ元20種で*Erigeron annuus*（CZ 77件、観察者47、100mセル74）、*Linaria vulgaris*（81/57/77）、*Geranium robertianum*（88/46/80）、*Artemisia vulgaris*（64/45/59）、*Plantago lanceolata*（61/47/59）、*Conyza canadensis*（73/32/67）も空間・観察者が複数あり、写真自体に測定可能な全草姿/花のROIが含まれるかを監査する価値がある。*Equisetum arvense*はCZ90件/55観察者で可用性ゲートを満たすが**胞子植物なので開花期研究に使用禁止**。*Vicia cracca*はCZ 41件ながらうちCCライセンス写真付きは28件で、著作権・全身視認可能性の別チェックが必要。

**大きな留保**：チェコ・スロベニアの2矩形を調べただけなので、要求された**全球の線路30m写真数は依然未確定**。最右列は世界に存在する写真付き観察機会であり、鉄道網に近い確率を世界へ外挿できない。次の全世界検証には公開iNaturalist bulk/DwC元観測（位置誤差と写真属性付き）と世界のOSM鉄道Polylineを統合する必要がある。

[v16.5 元API実件数＆近線路真ポリライン検証](../scripts/audit_20_original_rail_plants_inat_30m_two_regions_v165.py)・[v16.6 元40種×地理矩形のobserver/100mセル/対照写真数を再検証したコード](../scripts/audit_iNat_original20_SI_CZ_30m_source_eligibility_v166.py)・[v16.6真正集計CI](../.github/workflows/rail-iNat-20species-photo30m-source-gate-v166.yml)。v16.6の成果ZIPには40行の元国矩形別ソース集計CSVも含まれるが、元写真・観察者ID・厳密個別座標は共有しない。

## 事前に固定する写真解析への「次段階ゲート」（写真件数を見る前）

次段階で特定種の高さ・草姿・開花可視性を比較するための**最低のサンプリング機会**は、元観測日と写真を人為選別しない完全集計から、同じ種・同じ研究地域で線路30m以内の厳格な公開精度記録が**30件以上**、公開観察者10人以上、独立100mセル10地点以上、同じ国矩形内の100–1000m地理距離帯で**30件以上**とする。この基準は必要条件にすぎず、Dの種内適応を結論する十分条件ではない。現状のスクリプトは元年別分布と全草丈ROIの撮影可否・区域の真の管理履歴を検証していないため、このゲートを超えても直ちに遺伝的適応や一律の花期進化推定へ進まない。

**全草丈写真の尺度問題**：同一植物が写真に入っていても、撮影距離やカメラ投影が不明なら絶対草丈は原理的に復元できない。全景の高さ/茎節数/花の可視性などが測れるか別途人手で検証する。開花日を写真の撮影日・Flowers annotationだけで「初花日」とみなさない。*Equisetum arvense*については開花形質そのものが存在しないので、胞子嚢穂や肥沃茎に限定する。

**推定対象の境界**：基準を超えた候補だけで、地理地域・年・月・路線の現場景観や道路近傍をそろえた**可視的形質の種内差**を記述する。写真にない遺伝的適応・草刈り処理を推定しない。「もとの最適値θから遠い種ほどシフトが大きい」という仮説には、写真の選別とは独立のregional reference traitと、同じ努力で記録された非鉄道植物群のtrue site × species sampling frameが不可欠。

**CWMの計算仕様**は [純数式・無実データのturnover＋ITV関数](../src/tetsu/rail_turnover_itv_accounting.py) に実装済み。同一地域内の鉄道Rと対照Cで個体数/被度比pと実測値zがそろい、別途参照種平均z*が得られたときだけ、ΔCWM＝Σ(pR−pC)z*＋Σ[pR(zR−z*)−pC(zC−z*)] を厳密に分解する。**元サイトごとの植物種比・観測個体の実際の草丈/開花形質がない現段階で分解済みの生態効果は存在しない**。値が欠けた種は恣意的に補完しない。

## v16.5 30m以内の写真機会の正確な定義

- iNaturalist 2018-01-01〜2025-12-31、research-gradeかつ写真を持つ記録、元学名・承認された明示的な別名のみ。実ソースの全件をAPIの1〜200/ページで取得し、1秒に最大約1件、上限1万を超えるbulk要求はしない。
- 最初の空間範囲は**チェコ＋スロベニア2か国を覆う固定地理矩形**。必ずしも国境に沿ったpolygon集計ではない。両国の**源GISは Geofabrik元2023年CZと2024年SIの真のOSM Railway Way**。元PBF SHAを照合し、線路を形作るPolylineへ投影、観察点と**最寄り線路の線形**の距離が**30m以内**か判定。元の節点（OSM node）だけへの距離では、曲線や長い直線辺で大きな見落としがある。
- 必ず公開座標で、geoprivacy private/obscured・明示された栽培個体を除外。座標不確実性（positional_accuracy）が**0超30m以下**の実写真記録だけを厳格件数として数える。位置が確実に30m以内となるのは、線形までの距離＋位置不確実性半径 ≤30m のさらに狭いサブセット。30m以内の記録には、ホーム・駅周辺の公園・道路沿い・立ち入り許可区域外の写真も入り得る。
- 30–100mと100–1000mの同種記録を別に数える。観察者IDや公開座標の原本は成果に含めず、観察者の**ユニーク数・100mの小セル数の整数集計だけ**を公開。各種について月別に撮影日があっても、それだけで**開花開始日**ではない。写真をダウンロードせず、実際の全草丈が測れるかは別段のROI可用性監査が必要。
- 一種×国矩形のAPI研究品質写真観測が**8,000件を超える場合はAPIのbulk制限を尊重して「完全取得不能、30m件数は不明」**と報告し、0件とはしない。世界全域でのresearch-grade写真総数は**全球写真の観察機会数**として別取得可能だが、線路30m以内の世界総数ではない。
- **全球という本来の目標はまだ未達**。この第1波は原著2国の真正線路ベクターを使ったサンプル機会の確認。世界全体の30mをやるには世界鉄道ベクターとiNat大規模DwC/公開exportの空間結合が別途必要。最初に数が極端に少なければ、写真由来の種内モデルより同種共通栽培・文献統合を優先する。

[元source20種選択の再現](../scripts/select_two_original_rail_20_species_for_inat_v164.py)・[v16.4 verdaderamente元作者source CI](../.github/workflows/rail-original-20species-inat-feasibility-v164.yml)・[30m真正写真API×線路Polylineのコード](../scripts/audit_20_original_rail_plants_inat_30m_two_regions_v165.py)・[v16.5元GIS/写真source CI](../.github/workflows/rail-iNat-source-top20-photo30m-SI-CZ-v165.yml)。

研究はPR #1 Draft、main未マージ。写真・元座標は成果リポに複製しない。
