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
