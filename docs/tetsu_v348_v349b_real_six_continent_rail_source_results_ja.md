# v348–v349b：鉄道植物を**自前で**全球抽出できるか、6大陸の真正API検証（2026-11-24）

## 既存7研究が選ばれた理由と修正
以前のCZ/IT/HRBA/SI/SK/Kashmir/Australiaは「完全または部分的な種組成表・原著資料を見つけて再解析可能だった」研究で、**systematic global sampling frameではない**。その地域や符号を全世界へ外挿するのは無理がある。自分で地理的採択基準を定義し、公開鉄道路線と公開植物出現記録を全球で重ねられるか検証した。世界各地の鉄道に植物が発生しやすいか、温暖種が多いかという因果結果を直接示す段階ではない。

## どのようにsource抽出したか
1. 年**2025**に固定し、種出現結果を見ずに、東京／プラハ／ケープタウン／シドニー／サンフランシスコ／サンパウロを1大陸1都市ずつ予め固定。中心±0.012度の検索枠。都市部の鉄道ハブという便宜的な**技術パイロット**で、全球のrail length weighted確率標本ではない。
2. 植物を [iNaturalist API](https://api.inaturalist.org/v1/observations) のsource `taxon_id=47126`（Plantae）、`quality_grade=research`、`captive=false`、`geo=true`、2025年1月1日〜12月31日・WGS84 bboxで取得。元結果が3,600件以下なら200行/ページで**全ページ**取得、厳密に元総数と合うか照合。公式APIを約1req/secで扱い、sourceオリジナル観測記録・公開ユーザーID・座標・写真をGitHubに転載しない。6都市すべてsource Census完了、**170+536+96+89+133+17 = 1,041 source公開植物観測**。原著APIのsource成功は[CI #45328376103](https://github.com/zuizui0223/tetsu/actions/runs/45328376103)。
3. 元鉄道 geometryは最初OpenStreetMap公式Overpassを使用したが、クエリがタイムアウト/ネットワーク遮断となり、主OSM XML APIも5都市HTTP400、サンフランシスコのみ取得に成功。[v348 Original OSM source + actual iNat spatial check SUCCESS #45328203949](https://github.com/zuizui0223/tetsu/actions/runs/45328203949)。SFの原sourceはnear3件/2種、background38件/28種で事前支持条件を満たさない。他5都市は**鉄道が存在しないのではなくsource鉄道geometry取得不可**。
4. OSM全世界標準化の代替として、[Overture Maps 2026-11 transportation segment](https://docs.overturemaps.org/guides/transportation/)、`subtype=rail`、`class=standard_gauge/narrow_gauge/unknown`、トンネル/橋/工事中/廃止/覆道の各フラグを除いた元source lineworkをGeoParquetから取得。Source検証[v349 CI #45328700652](https://github.com/zuizui0223/tetsu/actions/runs/45328700652)では5都市取得、東京のみ最初のSTAC呼出HTTP429。これをリトライした**v349bでは6都市とも取得成功**した。植物検索窓の外側1.1kmまで線路元geometryを取得し、枠の縁の植物にとって線路が枠外にあるというバッファ誤判定を防止した。
5. GPS位置精度は元レコード `positional_accuracy <= 30m`、非公開・obscured、植栽（captive）、種でない同定、座標欠損を除外。最寄りrail中心線への距離 `d` と原測位誤差 `u` を使い、厳密near `d+u<=80m`、far `d-u>=300m AND d+u<=1000m`。その他はuncertainとして採用しない。支持はnear観測10/種5とfar観測30/種10を両方充足する場合のみ。これは**鉄道近傍の観測点であり線路バラスト内の実植物群集や非検出種を表していない**。

## 真正データ結果
[**v349b 公式Overture global source + iNat元全件・位置精度検証 SUCCESS #45328987737**](https://github.com/zuizui0223/tetsu/actions/runs/45328987737)。元iNat観測と世界鉄道source geometryはRunner内で突合、元の列座標・観測者・個別学名は共有しない。

| 大陸 | 都市 | 全2025元Plantae観測 | near観測／種 | far観測／種 | 二群の事前支持 |
|---|---|---:|---:|---:|---|
| アジア | 東京 | 170 | **0／0** | 82／48 | 不足 |
| 欧州 | プラハ | 536 | **60／41** | **150／88** | **充足** |
| アフリカ | ケープタウン | 96 | **0／0** | 32／21 | 不足 |
| オセアニア | シドニー | 89 | **15／12** | 18／15 | 不足 |
| 北米 | サンフランシスコ | 133 | **3／3** | 35／26 | 不足 |
| 南米 | サンパウロ | 17 | **0／0** | 1／1 | 不足 |
| **合計** | **6都市** | **1,041** | **78観測（重複種あり）** | **318観測（重複種あり）** | **1／6都市** |

元source Overture側は6都市全てで実rail系のジオメトリ候補が存在（候補features計55/208/173/167/60/70、合計733、**独立の733運行路線や群集ではない**）。近傍に観測植物が0という都市は植生の不在ではなく、特定の2025年iNaturalistユーザーがGPS精度などの条件を満たす写真をこの距離帯で記録した機会が足りないだけ。事前の二群支持ゲートはプラハのみ。ここから全球rail warm-affinity平均、差のp値、植物多様性・群集収斂や鉄道植物の実際の定着有無を推定してはいけない。

### 科学的示唆
**既存原著7研究に依存せず、自前の研究者選択で6大陸の鉄道位置×植物出現sourceの抽出が可能**と確認。一方、鉄道近傍のGPS確定記録は世界で不均等かつ少ない。都市中心の単年「記録密度」を種の生態学的優占、生活史戦略の因果効果、暖地性選別の証拠へ読み替えるのは危険。source観測面積と努力を固定した植物調査区ではないため、群集の本当の種数や未知の不在を導出できない。

**有効なのは次の二段階：**
1. 採択都市を恣意的に増やすのでなく、Overtureの約全球network全体をsampling frameとして、気候帯・大陸・都市化・運行形態と長さで層化し、**既定乱数seedで路線/空間ブロックを無作為抽出**。市街/農村双方を含める。記録の有無に基づく選択/欠測を全sampling strataで監査。大きいiNat/GBIF originalデータは公式bulk downloadで全件取得、source通常APIの10kページ制限を回避。
2. 真の `Track vs locally matched nonrail vegetation` と鉄道環境の気温・乾燥度を結ぶ独立の定面積植生調査で、鉄道での暖地植物選別を検証。グローバル気温形質はHarris et al 2022 EIDC `tmean`（全球分布のrealized niche）を検証後に全地区で共通に使う。EIDC元全件CSVは2026-11-24に公式ZIP転送中断で未取得のため、**今は温度差を算出していない**。鉄道近傍観測の時刻は2025年、線路geoデータは2026年snapshotである時点差にも留意。

## 研究用成果
- [固定デザイン](tetsu_v348_self_extracted_global_rail_proximity_frozen_design_ja.md)
- [確率的全球鉄道抽出へ移る正式計画](tetsu_v350_global_stratified_railway_sampling_plan_ja.md)
- [source圧縮済み6都市結果](../data/derived/v349b_self_extracted_six_continent_rail_proximity_source_counts.json)
- [v348 OSM元source / iNatコード](../scripts/pilot_self_extracted_six_continent_rail_flora_v348.py)
- [v349b Overture代替source / iNatコード](../scripts/pilot_Overture_rail_iNat_six_continent_v349b.py)
- [v349 Overture rail geometry品質コード](../scripts/audit_Overture_global_rail_segments_v349.py)
- [source iNat count v348b](../scripts/pilot_six_continent_plant_source_counts_v348b.py)
- CI 5/5 v349 source Overture tests、7/7 v348 rail geolocation / original public quality tests、2/2 v349b source rail extra-buffer tests。動作確認の全件スコープは実行ログ参照。

全て研究Draft [PR #56](https://github.com/zuizui0223/tetsu/pull/56)のみで、main変更なし。元データの著作権・プライバシー・観測者名・元のGPS位置の再公開を回避している。**全球観測検証と、本当に定面積で採った鉄道植物の全球生態学的推論を明確に区別**する。
