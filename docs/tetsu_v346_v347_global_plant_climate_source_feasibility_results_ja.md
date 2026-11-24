# v346–v347 全球種温度ニッチ × 鉄道フロラ：実行可能性・真正データ検証（2026-11-24）

## 判定（生物学的推定とデータ取得を区別）
ユーザー目標は「欧州用EIVEを南アジア/豪州に外挿せず、植物の全球分布に基づく共通温度ニッチで複数鉄道フロラを比較する」こと。

**データの存在/定義と全球一次ソースのファイル名は真正原本で確認済み。ただし全球植物ニッチCSV全件は取得完了できず、7鉄道研究との実結合効果は未計算。** 決してCSV最初の数KBから種の温度平均や世界的鉄道効果を作らない。

## 確認した全球植物形質
Harris, Brummitt, Cobbold & Reeve (2022) *Bioclimatic envelopes for global plant species*, NERC Environmental Information Data Centre, DOI [10.5285/ca339c86-3674-4030-b891-35326e71141e](https://doi.org/10.5285/ca339c86-3674-4030-b891-35326e71141e)、関連論文 [Harris et al 2022 GEB DOI 10.1111/geb.13564](https://doi.org/10.1111/geb.13564)。原論文ではGBIFの1億超の位置・時刻付植物出現記録（1901–2018）を入力し、記録時点と直前の気候から20万種超のglobal plant climate envelopeを作成。出現記録の地理的偏りに配慮し、EVI及び観測努力の補正を検討。**観測気候からのrealized nicheであって生理的限界/線路のバラスト温度/進化形質ではない。**

真正一次ソース/コード照合：
- 英国政府CKAN公式API `https://ckan.publishing.service.gov.uk/api/3/action/package_show?id=bioclimatic-envelopes-for-global-plant-species` の真正source descriptionとリンクを取得、EIDC公式公開データ情報 `https://data-package.ceh.ac.uk/data/ca339c86-3674-4030-b891-35326e71141e` のHTMLを実読。
- 正本の一次ZIPは `https://data-package.ceh.ac.uk/data/ca339c86-3674-4030-b891-35326e71141e.zip`。EIDC配布HTMLに実際にあるリンク。サポートZIPは `https://data-package.ceh.ac.uk/sd/ca339c86-3674-4030-b891-35326e71141e.zip`（SHA256 `7002c4d6a5e6e3bc6aa0011c94d555d58a973025dde65891232098ae6cff9f2c`；実ZIP7479バイト、3点の説明資料、**CSV無し**）。
- ソースのsupport資料 `metadata.md`, `ro-crate-metadata.json` と正本ZIPの初期25KBのローカルZIPヘッダーから、データファイルの名前を**`data/global_plant_climate_envelopes.csv`**と真正確認。データフォルダのパスは本物のarchive memberであり未確認の直接HTTPファイル公開URLではない。
- 原著の公式 Julia 元解析コード [boydorr/ClimatePref.jl examples/23_Output_initial.jl](https://github.com/boydorr/ClimatePref.jl/blob/main/examples/23_Output_initial.jl) に、出力ヘッダー `species`, `numrecords`, `tmean` と、気温 `tmean` が**ケルビン**、著者のEVI記録補正ロジックが明記されている。摂氏変換は `K−273.15`。コードの入力テーブルと未解決同名の扱いは実原CSV取得後に厳密再検査する。
- EIDC serverは一次ZIPに対する `Range: bytes=0-4095` と `Range: bytes=1048576-1052671` を無視し、HTTP200の動的なchunked ZIPを**毎回先頭から生成**（`Accept-Ranges` / `Content-Range`なし）。全件取得の最初の試行は受信途中に `http.client.IncompleteRead: IncompleteRead(391451 bytes read)`（約3分で切断）。この一次ZIPの完全なSHA256/CSVレコード数/生物学的効果は確認していない。
- 直接HTTPファイルの推定候補（公式 `catalogue.ceh.ac.uk/datastore/eidchub/uuid/filename`、`data-package.ceh.ac.uk/data/uuid/filename` と `/data/filename`）4経路＋indexは全て404。推測で存在しないURLから「ファイルダウンロード成功」と呼んでいない。
- 公式JSONとzip provenance: [v346a Actions #45309271134](https://github.com/zuizui0223/tetsu/actions/runs/45309271134) SUCCESS／[v346b #45309453248](https://github.com/zuizui0223/tetsu/actions/runs/45309453248) SUCCESS／[v346e #45311327352](https://github.com/zuizui0223/tetsu/actions/runs/45311327352) SUCCESS。HTTP Range [v346d #45311126545](https://github.com/zuizui0223/tetsu/actions/runs/45311126545) SUCCESS（取得不可を正しく報告）、公式CSV候補 [v346f #45311646529](https://github.com/zuizui0223/tetsu/actions/runs/45311646529) SUCCESS（各404）。v346c [#45309561328](https://github.com/zuizui0223/tetsu/actions/runs/45309561328) はデータ取得途中失敗であり**source climate calculations 0**。
- これらの成功は source access/metadata QA の成功であり、生物学的効果推定の成功ではない。

## 既存の本物鉄道フロラ：最大3大陸／7一次研究
`research/tetsu-v315-thermal-affinity-reference-20261122` および固定 `data/manifests/v314_matched_foreign_rail_road_sources.json` にすでに原著資料とSHA256があり、下記のリアル元sourceに同じ全球 `T_species` を当てることが可能：

| 一次研究/出版 | 本源鉄道標本支持 | 位置づけ |
|---|---|---|
| チェコCZ39 | 39駅現代フロラ、一出版研究 | Europe |
| イタリアIT28 | 28駅、1出版研究 | Europe |
| HRBA15 | クロアチア13 + ボスニア2、**同じ一出版研究** | Europe |
| スロベニアSI89 | 89線路内plant sites（ballast BS+BT）／15 sections、1出版研究 | Europe |
| スロバキアSK30 | 30駅の原著種別頻度（**地点ごとの完全種組成は未公開**） | Europe |
| カシミールIN | 2014年の同じ16駅 Track+Embankment vs Interior1+2（また31原著区間の軌道vs内陸区画あり） | Asia |
| 西豪州WA | 2008年鉄道用地の種出現、真の2,500m² 120調査区／本sourceから3,412種類×plotラベル、**稼働線路内とは限らない** | Oceania |

**全球**の植物種温度ニッチ情報を使うことで各大陸の種を同じ℃指標で比較できても、原著鉄道群集は大陸**3**（研究本数は最多7）、南北米・アフリカの独立鉄道植物調査を含まない。駅数・植物種数・同一論文の地点数を独立研究数に水増ししない。

### 本当に比較可能な対照
カシミールは2014年同駅の近隣内陸2区画が存在する（約150mまで。鉄道影響のない自然区画とは確証されない）。MIRENはCZ/IT/INの**道路**、SI/HRBAは都市植樹桝、SKは都市植生。どれも「同一地点・同面積・同時期・独立自然生息地の鉄道対照」ではないので、推定対象が揃わない7国の差を一本のメタ効果量としない。

## v347 Global temperature CWM 実装済み（数値結果はsource未入手）
`scripts/analyze_global_plant_temperature_railway_v347.py` に：
- 元EIDC CSVのexact source `species / numrecords / tmean` 名と単位、原著全CSV100,000行以上・種名二名法90,000種以上のhard source gate。
- 記録数20以上主、100以上感度、未知値を0扱いせず、原出現種≥15/既知種≥15/種数被覆≥50%/頻度量被覆≥50%を両群集に適用。
- 鉄道7研究の元種別出現頻度で加重した全球種の平均気温、各種等重み感度、Cynodon単独と暫定C4候補10種除外。
- 原著comparatorsの設計ラベルを区別した6研究7比較（CZ/IT/INのMIREN roads, SI/HRBA/SK都市資料、IN同駅内陸）と真の同駅16 pairs equal-station contrast。
- RAW original CSV / plant-by-plot matrices / GPS をrepoに再配布せずaggregatesだけ出す設計。
- `tests/test_global_plant_temperature_railway_v347.py` に合成数学・原著欠損/閾値の検査。**sourceを使う原生態学推定はまだ未実行**。

## 科学的発展の最小成立条件
一次EIDC original CSVのダウンロードが研究者環境（学術ネットワークやEIDC公式ページ）で完了すれば、SHA256を記録し上記parserで同じ7 source cohortへ結合できる。Global trait分布は全球だが真の鉄道選別・原因・ climate-limit release を語るには複数大陸で各々 local matched independent nonrail original botanical data と site macroclimate が不可欠。

**現時点の科学的 verdict: 全球の温度形質データは実在して定義も明瞭、7鉄道研究への統合コードはあるが、世界範囲の暖地性効果を新しく定量推定した事実はない。**

主張すべきものは **global-trait, three-continent feasibility design**。観察済みの真正温度差は以前のカシミール/欧州EIVE等の部分比較に限る。研究ブランチ `research/tetsu-v346-world-plant-climate-envelope-20261124`、[Draft PR #55](https://github.com/zuizui0223/tetsu/pull/55)、main未変更。
