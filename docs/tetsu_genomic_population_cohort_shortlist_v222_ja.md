# tetsu v22.2 — 地点・集団メタデータ付き遺伝型資源の一次資料候補（104候補から）

2026-11-16時点の外部公開資料（生物試料のサンプリングプロトコルと元データ公開ページ）に基づく**研究資源の在庫**。線路沿いの実sample数や鉄道に由来する選択の証明ではない。104候補の元リストは [v22.1](../data/derived/Czech_104_broad_rail_context_NCBI_genome_TSA_WGS_inventory_v221.csv)を参照。104は歴史を含むCzechVegの「RAIL文字列文脈」の探索候補であり、種の真の鉄道適応や独立再現を意味しない。

|元104中の種／順位|一次資料で確認したgenotypeの規模|集団位置の元メタデータ|型|集団鉄道比較への判定|
|---|---|---|---|---|
|Conyza canadensis (= Erigeron canadensis) / #14|**640個体、280集団**。VCF 11,501 SNP、19地域|Dryadの公開metadata.csvに Sample_ID/Population_ID/Latitude/Longitude/Human_Footprint/CWD/Biomass/Capitula あり。各列非欠測数は未算出|ddRADseq|**最優先：座標と種内形質を同時に監査**。FLC/EPSPS/CNV直接テストには匿名RAD SNPを参照座位へ再マップできるか別途評価|
|Bromus tectorum / #19|**307遺伝子型**（北米194、ユーラシア105、他8）|論文methodsは試料の収集地点GPSを利用と明記。元個票完全性・線路近接数は未監査|genome re-sequencing + common gardens|**有望：環境・祖先を統制した外挿検証**|
|Arabidopsis thaliana / #59|**1,135** 1001 Genomesアクセッション（公式）|第三者公開metadata.csvは1,135行、そのうち妥当なlat/lon 1,131行・異なる緯度経度829点（2026-11-16のmetadataスクリーニング、1001Genomes公式取得物と未照合）|WGS/VCF|**候補：FLC/FRI多型×線路距離**。座標精度/採集点同定必須|
|Capsella bursa-pastoris / #39|NCBI BioProject PRJNA1212728：**125 SRA実験／125 BioSample**|米国広域の採集と記載。公開BioSampleのGPSと同一実個体対応は未確認|再配列raw reads|**保留：独立個体・GPSの確認が先**。125は実験登録件数であって地点付き有効個体数ではない|
|Fraxinus excelsior / #17|2026年pangenomeの地理的に多様な50標本、既存の1,200超の遺伝型|英国seed sourceなど原資料に存在。個別の線路利用履歴はなし|pangenome + variation|ゲノム資源としては豊富だが、森林木本であり**高速生活環の機構仮説の主系から外す**|

**104の外の機構比較**：Arabidopsis arenosa（山地vs鉄道集団、CO/FLC対立遺伝子、2018 Baduel et al. DOI 10.1371/journal.pgen.1007510）、Ambrosia artemisiifolia（600超の現代/歴史WGS+高地理分解能の公開研究、Battlay et al. 2023 DOI 10.1038/s41467-023-37303-4）。これらを104に「鉄道候補」として無断追加しない。

一次資料：
- Conyza/Erigeron Dryad 2026 https://doi.org/10.5061/dryad.8w9ghx40k （metadata, VCF, phenotype）
- Conyza/Erigeron 2026 https://doi.org/10.1111/mec.70368
- Bromus 2025 https://doi.org/10.1038/s41467-025-64799-9
- Capsella BioProject https://www.ncbi.nlm.nih.gov/bioproject/1212728
- 1001 Genomes https://1001genomes.org/1135/datacenter
- 1135 metadata mirror https://github.com/hagax8/arabidopsis_viz/blob/master/data/dataframe_1001G.csv
- Fraxinus 2026 https://www.nature.com/articles/s41467-026-77194-9
- Arabidopsis arenosa 2018 https://pmc.ncbi.nlm.nih.gov/articles/PMC6049958/
- Ambrosia 2023 https://pmc.ncbi.nlm.nih.gov/articles/PMC10042993/

## 固定して次に確認する支持域

1. Conyza metadata.csv全行を原本取得し、Sample_ID→VCFとLatitude/Longitudeの型・重複・NA・不正座標・年代・測地誤差の有無をaudit。**640/280は公表データ全体の規模でありGPS完備と同義ではない。**
2. 年代に近い鉄道中心線を外部データから取得し、距離≤30/100/500/2000 m帯、近傍道路、都市域、気候、国・地域・祖先集団を整備。30mの線路近接を判定するにはその精度の元座標が必要。座標が粗ければ広帯域だけ評価し、鉄道区画試験とは呼ばない。
3. Rail距離変数はRail適応・強風・交通本数の代理ではない。先にsupport countsとmatchingを報告し、共通支持がなければ「選択圧」の議論を取り下げる。
4. Conyza ddRAD匿名SNPについてEPSPSやFLCに直接マップ可能かの検証を事前ゲートにする。FLC同祖性は主にBrassicaceae内で正確に解釈し、系統外の全種に一つのFLC遺伝子を仮定しない。
