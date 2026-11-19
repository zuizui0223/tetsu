# v29.8 — 真の系統群×植物の出現頻度の同時制約で草丈の空間構造を再検証

2026-11-19 JST。初回source実行の結果を見る前にsource study、同時制約、小群数と標本支持基準を固定する。新しい生物学的な係数は原資料ハッシュ・従来結果の完全再現・本帰無のSOURCE成功のあとに報告する。

## 問い
v29.6の欧州4真正鉄道独立原著では、V.PhyloMaker2元3系統樹のroot-to-tip中央値75%深度で切断した単系統群の内部のみでLEDA草丈・種子重のpaired vectorを置換しても、草丈CWMの地理的分化が4/4研究で正だった。v29.7の独立したMIREN山岳道路と元同一地点の内陸植生でも非ランダムな草丈分化は正。ただし元clade内部の植物種の出現頻度は置換では保存しなかった。頻出種と希少種の高さがclade内で体系的に違うことが既存結果の一部を作っている可能性がある。

今回の一次検証：**真の系統樹cladeと、元植物のsite出現頻度五分位を同時に固定しても、同じ方向の草丈地理的分化が残るか**。

## 元sourceと field 独立性
- 真のRail著者原著4件：チェコ39駅全域、イタリア28駅全域、HR+BA15駅（2国・1出版）、スロベニア89個100m道床地点。計171元siteは171独立研究ではない。元作者publisher PDF/XLSX、Oldenburg LEDA原著canopy heightとseed massのSHAを再照合し、v29.6の元study×3tree=12組の75%clade-only観測−帰無係数を1e-10以内でまず再現する。
- 真の非rail参照：MIREN Mountain ROADS v3 Zenodo DOI 10.5281/zenodo.15749536、元17,899,241 bytes、MD5 b640195cc44f252ccb6ddb28eb6fff13、167,804元植物観察行。元Plot1 Roadsideと同一site/year/100m²のPlot3 Interiorだけを採用。最初の適格年を使用。Primary（各zone元LEDA既知5種以上、元coverage25%以上）は元8 Region/22 Road/272地点、strict（10種50%以上）は5 Region/13 Road/192地点で凍結。v29.7の源3tree×2panel=6組のregion codes/物理Site数/75%clade-only道路際・内陸・その差の元効果を1e-10以内で再現できない場合は停止。
- 元LCVP/WP/TPLはJin & Qian (2022) V.PhyloMaker2原著Git commit 7af3fb5152f691af2e4ec9d5e2e467d1b50505e9の元RDataを正しいGit blob SHAで認証。元tree上のtipだけを使用し、欠落tipの推測graftingはしない。tipが元treeにあることは当該植物種のDNA配列直接測定とは違う。
- 元traitは両研究とも欧州鉄道3原著622植物種から事前固定したz(log canopy height, log seed mass)軸。MIRENの両形質が1096/5739 taxa（19.1%）にしか付かない偏りに注意。元site×taxa、species occupancy、元Roadside/Interior対応、原著地理群は動かさない。系統樹の3命名版は独立研究3件ではない。

## 同時制約の定義
1. 系統樹の実original root-to-tip枝長中央値75%の深さにおける実モノフィレティックcladeで全speciesを一意に分割。
2. 同じ研究の元植物種を元site出現頻度でソートし、タイブレーク順を固定して5つの等種数五分位binへ分割。MIRENではRoadside/Interior双方の元plot×species出現を合計して1つのspecies頻度を定義。
3. 本物のcladeと本物の頻度五分位の集合の交差を、互いに排他的なtrait-label交換群とする。1種のみの交差小群では元高さ＋種子重を固定。同一vectorの2つの成分を切り離さない。同じtaxonが道路と内陸両方にある場合は両habitatに同じ置換を適用。
4. original site×species presence、source group、元種の出現頻度、元種の共同出現、元植物相のsite richnessを全部変えない。

## 帰無実行前のsource-identifiability
元効果の符号を見ずに条件を固定：元交差小群の2種以上を含む群が最低5群、元site×species発生の25%以上がこの群で置換可能、Railの元著者地理区分が2群以上で各25%以上置換可能。MIRENではRoadsideとInterior両方で25%以上が置換可能かつ実Road2本以上で両zoneがそれぞれ25%以上交換可能。条件未達は NOT_IDENTIFIABLE_PHYLO_OCCUPANCY_JOINT_LABEL_NULL。計算不能を生態学的な効果ゼロと見なさない。解析後に有意性を求めて同時制約を緩めることは禁止。

各source適格研究について元2seed×199=398技術的なlabel assignment。元Rail4原著のうち3本以上source適格なら出版等重み記述平均とleave-one-publication-out、MIREN元Region3以上ならregion等重み記述平均とleave-one-region-out（ただし同じMIRENプログラム内の地域）。技術順位は野外field p値ではない。3本の系統樹は感度分析であり、12独立野外研究ではない。

## 事前の結果解釈
- Railと内陸の両方が正：真系統と種の頻度の双方を固定しても草丈の地理的構造が残り、鉄道専用適応という排他的な説明は弱い。
- Railだけ正、内陸ほぼゼロ：鉄道と内陸の観察分布は異なるが、調査粒度と系統・環境差があり、直ちに列車特有の因果選択とはならない。
- 両方のSOURCE支持不足：元系統＋頻度を同時に固定すると交換可能な植物taxaが足りないというデータ診断。以前のclade-only正の観察が取り消されたわけではない。
- Roadside−Interiorの符号が地域により混在：平均が負でも全球一般法則やメタ予測区間の非ゼロとはしない。

この検定では稼働中rail ballastと同じ調査面積・時期・元環境のnonrail open/gravel対照は未取得。列車種子散布、草刈り/風圧や現地形質の遺伝的適応はまだ測っていない。

## 実装
- scripts/phylogeny_occupancy_joint_label_null_v298.py
- scripts/analyze_4original_rail_joint_phylo75_occupancy_null_v298.py
- scripts/analyze_MIREN_roads_inland_joint_phylo75_occupancy_null_v298.py
- tests/test_rail_road_joint_phylo_occ_source_v298.py
- .github/workflows/tetsu-v298-original-RAIL-and-MIREN-road-inland-joint-clade-occupancy.yml

元plant×site row、元site GPS、元LEDA植物別値、原treeのtip個票は公開しない。出力は元出版および元Region単位のSOURCE支持数と派生効果だけ。