# v338（v337を見た後の事後検証）：種ごとのabundance分布・zoneを保持しても駅間内陸で暖地性種が優勢なのか

2026-11-23。元source原著 Dryad `Matrix_plot.xlsx` version99735（SHA256 `529f3b97521933018e12623b61f531f012c67241cae2602a6675dbd534192474`）、独立EIVE1.0温度選好SM08（MD5 `0f0ce19cd2a781eeed027d8b0acadbcb`）。Rashid et al. (2021)元著者受理原稿 [University of Antwerp source](https://repository.uantwerpen.be/docstore/d%3Airua%3A7552) 方法2.2による12個1m²方形区のspecies abundanceという表現を採用。ただし元XLSXセル値の合計/平均/換算の物理単位は断定しない。

## 既存v337の観測
2014年の緩いEIVE数値3種/30%（体系>=1）・地理と標高を合わせた駅間-駅14ペアで、`species presence weighted CWM`差+0.288560、`source abundance weighted CWM`差−0.037710、差（すなわち量的群集内形質対応の駅間-駅比較）**−0.326270**。v337では同じ野外plot内の種のabundance値を種間で交換し、探索的2側0.033。この帰無は種自体の生来の量的多さや生活史・成長形を保持しないため、種固有の多さが原因でも棄却され得る。**その制約を検証するのがv338**。

## 結果を確認する前に固定する再置換
1. 元v332による**2014年14ペア（G1 15km/75m、relaxed EIVE>=1/各区3種30%）**を主検査、厳格2014年5ペア（EIVE>=3/5種50%）・種数/形質被覆G3の7ペア、2017年EIVE>=3/数値1種以上12ペアを感度として保持。全マッチは元植物abundanceの結果を使わず固定。
2. **帰無S（species only）**：同じ調査年のKashmir元31site×4zoneの全真正観測で、speciesごとの全source正abundance値の**多重集合を固定**し、そのspeciesが実際に記録された元plot間でのみ値を入れ替える。各speciesの元sourceでの総量的値、全sourceでの出現種、各plotの種リスト、年・駅/駅間、ペア・EIVE種値を完全固定する。**元plotごとの総abundanceは固定されない**（再規格化CWMを必ず計算）。
3. **帰無SZ（species×habitat-zone）**：同じ調査年・同じ元`zone`（Track, Embankment, Interior Plot1, Interior Plot2）・同じspeciesの正のabundance値だけを、同じspeciesがそのzoneに出現した他の元siteへシャッフル。つまり種固有の一般量的優勢性、**Track/内陸というzone固有のabundance分布**も元源から保つ。地理・種出現頻度は固定するので、駅間/駅だけの局所的なabundance–温度選好の残差を検査。元zoneラベルをランダム割り付けた実験と呼ばない。
4. どちらの帰無も元sourceの全site・zoneから値を取得し、結果を求めた14ペアの4zoneだけを置換母集団にしない。こうして元sourceの種別abundance分布を使用する。**源に出現1回のspecies、同zoneで1回しかないspeciesは置換できない**。動かせた元source species群/zone群/セル数を監査し、帰無検定力が低ければ過大解釈を止める。
5. 各置換1999回、固定乱数338。各source site×zone×yearに元種名の異なるabundanceが割当されても、そのspeciesの数値値の多重集合は保持。評価する統計は前回v337の**abundance重みとpresence等重みの4群`BT−BI−ST＋SI`差**で、元v337と数値が完全一致することを停止条件とする。元source列に値がなければ欠測を0で埋めない。元site・元speciesのGPS/個票は成果物に出さない。
6. 報告はobserved、999→**1999置換**の技術的2.5/97.5百分位点、帰無平均、`observed−null_mean`の2側分率（帰無平均を中心とした中心化絶対差）、候補負方向の片側分率、2つの帰無S/SZを並べる。帰無平均が0とは限らず、元v337のplot内シャッフルと異なり**speciesの元abundance–温度指標相関を保持**する。
7. 2014年のSZ帰無も観測−帰無が負で明瞭なら、単純な種固有の高abundanceとTrack/内陸での種固有傾向だけでは説明しきれない駅/駅間局所的な種×abundance差の候補。超えなければ、v337の見かけの「内陸で暖地性種が優勢」は元からよくabundantな種・元生息区画別の普遍量的差で再現できる可能性が高く、温度特異的な新規機構を主張しない。

## 根本的な限界
種ごとの元source abundance値の量的単位（1m²方形区の平均/合計など）は未解決。元31地点は一つの鉄道路線の空間依存サイトで、2014/2017年は独立研究ではない。Species abundanceは微気候の測定ではない。アジアの植物に欧州EIVE指標を適用している。Abundanceのsite間置換は元plotの総量を固定せず、時系列・局所土壌・管理等の環境交絡を解決しない。多年生か一年生か、在来か外来かの分類も未結合。帰無は**結果を見た後の追加検証**で、多重検定補正済みの確証pではない。
