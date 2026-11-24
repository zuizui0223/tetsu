# v377 全球鉄道km×植物記録種数：観測者集中格子の除外による頑健性（2026-11-24、結果確認前に固定）

**v371正規化済み原世界100/200km M0+GERS rail-kmの元植物観測種数の大陸外予測改善と、v375気候帯条件付き偽鉄道地図との差を見てからの事後探索的反証**であり、事前登録された新仮説の確証的検定ではない。

- 元植物source：iNaturalist 2026-10-27 original open photo-linked source、全288,721,938行、2025年source精度・等級・canonical種条件後7,827,936真正植物観測。100/200km whole-landscape species richness/adjacent Jaccard元cell集計[CI #45375453375](https://github.com/zuizui0223/tetsu/actions/runs/45375453375)、source obs SHA256 `1c48a9646f3bf6d56ffdfdc91a1a176c173c98014830ad48fa162b2803c3b279`。
- 全世界GERS元source Overture2026-11-18.0全線路2,237,043.481km、100/200km実線路長地理分割[CI #45373985228](https://github.com/zuizui0223/tetsu/actions/runs/45373985228)。元WorldClim2.1 BIO1・Natural Earth検証済みsource。元非鉄道陸地も鉄道km=0として維持。
- **処理の順序**：v371の元観測>=50、distinct observers>=5、species>=20、landcentre・WorldClim macroclimateをすべて通した一貫した100km/200km植物格子集合を取得。ここで**植物種数（目的値）、植物相Jaccard、鉄道延長、源気候を見ず**、各大陸の元観測者数で降順ソートして上位 **5%** を除外（端数はceil、同数なら元cell indexを昇順にして決定）。さらに上位**10%**除外を第二感度。比較基準0%も含む。これらの閾値は元結果を確認する前に凍結。
- source隣接Jaccardは除外したcellに接するsource pairをそのまま外し、残った元隣接対のみ、原植物・鉄道源の数値を**一切再計算/補完しない**。種数とJaccardのM0（source記録件数/観測者数/BIO1/緯度、Jaccardはspecies richnessと差・方位も）とM1（+元railkm）を全く同じ元空間正規化ridge λ=.10、6大陸leave-one-entire-continent-out、1,000km spatial block equal weightedで比較。元v371の差は**0%で厳密再現されないと停止**。
- 主判定は source植物種数の0/5/10%における「平均M0−M1 heldout ΔMSE（大陸等重み）」、「改善した大陸数」、「各大陸外部foldごとの符号/推定対象数」。年平均気温とsource植物種の暖地性traitとは異なる。持続的な予測情報があっても鉄道自体の因果効果・列車種子運搬・植物の定着や真の群集豊富さは推定しない。
- 大陸によって残ったcell <100、隣接pair <50、訓練対<300等でfoldが減ったら、その欠損は発見の支持には使わず報告する。100/200kmは同一母集団の異なる空間粗視化であり独立生態研究に数えない。観測者重複・植物撮影アクセス・植栽/観測精度の生態バイアスはsourceフィールドに限界がある。ベースモデルに観測件数/観測者数は既に入っているので、この除外は局所的な高努力支持による感度に限る。
- 生成物は元公衆観測者ID、稀少植物・元species×緯経度、実GERS線路UUID/生路線点を公表せず、集計値・元source fingerprint・研究の限界のみ。
