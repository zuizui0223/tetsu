# v318 実原資料を読む前に固定：温帯オーストラリアで鉄道盛土と他攪乱環境の植生がどう違うか

2026-11-22 JST、v317のMell2024メタ原表（鉄道3効果、2原著）の抽出後、次段の**独立現場データの探索的再解析**として指定。v315/v316では同じKashmir回廊の植物温度親和性vs標高は揃わない結果だった。v318の目的は、真に別大陸の元植物×地点資料の個票を入手し、古典的な「人工生息地でのフィルター」を同じ研究内部の異なる人工環境間で診断すること。

## 原著と固定ソース
Everingham S.E., Hemmings F. & Moles A.T. (2019), *Austral Ecology* 44:702–712, DOI https://doi.org/10.1111/aec.12719。出版社公開の補足**Appendix S4: Raw data of species collected during fieldwork**、ファイル名 `aec12719-sup-0002-AppendixS4.xlsx`、約191.2kB。本研究はNSW5地域の道路芝生・道路肩と、sclerophyll森林の一地域で鉄道盛土・防火道・森林自然参照を含む攪乱度比較。原報告は在来種136種、合計調査種数とカバーデータ、造林・雑草の不在偽陽性に注意。

原元表は正式出版社の公開直接Supplementary URLのみから取得する。Wileyの公開パスが403、登録が必要な場合は未取得とし、権限を迂回しない。取得後はXLSX元バイト・SHA256・全シートとヘッダを記録し、元raw種×地点の行列・座標はリポへ再配布しない。

## 判定と厳格な文献再現
- Primary: sclerophyll forest地域の元**railway embankment**と**fire trail**の在来植物（native）の元植物種数とcover。Mell2024真正メタの`Everingham_2019_3` Hedges d **−1.63326905465261**（元24 rail vs24 fire trail）、同`_6` cover d **−0.516972129076364**（元24+24）を固定ソース確認用基準とする。元著者が添付 raw data を生物学的な本数24+24で示し、同じ原単位・被度測定を再構成できる場合だけdを再計算。Mellの元値は図を数値化した別種の推定値で、rawと完全一致しない可能性は前提、誤差を明示。
- Secondary: 元研究のその他の交通生息地（路面芝生、道路肩、防火道、鉄道盛土、森林）ごとの**単位面積・native richness、native cover、総植物種数、source各分類群の地点占有数**。同じ広域内で分類が確定する範囲だけ。原生植生の森林が同じsclerophyllの参照であれば元面積を合わせた比較を別に示す。
- Source site IDとhabitatカテゴリが実際にある場合、元1区画を生物学的独立候補とし、同じ調査地内でspatial cluster/pseudorepがあるなら地点ラベルを明示して再標本化する。単に種×地点セルを独立n扱いしない。空間単位がsourceで未同定なら有意差や95%CIの主張をしない。
- Hypothesis: 同じ調査地域の低攪乱参照に対し、鉄道でのnative種数・カバーデータは低い（元出版結果を知っているので**再発見ではなくsource replication**）。真の新規検定は、種の占有頻度を固定したとき鉄道と防火道で**同一種共有が変わるか**、また元trait/元lifeformが実在すれば、特定の機能型への選別が鉄道で過剰か。Source列にNative/Alien指定がない場合、種名から勝手に帰化状態を推測しない。
- 外部国間比較はrailway embankment vs ballast shoulder、fire trail vs natural nonrail、NSW vs Kashmirで設計が違うと明示し、等量のメタ効果に無条件で混ぜない。EIVE/Tichýはオーストラリア植物に数値ほぼつかないことをv315で確認済みなので、欠測種に温度指標を内挿して寒冷温暖の結論を作らない。

## Source gate
- 公開サイトの補足ファイル名とXLSX署名・byte/hash・original sheet名を確認。別の原著やHTMLエラーページをsourceとして読まない。
- フィールド正しいサンプル数、植生型、raw種×区画、native status、coverの列を特定できなければ、**生態学的元種共有に進まない**。
- 元sourceをinspect後に、推定量・行列グレインを本文/表で確定し、rawを理解せずdをfitして当てることはしない。
- PRは原著source監査+実行に必要なコード+限定的集計のみ、mainマージ・外部メール送信・公開でないsource突破はしない。
