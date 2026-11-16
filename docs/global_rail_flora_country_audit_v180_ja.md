# 鉄 v18.0：鉄道は地域種プールを上書きするか — 国別植物相一次資料の採録と検定設計

確認日: 2026-11-15。対象ブランチ: research/railway-seed-dispersal-v01。これは公開一次資料を追跡するメタデータ台帳であり、植物の全行を転記済み、効果量計算済み、あるいは統一した全球実験という意味ではない。

## 検定可能性の監査

公開原著に同定された植物名表が存在することを確認できた12か国: CZ, SI, HR, BA, IT, AT, SK, IN, PL, FR, TR, DE。独立した原著研究数や均質なバラスト調査数ではない。

そのうち、原著のsite-by-taxon名表を抽出する候補は8か国・7独立原著系列: CZ (39駅/761), SI (元92地点、既存QC89), HR+BA (13+2駅/359、同一論文), IT (Parma28駅/438), AT (東チロル6駅/398), SK (Košice 207とSpišská Nová Ves 204の歴史・現代別出現名表), IN (Dryad 31 transect、2調査年)。これらは同じ生息地/サンプリング面積ではなく、出現の欠如を比較可能な真のabsenceとみなす前に各原著の努力・面積・場所を再監査する。

追加の日本、英国、オランダ、ウクライナ、スイス、ブラジル、豪州等にも鉄道関連資料はあるが、公開表が部分的・historical selection・生息地相違・全文未取得などであり、上記12か国のフルネーム表確認数へ安易に足さない。

## 一次資料別台帳（A=サイト×種抽出候補、B=種リスト/頻度、C=追加原表要、EXCLUDE=別生息地）

|国|原著・対象|名表の範囲|階級|
|---|---|---|---|
|CZ|[Kutlvasr et al. (2024) 39 stations](https://doi.org/10.23855/preslia.2024.247)|761; station x taxon|A|
|SI|[Sabic et al. (2026) Slovenia railway microhabitats](https://doi.org/10.1007/s44473-026-00155-z)|312; site x taxon|A|
|HR;BA|[Jasprica et al. (2017) Mediterranean railway stations](https://doi.org/10.20302/NC.2017.26.21)|359; 15 station x taxon|A|
|IT|[Adorni and Ghillani (2020) Parma](https://www.museocivico.rovereto.tn.it/UploadDocs/21188_Art08_Adorni_C_agg_min.pdf)|438; 28 station x taxon|A|
|AT|[Stoehr and Brandes (2014) East Tyrol](https://www.zobodat.at/pdf/CAR_204_124_0631-0670.pdf)|398; 6 station x taxon|A|
|SK|[Majekova et al. (2016) East Slovakia stations](https://www.upjs.sk/app/uploads/sites/15/2022/10/173-188_Majekova_et_al-upr.pdf)|207 Košice; 204 Spišská Nová Ves; 2 station x taxon by historical or recent|A|
|IN|[Rashid et al. (2021) Kashmir railway](https://doi.org/10.5061/dryad.0p2ngf20h)|unverified; 31 transects x taxon in Dryad XLSX|A|
|IT|[Filibeck et al. (2012) Lazio embankment](https://hrcak.srce.hr/file/130139)|287; pooled taxon list|B|
|PL|[Warcholinska and Suwara-Szmigielska (2008) Lask](https://czasopisma.uni.lodz.pl/biologica/article/view/18835)|369; town-level pooled taxon list|B|
|PL|[Warcholinska and Suwara-Szmigielska (2009) Pabianice](https://czasopisma.uni.lodz.pl/biologica/article/view/7479)|382; town-level pooled taxon list|B|
|FR|[Hoff Strasbourg rail habitat Annex 12](https://www.societe-botanique-alsace.org/wp-content/uploads/2023/12/Annexe-a-la-flore-et-a-la-vegetation-de-Strasbourg.pdf)|unverified; rail habitat taxon names|B|
|TR|[Altay et al. (2015) Haydarpasa-Gebze line](https://www.jeb.co.in/journal_issues/201501_jan15_supp/paper_21.pdf)|194; published species appendices|B|
|DE|[Winkler (2026) Stuttgart railway track](https://www.mdpi.com/2673-4133/7/2/52)|68; taxon list with frequencies|B|
|AT|[Leonhartsberger and Poeltl (2019) Graz stations](https://www.zobodat.at/pdf/JoanBot_16_0031-0046.pdf)|177 vascular; taxon list with restricted station flags|B|
|IT;AT|[Bertolli et al. (2023) Verona-Innsbruck atlas](https://doi.org/10.53135/ANNMUSCIVROV20233903)|1177; reported 107 stations, 21301 archive records|C|
|NL|[Koster (1987) Netherlands railway flora](https://edepot.wur.nl/390615)|~1032; 292p historical regional survey|C|
|UA;PL|[Denisow et al. (2017) Lviv and Lublin](https://doi.org/10.5586/aa.1727)|509 all taxa combined; 70 invasive named; invasive subset with per-city frequency|C|
|UA|[Dziuba et al. (2022) Kyiv](https://doi.org/10.1007/s11756-021-00961-0)|unverified; 379 relevés aggregated in syntaxa tables|C|
|CH|[Tinner and Schumacher (2004) NE Switzerland](https://doi.org/10.1007/s00035-004-0703-6)|509; 4 stations complete list claimed in abstract|C|
|BR|[Mendes (2023) Sao Paulo railway flora](https://repositorio.ufscar.br/items/accded14-c5d2-4269-b357-a5874ad0d384)|unverified; thesis on four active sections and reference habitats|C|
|JP|[Abe et al. (2012) railway weed survey](https://www.jstage.jst.go.jp/article/weed/57/4/57_159/_pdf)|62 taxa reported; top 20 published; partial ranking from 109 sections|C|
|GB|[Sargent (1984) national railway land flora](https://hrcak.srce.hr/file/130139)|1632 secondarily reported; original checklist not retrieved|C|
|AU|[Western Australia Donnybrook-Katanning rail reserve](https://library.dbca.wa.gov.au/FullTextFiles/021122.pdf)|154 selected dominants; natural reserve floristic list|EXCLUDE|

## 研究仮説の変更

問い: Can repeated anthropogenic habitats override regional species pools? **バラスト/整備規格がどの国でも同一で選択圧も同一、は前提にしない**。盛土、バラスト間、バラスト肩、軌道脇、駅構内、貨物ヤード、休止線、トラム、鉄道保有林は別のsampling strata。生息地を混ぜれば同質化の過大評価・過小評価両方が起きる。

帰無集合は各地域で独立に記録された地域チェックリスト/標準植生調査（GIFTは地域粗度に注意）から、標本サイズ・分類粒度・対象時期・元データ検出努力を揃えて再標本化する。鉄道どうしの種名群集距離と、鉄道どうしの形質距離を別に計算する。鉄道距離−地域プール再標本化期待距離を観測量とする。負の差は地域プールの期待よりも均質化、正は相違の拡大。異なる大陸間の種Jaccardが0に飽和するときは属/系統類似を別応答にする。

同じ研究の複数地点、同一地域でのペア比較、HR+BAの2国を独立原著nとして扱わず、原著研究/国/地域に階層効果とブロック再標本化を設定。分類学的収束と機能収束の不一致を主要結論候補にする。外来/自生の層分け・生活史・散布形質・微生息地別感度・分類解像度の感度を必須とする。

機構候補は『自然のアナログ生息地の地域プール内の豊富さ』。乾性草地や砂礫河川敷が類似しているかは各地域の独立植生資料で確かめ、乾燥地域だから豊富だと決めつけない。輸送効果と管理選択はこれとは別推定対象であり、単一の群集均質化係数から列車種子散布/草刈り適応を主張しない。

## 次の実装ゲート

1. A原著7系列について元PDF/公開XLSXとサイト×種の行列を原本照合してCSV化。文献年代の違いと植生方形区のabsence意味を判定。学名をPOWO/COL/GBIFへ追跡し、元名保存。
2. site habitat class, active status, sampling years, season, survey area, station purpose, sample effort, species abundance/incidence, study identity, country, ecoregion, within-region independent pool の品質列を作る。
3. 可比な構造ごとに実残存数を先に集計。十分な気候・生物地理帯と独立原著が揃わなければ全球因果論をせず、群集集積の観測的証拠マップと機能収束の探索検定に下げる。
4. 地域種プール標準化で分類学的・機能的距離と不確実性、原著leave-one-out、地理層別、歴史/現代差、自然アナログ交互作用を評価する。まだ数値係数・有意性は存在しない。

関連する過去のtetsu原著検証はSI/CZで完成しているが、このv18.0は複数国の新規メタデータ採録のみ。mainの変更・PRのマージはしない。

CSV台帳: data/literature/global_rail_flora_source_registry_v180.csv
