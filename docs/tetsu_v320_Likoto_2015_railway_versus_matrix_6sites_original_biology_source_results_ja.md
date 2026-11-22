# v320 フランス・ベルギー *En marge* 2015 真正420頁報告：鉄道6地点の近接植物群集データを使えるか

2026-11-22 JST。原研究 *En marge... Paysage et biodiversité des accotements des grandes infrastructures de transport de l'eurométropole Lille–Kortrijk–Tournai*, décembre 2015。公的研究プログラムの[出版者掲載PDF](https://ittecop.fr/images/docman-files/docman-files/rapport-final/2008/RF_2008_En_marge.pdf)を新規に再取得。**96,783,781 bytes／420頁、SHA256 `25a1a38512156bf1264b697efd78984a1abf9217f1da21e6ef2c5623da509405`**。一部のLikotoサイトに置かれた別PDF `27a.pdf` は管理編25頁だけであり、全原著の植物調査区データと取り違えない。

## 結論

**この研究には、他国で足りなかった「同じ場所でのインフラ敷 vs 周辺のマトリクス」という野外植物区画設計が存在する。** 元報告の方法欄は24地点を選び、春・夏・秋にサンプリングを行い、**可能な場合は周辺空間と比較**すると明記。元地図の凡例には `placettes du botaniste (infrastructure / matrice)` があり、元植物調査担当はConservatoire botanique national de Bailleul。これは「単に周辺の山地道路を別国から借りた」v315の擬似対照より設計として有望である。しかし、**必ず24地点全てに真正な同面積・同年の近隣植物区画があるとは確認できない**。

原報告の目次は**24のインフラ地点**。鉄道については、TGV3＋通常鉄道3＝計**6候補**。元報告の現地記述章を420ページから場所名と構成で監査し、各鉄道候補における植物相の「La flore」節を検出したが、**このソースPDFだけで6地点×鉄道/非鉄道×各種の完全出現行列、植物被覆・調査面積の等価性、真正対の種数はまだ検証できない**。したがって鉄道の種数・温度親和性・Jaccard等の**新しい実測差は0件承認**とする。0認証は、現地に本当の対照区が0個あったという意味ではない。

| インフラ候補（source site ID） | 原報告書の鉄道線 | 主体の型 | 原PDFの現地記述頁 | 原文「La flore」の頁 | ボタニスト地点図のインフラ/マトリクス凡例を文字層で検出 |
|---|---|---|---:|---:|---|
| 17 | Paris–Bruxelles (Esplechin) | 高速鉄道 | 263–268 | 268 | 263 |
| 18 | Paris–Lille (Fretin) | 高速鉄道 | 269–275 | 274 | 269、275 |
| 19 | Lille–Londres (La Chapelle d'Armentières) | 高速鉄道 | 276–283 | 282 | **未検出（区画不存在とは解釈しない）** |
| 20 | Lille–Douai (Phalempin) | 通常鉄道 | 284–291 | 290 | 285 |
| 21 | Lille–Kortrijk (Lauwe) | 通常鉄道 | 292–299 | 292 | 293 |
| 24 | Lille–Dunkerque (Lompret) | 通常鉄道 | 324–344 | 324、334、336–338、342 | 325 |

※頁番号は表紙から数えた**PDFの1-based**。原報告は約420頁の長文で、特にsite24の範囲は最後の現地章・総括の境界近く。文字層に名称・凡例を検出しただけで、地図上の対照**区画の座標や種の完全記録**が得られたわけではない。「La flore」の1～数パラグラフは元現地での**生態学的叙述**で、無条件に全種名リストへ変換できない。

## 方法の厳格な意味

1. 元24地点は1研究の現地プログラムに従属し、24独立研究、6独立国、6気候帯ではない。原報告は水路・道路・高速道路・高速鉄道・通常鉄道を扱う。インフラ種類と周囲の植生や地形の差が交絡するので**道路・水路の調査地点を自然植生の同一地点対照とすることは禁止**。
2. 元論文の「placettes du botaniste (infrastructure / matrice)」は試験設計の**存在を示す**が、対応する元学名表が報告書に完全収録されることを保証しない。サイト17–21/24各区画の`PlotID / original infra or matrix / survey season / area / all species and cover`が必要。
3. 生の原資料を入手できたら、同一年/季節/面積/斜面/線路管理の明示IDで**鉄道側−マトリクス側の種数**、**在来外来組成**、元source国内に妥当な種の**温度親和性指標**を対応で評価する。同じ原論文内でも鉄道の現地区画が不成立なら`not_identified`にする。種の場当たりな出現を作らない。
4. 既存v315の欧州参照差（特にHRBA +1.33、チェコ+1.11）は別地の道路/植樹枡を対照としていた。本研究が有用なら、それらの差を「local matched pairsによる検定」に置き換えても、鉄道内の温暖種親和性があるか**より適切に反証**できる。
5. 単一研究の鉄道6候補では荷物列車頻度・交通OD・種子活力は原取得されておらず、輸送と環境選別の分解はまだ不可能。高速/通常鉄道は元インフラタイプの名称だけであり、旅客/貨物本数の直接処置ではない。

## 今回までの原データ取得と不足を明確化

- **原報告取得は成功。** [GitHub Actions #45013059824](https://github.com/zuizui0223/tetsu/actions/runs/45013059824) にて元原PDF SHA/420ページ/24項目/鉄道6現地章/La flore章/5地点の植物区画地図凡例を再検証。元報告PDFや生地点GPS・種×区画画像をリポへ再配布していない。[QA集計アーティファクト #14750288502](https://github.com/zuizui0223/tetsu/actions/runs/45013059824/artifacts/14750288502)。
- **入手できていないもの：** 元6地点のインフラ/マトリクス調査の**植物種×区画（在来/外来および被覆を含む）**、各調査季節（元3回）、PlotIDと原区画面積・地理。元報告の生物調査を担当した **Conservatoire botanique national de Bailleul / CBNBL** の生調査票が存在するか、責任著者・生物学担当者へ原資料と再利用条件を確認するのが筋。**著者への問い合わせはまだ送っていない**。
- **元植物に関する叙述**は各現地章に掲載されるが、全元種名・不在まで含む完全な統一フォーマットであることを確証していない。文字列上の属種らしい組を種出現事実として推定しない。
- **他独立候補**：Prober 1996 NSWの鉄道沿いWhite-box woodlandの広域種組成; China–Laos 2026同駅近接植物区画; Czech Winkler et al. 2026 five habitat at four localities。いずれも元science plot×speciesを再分析するには元原個票/正しい参照対を要する。

## 研究判断

**本プロジェクトで最も重要な設計上の進展は、「別々の山地道路や都市植樹枡で代用するしかなかった非鉄道対照」を、実際にインフラ敷と周辺を同時調査した原研究で置換できる見込みを具体化したこと。** しかし現段階では原個票が手元にないので、肯定・否定の新しい生態効果や普遍性の統合ではない。

本成果は**証拠探索のソースゲート**であり、現地原種×区画データを取得してはじめて真の鉄道特異的な局所生態学効果に進める。安易に205レビュー効果量や24site数を用いてメタ解析・CIを発行しない。研究リポジトリ上では既存v317のmeta／スペインv319のRCTと異なる独立出版物として識別する。

### 出典と再現

- [元2015研究全報告PDF（420頁）](https://ittecop.fr/images/docman-files/docman-files/rapport-final/2008/RF_2008_En_marge.pdf)
- [原著全文のsource source/ページ/地理点を監査するPython](../scripts/audit_Delbaere_2015_FR_BE_original_24site_railway_flora_report_v320.py)
- [6鉄道候補現地章からの配置/植物節確認](../scripts/audit_Delbaere_2015_six_railway_site_chapter_flora_source_v320.py)
- [6現地章の「La flore」および元マトリクス区画凡例の存在監査](../scripts/audit_Delbaere_2015_railway_six_flora_narrative_source_v320.py)
- [原PDF全再取得Actions #45013059824](https://github.com/zuizui0223/tetsu/actions/runs/45013059824)
