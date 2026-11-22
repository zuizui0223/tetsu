# v321 解析前のソース確認仕様：En marge (2015) Rapport final の「Annexes RF」に原植物区画表が存在するか

2026-11-22 JST。前回v320では真正420頁のフランス・ベルギー En marge 報告（原著24インフラ地点、うち鉄道6地点）を全文監査したが、インフラ敷vs周辺マトリクスの **完全な種×プロット**を検証できなかった。今回公式の[原研究プロジェクトページ](https://www.ittecop.fr/fr/tous-les-projets/recherches-2008-2009/item/33-en-marge)に「Rapport final」と**別リンク「Annexes RF」**があることを知った。これを元未取得の付録/データリンクとして検証する。

## 最初にすべきこと（結果や植物種名を検索する前に固定）
1. 出版者公開のプロジェクトページをGETしHTML（またはJSON-LD文書）を検証。元の見出し`En marge`、`Annexes RF`、`Rapport final`の存在を確認し、両者のhrefを**相対URLから真正同一ドメインの完全URL**として復元。出典ページのSHAと発行媒体、HTTP応答を記録。別のブログや非公開文書を真正付録として使わない。
2. `Annexes RF` からリンクされる公開ファイルを同じ著者/原研究IDで照合し、PDF/ZIP/XLSX/CSV形式、バイト・SHA、資料内見出し・页数・目次を記録。HTTP403/ログイン要求なら無断突破せず`not_accessible`。HTML表示ページであれば、公式から辿れるリンクだけをたどる。
3. 植物rawデータの真の支持条件は元 **site17–21/24×(infrastructure,matrice)×plot ID×season×taxon×cover** が揃うこと、元学名とplot unitが全24siteを共通で解釈できること。原付録に種数の文章やラテン学名の部分リストがあるだけなら完全な原出現行列と解釈しない。区画ID・原面積・季節に欠損があれば厳格pair-onlyに限定。
4. 真の対応で**≥3鉄道地点、各最低1対応マトリクス対照＋全植物種出現**がソース確認できた場合のみ、鉄道−近接対照のサイト等重みの種数差、外来割合、温度親和性差を別実装で計算する。元原資料の全現地観測数、季節、マトリクスと盛土の物理位置が照合できなければ自動で推定禁止。
5. 初回は原資料のソースゲートと妥当な個票支持量を算出。原PDF・著作物や種×plot rawやGPSはリポへ再配布せず、公開URL＋SHA、シート・行列寸法、匿名集計のみ。鉄道6候補は1研究1都市圏の一部、6独立研究・6気候帯ではない。
6. 原付録未取得でも、責任者 **Denis Delbaere** と植物部門 **CBNBL / Conservatoire botanique national des Hauts-de-France, Digitale2** が元記録の対象管理者。プロジェクトの進展としては、原付録の真正公開リンク探索・機関DBの問い合わせ仕様までを元GitHub Issue #50に追加し、メールの実送信はユーザーの明示了承まで行わない。

## 外部一次資料の優先順位
- プロジェクト公式ページ「Annexes RF」、原真正報告「Rapport final」（v320 verified PDF SHA256 `25a1a38512156bf1264b697efd78984a1abf9217f1da21e6ef2c5623da509405`）。
- 同所の文書・付録、ITTECOP公開ダウンロード。
- CBN des Hauts-de-France [Digitale2](https://jardins.cbnbl.org/ressources/) は植物観測の大規模DBだが、一般市町村の種観測をEn margeの厳密プロット対と**同一資料だと偽装しない**。

**計算は真の原表取得を満たしたときだけ**。進捗の目的は、ソースゲートを増やすことではなく、厳密な対応植生効果が原資料で同定できるかを確定させること。
