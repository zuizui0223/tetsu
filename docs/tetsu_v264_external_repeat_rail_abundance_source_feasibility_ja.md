# 鉄 v26.4 — 初期量的記録→後年再記録の独立研究拡張：公開資料監査

確認日 2026-11-18 JST。v26.1–v26.3は1つのKashmir鉄道調査（2014/2017）を使った3感度系列であり、独立野外研究3本ではない。複数国のマクロ解析へ展開するため、同じ現地地点・同じ区画で初回の種別abundanceと後年の再記録を測った独立データを探した。**以下の未解析資料を鉄道の新しい肯定結果に加算しない。**

## 真の採用ゲート

- 物理的な線路/道床として現地同定でき、同じ地点・区画を少なくとも2時点で再調査した。
- 初回の種×plotの数値abundanceと後年の同種再記録を復元でき、面積・調査季節・観測努力・再配置履歴が明記される。
- 元ファイルの利用許諾、SHA/MD5、site ID、taxon名と測定単位を検証する。検出されない種を確実な不在/死亡としない。
- 同じ研究の地点・年月を独立国や独立論文へ水増ししない。leave-study-outに使える独立研究が複数揃った場合だけ一般化を検討する。

## 一次ソースごとの判断

| ソース | 真に公表されていること | 現在の分類 |
|---|---|---|
| [Kashmir Rashid et al. 2021 Dryad](https://doi.org/10.5061/dryad.0p2ngf20h) | 31地点、2014/17年、4種別区、作者真正122,187B XLSXから2014量的順位→2017再記録を27地点で検証成功 | **確認済み鉄道反復：1原著** |
| [Slovenia Šabić et al. 2026](https://doi.org/10.1007/s44473-026-00155-z) | 89地点・173同側道床肩/法尻ペア、2022/23年、原種リスト | 同地点・同zoneの再訪とは未確認。横断的鉄道微環境研究として扱う |
| [ReSurveyGermany Jandt et al. 2022](https://doi.org/10.1038/s41597-022-01688-6) / [CC BY公開データ](https://doi.org/10.25829/idiv.3514-0qsq70) | 7,738 semi/permanent plots、23,641 plot observations、458,311種別cover記録、1,794維管束植物、1927–2020、EUNIS生息地分類・plotメタデータ。ソース表名：Header_ReSurveyGermany.csv、ReSurveyGermany.csv | **最優先の公開メタデータ監査対象**。物理的railway plotがあるか不明。一般植生対照と鉄道サンプルを混同しない |
| [ReSurveyEurope 2024](https://doi.org/10.1111/jvs.13235) / [公式利用ルール](https://euroveg.org/resurvey/) | 79,190 plots、283,135観測、449再調査project、約50,720 plotsに量的記録 | 共同研究・提供者の許可制。全件オープンデータではない。railway対象の原plotがあるか未検証 |
| [UK Countryside Survey 2007](https://catalogue.ceh.ac.uk/documents/4062f157-8990-46de-8f50-565eaa9ef847) | 1978/1990/1998/2007の植生再調査で種別presence/cover（一部）を記録 | 一般植生の反復対照候補であり、現場railwayという証拠なし |
| [Stuttgart station Winkler 2026](https://doi.org/10.3390/ecologies7020052) | 2024–2025、単一大駅・68 taxa | 元site×taxon×abundanceの自由公開は未確認。独立多地点とは数えない |
| [China–Laos rail Li et al. 2026](https://doi.org/10.1016/j.jenvman.2026.130380) | 14駅・84植生/土壌/菌類plot、物理的station対未攪乱地の横断比較 | 植物の元84plot完全行列と反復年情報の公開は未確認。空間対照には有望、再訪研究とは別 |
| [France Rakotoson et al. 2025](https://doi.org/10.1111/wre.70025) | 2021年、同一研究の2m²鉄道128・法面266plot | 原plot出現データは著者申請制。別国の横断的生息地対照として重要 |

## 優先する次の真正ソース検証

1. ReSurveyGermanyのオープン Header メタデータを現実に取得し、鉄道・道路・荒地などの明示的source habitat labelがあるか、そのサンプル数・同一plot反復年とabundance completenessを監査する。EUNIS general anthropogenic は rail現地同定と同義ではない。
2. ReSurveyEuropeへ承認制の正式研究提案を検討する。railwayのみの該当resurvey plot/supportがあるか管理者に確認し、許可なしにrestricted raw dataを入手しない。
3. 少なくとも2件目の独立真正鉄道再訪が得られるまでは、Kashmirで得た種×habitat条件付きのrank予測改善を『全球的な鉄道群集形成の普遍則』と呼ばない。

統合論文は時間軸別の2層に分離する。縦断は2014→2017に相当する再記録・量的初期値、横断はCZ/IT/HR/BA/SI/SK等の駅間の種置換と生活史・機能形質。異なる観測単位の値を同じメタ回帰効果量に無理にプールしない。
