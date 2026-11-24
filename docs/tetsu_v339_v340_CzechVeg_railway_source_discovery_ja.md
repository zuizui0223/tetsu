# v339–v340：CzechVeg-Open 2026 年版 117,739区画から鉄道植生を抽出した原ソース監査

**真正新規大規模source監査 SUCCESS、2026-11-24。**
- 独立公開データ [Czech Vegetation Database v1.1, Zenodo DOI:10.5281/zenodo.19231152](https://zenodo.org/records/19231152)。117,739元野外plotを収録、原ヘッダー`CzechVeg-Open_headers_v1-1_2026-03-26.csv`のMD5 `0bc920462d1f9eac473e7b2a231e6782`を実際にActions上で照合。元仕様に`HABITAT, LOCALITY, REFERENCE, YEAR, SURF_AREA, DEG_LAT/LON, PRECISION`等がある。
- [v339全sourceテキスト分類CI #45372085741](https://github.com/zuizui0223/tetsu/actions/runs/45372085741)：合成テスト8/8成功、全117,739元headerを解析。
- [v340厳密認定と元対照支持CI #45372596136](https://github.com/zuizui0223/tetsu/actions/runs/45372596136)：合成テスト9/9成功。生息環境記述を厳密分類し、実際の元緯経度・調査年・面積・位置精度・引用元を使い、近隣植生区画の**存在可能性のみ**を先に確認。
- 原著全plotの元GPS・species×plot情報・元個票は再配布しない。source plot単位の気候変化や生態効果は**この2段階では推定していない**。

## 1．単純な「鉄道という文字がある」だけでは使えない

| 元sourceヘッダーの鉄道キーワード検索状況 | 元vegetation plots |
|---|---:|
| `HABITAT`（**実際の環境の文章**）に鉄道関連語があるplot | **239** |
| `LOCALITY`（**地理的目印の文章**）に鉄道関連語があるだけのplot | **3,573** |
| どちらにも鉄道関連語がないplot | **113,927** |
| **合計** | **117,739** |

LOCALITYに「線路・駅の近く」とあるのみでは、そもそも何を植生調査したか確認できないので**鉄道生息地と認定しない**。

さらに`HABITAT`の239件を実際の文章で見ると、厳密分類は：
| 元source`HABITAT`が示すもの（あくまで文面の候補） | 件数 | 本次解析での位置づけ |
|---|---:|---|
| 線路構内・鉄道ヤード・鉄道構造物上 | **61** | **railway候補**、ただし一部の石積み橋壁等はバラストと別の微生息地 |
| 鉄道の実際の盛土・側溝・切土 | **23** | **railway候補**、線路構内とは別クラス |
| 線路の上方・下方・200m離れ・近くにある湿地/がれ場等 | 62 | 直接鉄道採録ではない可能性が高く除外 |
| 駅構内の敷地・駅周辺とだけ書かれた区画 | 45 | 線路上/線路外が不確実なので除外 |
| 鉄道を含むが微生息地が不明な文章 | 44 | 保留 |
| 使われていない鉄道橋・未敷設の線路構造物 | 4 | **現役線路とは別**で除外 |
| **合計** | **239** | **84候補＝61＋23** |

この分類は自動的なテキスト抽出＋保守的認定であり、**まだ原研究を個別に読んで実際の線路中心からの距離・列車の有無・路床断面を確証したものではない**。

## 2．確実に出てきた新しい母集団とメタ解析の難点

厳密認定した**84plotは13元出典（REFERENCE）**、実調査年1960〜2024年、plot面積の中央値10m²。84区画はすべて面積・年が判明、83はチェコ国内に妥当な小数度座標を持ち、位置精度が100m以下なのは**22**、1km以下は**56**しかなかった。したがって文献数が増えたといっても、**異なる研究の設計・年代や座標精度が著しく違う**。

元sourceの非鉄道候補は、`HABITAT`にも`LOCALITY`にも鉄道関連語が書かれていない植生plot（**記載上はrailway未指定**）を採用。これは**必ず鉄道から無影響の真の非鉄道の同一サンプリング対照とは限らない**。結果を示すときは`nearby unlabelled reference vegetation`と呼ぶ。

| 原地理・年・面積・位置精度ゲート | 元84rail候補のうち潜在的近隣区画あり | 元sourceの独立文献支持 |
|---|---:|---:|
| **G1 同じREFERENCE、距離0.5km、調査±3年、面積比<=2、精度<=100m** | **2plot** | **1元出典** |
| **G2 同じREFERENCE、距離2km、調査±7年、面積比<=4、精度<=1000m** | **41plot** | **5元出典** |
| G3 異なるREFERENCE許可、G1と同条件 | 2plot | 出典混合 |
| G4 異なるREFERENCE許可、G2と同条件 | 51plot | 出典混合 |

元v340の`G2`で41plotに近隣候補があるが、同じREFERENCEで何区画が実際に**重複なし1対1**で組み合わさるかは別の問題。次v341で元species×plotを実際に結合し、形質被覆率と元species/plot数の基準を実行して**真に残る解析群を初めて確定する**。

## 3．論文の扱いと最重要の留保

CzechVeg v1.1 は**新しいフィールド原著ソースではなく、複数の元野外研究を収録した統合データベース**である。13参考文献REFERENCEは「13独立同条件のrail/nonrail比較実験」を意味せず、G2を満たし得る研究5件も対照の確実な非鉄道・面積・habitat matchingは保証されない。散布・バラスト微気候・植物の在来外来・繁殖成功・実際の列車頻度はこれらヘッダーだけからは推定できない。

もしv341でEIVE温度中心の差が正でも、それはあくまで`HABITAT`指定されたrailと近隣の未指定植生の**観察的比較**に留める。元の「寒冷地域ほど鉄道による気候制約緩和が強まる」仮説Aの直接検定には、独立気候勾配×真の鉄道/周囲非鉄道対照がまだ必要。

## 4．真正source識別と再解析用コード

`scripts/audit_CzechVeg_railway_headers_source_v339.py`
`scripts/audit_CzechVeg_actual_railways_nearby_controls_v340.py`
`tests/test_CzechVeg_railway_headers_source_v339.py`
`tests/test_CzechVeg_actual_railways_nearby_controls_v340.py`
`.github/workflows/tetsu-v339-CzechVeg-original-railway-label-audit.yml`
`.github/workflows/tetsu-v340-CzechVeg-real-railway-control-feasibility.yml`

次のv341は事前固定した[元source species/EIVE 対応解析ゲート](tetsu_v341_CzechVeg_rail_vs_reference_EIVE_source_frozen_protocol_ja.md)と対応。元ヘッダー・元source種×plotの著者生ファイルをGitHubに再掲せず、Zenodoと全ファイルMD5で照合して解析を再現する。
