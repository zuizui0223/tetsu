# v342 — 全球鉄道暖地化検証の原データ取得ゲート（2026-11-24）

## 今回の検証結果
GitHub Actions [#45239947336](https://github.com/zuizui0223/tetsu/actions/runs/45239947336) で原著補足資料・GBIF分類API・WorldClim実データのアクセスを試行。5/5合成QAテスト PASS、原物に対するデータ出所／整合性は以下の通り。

| 原資料 | source正本の実取得 | 認証・検証 | 生物学的効果 |
|---|---|---|---|
| Everingham et al. (2019), DOI 10.1111/aec.12719, Appendix S4 (191.2 KBと出版社記載) | **未取得**: Wiley 公式2経路いずれもHTTP403 | 0件の真正source XLSX rowを閲覧 | 推定禁止。データが無いことは効果0ではない |
| GBIF v1 API, Cynodon dactylon | EXACT accepted species match, usageKey 6109637, geospatial-quality coordinate-positive counts **154,547** | API成功。ただし実際に全分布点を取得した訳ではない | 温度ニッチ未計算 |
| GBIF v1 API, Poa annua | EXACT accepted species match, usageKey 2704179, coordinate-positive count **952,090** | API成功。ただし検索ページ最大300、オフセット+件数100,000制限 | 温度ニッチ未計算 |
| WorldClim v2.1, 10m bioclimatic archive | **実取得成功 49,869,449 bytes** | ZIP SHA256 `00513224583665ec0f2f955a4ec252730c4deb2004cce9e793492a3f26df4dcf`、**19気候TIFF**、`wc2.1_10m_bio_1.tif`を検証 | 気温–植物種の結合未実行 |

WorldClim 2.1の過去の基準は1970–2000年で、微気候のバラスト温度を直接測るものではない。現代の出現記録を無条件に過去の気候平均へ対応づけても真の種温度ニッチにならない。GBIFの分布レコードに非在来の侵入・採集努力・位置誤差が含まれるため、出現点単純平均は真の至適温度ではない。正式な全球種気候ニッチにはGBIFの一括DOI取得、対象年代・native-rangeや公開地理的偏りの扱い、空間間引き・生育地を制御した気候空間の補正などが必要。

EveringhamのAppendix S4は原著ウェブ公開欄に存在を確認したが、取得した実ファイルは**0**。真のsource種×調査区データがないまま、本研究との「オーストラリア再現効果」やKashmirとの差を計算・宣言しない。元著者研究にはrail embankments, fire trails, roadside lawnなど撹乱勾配が含まれるが、original rail-sclerophyll plot と非鉄道のstrict同地点対応を元原票で検証できていない。

## 独立地域の一般化について
真正sourceで同一地点・類似面積・同方法のrail/非rail比較が複数独立気候地域から得られないため、**元仮説A（より寒い広域気候ほど鉄道群集が暖地性に偏る）は未検証**。スロベニアのBS/BTは共に線路内、MIRENは道路、カシミールのInteriorは線路付近約150mで自然未影響区画とは確約されない。この三者の効果を単一のrail-versus-nonrail meta効果として統合しない。

## 正しく次段へ進む手順
- CzechVegの既取得データ体系を用いた別検定 v343（同一出版論文×植物地理区×年代×面積帯の歴史的RAIL-text/ROAD-text）。こちらも同一地点の自然対照ではなくCzechという**一つの国の研究集積**なので、方向の外部点検だけが目的。v343の原source実行成功と実際の支持数を確認してから結果を記載する。
- 豪州本物のAppendix S4取得は出版社/著者の正式公開経路または著者同意が必要。元原票の存在だけでは温度指標との結合は不可。
- DOIとデータの出典、未取得制約は推定値とは区別し続ける。未計算の気候効果や擬似種リストを再現された知見と呼ばない。

## 再現
実装 `scripts/audit_external_railway_global_thermal_sources_v342.py`、テスト `tests/test_external_railway_global_thermal_sources_v342.py`、実行 `.github/workflows/tetsu-v342-railway-global-thermal-source-gates.yml`。初回は合成テストのシグネチャバグで停止し、修正後のsourceゲートCI #45239947336が真正アクセス検証を成功。原著元plot列/元GPSを転載しない。研究枝`research/tetsu-v342-independent-thermal-reference-preflight-20261124`、Draft PR #54；main未変更。
