# v29.9 — 山岳道路・内陸植生の草丈の地理的分化を標高勾配で直接検証する

2026-11-20 JST、実数値を見る前の原sourceに基づく解析仕様。元のMIREN山岳**道路**データと内陸対照を用いる。鉄道の新しい現地データや、植物個体の適応進化を測ったものではない。

## 既存の確定結果と今回の問い

v29.8は本物の欧州鉄道4原著と独立したMIREN道路／内陸のデータに対し、植物種の配置を完全固定し、元公開V.PhyloMaker2系統樹のroot-to-tip深度75%の真単系統群内、かつ**元species出現頻度の五分位内**だけでLEDA高さ＋種子重をペアで交換しても、草丈CWM地点間距離の観測−帰無期待値が多数の地域・原著で正と確認した。これは植物の形質と空間的な種配置の規則性であって、鉄道特異的な選択・列車散布・草刈り/車両風・同種内適応を証明していない。

次の検証は**地理構造が実環境勾配（元現地Elevation）に対応するか**。元のMIREN観測CSVには各PlotのElevation(m)が収録されている。標高は実測環境軸だが、気温・降水・森林開放度と共変しうるため「気候の因果媒介そのもの」と呼ばない。

## Source-only SUPPORT gateを既に独立監査

[真正原著標高メタデータQA 成功 #44656224041](https://github.com/zuizui0223/tetsu/actions/runs/44656224041)、元Zenodo MIREN v3 CSV：17,899,241 bytes、MD5 b640195cc44f252ccb6ddb28eb6fff13、元植物観測167804行。LEDA原著height/seed元source2 SHAも照合。原著元6446のplot elevationメタデータが完全かつ一貫して使えた。

原著のRegion×Road×Site×YearについてRoadside Plot1とInterior Plot3を同年・同地点で組とし、両PlotのElev差が100mを超えないものだけ平均してサイトの標高を定義する。元各Plot内の重複植物種行を独立した標高サンプルとして数えない。元同一Road×Siteの最初の適格調査年1件だけ使用する。

最初からSOURCE固定された形質条件はprimary（各Plot両LEDA既知種5以上、元植物相25%以上、8 Region 22 Road 272 Site）とstrict（10種50%以上、5 Region13 Road192 Site）。v29.9環境のみの支持域監査で**元8 Region中7、strict5 Region中5**が「1 Roadにつき真4Site以上かつ標高範囲150m以上、Regionに真Road2本以上、元Siteの80%以上を保持」という事前条件を通過した。ただしこの比率はphylogenetic matching前であり、最終source解析対象がさらに減る可能性を明示する。

## 真正source×phylogeny×taxon-occupancyでの条件付き環境検定

3本の元公開V.PhyloMaker2系統樹LCVP/WP/TPLは、authorの固定commit 7af3fb5152f691af2e4ec9d5e2e467d1b50505e9の元R-dataをSHAで認証。直接tip照合した植物だけを使い、欠落種のGraftをしない。元LEDA形質は欧州鉄道3原著622種で凍結した共通z(log height, log seed mass)軸を適用する。**元source植物site×taxonと元Elevationは固定**する。

新しい環境アウトカムを2つ事前定義：

1. **Primary：同じRoad内の草丈CWMの標高への符号付き傾き**。元作者の真サイト標高(m)をkm単位とし、同Roadのsite≥4でOLS傾きを推定。単位は z(log height) / 標高1000m。正は高標高ほど高い植物、負は高標高ほど低い植物。Road単位の傾きをRegion内等重み。
2. **Secondary：同じRoad内の異なるSite間の草丈CWM絶対差に対するSite標高差の絶対値の傾き**。単位は z(log height) 距離 / 標高差1000m。同じRoad内の全元Site対を用いるが**サイト対は独立野外研究ではない**ためRoadごとに先に集約し、Region内等重みとする。

同じ元Regionの植物種の元サイトへの割付・元species出現頻度・元Siteの標高・Year、元LEDAのheight+seed mass2D vectorは保存。元phylogenetic clade75×species original occupancy 5分位の**交差群内部**のみで2D trait vectorの種名割当てをシャッフル、二固定seed×199＝398個の技術ラベル割当て。Roadside/Interiorに同じ置換されたspecies trait vectorを使用する。各Roadに最低4真Site、標高差150m以上、Regionに最低2真Road、元direct-tree-tip合法Siteを80%以上保持することを結果符号を見ず固定。全てを満たせないRegionは**NOT IDENTIFIABLE**で、新しい効果量0とみなさない。

草丈CWM環境傾斜の実観測と元conditional trait-label null平均との差を、Roadside, Inland, Roadside−Inlandで計算。元3 source Regions以上が適格ならregion等重みの記述平均とleave-one-region-outを併記するが、元Regionは**MIREN単一原著プログラム内部の反復地域**であり、3本の系統樹も別野外研究ではない。これらの技術シャッフル順位は野外field p値でない。

## 時間交絡の監査・感度

元MIRENの「最初の適格年」は同じRoadでもSiteごとに異なりうる。前source auditではPrimaryのCLS・IND・OREなどで元Road内複数survey yearsが存在した。そこで元のRoadの**有効site全てが同じ調査年**のRoadだけを感度解析に残す。そこでも2本以上のRoad/Region、元direct-tip Siteの60%以上を保持した場合だけ結果を算出。支持不足なら感度解析はNOT IDENTIFIABLEであり、基準解析の「有意/非有意」を代入しない。

## 解釈の判定

標高に沿う符号付きCWM傾きがRoadside/Interior双方で負→高標高で草丈の低い群集に置換される候補。元clade/頻度を条件付けても負の超過が残れば、地域系統種プールの地理分布だけでは説明しきれない環境関連を直接確認したことになる。ただし観測された標高と共変する温度・開放度・降雨/人間の管理を分解したものではない。0近傍や地域で符号不一致なら、これまでの草丈の非ランダムな空間配置が単純な一方向標高勾配と同一ではない可能性がある。木本/草本の比率・trait gapによる植生組成依存にも注意する。

元Rail4出版は現時点で同様の元site Elevation+survey grain照合をしていないので、今回のMIREN標高係数を**鉄道の係数として流用しない**。新規森林・開放砂礫地対照や元気候・管理データを別途調べる必要がある。

## コード・真source

- scripts/audit_MIREN_original_native_elevation_source_v299.py
- scripts/test_MIREN_source_elevation_phylo_occ_stature_v299.py
- tests/test_MIREN_original_Elevation_source_v299.py
- tests/test_MIREN_native_elevation_phylo_occ_mechanism_v299.py

新しい生態学的効果をSource Actionsで検証するまで、ここには実係数を記載しない。元の個別Site×species、Elevation・GPS値、元LEDA個別trait値、元phylogenetic tipラベルはGitHubに転載せず、研究/Region/Roadの支持数と推定傾斜の集約値のみを公開する。