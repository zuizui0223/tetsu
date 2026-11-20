# v29.9 SOURCE VERIFIED — 標高に沿う草丈の群集交代は、系統と種の出現頻度を固定しても残る

**元原著と数値の認証日：2026-11-20 JST。** [GitHub Actions #44659330637 SOURCE SUCCESS](https://github.com/zuizui0223/tetsu/actions/runs/44659330637)、[元原資料由来の効果・独立QA artifact #14594079180](https://github.com/zuizui0223/tetsu/actions/runs/44659330637/artifacts/14594079180) に実コード・原著入力SHA・元地域別効果値の正本がある。

## 0. 何を「原著データで測定した」のか

MIREN Global Mountain ROADSの元現地調査（2007–2023年、植物出現167,804件）、同じ元Site×YearのRoadside Plot1とInterior Plot3の各100 m²で測定・記録された **Elevation (m)** と、元出版LEDAのspecies median canopy heightとseed mass、真のV.PhyloMaker2（LCVP/WP/TPL）の枝長付き元published treeを真正データとして使用した。

先行v29.8で確認した本物の原著cohort、実元Region、Road、Site、Year、種×siteの0/1出現、元の標準化した z(log canopy height, log seed mass) 2形質ペアを動かさない。元tree root-to-tip枝長中央値の75%で切断した実単系統cladeと、元sourceの植物種出現頻度五分位を**同時に固定**する交差群の内部のみでspeciesの2形質ベクトルをシャッフルし、元植物の標高勾配との対応がこの帰無に対して何を上回るかを測る。

**新しい主要量**：同じ元Road上の物理Site間で、元plotの高さCWM（z(log height)）を現地Elevation(km)に回帰した **方向付きの傾き**。元Roadごとに平均、次に元RegionでRoad等重み、最後に有効Regionを等重み。元標高に実在する地点だけ、同roadに4 site以上・標高幅150m以上、同Regionに2 roads以上、前のtrue-tip siteを80%以上保持する条件。補助量は同Road上で元site間の標高差の絶対値に対する CWM高さの絶対差の傾き（Siteペアを別の独立研究とはカウントせずRoadごと集約）。

**SOURCE QA**：元plotのElevation情報6,446件が本物の現地plot metadataとして数値/IDを検証通過。先行のv29.8 3 original tree×primary/strict 6通りの元source有効Region数とRoadside/Inlandの真joint phylo+occupancy帰無height beta効果を1e−10以内で完全再現してからのみ、新しい標高帰無を走らせた。元site×taxon/plot Elevation/GPSはartifactに再配布していない。**実作者データから2,256個の数値差分・符号・LOO・地理/調査年に関する独立QA checksが通過**。数学テストは元plot metadata 7件と環境傾斜の数理10件に成功。

## 1. NEW 主結果：同じ道路内で標高が高いほど草丈CWMは低い

Primary＝元両LEDA trait既知植物5種以上・元区画種被覆25%以上。元MIREN8 Region/22 Road/272 Siteの初期cohortから、元系統樹tip・元標高ソースの事前基準でLCVP5 Region/14 Road/210 Site、WP5/14/210、TPL6/16/224が対象。以下は元Region等重み記述平均。単位は**z(log canopy height) / Elevation 1,000m**。

| 真元mega-tree | 道路際：実際の高さCWM標高傾斜 | 道路際：観測−joint phylo×occupancy null | 内陸：実際の標高傾斜 | 内陸：観測−joint null | 道路際−内陸の過剰傾斜 |
|---|---:|---:|---:|---:|---:|
| **LCVP** | −0.408509 | **−0.266551** | −0.609980 | **−0.328002** | +0.061452 |
| **WP** | −0.428897 | **−0.296026** | −0.656502 | **−0.350704** | +0.054678 |
| **TPL** | −0.396809 | **−0.212749** | −0.610287 | **−0.364725** | +0.151976 |

対象の元source全Regionで、Roadside高さ傾斜の観測−帰無差は負（LCVP5/5、WP5/5、TPL6/6）。同じ現地Roadの内陸でも5/5、5/5、6/6が負。Regionひとつずつ除いた平均の負の範囲も、元樹とzoneの全組合せでゼロ未満。これは **実測した標高勾配に沿う群集の草丈の低下が、同phylogenetic clade・似たsite頻度の種の間でtrait-labelをランダムに付けた場合の期待より大きい** という観測的・機構候補を示す。

## 2. 補助結果：Siteの標高差が大きいほど草丈差も増える

同じ原Road上の異なるSite間、CWM草丈の絶対距離を元サイトのAltitude差の絶対値で回帰した傾き（原Road内集計→等Road→等Region）も、観測−nullが正。

| 元tree | 道路際：距離増加の超過傾斜 | 内陸：距離増加の超過傾斜 |
|---|---:|---:|
| LCVP | +0.138025 | +0.286884 |
| WP | +0.149583 | +0.321822 |
| TPL | +0.105257 | +0.126085 |

「元植物の地理的形質構造」という以前の記述のうち、少なくとも一部は**実測した環境差で直接予測できる方向**にある。今回の元sourceは道路沿いの単一プログラムなので、独立した鉄道原著4本に同じ値を流用しない。

## 3. 厳格な元trait source supportと調査年感度

形質既知植物10種以上/元plot flora50%以上のStrictでもsource真treeが有効なLCVP5 Region/13 Road/185 Site、WP4/10/129、TPL5/13/186で観測−joint nullの符号は以下のとおり：

| original mega-tree（Strict） | 道路際 高さ×標高の追加傾斜 | 内陸 高さ×標高の追加傾斜 | 道路際−内陸 |
|---|---:|---:|---:|
| LCVP | −0.272918 | −0.342623 | +0.069705 |
| WP | −0.234098 | −0.290858 | +0.056761 |
| TPL | −0.231654 | −0.286245 | +0.054590 |

**元Road内の調査年が混在しないRoadのみ**残す事前の感度解析も実施し、PrimaryではLCVP5 Region、WP5 Region、TPL5 RegionでRoadsideの観測−null高さ標高傾斜は順に **−0.282348、−0.297200、−0.233277**、Interiorは **−0.344301、−0.352852、−0.287277**。StrictでもLCVP4 Region、WP3 Region、TPL4 Regionに減るがRoadside−0.307828/−0.282416/−0.257785、Inland−0.420898/−0.384069/−0.342624で負方向を維持。元同Roadで異なる調査年の混合が主効果の唯一の原因という説明には整合しにくい。

## 4. 道路際と内陸の違い／誤った因果を避ける

道路際の高さ×標高の負傾斜は、内陸より**負方向が弱い**（Roadside−Inlandの観測−帰無差は3元treeで平均+0.055〜+0.152）。ただし個別元Regionの差の符号は一様ではなく、LCVP基本1/5、WP2/5、TPL2/6 Regionで道路際−内陸の超過傾斜は負である。平均の一貫性が全球での一般法則や道路という処理の因果効果を保証するわけではない。

標高は気温、水分、土壌、草原・森林の開放性、除草・道路管理と共変しうる。今回のモデルはElevationを元の標高として**実測環境座標と植物形質地理の条件付き対応**を測るもので、標高が直接の選択圧である証明や、気候と草刈りの機構分解ではない。元LEDA traitは種中央値であり、現地個体の形態可塑性や遺伝的適応を測っていない。

MIRENは**独立野外調査プログラム1件**。元source Region5/6は論文5/6本ではなく、3 megatreesも異なる野外独立実験ではない。Rail v28.5–29.8の4原著に元siteごとの同じ品質の標高観測を確認したわけではないため、このMIRENの勾配を「鉄道で証明した環境機構」と誤認しない。現実の稼働中rail ballast vs同面積・同年・開放非rail砂礫地をさらに検証する必要がある。

## 5. 新しい一般生態学の含意

**実際の形質−地域配置の結び付きは、地域植物相に属す系統とspecies prevalenceだけで説明されない。少なくとも山岳道路沿いと隣接する通常植生では、環境勾配の実測値が残余の高さ形質構造と系統的に対応する。**

これは複数の異なる元群集形成段階を分ける次の理論仮説の候補になる：
1. 地域種プールの到達可能性と系統構成（以前のgenus/family/true phylo null）
2. 同系統群内の種の出現頻度（v29.8 source original）
3. 元現地Elevationに沿った追加の高さCWM sorting（今回v29.9）
ただし各段階のnull差は加法的な分散寄与率ではなく、ここから割合を合計したり「自然選択が何%」とは言わない。

## Source codeと正本

GitHub Actions： https://github.com/zuizui0223/tetsu/actions/runs/44659330637

真正元SOURCE artifact： https://github.com/zuizui0223/tetsu/actions/runs/44659330637/artifacts/14594079180

コード：
- scripts/audit_MIREN_original_native_elevation_source_v299.py
- scripts/test_MIREN_source_elevation_phylo_occ_stature_v299.py
- scripts/validate_MIREN_native_Elevation_environmental_source_v299.py
- tests/test_MIREN_original_Elevation_source_v299.py
- tests/test_MIREN_native_elevation_phylo_occ_mechanism_v299.py
- .github/workflows/tetsu-v299-source-MIREN-native-Elevation-vs-true-phylo-occupancy-trait-null.yml

元original source species×site行列、GPS、個別Elev値、LEDA植物別値やmega-treeの個別tipはGitHubに再配布していない。元のsource study/region/road/site件数と環境傾斜・帰無推定を元研究単位で保存する。
## v29.9 最終Source公表数値照合（2026-11-20）

**[最終原著再実行・元source集計独立検証 SUCCESS #44659330637](https://github.com/zuizui0223/tetsu/actions/runs/44659330637)**。真のMIREN元CSV MD5、LEDA元2表SHA、V.PhyloMaker2元3樹Git blob SHA、数学17テストに加え、前source v29.8の全6組の元対応群集baselineを再現し、v29.9元地域別の草丈・標高観測／条件付き帰無／差分・時間感度と地域除外の **2,256項目** を照合した。その上でGitHubブランチに凍結した読みやすい集計JSONと原sourceの**54項目**を独立照合し、すべてPASS。新しいSource効果量や閾値は前の成功runから変更していない。真正元source派生アーカイブ **[artifact #14594079180](https://github.com/zuizui0223/tetsu/actions/runs/44659330637/artifacts/14594079180)** が正本。元site個別Elevと元species/taxonは非公開。
