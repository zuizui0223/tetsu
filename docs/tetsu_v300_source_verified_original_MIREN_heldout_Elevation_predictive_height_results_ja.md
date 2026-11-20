# tetsu v30.0 SOURCE VERIFIED — 実測標高で未使用地点の植物群集草丈を予測できるか

更新：2026-11-20 JST。**[真正元原著ソースのsite held-out CV source workflow SUCCESS #44661986353](https://github.com/zuizui0223/tetsu/actions/runs/44661986353)**。原source成果物 [#14596295975](https://github.com/zuizui0223/tetsu/actions/runs/44661986353/artifacts/14596295975) は原元site×taxon観察行、individual Elevation/GPS、LEDA種別trait、系統樹tip名を再公開せず、元Region/無記名Road/site採用数とsource CV結果のみ保存。

## 0. 生物学的研究の独立性とデータ同一性

MIREN Global Mountain ROADSの著者原著（2007–2023、植生観察167,804件、元真Plot Elevation field 6,446 plot、各Plot100m²、同Site×Year Roadside Plot1 / Inland Plot3）とOldenburg LEDAの種中央値canopy height・seed mass。実V.PhyloMaker2公開LCVP/WP/TPLの同一source phylogenies root-to-tip中央値の75%枝で切断した単系統clade×species site occupancy五分位の**交差**小群内のみ、paired height + seed massをspecies間に割当変更。元speciesと元Siteの対応・本当のElevationとsurvey Yearは全部元実データから固定。

原著4独立欧州Rail studyのsourceはここに含まない。MIRENは独立植生野外プログラム**1件**、元5〜6個Regionは1元field programmeの空間反復、3樹は同元現場植物群集の異名標準化感度。4/4鉄道原著をこの新しい実測標高CVの実験と勘違いしない。

本v30.0のsource scriptは先行v29.9**元tree×primary/strict計6組**の元原著Region数、実Road/実paired Site採用数、道路際・内陸のCWM高さ×標高signed slope observed−元phylo×occupancy label nullの既出版係数を1e-10以内で**原著から完全再現**するまでは新しいCVを報告しない。元Zenodo CSV original MD5、LEDA2表SHA256、元3実RData系統樹Git blob SHAも原資料再取得で成功。

## 1. 新しい真の未使用地点予測テスト

各元Roadに真の元paired physical Siteが6地点以上、元Elevation幅150m以上ある場合に限り、一つずつ本物のSiteをheld-outし、(A) **訓練地点のみ**のCWM高さ平均で予測、(B)同じ訓練地点だけから標高に沿う直線をfitしてpredict。真の未使用Siteの高さCWMの誤差を各元RoadでRMSE計算し、**ΔRMSE = RMSE(元平均baseline) − RMSE(元標高直線)**を主要推定量とした。原サイト自身の高さはそのモデルの学習に一切使わない。数理テストでは全観測Siteを一回ずつ除いた真の直接refitとPRESS/hat matrix式が小数点12桁で一致。

ΔRMSEが正なら、同一Roadの真の未使用地点では標高を使う方が平均値のみより予測誤差が小さい。SOURCEの追加期待と独立に「実際のΔRMSEが正か」を必ず表示。同時に、元植物のphylogeny75% clade×元taxon occupancy quintileの交差内のみheight+seed traitsをシャッフルした**398技術的ラベル割当**下で同じLOOCVを再計算し、実ΔRMSE−帰無期待ΔRMSEを出す。これは生態学的p値や地域間での予測区間ではなく、元種・標高と頻度の対応からの条件付き超過情報の診断。

## 2. Source 主結果：元basic flora support、共通5 Region / 14 Road / 210 paired physical Sites

| 元真正mega-tree | Roadside 実観測LOOCVでのRMSE減少 | Roadside 実観測−元系統×頻度label nullの予測利得 | Inland 実観測−元条件付き帰無の予測利得 | Roadside−Inlandの追加利得 |
|---|---:|---:|---:|---:|
| LCVP | **+0.035156** | **+0.028757** | **+0.040482** | −0.011725 |
| WP | **+0.037467** | **+0.030722** | **+0.044512** | −0.013790 |
| TPL | **+0.034798** | **+0.027366** | **+0.036327** | −0.008961 |

観測予測のRMSE減少の単位は元の欧州Rail3publicationで校正固定した z(log LEDA canopy height)。元Roadを等重み→Region内→Region等重みで記述。LCVP/WP/TPL**すべて同一のAUT、CHE、IND、ITA、NOR**に揃い、実Sourceの元Road14本・本物の物理Site210ペアという同一環境field supportで対応。新しく平均値の符号を揃えるために地域を削除したのではなく、**6か所以上の元未使用Site推定可能性という先のsupport gate**からTPLにも5Regionだけ残った。

**原Regionの異質性は重要**。Source Region単位ではRoadside実ΔRMSEが正3/5 Region、元traitラベル帰無を上回る利得は正4/5 Region。LCVP道路際の実ΔRMSEはITA+0.1014、NOR+0.0793、CHE+0.0121、IND−0.0032、AUT−0.0138。つまり元の予測性能は強くイタリア/ノルウェーに支えられ、一般化に慎重であるべき。元内陸の実改善は3/5 Region、帰無超過の符号は4/5 Region。LCVP/WP/TPLでの元regionを一つずつ除いたRoadside帰無超過利得平均の範囲は順に LCVP+0.01847〜+0.03843、WP+0.01727〜+0.04073、TPL+0.01547〜+0.03635 で正だが、これはsource Region5地点集合の記述的感度であり、独立研究をランダム抽出した95%CIではない。

## 3. Strict flora source supportと同一年Road感度

厳格な元Botanical floristic trait sample（各Plot元LEDA高度・種子重量両形質植物10以上、元floraの50%以上）と本sourceの6サイト／標高幅150m／2Road／70%元Altitude site保持条件にすると、MIREN元5Region/13Road/192Sitesの後、**LCVP4 Region/11Road/174 Site、WP3 Region/8Road/118 Site、TPL4 Region/11Road/175 Site**となった。

| 元tree strict | Roadside実ΔRMSE | Roadsideの帰無超過改善 | Inlandの帰無超過改善 |
|---|---:|---:|---:|
| LCVP | +0.047063 | +0.029429 | +0.054653 |
| WP | +0.041080 | +0.027841 | +0.039189 |
| TPL | +0.047696 | +0.027406 | +0.046084 |

Strictのsource Roadside実予測利得が正はLCVP3/4、WP2/3、TPL3/4。WP strictは元Regional delete-one感度のsource帰無超過改善範囲が**−0.00108〜+0.04641**で0をまたぐ。したがって「どの地域でも標高の追加予測が頑丈」と言わない。

元同Roadで異なる調査年が混ざる場合、**1つのsurvey YearのみのRoad**に限り、源original年度を前から固定した感度も実行した。basic source5RegionではRoadsideの帰無超過改善 LCVP+0.028159、WP+0.030489、TPL+0.027018、Inland+0.039092/+0.043198/+0.034250。strictの残存4/3/4RegionでもRoadside+0.027167/+0.029064/+0.025764と、元field年の混在だけが主平均利得を作ったとは考えにくい。

## 4. 二つの大事な区別

**A: 本物のheld-out site prediction と trait label null の超過は違う。** たとえばINDはRoadsideの実ΔRMSE ≈−0.0032とほぼ0/負だが、元系統・元頻度条件付き帰無を上回る差は+0.02であった。これは「帰無より当たりやすい」≠「標高予測が元平均モデルより優秀」を意味する。したがって真の環境予測としては**実ΔRMSEの符号を主に報告**、条件付きnullは空間的な種形質配置の追加診断とする。

**B: Roadside−Inlandの追加予測利得は平均では負だが、元Regionで符号が混在する。** basic LCVP、WP、TPLの source Roadside−Inland 追加超過改善は−0.0117/−0.0138/−0.0090で、Region削除平均範囲はそれぞれ −0.0267〜+0.0049、−0.0303〜+0.0020、−0.0226〜+0.0090 と**ゼロをまたぐ**。これは「交通回廊が環境勾配への依存を普遍的に弱める」という因果的な法則にはならない。

元植物相源source自体のLEDA2形質が一致したのはMIREN5,739植物のうち1,096 taxa（19.1%）だけで、気温・降雨・森林開放度・標高の関係は分離していない。真系統treeのtip掲載は当該元species自身の直接分子配列がある保証でもない。未観察の別Road、別Region、異なる植生調査protocolへ**持ち出し可能なtransferabilityは未検証**。元Rail4publicationにも同程度に認証されたElevation原Site列はまだ対応していないので、今回のΔRMSEを「鉄道独自の環境適応」に使わない。

## 5. 科学的な位置付け

v29.9は「高さCWMと標高の元現場内相関・傾きが、source clade×plant occupancy trait-label nullより強い」ことを示した。v30.0は一段強く、「同じ元Roadの未使用physical Siteで標高線形モデルがSite平均モデルより**平均的に**高い予測精度を持ち、元Source植物相の系統・頻度だけに予測性を帰せない」ことを示した。一方、source地域により実CV利得が負になるため**普遍性・機構の局地条件付き性**が新たな焦点になった。

次に科学的に反証力が高いのは、同じregionにあっても元Road間で標高への草丈の応答の**傾き**がtransfer可能かどうか、元Roadを伏せた予測で検定すること（Road interceptをどう推定するかは独立に設計する必要）。植物群集の環境傾斜は種中央値の物理的置換を反映し、同種内草丈可塑性/鉄道風圧/列車種子運搬を直接見ていない。

## SOURCEと再現

- 真原データ再実行SUCCESS https://github.com/zuizui0223/tetsu/actions/runs/44661986353
- 元source数値正本 https://github.com/zuizui0223/tetsu/actions/runs/44661986353/artifacts/14596295975
- 追加源site RData SHA/LEDA2表/source 17万observations元QAは別源出力として同原runへarchive
- scripts/evaluate_MIREN_native_elevation_heldout_stature_predictability_v300.py
- tests/test_MIREN_true_site_elevation_LOOCV_v300.py
- .github/workflows/tetsu-v300-original-MIREN-LOOsite-Elevation-phylo-occ-source.yml
- data/derived/tetsu_v300_original_MIREN_heldout_Elevation_predictability_summary.json

元individual botanist site×plant、元Elevm数値の個別値、GPS/author plot names、種のLEDA medianは本リポ/artifactへ配布していない。
