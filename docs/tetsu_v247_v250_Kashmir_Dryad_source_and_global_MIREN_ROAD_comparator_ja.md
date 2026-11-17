# 鉄 v24.7〜v25.0：カシミール鉄道真正原データ取得監査とMIREN全球道路植生の外部比較

**2026年11月17日、公開原著での新しい広域研究検証。** tetsu既存研究ブランチ `research/tetsu-slovenia-identity-guild-v210-20261116` のみに実装、main未変更。今回のMIREN結果は**道路であって鉄道ではない**ことをすべてのコード・報告で明記し、既存の欧州5か国・4原著を5つ/6つの鉄道原著へ水増ししない。

## I. インド・カシミール鉄道：31地点×4区×2014/2017年の独立原著は実在、ただしデータbyteに認証ゲート

**原著** Rashid, Haq, Lembrechts, Khuroo, Pauchard & Dukes (2021) *Journal of Applied Ecology* 58:1967–1980, DOI https://doi.org/10.1111/1365-2664.13961。元研究者はカシミール鉄道の**31のT型空間地点（16駅・15駅間）**で**2014年・2017年**に同地点を再調査。各地点で**線路側、法面、内陸1、内陸2の4つの50×2m方形区**を配置、原研究は植物123種を記録、内陸は鉄道から最大約150m。**同一地点で鉄道からの距離効果と時間変化が直接分けられる設計**。124区画/調査年であり、「248の独立した駅」ではなく同じ124区画の再調査である。

2026-11-17真正Dryad RESTでタイトル・DOI・version・sourceファイルdigestを確認：
- 真のDOI [10.5061/dryad.0p2ngf20h](https://doi.org/10.5061/dryad.0p2ngf20h)（原著の引用の一部に`0pngf20h`という脱落表記があるが、公式RESTが真正と確認したのは`0p2ngf20h`）
- Dataset ID 60988、正式version ID **99735**、file ID **546003**、ファイル名 **`Matrix_plot.xlsx`**、元メタデータ表示サイズ**122,187B**、元source SHA-256 `529f3b97521933018e12623b61f531f012c67241cae2602a6675dbd534192474`。
- [真正出典確認v24.7成功 #44012736064](https://github.com/zuizui0223/tetsu/actions/runs/44012736064)（first-party APIの正式versionへのリンクからfile IDとSHAを取得）
- **[v24.8実ファイル読取 #44012873510](https://github.com/zuizui0223/tetsu/actions/runs/44012873510)は原バイト取得に失敗**：`https://datadryad.org/api/v2/files/546003/download`が同意済み認証なしGitHub Actionsから**HTTP 401**。Dryadデータセットのメタデータが公開されていることと匿名APIファイルダウンロードは2026時点で別。**元Excelのシート名・site×taxon構造や種別効果は未検証であり、偽のIndia 31×123の係数を作っていない**。
- ユーザーが公開ブラウザから正規ファイル`Matrix_plot.xlsx`をダウンロードし、この会話に添付すれば、元SHAと元field/plot/yearラベルを検証して本研究の直接のrail→interior 2014/2017の科学的検証を追加できる。公開上のダウンロード障壁をすり抜ける偽装やAPIへの未許可アクセスはしない。

**中国–ラオス鉄道** Li, Zheng & van Kleunen (2026), *Journal of Environmental Management* DOI https://doi.org/10.1016/j.jenvman.2026.130380、**14駅・84vegetation/soil/fungi plots**に同一調査計画の disturbed/un-disturbed 対照があり、真正種×84plot行列が取得できれば最重要独立の熱帯・亜熱帯鉄道検証。ただし原完全行列が配布可能なことは現状未確認。**フランス** Rakotoson et al. (2025), *Weed Research* DOI https://doi.org/10.1111/wre.70025 は**2m²の394真正植生plot（線路128、法面266、255植物種）**を2021年に同一研究で調査しており極めて強い同一努力比較だが、本文Data Availability Statementは**著者へのデータ提供依頼制**と明記。これらを勝手に「解析済み観測」へ足さない。

## II. 別の真正オープン研究：MIREN全球道路の同方法100m² road-vs-interior対照

**出典** Seipel, Haider & MIREN Consortium. [MIREN survey of plant species in mountains v3](https://doi.org/10.5281/zenodo.15749536) 、元CSV`MIREN_road_survey_(2007-2023)_zen_v3_002507.csv`、**実17,899,241bytes**、元MD5 **`b640195cc44f252ccb6ddb28eb6fff13`**。対応出版 Haider et al. (2026) *Global Ecology and Biogeography* DOI https://doi.org/10.1111/geb.70284 。これは**道路**のデータで、源`Plot 1`は道路沿い、`Plot 2`は道路から中間、`Plot 3`は内陸、全100m²（2×50m）。元観測年2007〜2023のサイトと再調査があり、個々のsource Region×Road×Site×Year×Plotを元sourceラベルから復元。**どのMIRENサイトも鉄道調査地点に置換していない**。

**[v24.9 SOURCE SUCCESS #44013341079](https://github.com/zuizui0223/tetsu/actions/runs/44013341079)**：
- 実元CSV **167,804**種観測行、原ソース`Region`コード**26種類**。出版元本文では**25地域**と称している。25という説明とCSVコード26の相違を黙って潰さず保持し、地域コードを独立国家・独立論文数として扱わない。
- 元二名法に解析できた植物**5,739名**中、元大学LEDAの**葉群高+種子重とも対応できたのは1,096名（19.10%）**。地域・系統での形質非対応が非常に大きく、いわゆる全球効果に一般化できる範囲を強く制限する。
- 元Plot 1と3が同じ空間Source Site×Yearで実在する組だけを使用。地域×道路×地点ごとに**形質要件を満たす最初の実調査年のみ**使用（同じ地点の他年を独立としない）。両元plotの形質既知種**各5種以上**・元分類群種数の**25%以上**、さらに1地域内で少なくとも**2本の真正独立道路×各4地点**を要求した支持域は**8 Regionコード・22道路・272空間地点**。25地域全体の代表と呼べない。
- 2D CWM はlog LE DA 2形質を元MIRENの真の生物種1,096種で同一z軸に変換、種の地点内等重みのCWM2D間ユークリッド距離とした。**同じ道路内・同じ本来の道路地点ペア**で、roadside Plot1同士とinterior Plot3同士の距離差`(roadside1–interior3)`を作り、道路内平均→地域平均。個々の区画ペア千件を独立の効果量として誤用しない。
- 地域符号別の機能距離差：AUT −0.13312、CHE −0.02852、CLA −0.16924、CLS −0.09107、IND **+0.08902**、ITA −0.23999、NOR −0.10915、ORE −0.14841。**8地域のうち7/8でroadsideの方が機能的に収束した方向**、原地域等重みの差 **−0.103811**。
- 同一Siteで形質既知分類群数を`min(roadside,interior)`に等化する59技術希薄化では、元8地域等重み距離差**−0.056138**、図上の技術2.5–97.5% **−0.07640〜−0.03541**。この値は本当の独立地域間の95%CIではない。元税onomic Jaccard差は調整前 **+0.053171**、種数等化後**−0.022057**と符号が変わる。**種名の均質化は豊富さの違いに脆弱**、機能的重心距離は同じ方向に留まる。

## III. 最重要：LEDA missingnessを厳格化したら残るのか

**[v25.0 SOURCE SUCCESS #44013611739](https://github.com/zuizui0223/tetsu/actions/runs/44013611739)**。同じMIRENデータ/同じ事前の測定元sourceを使用し、元本来の路線×地点×年の区画の本物の形質カバーを5/10種、25/50%で多段階に厳格化した。これらは**観察研究の感度解析であり、独立した国の研究を増やしたわけではない**。原source地域を1個除いたときの平均方向も計算。

|解析要件（各Plot1と3で元LEDA両形質）|元解析で残るsource地域コード|道路数|空間地点|元source地域等重みのroad−interior CWM距離差|負の地域数|
|---|---:|---:|---:|---:|---:|
|**5種以上・25%以上**（v24.9基準）|8|22|272|**−0.103811**|7/8|
|**5種以上・50%以上**|6|15|207|**−0.117164**|6/6|
|**10種以上・25%以上**|5|13|199|**−0.113090**|5/5|
|**10種以上・50%以上**|5|13|192|**−0.106407**|5/5|
|**2015年以降、5種以上・25%以上**|7|20|244|**−0.105234**|6/7|

厳格10種+50%で残る地域コードは**AUT, CHE, IND, ITA, NOR**。熱帯・南半球地域の多くは消滅。INDは**最初の緩い基準では正+0.089**だったのに、厳格10+50%では**−0.0683**へ変わる。これは純粋に標本構成と形質欠測の品質フィルタ変更で符号が動いたものであり、道路がインドで生態学的方向を転換した証拠ではない。厳格化が推定対象の道路や地点を変えているので、`7/8`と`5/5`を独立反復と扱わない。

この結果は**全球道路の標準化したtransport corridor内部の形質的偏り**という仮説を探索的に支持する一方、広域の植物相・気候の変化に対する一般因果ではない。全球original 25地域中実際にテスト可能だった8地域、真に厳格な5地域だけで、**MIREN roads are not railways**。LEDAのヨーロッパ偏りと地形・樹木/草本の原地域プール違い、source MIRENのPlot1/Plot3の周囲微環境違いも残る。

## IV. 鉄道と道路をどう統合すべきか／何を混ぜないか

既存鉄道証拠は元4独立原著・5か国：チェコ39駅・イタリア28駅・クロアチア/ボスニア15駅・スロベニア89か所100m線路区間で、元鉄道内部の地理区・種数・種出現頻度固定帰無からのCWM地点間距離**正**偏差（CZ +0.0441, IT +0.0387, HR/BA +0.0747, SI +0.0650）を確認。ここでの分母は**既に鉄道にある種プールの再配置**。チェコ14真駅vs周辺10km内の別研究CzechVegの非交通文脈5plot合成poolの元距離差 **−0.4267** はまた別の分母。MIREN道路では同じ地点×同じ100m²のPlot1 vs Plot3で **−0.106（厳格）**。各研究でLEDA z軸を作る元species universeと調査areaが違うため、**−0.4267 vs −0.106を効果量として直接比較・統合しない**。内部帰無に比べた`正`と独立外部参照に比べた`負`が同居することも数学的に矛盾しない。

今後の共通解析に必要な比較単位：
`D_{region} = [beta_functional(infrastructure sites)] - [beta_functional(matched non-infrastructure sites)]`、元同調査年・同source調査面積/努力・同近傍地域・同形質z軸で整える。**本物のスロベニアBS/BT微生息地対比とMIREN道路Plot1/Plot3は物理的に同じ植生環境でない**。因果的な世界一般則は`rail/road×region×source study`階層効果と同一方法での本物の鉄道/非鉄道調査を要求する。列車による種子輸送・草刈りや風への個体適応を群集出現から推定しない。

**検証した次段として最優先**：著者公開の真正カシミール`Matrix_plot.xlsx`（Dryad ID546003、元source SHA256上記）を正規のアクセス許可のある環境でダウンロードし元と一致するか照合。31地点×4区×2014/17年の種×区行列が正確に復元できたら、(1)線路側vs内陸3区の**同法100m²ペア**、(2)2014→2017の追跡と地点別距離変化、(3)STN/BTN（16駅 vs 15駅間）という鉄道本来の連結位置に対する層を実解析し、MIREN全球道路の同法対比と比較できる。中国雲南84plot・フランス2025年394plotの元出現行列は著者問い合わせ待ちの候補であり、**配布されていると捏造しない**。

本研究ログ：`scripts/audit_India_Kashmir_railway_Dryad_original_files_v247.py`, `scripts/audit_India_Kashmir_Dryad_2014_2017_xlsx_source_structure_v248.py`, `scripts/analyze_global_MIREN_road_interior_LEDA_2D_functional_beta_v249.py`, `scripts/validate_MIREN_global_ROAD_2D_LEDA_coverage_sensitivity_v250.py`, GitHub Actions本文・集約jsonで再現できる。真正生物種×地点原行列や元GPSを再配布せず、主要統計と元MD5/SHAの検証だけを保存する。

MIRENデータ論文は研究計画の共有とMIREN Steering Committeeへの予定分析連絡を推奨している。公開data元利用時は原DOI/元region contributorsを正確に引用し、投稿の際は元共同著者方針を考慮すること。
