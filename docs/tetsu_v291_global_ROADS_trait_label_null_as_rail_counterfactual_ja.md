# tetsu v29.1 — 鉄道4研究の「草丈の空間的形質配置」は山岳道路でも生じるか

2026-11-19 JST、原著source CI実測完走前は新たな正負効果を主張しない。
作業ブランチ \`research/tetsu-MIREN-road-trait-shuffle-counterfactual-v291-20261119\`；\`main\`は変更しない。

## 研究の主張を機構と別に固定

鉄道v28.5では、欧州の**4独立出版原著**（CZ39, Italy28, HR+BA15, SI89 100m道床site）すべてで、現実のsite×species出現構造を固定して、同じ元publicationの植物種に付く Oldenburg LEDA 草丈・種子重の **paired 2D trait vector** だけをシャッフルしても、元source同地理群内の草丈CWM地点間距離 **観測−trait-label帰無** が正だった。これは鉄道に固有の選択圧か、より一般的な地域種プール×環境×種の形質構造なのか分からない。

MIREN山岳道路は独立の**実際の交通環境の参照研究プログラム**だが、同等の鉄道野外研究としてカウントしない。roads \`Plot 1\` と同一site・同年・同100m²の interior \`Plot 3\` を現実出現行列で比較する。同一MIREN原著・地域・種プールと気候差を持つ、真の道路環境内外の「形質種名対応」帰無を鉄道と近い定義で同じように使う。

### 元データと真正セキュリティ

- Seipel, Haider, MIREN Consortium, Zenodo v3 https://doi.org/10.5281/zenodo.15749536、\`MIREN_road_survey_(2007-2023)_zen_v3_002507.csv\`、原真正bytes 17,899,241、MD5 \`b640195cc44f252ccb6ddb28eb6fff13\`、真source 167,804種×plot出現レコード、source 26 Region codes（出版は25地域）。
- Oldenburg LEDA source canopy_height(m) SHA256 \`cbf5e761390d78a843512dfc1139acfc1154d46d6c1995a827354840a9f78a92\`、seed_mass(mg) SHA256 \`34aa12c333fc1ba6de2e560a495ac5d693aa0722a09ecf68418ad94ac1c19cef\`。各植物の二名法厳格一致、形質0/不明を補完しない。
- 先行v24.9/v25.0のsource受理条件を*変更しない*。基本>=5両形質既知taxaかつsource元speciesの25%以上をroadside/interiorの両方で満たす、1 original source roadに>=4真site、1 Regionに>=2本のsource road、各本来のRegion×Road×Siteで初回調査年のみ。**8 source Region codes,22 real roads,272 unique original sites**。厳格10種かつ50%以上は**5 codes,13 roads,192 sites**。元source v24.9の主指標 roadside−interior CWMβ = **−0.10381093795211047** をまず再現できなければ新帰無解析を停止。各感度パネルは同じ元研究の重複で、2つの独立publicationではない。
- z(log height, log seed)軸は v28.5と同じ**欧州CZ/IT/HRBA元622種で予め固定した基準**、mean (-0.5802365754696401,+0.2936113186962433), sd (1.3801000545693773,2.1657660559995824)。以前のMIREN自身で作ったz軸の差とは同じ係数としてプールしない。

### 帰無モデルの真の検定対象

同じRegion内の真source \`(Road, Site, first Year, Plot1 / Plot3, binomial species)\` の全incidence・species prevalence・co-occurrence・richness・現地地点を完全に固定。比較するのは**元Road単位で同じ道路上の異なるsource site間**の植物trait CWM2D/height/seed distance（元各roadの平均→Region内road等重み）。

処理は2モード×2seed×各99 technical assignments:
- A. 各Regionの両Plotに現れた本物のtrait-known species集合で、元speciesに割り当てる\`(zlog height,zlog seed)\`ペアをランダムに再配置。MIREN各地のtrait分布を保存し、植物の真地理位置・同所的共存をすべて固定。
- B. さらにRegion内の元source taxon occupancyを5分位に固定し、各quintile内でのみペアvectorを交換。元source種出現頻度×trait平均の粗い関係を制御する。
- 2技術seedは2つの野外調査ではない。順位は形質帰無下の診断で、野外field inference p値ではない。植物群集データの既知空間自己相関・種間系統・trait missingness・元survey year・本来のsample abundanceを完全制御するわけではない。
- outcome: \`observed-minus-trait-label-null\`（2D草丈＋種子重、height 1D、seed mass 1D）を、同じ元source sitesのroadsideおよびinteriorで別々に出力する。さらに両者の差分を出すことで、元roadsideの形質構造が単に地域の植生全般を反映するだけなのかを検定しようとする。
- どのRegionも独立したpublicationやcountryとは呼ばない。Region等重み・leave-one-region-outは記述的な堅牢性であり、全球的予測区間ではない。

### 事前に固定した科学的解釈

| 元検定 outcome | 対立する含意 |
|---|---|
| Roadside 正 / Interior も同程度に正 | **通常の地理的群集形成**でも鉄道v28.5と同様の正の高さ効果が生じる。鉄道固有の選択だという仮説は弱まる。 |
| Roadside 正 / Interior ≈ 0 | 道路沿い環境でextra functional sorting候補。**鉄道限定ではない**。 |
| Roadside 負 / Interior 正 | 道路が種プールによる通常の形質分化を相対的に抑制する可能性。 |
| Roadside とInteriorとも正、だが道路の実CWMβは内陸より小 | 「固定真種群集内の非ランダム形質配置」と「外部対照に対する相対的収斂」は**両立**する。これはv28.5対v28.6の概念的緊張と同じ二階層構造。 |

単なる草刈り・高速車の風圧・列車による散布・植物遺伝適応の実観測はないので、肯定結果は形質と群集の空間構造という観察結果に限定する。MIRENでpositiveになっても「鉄道と道路の同一的環境フィルター」とは断定できない。道路plot1とplot3は物理的に異なり、experimental random assignmentではない。

### 現行原著ソース入力・出力

\`scripts/validate_MIREN_ROAD_original_site_frozen_trait_label_counterfactual_v291.py\`
\`tests/test_MIREN_original_ROAD_trait_label_source_placebo_v291.py\`
\`.github/workflows/tetsu-v291-original-MIREN-road-plot-frozen-trait-label-null.yml\`

公開結果はRegionコードごとの数値と真サンプル数だけ。真元species×plot taxon identity、元lat/lon、LEDA rawテーブル、元survey CSVをGitHub artifactへ再配布しない。


## v29.1 新しい原著source実測完走結果（2026-11-19）

**[SOURCE SUCCESS #44462368280](https://github.com/zuizui0223/tetsu/actions/runs/44462368280)**、[公開派生source artifact #14500352634](https://github.com/zuizui0223/tetsu/actions/runs/44462368280/artifacts/14500352634)。元Zenodo CSV 17,899,241 bytes MD5、真正167,804出現行、Oldenburg LEDA height+seed SHA、数学テスト8/8完走。旧v24.9の実元MIREN道路Plot1−Plot3 2D CWMβ **−0.10381093795211047** をコードから元source再計算して1e−9以内で照合後にのみ新trait-label shuffleを実行した。既存原植生は一切移動しない。

主解析は事前定義されたregion内の本物のsiteでpaired(height,seed) traitラベルを同Region出現頻度5分位内で入替え、**2固定seed×299=598 technical null assignments**。帰無はRegion内の元植物種プール・site×species・各種の源頻度・元site richnessを変えない。標準化座標は以前欧州3原著622 taxaで固定したもの。元空間field unitは依然として基本**8 region codes/22 roads/272 sites**、厳格**5/13/192**、実在資料の複数region・道路を数十独立出版とはカウントしない。

### 1．原source真の植物出現配置を固定した「観測−trait-label帰無」の新しい効果量

| 元MIREN道路、源region等重み、2D(草丈+種子重) trait-label occupancy-5分位shuffle | 道路際Plot 1 | 内陸Plot 3 | **道路際−内陸** |
|---|---:|---:|---:|
| 基本 5種25%、8 regions 22 roads 272 sites | **+0.097816**（7/8地域正） | **+0.146879**（7/8正） | **−0.049062**（6/8地域で負） |
| 厳格 10種50%、5 regions 13 roads 192 sites | **+0.137866**（5/5正） | **+0.201423**（5/5正） | **−0.063557**（4/5で負） |

| 草丈のみ、同じ元Europe-fixed z(log canopy height) | 道路際Plot 1 | 内陸Plot 3 | **道路際−内陸** |
|---|---:|---:|---:|
| 基本5種25% | **+0.073089**（6/8地域正） | **+0.136309**（7/8正） | **−0.063220**（7/8で負） |
| 厳格10種50% | **+0.116046**（5/5正） | **+0.167266**（5/5正） | **−0.051220**（4/5で負） |

主解析の**同じsource真サイトと同じEurope z軸**で草丈の生観測Roadside−Interior β距離は **−0.0931493207367427**。元species分布と元occupancyを保つ形質帰無で期待されたRoadside−Interior距離差は **−0.02992917507743902**。したがって元真の差−帰無期待は**−0.06322014565930367**。これは単なる「道路際の方が元植物種数が違うので地理的に類似して見える」現象のみではなく、植物species→trait値の特定の対応による追加成分がある、という観察学的結果。

**Region元独立leave-one-region-out descriptive mean** の草丈超過差は基本 **−0.085625〜−0.049029**、厳格 **−0.064927〜−0.029963** と全て負。これはMIRENの元地域コード分布に依存する同source繰り返し感度で、95% field CIではない。基本sourceのIND地域のみ道路際−内陸の高さtrait-null差が **+0.093616** と逆向き（他7地域負）。厳格ではIND **+0.003606** とほぼゼロ、他4負。一般法則の符号が例外なく全地域に成立したとはしない。

Unrestricted nullでも草丈のRoadside−Interior超過差は基本 **−0.060200**、厳格 **−0.053264**、いずれも負。帰無指定一つだけの恣意的な符号ではない。新たな技術draw598回の上側順位を真の野外field p値に見せない。

**種子重のみは基本panelでは不均一**：草丈とは違い、基本Roadside−Interior excess **−0.007312** は負4/8のみで、地域leave-one-outの平均値範囲 **−0.018995〜+0.005921** とゼロをまたぐ。厳格サブセットでは **−0.038453**、4/5地域負だが、被覆率で地域sample setが変化している。種子重より草丈が特に再現的という現時点の探索的メッセージは、源sourceと地域制限を必ず添える。

### 2．重要な生態学的解釈

- MIREN山岳道路の**道路際も内陸植生も、元種出現配置を完全に保持したtrait-label nullより草丈の空間的分化が強い**。したがって鉄道v28.5の欧州4原著「草丈が地理に沿い非ランダムに配置する」という正方向の観察は、**rail substrate固有の適応証拠にならない**。地域種プール・地形・環境傾度でも同じ数理的構造が現れる可能性が、真の道路＋内陸別データで示された。
- さらに「内陸も正方向だが、道路際の追加的な高さの空間分化は内陸より**小さい**」。これは**非ランダムな形質の地域配置と、インフラ沿いの相対的な機能収斂が両立する**ことを、同じ面積、同じ年、同じ現地source道路site集合で示す。もちろん道路沿い/内陸に無作為処理したわけでなく、地形/種プール/森林の違いが残るので因果選択ではない。
- 道路の8 Regionコードは真の25出版地域を代表する無作為分布ではない。LEDA二形質既知は元MIREN 5,739名中1,096名（19.1%）で欧州・系統偏りが大きい。厳格5地域は主に温帯・北半球。ここから「全球どこでも同じ大きさで道路特有の選択が働く」と一般化しない。
- 欧州4鉄道原著のv28.5 observed−label-shuffle trait excessを、MIREN1道路原sourceのデータと合わせて**5つの独立鉄道研究としてメタ平均しない**。これらは近い概念のcounterfactualだが、国/地域構成・本来のfield grain・地区の組み方が違う。次の機構推論には同source同面積で物理条件を確認できる**実railバラストとopen/gravel非rail植生**が少なくとも2以上の独立野外調査で必要。

### 3．科学的結論

**「人工回廊では、地域由来のtrait–taxon空間構造が残ったまま、その分化の大きさが周辺植生に比べて縮小しうる」** を、山岳道路という別の野外基盤から支持した。ただし鉄道の全球一般性の肯定ではない。

GitHubの派生再利用用の固定値は \`data/derived/tetsu_v291_MIREN_road_vs_interior_trait_label_null_aggregates.json\` に保存。**元source artifactを正本**とし、本文表は表示用丸めである。
