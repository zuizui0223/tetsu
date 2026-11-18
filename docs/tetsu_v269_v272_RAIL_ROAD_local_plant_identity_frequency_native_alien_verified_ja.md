# tetsu v26.9–v27.2：鉄道・道路の局所植物相と種の置換、在来・外来の条件付き比較

更新 2026-11-18 JST。鉄道・道路の真正原データ2研究から新たに再計算した探索的マクロ解析。研究ブランチのみで、mainは未変更。原種×区画行列、原観測座標、個別学名の出現記録を再配布しない。

## 1. 一次資料とソース保証

Kashmir RAIL: Rashid et al. (2021) JAE DOI https://doi.org/10.1111/1365-2664.13961, Dryad https://doi.org/10.5061/dryad.0p2ngf20h version 99735。作者真正 Matrix_plot.xlsx 122,187 bytes / SHA256 529f3b97521933018e12623b61f531f012c67241cae2602a6675dbd534192474。元31地点・2014/2017年・4×100㎡ zone。欠損4セルは不在0へ置き換えず、全zone×両年で完了する駅15・駅間12地点のみ採用。

MIREN mountain ROADS: 元MIREN consortium公開 https://doi.org/10.5281/zenodo.15749536 と2026 GEB data paper https://doi.org/10.1111/geb.70284。原CSV v3 17,899,241 bytes、167,804出現行、MD5 b640195cc44f252ccb6ddb28eb6fff13。原Region/Road/Site/Year/Plotコードを真の調査区単位とする。同年のPlot1道路際 vs Plot3内陸が完全なsiteを、初回年だけ採用。初期2007/08年に外来種中心だった原プロトコルを考慮し、主解析は初回2012年以降の16Regionコード・48道路・669地点。各plot 100㎡。元Regionコードを独立した16か国・16論文と数えない。

成功した作者原データ照合：
- [v26.9 Actions #44247105966](https://github.com/zuizui0223/tetsu/actions/runs/44247105966) - 数式テスト7件成功
- [v27.0 Actions #44247682744](https://github.com/zuizui0223/tetsu/actions/runs/44247682744) - source fixed-frequency nullと数学テスト6件成功
- [v27.1 Actions #44248209721](https://github.com/zuizui0223/tetsu/actions/runs/44248209721) - 調査年固定・独立乱数系列検証4件成功
- [v27.2 Actions #44248849286](https://github.com/zuizui0223/tetsu/actions/runs/44248849286) - Region別在来外来・同一site支持域・QA5件成功

## 2. v26.9：元現場同地点で「単なる脱落」より種の置換が大きい

同site・同調査年のcorridor集合とreference集合に共通a、corridor-only b、reference-only cを数える。原Jaccard非類似度 βjac=(b+c)/(a+b+c)、置換成分 βjtu=2min(b,c)/(a+2min(b,c))、nestedness-resultant βjne=βjac−βjtu。これは同時点の分類群差の数式成分であり、時間的な実置換や種子移動の実測ではない。

| 同一site 100㎡比較 | 原地点数 | βjac | βjtu | βjne | corridor-only taxa/site |
|---|---:|---:|---:|---:|---:|
| MIREN ROAD 2012+ Plot1 vs Plot3 | 16地域・48道路・669地点 | 0.71497 | 0.59877 | 0.11620 | 17.32 |
| MIREN 2017+ first year subset | 8地域・25道路・372地点 | 0.72290 | 0.62455 | 0.09835 | 14.25 |
| Kashmir RAIL Station 2014 Track vs Interior2 | 15地点 | 0.92370 | 0.89752 | 0.02618 | 6.33 |
| Kashmir 同駅2017 | 同15地点 | 0.94286 | 0.84000 | 0.10286 | 2.47 |
| Kashmir Between 2014 Track vs Interior2 | 12地点 | 0.91053 | 0.79716 | 0.11337 | 3.00 |

原MIRENの全16地域でβjtu>βjne、Kashmir駅でもβjtuが大きい。局所環境の植物は周辺の単なる完全な部分集合ではない。

数学的な停止条件：同じ元siteの両plotの種数を意図的に等化するとexclusive数b=cとなり、必ずβjne=0となる。79回のsource species rarefactionで「置換100%」となったのは定義上の自動結果なので、**生態の100%置換を示す証拠として使わない**。種数差への調整には元site richnessとtaxon出現頻度を保持する別帰無が必要。

## 3. v27.0：siteごとの種数と各taxonのROAD内出現頻度を厳密に固定

真のcorridor側の種集合は固定。reference側のsite×species incidenceだけを2×2 checkerboard swapsでランダム化する。各元siteのreference種数を固定し、各元Roadで各種のreference頻度を厳密に固定したままsite間の対応だけを崩す。MIRENは元Road内、Kashmirはstation/between層別×年内で同じ方法を適用。希薄化と違いβjneを機械的にゼロにしない。観測した同地点隣接植生との種共有が、元Roadの種頻度だけで生まれるかを問う。Markov連鎖drawは相関し、順位は独立野外試験p値ではない。

| MIREN 2012+、16地域・48道路・669地点 | 観測 | 元source頻度固定帰無 | 観測−帰無 |
|---|---:|---:|---:|
| βjac | 0.71497 | 0.83224 | **−0.11727** |
| 置換βjtu | 0.59877 | 0.77459 | **−0.17582** |
| nestednessβjne | 0.11620 | 0.05765 | +0.05855 |
| 同site共有taxa | 10.6259 | 6.7145 | **+3.9114** |
| corridor-only taxa | 17.3214 | 21.2328 | −3.9114 |

**全16Regionコードで、同siteのβjacは元source species頻度のみから組んだ帰無より低い**。したがって、生態解釈上「回廊に固有の種が多い」と「地元周辺植生との共有も予想以上に強い」は両立する。MIRENの同一siteの微地形・気候・種供給・管理と地理距離が共変しているので、局所種子供給の因果効果と確定しない。

Kashmir RAILのstation Track vs Interior2 βjac差（観測−帰無）は2014年 −0.01229、2017年 −0.00379。between層では2014年 −0.02546、2017年 +0.00248。RAILでは統計的な強い普遍的なlocal couplingとは扱わない。**MIRENのnullは元Road内、Kashmirは広いstation/between層別**で、空間条件の等価性はない。異研究間の差から鉄道vs道路の因果効果は計算しない。

## 4. v27.1：独立乱数系列と本物の初回調査年を固定した感度

MIRENの原48Roadでは採用元Siteの初回調査年が各Road内でそろっていた。したがって元ROAD×初回年の層を作っても対象は16Region/48Road/669Site。v27.0とは独立のseedおよびMCMC stepを用いた帰無βjacは0.83210（元v27.0は0.83224）、初回調査年で層別した帰無は0.83194。観測βjacは0.71497で不変。**16/16地域が観測<頻度帰無**を維持した。48元Road中1Roadではcheckerboard swapsを受理できず、出現頻度条件付き交換可能性に情報がないことを明記する。

## 5. v27.2：在来と外来の両方が、局所種共有と回廊固有組立てに寄与

MIRENの元Status Native/Alien/Unknownを作者原Region×taxonの分類に限って使用。原sourceの全167,804行を再走査して4つのRegion×taxon組でStatus区分の不一致を発見（AlienとNativeが直接食い違う1組、およびUnknown vs既知の3組）。全4つをUNKNOWNに保守的処理し、何も外来・在来へ捏造しない。以前のv9.9の「2つの矛盾」という説明は分類条件が異なり、今回の定義ではUnknown対既知の不一致も算入するため4組。Statusは地域別在来/外来区分であり、個々の外来植物の侵略性や具体的な車両輸送を表さない。

| 本当のStatus source plotペア | 元Region/Road/Site | βjac観測 | βjac源頻度固定帰無 | 観測−帰無 | 正しい方向のRegion |
|---|---|---:|---:|---:|---:|
| Nativeのみ | 15地域・43道路・631地点 | 0.72938 | 0.85789 | **−0.12851** | 15/15 |
| Alienのみ | 11地域・33道路・314地点 | 0.70005 | 0.78194 | **−0.08189** | 11/11 |
| 両Statusが同siteの両plotにある**同じ277地点** Native | 10地域・28道路・277地点 | 0.72628 | 0.85187 | **−0.12559** | 10/10 |
| 同じ277地点 Alien | 同上 | 0.72754 | 0.81575 | **−0.08821** | 10/10 |

同じ277地点・10Regionの同じsource sampling supportで、Nativeの観測−nullはAlienより負の方向が強いRegionが9/10。地域等重みのnative-minus-alienの超過非類似度差 = −0.03738。しかしNative/Alien間でtaxon母数・richness・頻度・生育場所が違うためこの差を「本当の選択確率の差」としない。

Source corridor-only taxa vs interior-only taxaに対する元Region statusから、両方の外来割合が計算できる**16Region・46ROAD・628site**に限定し、site→road→regionの順に等重みで算出。外来割合は corridor-only **36.33%**、interior-only **15.64%**、差 **+20.68ポイント**（15/16地域でcorridor-onlyの方が外来割合が高い）。これは元長期MIRENのalien-roadside enrichmentと整合するが、車両による種子運搬の直接検証ではない。

**総合的な説明候補**：交通回廊では、外来種も含む独自の植物群が増えるのと同時に、同じ場所の在来植物相との地域固有の結びつきも維持される。これは同じ人工基質が種同一性を全球一律に決定しないことと整合する。単純な環境フィルターのみ、単純な種子流入のみ、単純な生息地周辺の部分集合のみ、のどれか一つに全群集結果を還元しない。

## 6. 論文の射程と識別限界

この新しい結論の大部分は**MIRENの1全球道路データセットの内部地域反復**である。ROADの16Regionは16独立論文ではなく、RAILのKashmirは独立鉄道野外研究1件の縦断資料。今回のzone距離・地域植物相・現地設計・調査年と面積の一部が揃っていても全体を同一の因果処理にしてはならない。列車の運行頻度、主線/支線、列車の種子搬送、現地剪定、風圧ストレス、個体の生残・適応遺伝子・進化は調べていない。

**次の真のマクロ研究に必要なデータ**：同じ季節・近傍のrailway/roadside/荒地の元plant × site incidence、同面積の野外調査、独立した複数地域・研究、路線接続性とactual traffic、最低限site-level native/alien status。先行研究の地理・資料粒度を保持してleave-study-outを試し、同じ源が多い/少ないと種の置換と機能的収斂がどう対応するかを検証する。

再現コードと原本QA:
- v26.9 scripts/analyze_original_RAIL_vs_global_ROAD_taxonomic_microhabitat_beta_v269.py / tests/test_original_RAIL_ROAD_taxonomic_beta_v269.py
- v27.0 scripts/audit_RAIL_ROAD_source_frequency_fixed_reference_null_v270.py / tests/test_RAIL_ROAD_fixed_source_frequency_null_v270.py
- v27.1 scripts/audit_MIREN_year_stratified_frequency_null_v271.py / tests/test_MIREN_year_stratified_reference_null_v271.py
- v27.2 scripts/analyze_MIREN_native_alien_road_local_reference_frequency_null_v272.py / tests/test_MIREN_native_alien_original_source_frequency_v272.py

Corresponding GitHub Actions workflows tetsu-v269/v270/v271/v272 all passed after a source QA fix in v27.1 that counted number of year strata, not number of Roads. MIREN's source data users should respect the consortium's study-planning and contributor involvement guidelines before publication.
