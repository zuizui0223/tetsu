# 鉄道植物：中心主張と実証証拠の対応表

**2026-11-20の更新：** 以下は11月15日時点の輸送・微生息地主軸の整理。現在の全結果の扱い・論文化方針は[全結果統合案](tetsu_all_results_publication_synthesis_ja.md)へ移した。過去の「固定ルール」は当時の編集範囲を示す履歴であり、後続成果の検討を禁止するものではない。

**主張固定日：2026-11-15｜研究ブランチの既存 Draft PR #1 に追記。新しい v 番号・別主張・未取得の結果は増やさない。**

## 本文の中心命題（支持範囲を守る表現）

> **貨物鉄道は作物由来種子の物理的供給経路となり得る。一方、既存の野生植物群集データでは、鉄道経路距離に固有の種共有予測力よりも、線路内の微小生育地に伴う群集組成・生活史形質の違いのほうが明瞭に確認できる。**

「列車は野生植物を運ばない」「環境フィルターが直接の原因だ」とは断定しない。後者は **delivery → viability → establishment → community filtering** という異なる段階の未接続な証拠である。

## 証拠対応表（本文の行と原本を1対1で追跡）

| 束／対象段階 | 支える命題 | 元の観測・再解析で確認した数値 | データの単位・独立性 | 結論の上限／未同定点 | 根拠・原本への導線 |
|---|---|---|---|---|---|
| **1-A カナダ：貨物→軌道への入力** | 列車から作物由来穀粒が実際に落ちる | **134 km・19地点**で穀粒落下を捕捉。計1,979採集器×時点サンプル。**年間約110 tは著者のモデル外挿値** | 1原研究・19地点。1,979は独立鉄道区間数ではない | 種子の発芽率、同一種子からの現地定着、野生植物の輸送は未測定 | [Canada v13.0](Canada_direct_freight_train_grain_deposition_v130.md)、Gangadharan et al. 2017 DOI [10.1111/acv.12336](https://doi.org/10.1111/acv.12336) |
| **1-B スイス：貨物内混入種子→発芽** | 作物貨物に混ざるアブラナ種子には発芽能力がある | 輸入小麦に混入した種子**300播種→75実生→40 GM実生** | **2製粉所の2試料群**。300は独立輸送便ではない | 実験室発芽は線路への落下、定着、繁殖率そのものではない。GMと非GMの相対発芽優位も不明 | [Swiss v13.4–13.5](Swiss_wheat_seed_viability_and_route_monitoring_v134_v135.md)、Schulze et al. 2015 DOI [10.1007/s11356-015-4903-y](https://doi.org/10.1007/s11356-015-4903-y) |
| **1-C スイス：輸送結節点→野生化個体** | 港湾・駅で野生化作物が現存し、貨物供給と整合する | 鉄道関連**79区域**中**4区域**にGMアブラナ陽性。遺伝子検査は鉄道関連の**50/2,398個体** | 1国内調査。駅と港を含み、4区域と50個体は独立の貨物到着を表さない | 同一輸入種子→落下→実生→成熟株を追跡していない。港は船・トラック・列車の結節点 | [Swiss v13.1–13.3](Swiss_rail_freight_GM_feral_transition_2011_2020_v131_v133.md)、Schoenenberger & D'Andrea 2012 DOI [10.1186/2190-4715-24-23](https://doi.org/10.1186/2190-4715-24-23) |
| **1-D バーゼル：供給停止と継続** | 貨物の有無と野生化GM個体の動態に観察的整合性がある | 種子貨物取扱停止後のSt. Johann駅では**2016年以降GM非検出**、稼働するライン港ではGMを継続検出。**2013–2020年、2地点** | 2地点×8年。上記のスイス試料・現地調査と一部重複する | 除草・検査努力・既存埋土種子・港の船舶入力が交絡。停止の因果効果は未識別 | [Swiss v13.1–13.3](Swiss_rail_freight_GM_feral_transition_2011_2020_v131_v133.md)、Basel-Stadt official report 2020 |
| **2-A スロベニア：微小生育地→機能群** | 同一鉄道内で生活史・散布戦略と位置が対応する | **89地点**の同地点・同側の道床肩（BS）−法尻（BT）：厳密一年草比率 **+0.054703**、風散布戦略比率 **−0.050740**、植生被度差 **−35.18元資料単位** | **173左右側ペア→89地点→15原著鉄道区間**。89を独立国／173を独立実験としない | 実際の種子到着、刈取り・風速・現地の選択圧は測っていない。BS/BTの実測面積が同等か不明 | [Slovenia v10.3](Slovenia_real_RAIL_microhabitat_plant_functional_v102_v103.md)、[v10.4](Slovenia_real_RAIL_diaspore_wind_vs_annual_v104.md)、Šabić et al. 2026 DOI [10.1007/s44473-026-00155-z](https://doi.org/10.1007/s44473-026-00155-z) |
| **2-B フランス：線路内の異なる環境→種数差** | 線路と法面の固定面積区で植物群集が違う | RATPの**394区画×各2 m²**。線路上**2.4種/区**、法面**7.0種/区**。差は法面−線路 **+4.6種/区** | 線路**128区**、法面**266区**。**同じ地点の真の対応ペアではない**。元区画×種行列は未取得 | 密度差を貨物散布の欠如、刈取り因果、成長抑制の適応と読まない。元分散がないため効果量SEも未推定 | [France RATP v16.1](RATP2025_original_track_vs_embankment_staff_ALERT_evidence_v161_ja.md)、Rakotoson et al. 2025 DOI [10.1111/wre.70025](https://doi.org/10.1111/wre.70025) |
| **2-C 日本：微小生育地の種特異性** | 軌道・法面で報告種の構成が異なる | Abe et al. 2012の**109鉄道区間**の管理者報告。スギナは**軌道24/70報告（34.3%）**、法面**22/642報告** | 管理者による種×区分の報告頻度。固定面積の現地植生区画ではない。上位10種＋残余で打切りあり | スギナの生育適応度・種子散布・群集全体の線路選好率は推定不能。仏SIと共通のHedges gは作れない | [Japan Abe 2012](japan_rail_microhabitat_reanalysis_v10.md)、DOI [10.3719/weed.57.159](https://doi.org/10.3719/weed.57.159) |
| **反証-N1 スロベニア：野生群集の経路距離** | 道路距離・地理・形質を越える「線路が短いほど野生種が共有される」予測利得は確認できない | **89地点・15区間**、道路のみから鉄道経路比を追加した区間外MSE改善率 **−0.5371%（悪化）**。24感度仕様すべて負方向 | 1原著、元地点ペアは依存する。実OSM物理経路であり、貨物実運行・到着種子は未観測 | **検証したモデルの非支持**であり、野生植物の列車散布がゼロだという証明ではない | [SI v15.0–15.2](Slovenia_2024_original_OSM_rail_road_control_v150_v151_ja.md) |
| **反証-N2 チェコ：駅単位の経路距離** | 異なる野外調査単位でも、経路距離の追加予測は頑健でない | **39原著駅、代表点厳密照合35駅**。同元5路線区間留保の全植物相で鉄道経路の頑健な追加利得を支持せず。新帰化種の小さな改善は調整法で消える | 1原著の駅全体群集。植物実調査区GPSでなく後年の駅代表点、35駅間のペアは独立実験ではない | SIとの違いは**観測単位の非互換性**。チェコの失敗を微小生育地フィルターの反証とも再現性確認とも呼ばない | [CZ v15.3–15.4](Czech_original_39station_SR70_geography_falsification_v153_v154_ja.md)、[v16.0](Czech_original_arch_neo_real_2023_railroad_5stretch_CV_v160_ja.md) |

## 図と論文の主張を対応させる

**Figure 1 — 供給の段階別証拠（輸送経路、発芽、野生化、長期検出）**  
Canada と Swiss／Basel の著者原研究を「1研究の複数独立実験」に水増しせず、原観測単位を明記する。これらは段階別の証拠であって、連結された同一種子の生存率推定ではない。

**Figure 2 — 線路内フィルター（独立国別結果を並べて比較）**  
SI：同一地点BS−BTの一年草・風散布型の構成差。France：同面積2m²区の線路 vs 法面の種数差（非対応区）。Japan：管理報告内の軌道・法面の種別頻度。*同じ生態学的質問でも観測単位と効果量は現時点で異なるため、1つの統合Hedges gや「3研究で同方向」とは言わない*。仏の元plot分散／plot×種行列、日本の元区間別種×区分表の取得までは設計別に並置する。

**Figure 3 — 経路散布の予測的反証と作物／野生種の署名**  
SI/CZで「道路・地理を越える物理鉄道経路距離の独立予測利得が支持されなかった」という**限定的否定結果**を示す。FAF6の2022作物貨物量とトラック対照を同じ地理単位で野生化作物 vs 野生種に適用する比較は今後の決定的テスト候補。現行US FAF6 v8.9は鉄道正相関を確認していない（[v8.9 source result](US_FAF6_directional_state_FE_node_label_null_v89.md)）。**未読・未成功の検定をFigure 3の肯定結果に追加しない。**

## 本文の3段落に対応する文章（主張の上限）

1. **Delivery is demonstrable for cargo-associated crop propagules.** Real train-kernel fallout in Canada and viable contaminant oilseed-rape seeds in Swiss imported wheat establish distinct physical and germination stages. GM feral plants at freight hubs and site-specific monitoring dynamics are compatible with repeated cargo-associated introduction, without tracing each seed through establishment and reproduction.
2. **Railway vegetation is structured at the microhabitat scale.** Slovenia's directly paired ballast shoulder/toe contrasts show higher annual fractions but lower measured wind-diaspore fractions close to tracks; French fixed-area plots show lower richness on track than embankment; Japanese railway management reports independently show taxon-specific occupation of track versus slope. These results indicate within-corridor compositional sorting, but do not identify a common causal management or adaptive mechanism.
3. **Physical rail-network proximity is not a sufficient proxy for seed delivery or establishment.** In the available SI and CZ original floras, adding rail route length over geographic/road controls does not provide a robust extra prediction of shared flora. These are source- and unit-specific negative predictive results, not proof of absence of train-mediated dispersal by wild plants.

## 論文の固定ルール

- 中心は **貨物由来作物種子の散布 × 生育地フィルターの段階分解**。種内適応／可塑性（D）、写真での形質進化は本文の直接検証から外す。
- 広域GBIFマクロ、沿線リモートセンシング、FAA/FAF6の未確証の植物到達効果は補遺／将来の検定候補。FAF6貨物量＝生きた種子の輸送数ではない。
- フランスRATPの非対応2m²区、スロベニア同地点のBS−BT、日本の区間管理報告を、統一されたpaired plotの「メタ解析済み効果量」と呼ばない。
- 2012年フランスPenone論文・博士論文は別原著として加算しない。スイス同地点・同年の遺伝子研究と行政監視も独立反復として合算しない。
- **実際に書ける結論：輸送の直接観測は貨物作物、群集の実測構造は線路内微小生育地。まだ書けない結論：野生植物は列車で運ばれない／フィルターが唯一の原因／鉄道で適応進化が起きた。**
- 新しい解析版番号や新たな仮説探索を本文の固定主張に持ち込まない。続く本文作業は「対応表に列挙した原著を使った設計別メタ解析可能性の監査」のみ。
