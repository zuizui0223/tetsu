# 初期研究系列の全件索引と論文化上の扱い（v0–v17.6）

確認日：2026-11-20。対象は現在のブランチに存在する初期報告・設計・監査文書。ここでは原データを新たに再解析していない。各文書の本文を読み込んで見出し・結果表・結論を照合し、主要系列では後続の反証・感度結果まで追跡した。全数値を原データから再現した監査ではない。下記の一覧には設計だけの文書も残す。実証結果、データ取得成功、数学的実装、既刊文献の再集計は異なる成果である。

## 全結果をまとめるための分類

初期系列だけでも、(1)輸送網の形状・輸送量、(2)観測された植物組成、(3)局所の生育地選択、(4)種・形質による予測、(5)形態・開花・進化、(6)測定・観測過程、の六段階がある。同じ「肯定」で一列に並べない。複数のバージョンは追加の独立研究ではなく、同じ問いの更新として一行の証拠系列にまとめる。先行研究の種子実験や進化実験を本リポの新発見に含めない。

| 結果系列 | 残す結果と単位 | 反証・条件・更新 | 論文での扱い |
|---|---|---|---|
| 日本写真の初期監査 | 8分類群×3都市の注釈189件中、非開花・非結実の明示は4件のみ。写真・位置品質を満たす132記録 | 注釈欠測を非開花に変えられず、写真の尺度・全個体像も不足 | データ不足の補足表。開花効果の結果にしない。 [初期監査](real_photo_feasibility_2026-11-07.md)、[ROI](real_photo_roi_review_2026-11-07.md) |
| 東京鉄道–道路写真比較 | 地理的な写真記録分布の探索的比較を実行 | 初回重みの誤りによる0.1159は訂正済み旧値。観測機会と植物生態を分離 | 訂正履歴・測定監査。訂正前後を2結果として数えない。 [結果](real_rail_vs_road_tokyo_2026-11-07.md) |
| USA-NPN開花 | MarylandのCercisは近鉄道で初回YES区間中点が−4.75、−7日、ArizonaのProsopisは+4、+0.75、+7.75日 | 各年近側2地点、区間幅に比べ効果が小さい。植物個体の鉄道敷内生育ではなく観測地点中心の距離。GDD差も年で変動 | 相反する探索記述を同じ表に。統合した「早期開花効果」は不可。 [中点](npn_observed_onset_midpoint_v08.md)、[気温](npn_daymet_thermal_control_v09.md) |
| 日本の既刊鉄道微生育地 | 温暖・寒冷地域の報告頻度合計1,334。スギナの道床報告割合24/70、法面22/642、割合比10.01 | 調査区内実個体数や同面積出現率ではなく、既刊の種×生育地報告頻度 | 局所的生育地差の先行証拠／原表二次集計。 [原表再集計](japan_rail_microhabitat_reanalysis_v10.md)、[管理史](japan_rail_mowing_historical_baseline.md) |
| 宇都宮LRT衛星 | 新規開業の年代と衛星計測を確認、近傍と対照のNDVI差に負方向の例 | 校正・季節・都市造成・土地被覆・道路混在。20mオフセットでも68点中60点が道路25m内、厳格対照1–2組 | 生物効果ではなく空間支持域の監査。 [NDVI](utsunomiya_s2_measured_ndvi_v12.md)、[20m](lrt_true_rail_verge_20m_v14.md) |
| 通常鉄道・道路独立対照 | 370候補で道路から離すほど樹林割合39.7%→51.0%→73.7%→90.0% | 道路交絡を減らす操作で対象生育地が変わる。10m/20m画素は鉄道植生の確証ではない | 比較設計による対象集団の変化の具体例。 [370地点](japan_worldcover_370_site_posthoc_discovery_v16.md) |
| 初期10都市・GBIF近接 | 10都市で51万件超の植物記録と鉄道・道路形状を取得。後期外来種の接続差は5都市中3正、中央値+0.000328 | 両地点に鉄道がある対照へ限定すると2都市しか残らず両方負。国別外来リスト不在・不足あり | 取得量は実証規模と分ける。都市結果の更新前履歴。 [v1.9](overture_gbif_10country_macro_v19.md)、[v2.5](rail_griis_macro_results_v25.md) |
| 50km格子のネットワーク | 10ブロック8国160格子。鉄道のある格子同士で接続／非接続の対照が消失 | 規模拡大しても二値接続は識別できず、連続距離へ移行 | 無効なゼロ効果として扱わず、支持域不足の結果。 [v3.0](50km_macro_rail_opportunity_v30.md) |
| 6地域の鉄道–道路相対距離 | 6国各1地域、後期の未標準化Jaccardでは仮説方向5/6、K20記録希釈では3/6、中央値+0.021259 | K80は事後選択で5/6になったが年代・座標アンカー・提供元除外で不安定。フランスK80は提供元除外−0.059143→+0.101852 | 仮説に合った初期結果と最終的な不安定性を併記。全球の鉄道均質化を主張しない。 [6国](6country_macro_rail_circuit_flora_results_v46.md)、[提供元](GBIF_GRIIS_dataset_leaveout_macro_results_v57.md) |
| 回路理論・迂回経路 | 実鉄道グラフの代替経路構造を定量化、道路対照・合成グラフ検証 | 構造的冗長性は実運行・種子流量でない。植物効果は一貫しない | 必要なら方法補足。回路計算の成功を生態的肯定にしない。 [v3.8](circuit_macro_results_v38.md)、[対照](circuit_rail_vs_road_source_validation_v41.md) |
| 木本／草本の広域差 | 初期の木本相対距離相関は5/5仮説方向、中央値−0.144935 | 同じ4地域・格子で最大提供元を除くと木本4/4→2/4。草本も変化 | 機能群でも観測バイアスは消えないという結果。 [v7.2](GIFT_growthform_provider_leaveout_macro_v72.md) |
| フランス原イベントによる確認 | 1地域の原イベント・座標誤差≤5km・木本K20でrho−0.208302、16格子除外すべて負 | 同じ標本のraw Jaccardは+0.049118。単一提供元・完全な地域調査でない | 単地域の再計算成功を国際再現と数えない。 [v7.3](FR_original_GBIF_events_growthform_rail_macro_v73.md) |
| GBIFのturnover／nestedness | 木本のraw betaは提供元除外後4/4で仮説方向、中央値+0.141035 | 同条件K20 betaは2/4、中央値+0.024878。turnoverとnestednessの寄与も不安定 | rawの肯定を保持し、指標依存性を本文か補足で明示。 [v7.4](sixcountry_beta_turnover_nestedness_rail_v74.md) |
| 米国列車本数と木本 | Midwestの正の列車本数がある13格子で木本K20 rho−0.225759 | SoutheastとPacificは厳格支持不足、廃止踏切割合や観測量も関係。1地域の条件付き相関 | 列車圧の因果効果として使わず、運行指標の探索系列。 [v7.5–7.7](US_FRA_train_traffic_vs_GRIIS_GIFT_macro_v75_v77.md)、[対照](US_FRA_woody_observation_negative_controls_v78.md) |
| 米国貨物網の形状 | 農産物の州外輸送割合は鉄道63.43%、トラック17.76%。共通32出発州の29州で鉄道の有効到着州数が少ない | 経済ODの推定トン数であり、通過線路や散布種子ではない | 広域輸送機会の実証。種供給の仮説背景として残す。 [v8.0](FAF6_agricultural_OD_rail_vs_truck_network_v80.md) |
| 貨物ODと植物類似性 | 正の鉄道ODに限定すると小さな正相関、作物・資源両方正の41ペアでは作物rho+0.284174 | 全州組では負で49/49州除外でも負。raw応答・貨物限定で変化。州ペアは独立でない | 有利な限定集合だけで輸送効果を主張しない。 [v8.0–8.3](US_state_FAF6_crop_OD_vs_GRIIS_flora_v80_v83.md)、[資源対照](US_FAF6_crop_resource_rail_floristic_specificity_v84.md) |
| 貨物方向と後期のみの記録 | 前期起点種が後期終点に新たに記録される割合を検証 | 正鉄道148方向K20 rho−0.024113、両交通正123方向−0.091896。州固定効果・地域制約帰無も一般的鉄道肯定なし | 真の初侵入ではなく「新規記録」。時系列の反証として採用。 [v8.5–8.8](US_FAF6_directional_temporal_GRIIS_results_v85_v88.md)、[固定効果](US_FAF6_directional_state_FE_node_label_null_v89.md) |
| 観測過程の広域構造 | 同じ観察者×50km地点の再出現率と後期のみの種の割合は、48州除外すべて負、rho−0.456612～−0.347959 | 絶対新規記録数ならrho+0.072024。観察者交代が生物侵入を起こしたとはいえない | 観測過程の方法論的結果。輸送生態とは別の副論文候補。 [v9.5–9.7](US48_GBIF_iNaturalist_observer_site_recording_process_v95_v97.md) |
| MIREN道路沿い外来種数 | 同年同地点100m²道路際–内陸で24地域コード中23正。20地点以上22コードの中央値+2.155953種 | 26コードと著者25山地地域は同一でない。地域コード・地点・年を独立研究数にしない | 広域で再現する道路側の外来種富化。既刊研究と本再集計の新規性を区別。 [v9.8–9.9](MIREN_real_roadside_alien_flora_macro_v98_v99.md) |
| MIREN外来種交代 | 両側外来種あり583地点でturnoverがJaccard非類似度合計の62.67%。10地点以上17地域中15でturnover優勢 | 1,005地点中422は両側陽性でなく分解対象外。相対新規増加は道路−内陸−0.036670、11地域中4正 | 単なる内陸の部分集合でないことと、増加率の非普遍性を並記。 [v10.0–10.1](MIREN_standardized_roadside_alien_species_replacement_v100_v101.md) |
| Slovenia近線路微生育地 | 89地点173左右対応、15区間。BS−BTの一年草割合+5.47pp、被度尺度−35.18、風散布割合は低い | ゾーン面積同一の保証なし。外来割合も−1.90pp。単純な「過酷な場所ほど風散布・外来」は支持せず | 直接の群集内対照として重要。管理や列車風の原因同定は別。 [v10.2–10.3](Slovenia_real_RAIL_microhabitat_plant_functional_v102_v103.md)、[分解](Slovenia_rail_annual_mix_vs_within_lifespan_wind_v106.md) |
| Slovenia頻度による符号反転 | 同じ261種で風散布群のBS選好：種等重み+0.189388、出現機会重み−0.098125。15区間個別除外で両符号維持 | 加重側の区間・属ブートストラップ区間は0を含む。差の区間反転検定はp=.002075、5比較補正.010376だが各効果の有意性とは違う | 独立した推定対象の違いとして論文の有力な柱。普遍的な風散布選択ではない。 [v10.9–11.0](Slovenia_species_frequency_weighting_reversal_v109_v110.md)、[帰無](Slovenia_RAIL_section_flip_and_strategy_placebos_v115_v116.md) |
| CzechVeg独立対照 | 同論文・地域・年代・面積・植生クラス105鉄道文脈／167道路文脈、18原研究で種等重み−0.002921、頻度重み−0.023405 | 重み差区間[−0.002863,+0.045664]。厳格HABITATでは支持不足。SloveniaのBS/BTとは別比較 | 独立反証として残す。再現しない理由を「国差」だけに帰属しない。 [v12.7](CzechVeg_vegetation_class_same_study_wind_results_v127.md) |
| 実鉄道網と種共有 | SI同区間Jaccard差+.052644、CZ+.02394。地理条件で大きく減衰 | 道路を加えたモデルへの鉄道距離追加はSI予測改善−.537%、CZ−2.633%（悪化）。機能群・帰化年代別の微小肯定も種数調整で反転 | 「同一路線なら似る」と「鉄道距離が未知区間を予測する」を分ける。 [SI](Slovenia_2024_original_OSM_rail_road_control_v150_v151_ja.md)、[CZ](Czech_2023_true_rail_road_source_design_v155_v156_ja.md)、[4群](Two_country_original_four_guild_source_robustness_v163_ja.md) |
| 種名か粗形質かの未知路線予測 | SI15区間留保で4形質群のみより種固有履歴を加えると対数損失4.695%改善、14/15区間有利。上位10種除外でも1.591%、12/15 | 1原研究、既知種の別区間予測。全く未知の種・国への転移とは違う | 草丈だけに狭めない中心候補。 [v17.2](Slovenia_original_species_vs_guild_15route_microhabitat_predict_v172_ja.md)、[v17.3](Slovenia_original_common_species_leaveout_shrinkage_15route_v173_ja.md) |
| Czech鉄道文脈種・写真 | 36原研究の44共通ブロックで104種の鉄道−道路文脈頻度差が正。Echiumは+.06909、該当研究除外最小+.04558 | 現代・厳格HABITATに限定した支持は弱い。4,176写真の人手判読は未完。駅39/39と文脈差は別比較 | 種選定の根拠として保持。写真の表現型効果は未測定。 [v17.4–17.5](CzechVeg_original_39station_rail_context_to_all_Echium_photo_annotations_v174_v175_ja.md)、[台帳](Czech_Echium_public_iNaturalist_all4176_photo_label_REVIEW_queue_v176_ja.md) |
| 物理供給・発芽・野生化の既刊証拠 | カナダ19地点1,979試料で穀粒落下、スイス貨物混入300種子→75実生、オーストリア60地点2年の野生化植物、スイス鉄道79調査地中58でアブラナ | 同一コホートの輸送→発芽→定着連鎖ではない。年約110tはカナダ著者のモデル外挿。発芽は2製粉所試料 | 背景の機構証拠として引用。本リポの独立野外結果として合算しない。 [カナダ](Canada_direct_freight_train_grain_deposition_v130.md)、[スイス](Swiss_wheat_seed_viability_and_route_monitoring_v134_v135.md)、[野生化](Swiss_rail_freight_GM_feral_transition_2011_2020_v131_v133.md) |
| 鉄道進化・既存メタ解析 | Arabidopsis arenosaの共通栽培・交配・ゲノムから早期開花等の適応証拠は既存。Mellレビュー316研究・メタ47研究205事例 | 本リポの写真から進化を推定したわけでない。メタから抜粋した仏正・豪負のgを新メタに見せない | 新規性の上限を明示する必須先行文献。 [既刊の照合](Mell2024_true_railway_meta_and_Arenosa_evolution_v147_v149_ja.md) |
| フランス線路–法面・管理 | 394個2m²区画で線路より法面の種数が多い、線路は一年草・ruderal、法面は多年草・competitive | 170件は管理要請で実際の草刈り処理でない。2019政策の前後比較もない | 生育地による群集差の独立文献、管理因果の結果にしない。 [RATP](RATP2025_original_track_vs_embankment_staff_ALERT_evidence_v161_ja.md) |
| 数学・測定仕様 | CWMの種交代／種内変化／交互作用／非対応種の分解、画像・刈取りイベント解析の実装 | 合成テストは生物学的肯定ではない。周辺固定区画・同種実形質が未完備 | 方法付録か将来計画。 [分解](Rail_CWM_turnover_true_paired_ITV_interaction_unpaired_v169_ja.md)、[刈取り](mowing_roi_event_study.md) |

## 初期系列から引き継ぐ論文構成上の判断

最も広い一つの問題は「鉄道による供給機会が、どの空間スケール・どの生物学的単位で観測群集の再現性へ変換されるか」である。この問いなら貨物、ネットワーク、微生育地、種固有予測、形質、観測バイアスを捨てずに扱える。ただし異なるデータで測った段階を一本の検証済み因果経路にはできない。

本文候補は、実群集での局所的再現性・種固有情報・形質配置と、遠距離輸送代理指標の追加予測が弱いことの対照。GBIFの豊富な探索結果は全数を補足証拠台帳で保持し、提供元除外・応答定義による更新を1系列として示す。衛星・写真の未識別結果は「方法の失敗」と隠すのでなく、観測スケールが比較可能な単位を決める証拠として配置する。全部を同じ効果量に統合する必要はない。

初期データだけで別論文にできる候補は二つある。第一は種等重みと出現機会重みの生態学的問いの違い（Sloveniaの反転、Czechでの非再現、MIRENでは両方正）。第二は生物記録データで輸送網の効果を読む際の観測過程（希釈深度、提供元、観察者×地点継続）。どちらも、既存の方法論との重複と複数独立資料での追加検証を文献レビューで確認してから投稿単位を決める。

## 対象文書の全件索引

以下はファイルの全件索引であり、独立した研究や肯定結果の件数ではない。設計、取得監査、訂正、更新報告も含む。読み取りは文書本文と結果表の監査であり、原資料全行の再計算ではない。

対象：114文書。

- [50km_macro_rail_opportunity_v30.md](50km_macro_rail_opportunity_v30.md) — v3.0 — From ten small urban windows to independent 50km equal-area macro source regions
- [50km_macro_rail_road_plant_results_v36.md](50km_macro_rail_road_plant_results_v36.md) — tetsu v3.6 — 50km macro railway vs road *continuous* accessibility and observed introduced plants
- [6country_macro_rail_circuit_flora_results_v46.md](6country_macro_rail_circuit_flora_results_v46.md) — tetsu v4.6 — Six-country 50km equal-area macro study: physical rail, road, circuits, and country-listed introduced plants
- [Austria_feral_plants_vs_roadseed_airflow_mechanism_v128_v129.md](Austria_feral_plants_vs_roadseed_airflow_mechanism_v128_v129.md) — tetsu v12.8–v12.9: separate the mechanism of arriving seeds from sampled feral establishment and source railway functional filters
- [Canada_direct_freight_train_grain_deposition_v130.md](Canada_direct_freight_train_grain_deposition_v130.md) — tetsu v13.0 — DIRECT measured grain/seed-kernel deposition from moving freight TRAINS (not just observed feral plants)
- [circuit_macro_results_v38.md](circuit_macro_results_v38.md) — tetsu v3.8 — Circuit-theory railway network: source novelty versus biological test
- [circuit_rail_vs_road_source_validation_v41.md](circuit_rail_vs_road_source_validation_v41.md) — tetsu v4.0–v4.2 — independent circuit-network controls
- [circuitscape_network_design_v37.md](circuitscape_network_design_v37.md) — v3.7 — Circuitscape NETWORK mode for railway-driven plant redistribution
- [claim_evidence_crosswalk.md](claim_evidence_crosswalk.md) — 鉄道植物：中心主張と実証証拠の対応表
- [Czech39_real_railway_functional_vs_taxonomic_source_v143.md](Czech39_real_railway_functional_vs_taxonomic_source_v143.md) — v14.3 — Independent Czech 39-station field flora and Pladias four-guild functional–taxonomic identity diagnostic
- [Czech_2023_true_rail_road_source_design_v155_v156_ja.md](Czech_2023_true_rail_road_source_design_v155_v156_ja.md) — 鉄 v15.5–v15.6：チェコ実駅×真の歴史的鉄道・道路網の独立反証（設計凍結）
- [Czech_archaeophyte_vs_neophyte_paired_source_contrast_v158_ja.md](Czech_archaeophyte_vs_neophyte_paired_source_contrast_v158_ja.md) — 鉄 v15.8：考古帰化植物と新帰化植物の鉄道区間共有差を同じ比較ペアで検証
- [Czech_Echium_public_iNaturalist_all4176_photo_label_REVIEW_queue_v176_ja.md](Czech_Echium_public_iNaturalist_all4176_photo_label_REVIEW_queue_v176_ja.md) — 鉄 v17.6：画像を再配布せず4,176個の公開iNaturalist観察を実測ラベル化する台帳
- [Czech_original_39station_SR70_geography_falsification_v153_v154_ja.md](Czech_original_39station_SR70_geography_falsification_v153_v154_ja.md) — 鉄 v15.3–v15.4：チェコ鉄道39駅の真正位置対応と地理的種共有の独立原著反証（2026-11-15）
- [Czech_original_arch_neo_real_2023_railroad_5stretch_CV_v160_ja.md](Czech_original_arch_neo_real_2023_railroad_5stretch_CV_v160_ja.md) — 鉄 v16.0：チェコ真正鉄道・道路網で考古帰化植物と新帰化植物の種共有予測を直接比較（2026-11-15）
- [Czech_original_archaeophyte_neophyte_geographic_identity_v157_ja.md](Czech_original_archaeophyte_neophyte_geographic_identity_v157_ja.md) — 鉄 v15.7：考古帰化植物・新帰化植物別の鉄道区間内種同一性――元チェコ原著の分類を再使用（2026-11-15）
- [Czech_original_flower_photo_same_month_local_control_v167_v168_ja.md](Czech_original_flower_photo_same_month_local_control_v167_v168_ja.md) — 鉄 v16.7–v16.8：線路30m写真で「同種の花・草姿」まで比較できるか
- [Czech_real_railway_stations_source_audit_and_area_alien_sensitivity_v107_v108.md](Czech_real_railway_stations_source_audit_and_area_alien_sensitivity_v107_v108.md) — tetsu v10.7–v10.8: Original independent Czech railway stations – documented non-significant station-area proxy result
- [Czech_Slovenia_original_road_rail_CV_robustness_v159_ja.md](Czech_Slovenia_original_road_rail_CV_robustness_v159_ja.md) — 鉄 v15.9：チェコ・スロベニア原著2研究の鉄道距離追加予測、区間影響と重み付け感度（2026-11-15）
- [CzechVeg_original_39station_rail_context_to_all_Echium_photo_annotations_v174_v175_ja.md](CzechVeg_original_39station_rail_context_to_all_Echium_photo_annotations_v174_v175_ja.md) — 鉄 v17.4–v17.5：独立CzechVegの植物種別鉄道文脈指標から、全写真の開花可視性ラベル可用性へ
- [CzechVeg_Pladias_independent_rail_source_v117_v118.md](CzechVeg_Pladias_independent_rail_source_v117_v118.md) — tetsu v11.7–v11.8: Independent Czech vegetation-plot macro source and exactly homologous Pladias diaspore strategies (original-source verification)
- [CzechVeg_Pladias_real_railroad_diaspore_2026_v119_v120.md](CzechVeg_Pladias_real_railroad_diaspore_2026_v119_v120.md) — tetsu v11.9–v12.0: authentic CzechVeg + Pladias plant trait join, independent railway-text versus road-text habitat diagnostic
- [CzechVeg_Pladias_same_region_year_area_study_gates_v121_v122.md](CzechVeg_Pladias_same_region_year_area_study_gates_v121_v122.md) — tetsu v12.1–v12.2 — CzechVeg railway-vs-road source comparison after original region, date, area and original publication controls
- [CzechVeg_Pladias_same_study_railroad_traits_v123.md](CzechVeg_Pladias_same_study_railroad_traits_v123.md) — tetsu v12.3 — Same ORIGINAL source publication reveals that RAIL/ROAD wind-trait sign is study-selection sensitive
- [CzechVeg_syntaxon_exact_to_verified_class_v124_v126.md](CzechVeg_syntaxon_exact_to_verified_class_v124_v126.md) — tetsu v12.4–v12.6: original CzechVeg vegetation SYNTAXON matching reveals source selection bottleneck; class hierarchy permits wider diagnostic
- [CzechVeg_vegetation_class_same_study_wind_results_v127.md](CzechVeg_vegetation_class_same_study_wind_results_v127.md) — tetsu v12.7 — Exact original Czech vegetation CLASS-matched wind-diaspore plant comparison fails to reproduce a robust Slovenian-style frequency sign inversion
- [Echium_29pair_original_CC_photo_pixel_blind_review_v170_ja.md](Echium_29pair_original_CC_photo_pixel_blind_review_v170_ja.md) — 鉄 v17.0：チェコ原著Echium 29組の真正iNaturalist実写真をライセンス確認してブラインド化する
- [Echium_CC0BY_rights_first_rematch_real_photo_review_v171_ja.md](Echium_CC0BY_rights_first_rematch_real_photo_review_v171_ja.md) — 鉄 v17.1：CC0・CC BY写真を先に選んだEchium再対応と、実画像3組のブラインド暫定レビュー
- [empirical_feasibility_2026-11-07.md](empirical_feasibility_2026-11-07.md) — Empirical data feasibility gate — 2026-11-07
- [FAF6_2022_crop_rail_vs_truck_macro_source_v79.md](FAF6_2022_crop_rail_vs_truck_macro_source_v79.md) — tetsu v7.9 — Genuine 2022 FAF6 agricultural commodity RAIL vs TRUCK: macro freight opportunity, NOT observed seed dispersal
- [FAF6_agricultural_OD_rail_vs_truck_network_v80.md](FAF6_agricultural_OD_rail_vs_truck_network_v80.md) — tetsu v8.0 — Actual FAF6 state OD freight network: agricultural RAIL versus TRUCK, no seed-transport causality
- [FR_original_GBIF_events_growthform_rail_macro_v73.md](FR_original_GBIF_events_growthform_rail_macro_v73.md) — tetsu v7.3 — real original GBIF events versus legacy facets: French GRIIS introduced flora × GIFT growth form × physical railway/road paths
- [France_Penone2012_original_thesis_occupancy_v145.md](France_Penone2012_original_thesis_occupancy_v145.md) — v14.5 — フランス71鉄道調査地点、原著博士論文Appendix Bの出現頻度を発見
- [France_Penone2012_source_marginal_nonidentification_v146.md](France_Penone2012_source_marginal_nonidentification_v146.md) — v14.6 — 同じ原著種頻度から正負の「路線内−路線間種共有」をどちらも構成できる
- [gbif_10country_source_v18.md](gbif_10country_source_v18.md) — v1.8 — Ten-country GBIF Plantae record source census (real API)
- [GBIF_dataset_recording_bias_and_K80_v52.md](GBIF_dataset_recording_bias_and_K80_v52.md) — tetsu v5.2 — Recording-source concentration, dataset overlap and fixed-cohort rarefaction
- [GBIF_GRIIS_dataset_leaveout_macro_results_v57.md](GBIF_GRIIS_dataset_leaveout_macro_results_v57.md) — tetsu v5.7 — Full introduced-plant source vs leaving out the single dominant GBIF dataset
- [GBIF_original_DWCA_event_ingestion_v58.md](GBIF_original_DWCA_event_ingestion_v58.md) — tetsu v5.8 — Original GBIF DWCA event ingestion is implemented, but authentic occurrence downloads have NOT been obtained
- [GBIF_public_FR_Burgundy_original_events_v66_protocol.md](GBIF_public_FR_Burgundy_original_events_v66_protocol.md) — v6.6 frozen public original-GBIF France Burgundy 50km source-event quality protocol
- [GBIF_rarefaction_and_status_contract_v49.md](GBIF_rarefaction_and_status_contract_v49.md) — tetsu v4.7–v4.9: rarefaction-depth sensitivity and GBIF true detection/status source audit
- [GIFT_growthform_provider_leaveout_macro_v72.md](GIFT_growthform_provider_leaveout_macro_v72.md) — tetsu v7.2 — The woody-flora railway signal does not survive GBIF source-provider exclusion
- [holdout_us_southeast_v53_protocol.md](holdout_us_southeast_v53_protocol.md) — tetsu v5.3 — Locked out-of-analysis-region US Southeast replication check
- [horsetail_spore_roi_inputs_v10.md](horsetail_spore_roi_inputs_v10.md) — v1.0 horsetail spore-cone ROI inputs (not flower ROI)
- [inat_rail_phenophase_pilot.md](inat_rail_phenophase_pilot.md) — iNaturalist railway phenophase feasibility pilot (v0.4)
- [japan_ordinary_rail_road_separation_v15.md](japan_ordinary_rail_road_separation_v15.md) — v1.5 — Japan ordinary railway sites with independent road-only GIS controls
- [japan_ordinary_rail_s2_current_v15.md](japan_ordinary_rail_s2_current_v15.md) — v1.5 — real 2025 Sentinel-2 vegetation verification of Japan rail/road separated candidates
- [japan_ordinary_rail_worldcover_2021_v16.md](japan_ordinary_rail_worldcover_2021_v16.md) — v1.6 — ESA WorldCover 2021 reveals tree-covered 2025-SCL4 railway candidates
- [japan_rail_microhabitat_reanalysis_v10.md](japan_rail_microhabitat_reanalysis_v10.md) — v1.0 — Japanese railway vegetation: original reanalysis of published microhabitat counts
- [japan_rail_mowing_historical_baseline.md](japan_rail_mowing_historical_baseline.md) — Railway maintenance is a separate selection gate — documented Japanese baseline
- [japan_rail_weed_photo_feasibility_v10.md](japan_rail_weed_photo_feasibility_v10.md) — v1.0 modern Japanese railway-weed photo source audit (2018–2025)
- [japan_two_grass_sites_S2_GSI_v16.md](japan_two_grass_sites_S2_GSI_v16.md) — v1.6 — Two mapped grassland rail candidates: actual S2 2025 and GSI aerial-source cross-check
- [japan_worldcover_370_site_posthoc_discovery_v16.md](japan_worldcover_370_site_posthoc_discovery_v16.md) — v1.6 — a real sampling confound: road-isolated railway pixels become overwhelmingly forest
- [lrt_n13_road_prelandcover_gates_v13.md](lrt_n13_road_prelandcover_gates_v13.md) — v1.3 — Independent road and baseline land-class validity gates, one Utsunomiya LRT
- [lrt_true_rail_verge_20m_v14.md](lrt_true_rail_verge_20m_v14.md) — v1.4 — Are any 10m Sentinel pixels actually suitable for railway-verge vegetation?
- [macroecology_rail_network_v17.md](macroecology_rail_network_v17.md) — tetsu v1.7: macroecology of railway networks and plant filtering
- [Mell2024_true_railway_meta_and_Arenosa_evolution_v147_v149_ja.md](Mell2024_true_railway_meta_and_Arenosa_evolution_v147_v149_ja.md) — 鉄 v14.7–v14.9：世界の鉄道植物メタ解析可能性と遺伝的適応の原著証拠（2026-11-15）
- [MIREN_real_roadside_alien_flora_macro_v98_v99.md](MIREN_real_roadside_alien_flora_macro_v98_v99.md) — tetsu v9.8–v9.9 — Independent standardized road-verge vegetation survey reveals near-universal roadside alien-species enrichment (not yet a RAIL effect)
- [MIREN_standardized_roadside_alien_species_replacement_v100_v101.md](MIREN_standardized_roadside_alien_species_replacement_v100_v101.md) — tetsu v10.0–v10.1: Roadside alien plants show higher gains and losses, plus species replacement rather than simple nested additions
- [mowing_roi_event_study.md](mowing_roi_event_study.md) — Fixed-plot grass cutting × photograph ROI pilot
- [npn_active_rail_2024_results_v07.md](npn_active_rail_2024_results_v07.md) — v0.7 verified 2024 USA-NPN rail-near site census — results and identification gate
- [npn_active_rail_eligibility_v07.md](npn_active_rail_eligibility_v07.md) — v0.7 — Actual 2024 NPN plant phenology activity at registered rail-near sites
- [npn_candidate_status_verification_v08.md](npn_candidate_status_verification_v08.md) — v0.8 — Real individual 0/1 flowering eligibility after species/data and geography screening
- [npn_daymet_thermal_control_v09.md](npn_daymet_thermal_control_v09.md) — v0.9: Real 1-km Daymet weather at matched NPN rail-near vs rail-far flowering sites
- [npn_observed_onset_midpoint_v08.md](npn_observed_onset_midpoint_v08.md) — Observed Open flowers interval midpoint by site — real-data v0.8 diagnostic
- [npn_site_pheno_2023_25_v08.md](npn_site_pheno_2023_25_v08.md) — v0.8 — NPN national Site Phenometrics joined to rail geography (2023–2025)
- [overture_gbif_10country_macro_v19.md](overture_gbif_10country_macro_v19.md) — v1.9 — International rail graph × plant occurrence spatial-source audit
- [pairwise_inputs.md](pairwise_inputs.md) — Descriptive pilot CSV contract
- [photo_and_phenology_inputs.md](photo_and_phenology_inputs.md) — Photo ROI and flowering-phenology input contract
- [protocol.md](protocol.md) — v0.1 Study design — trains as a possible seed dispersal network
- [Rail_CWM_turnover_true_paired_ITV_interaction_unpaired_v169_ja.md](Rail_CWM_turnover_true_paired_ITV_interaction_unpaired_v169_ja.md) — 鉄 v16.9：種間の入れ替わりと真の同種形質シフトを分けるCWMの厳密分解
- [rail_delivery_vs_functional_filter_macro_v71.md](rail_delivery_vs_functional_filter_macro_v71.md) — v7.0–v7.1 — Railway delivery × functional filtering macroecology: measured first-stage
- [rail_environment_traits_protocol.md](rail_environment_traits_protocol.md) — Railway life-history and phenotypic response workstream (v0.3)
- [rail_griis_macro_results_v25.md](rail_griis_macro_results_v25.md) — tetsu v2.5 — GRIIS-introduced plant composition vs mapped railway connectivity
- [rail_macro_source_readiness_v55.md](rail_macro_source_readiness_v55.md) — tetsu v5.5 — Macro railway–plant analysis: source reliability and next data acquisition
- [rail_road_control_pilot.md](rail_road_control_pilot.md) — v0.5 fixed rail-vs-road observation contrast (2026-11-07)
- [rail_vs_global_road_weighting_frequency_v113_v114.md](rail_vs_global_road_weighting_frequency_v113_v114.md) — tetsu v11.3–v11.4: frequency-sensitive RAIL dispersal guild vs frequency-robust global ROAD alien enrichment
- [railway_history_n05_v11.md](railway_history_n05_v11.md) — v1.1 — Japan official N05 rail history (1950–2025) and independently validated event design
- [RATP2025_original_track_vs_embankment_staff_ALERT_evidence_v161_ja.md](RATP2025_original_track_vs_embankment_staff_ALERT_evidence_v161_ja.md) — 鉄 v16.1：フランスRATP原著から「草刈り要請」と「鉄道植物の環境フィルター」を分離（2026-11-15）
- [real_photo_feasibility_2026-11-07.md](real_photo_feasibility_2026-11-07.md) — First real iNaturalist phenophase × rail proximity pilot — 2026-11-07
- [real_photo_roi_review_2026-11-07.md](real_photo_roi_review_2026-11-07.md) — First real-photo ROI feasibility review — 2026-11-07
- [real_rail_vs_road_tokyo_2026-11-07.md](real_rail_vs_road_tokyo_2026-11-07.md) — Tokyo real-data railway vs road photographic-occurrence audit — v0.5
- [road_adjusted_rail_alien_macro_results_v27.md](road_adjusted_rail_alien_macro_results_v27.md) — v2.7 — Ten-country macro pilot: introduced plant composition, mapped railway connectivity, and actual road-density matching
- [sixcountry_beta_turnover_nestedness_rail_v74.md](sixcountry_beta_turnover_nestedness_rail_v74.md) — v7.4 — REAL six-country railway macro: species turnover versus nestedness
- [Slovenia_2024_original_OSM_rail_road_control_v150_v151_ja.md](Slovenia_2024_original_OSM_rail_road_control_v150_v151_ja.md) — 鉄 v15.0–v15.1：元のスロベニア89地点に現実の鉄道・道路ネットワークを接続して反証する
- [Slovenia_functional_vs_taxonomic_spatial_source_v141_v142.md](Slovenia_functional_vs_taxonomic_spatial_source_v141_v142.md) — tetsu v14.1–v14.2: REAL GPS-matched taxonomic identity versus functional profiles in Slovenian railway communities — geographic confounding changes the interpretation
- [Slovenia_original_common_species_leaveout_shrinkage_15route_v173_ja.md](Slovenia_original_common_species_leaveout_shrinkage_15route_v173_ja.md) — 鉄 v17.3：頻出植物に依存しないかを再学習して確かめた鉄道内フィルターの種固有性
- [Slovenia_original_species_vs_guild_15route_microhabitat_predict_v172_ja.md](Slovenia_original_species_vs_guild_15route_microhabitat_predict_v172_ja.md) — 鉄 v17.2：鉄道微生育地のフィルターは形質群だけか、植物種自身の個性も別路線で再現するか
- [Slovenia_rail_annual_mix_vs_within_lifespan_wind_v106.md](Slovenia_rail_annual_mix_vs_within_lifespan_wind_v106.md) — tetsu v10.6: original RAIL wind-diaspore deficit persists WITHIN source lifespan strata
- [Slovenia_rail_section_cluster_and_missing_trait_sensitivity_v105.md](Slovenia_rail_section_cluster_and_missing_trait_sensitivity_v105.md) — tetsu v10.5 — original Slovenia railway wind-vs-annual trait result survives section clustering and unknown-trait worst-case bounds
- [Slovenia_RAIL_section_flip_and_strategy_placebos_v115_v116.md](Slovenia_RAIL_section_flip_and_strategy_placebos_v115_v116.md) — tetsu v11.5–v11.6: source-verified side aggregation, exact section-flip and dispersal-class placebo sensitivity
- [Slovenia_real_RAIL_diaspore_wind_vs_annual_v104.md](Slovenia_real_RAIL_diaspore_wind_vs_annual_v104.md) — tetsu v10.4: Real railway microhabitat annual-life-history enrichment does NOT imply wind-dispersed propagule enrichment
- [Slovenia_real_RAIL_microhabitat_plant_functional_v102_v103.md](Slovenia_real_RAIL_microhabitat_plant_functional_v102_v103.md) — tetsu v10.2–v10.3: Actual Slovenian railway plant microhabitats and life-history filtering proxies
- [Slovenia_species_frequency_weighting_reversal_v109_v110.md](Slovenia_species_frequency_weighting_reversal_v109_v110.md) — tetsu v10.9–v11.0 REAL Slovenian railway source: same species, opposing wind trait contrasts when abundance opportunity changes weights
- [sources.md](sources.md) — Sources verified 2026-11-07
- [Swiss_rail_freight_GM_feral_transition_2011_2020_v131_v133.md](Swiss_rail_freight_GM_feral_transition_2011_2020_v131_v133.md) — tetsu v13.1–v13.3 — genetically traced feral oilseed rape recruitment along Swiss freight railway and port interfaces, with full source-unit and 8-year monitoring audit
- [Swiss_wheat_seed_viability_and_route_monitoring_v134_v135.md](Swiss_wheat_seed_viability_and_route_monitoring_v134_v135.md) — tetsu v13.4–v13.5：貨物中の発芽可能な種子と、沿線での非検出を分ける
- [tetsu_two_track_real_field_evidence_v141_v144_ja.md](tetsu_two_track_real_field_evidence_v141_v144_ja.md) — 鉄 v14.1–v14.6：二本立てマクロ解析の原表再解析と段階別証拠マップ
- [Top20_original_rail_flora_iNaturalist_real_30m_photo_feasibility_v164_v165_ja.md](Top20_original_rail_flora_iNaturalist_real_30m_photo_feasibility_v164_v165_ja.md) — 鉄 v16.4–v16.5：実原著から候補20種を固定し、iNaturalist写真×線路30mの可用性を監査
- [Two_country_original_four_guild_source_robustness_v163_ja.md](Two_country_original_four_guild_source_robustness_v163_ja.md) — 鉄 v16.3：真のスロベニア・チェコ原著2研究、風散布×生活史4群の追加鉄道経路効果の頑健性監査
- [Two_country_original_railroad_lifeform_diaspore_guilds_v162_ja.md](Two_country_original_railroad_lifeform_diaspore_guilds_v162_ja.md) — 鉄 v16.2：風散布形質×生活史の4群で「鉄道経路種共有」と「環境フィルター」を区別する原著2国再解析
- [US48_GBIF_iNaturalist_observer_site_recording_process_v95_v97.md](US48_GBIF_iNaturalist_observer_site_recording_process_v95_v97.md) — tetsu v9.5–v9.7 — Original iNaturalist observer/site continuity and source "first-recorded" botanical categories
- [US_FAF6_crop_resource_rail_floristic_specificity_v84.md](US_FAF6_crop_resource_rail_floristic_specificity_v84.md) — tetsu v8.4 — Does the US cargo–floristic relationship distinguish crop freight from resource freight?
- [US_FAF6_directional_state_FE_node_label_null_v89.md](US_FAF6_directional_state_FE_node_label_null_v89.md) — tetsu v8.9 — Directional freight–plant-record associations with origin/destination state fixed effects and constrained node-label nulls
- [US_FAF6_directional_temporal_GRIIS_results_v85_v88.md](US_FAF6_directional_temporal_GRIIS_results_v85_v88.md) — tetsu v8.5–v8.8 — Genuine 2022 U.S. crop-freight directions versus earlier-origin and later-destination GRIIS-listed plant records
- [US_FRA_train_traffic_vs_GRIIS_GIFT_macro_v75_v77.md](US_FRA_train_traffic_vs_GRIIS_GIFT_macro_v75_v77.md) — tetsu v7.5–v7.7 — railway TRAFFIC observation versus functional plant source patterns
- [US_FRA_woody_observation_negative_controls_v78.md](US_FRA_woody_observation_negative_controls_v78.md) — tetsu v7.8 — FRA train reports × plants: all-Plantae response comparator and closed-crossing infrastructure comparator
- [US_GBIF_2018_21_2023_25_reporting_effort_v91.md](US_GBIF_2018_21_2023_25_reporting_effort_v91.md) — tetsu v9.1 — New recorded plant taxa coincide with a nationwide GBIF reporting shift
- [US_GBIF_iNaturalist_dominance_and_FAF_leaveout_v92_v94.md](US_GBIF_iNaturalist_dominance_and_FAF_leaveout_v92_v94.md) — tetsu v9.2–v9.4: A single citizen-science GBIF dataset dominates US recorded plant source assemblages, but removing it does not recover a general rail signal
- [US_GBIF_original_coordinate_2025_Census_border_QA_v90.md](US_GBIF_original_coordinate_2025_Census_border_QA_v90.md) — tetsu v9.0 — Original GBIF stateProvince source geography QA against independent Census 2025 state boundaries
- [us_southeast_heldout_K80_source_feasibility_v53.md](us_southeast_heldout_K80_source_feasibility_v53.md) — tetsu v5.3 — Frozen US Southeast out-of-analysis-region validation: NOT ESTIMABLE
- [US_state_FAF6_crop_OD_vs_GRIIS_flora_v80_v83.md](US_state_FAF6_crop_OD_vs_GRIIS_flora_v80_v83.md) — tetsu v8.0–v8.3 — Freight network opportunity and recorded introduced-flora geography are decoupled in state-scale U.S. source data
- [usa_npn_real_phenology_rail_v06.md](usa_npn_real_phenology_rail_v06.md) — v0.6 — Real USA-NPN phenophase × federal rail GIS feasibility (2026-11-07)
- [utsunomiya_lrt_photo_sentinel_v11.md](utsunomiya_lrt_photo_sentinel_v11.md) — v1.1 — New 2023 Utsunomiya LRT passenger service: botanical photo and satellite source checks
- [utsunomiya_s2_measured_ndvi_v12.md](utsunomiya_s2_measured_ndvi_v12.md) — v1.2 — Utsunomiya LRT actual Sentinel-2 pixels: cloud-masked NDVI, not railway causal effect
