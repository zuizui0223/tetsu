# tetsu v28.9 — Kashmir authentic 27-site strict annual sign identification after Kew-curated perennial sources

検証 2026-11-18 JST。研究 branch のみ（main は変更していない）。一次文献により厳密に採録した多年草の原学名2件と、2023全球植物life-cycleデータを同一原本Dryad鉄道サイトで結合。**統計的な全球普遍性や因果的な鉄道適応は主張しない**。

## 原データと出典
- Rashid et al. (2021) Journal of Applied Ecology, DOI https://doi.org/10.1111/1365-2664.13961、Dryad DOI https://doi.org/10.5061/dryad.0p2ngf20h original version 99735、`Matrix_plot.xlsx` original SHA256 `529f3b97521933018e12623b61f531f012c67241cae2602a6675dbd534192474`、31独立空間地点（駅16、駅間15）、2014/2017の四方形区100m²、123原学名taxa列、元種記録欠損4cell。4区×両年が真正揃う元27地点（駅15、駅間12）を固定。
- Poppenwimer, Mayrose & DeMalach (2023) Nature DOI https://doi.org/10.1038/s41586-023-06644-x / source CSV Figshare ID 20764777, collection https://doi.org/10.6084/m9.figshare.c.6176239、原 `LifeForm_Trait_Database.csv` 235,979 rows、original SHA256 `83b23e7874d73387a0c24aeae3e2931c50ec7c78cad52ff26cdfc99d1b6623b1`。元accepted namesはlowercaseだが原二名法を厳格に復元して照合。104/123 original taxaと一致。Global trait database is **not** another railway field study.
- Kew POWO original-name evidence-based manual curations (2026-11-18): `Cynodon dactylon (L.) Pers.` perennial, https://powo.science.kew.org/taxon/urn%3Alsid%3Aipni.org%3Anames%3A397065-1/general-information ; `Breea arvensis (L.) Less.` Kew homotypic synonym of accepted `Cirsium arvense (L.) Scop.` https://powo.science.kew.org/taxon/urn%3Alsid%3Aipni.org%3Anames%3A186181-1 ; accepted taxon perennial https://powo.science.kew.org/taxon/urn%3Alsid%3Aipni.org%3Anames%3A195034-1/general-information 。元original nameを曲げず文献を明記したCSV `data/literature/tetsu_v289_original_name_Kew_provenance_lifespan_overrides.csv` に2種のみ追加。分類済み106/123、未分類17 taxa。

## 主結果：original all 27 site×2年×元Zoneを固定した部分識別
主対比 D = original spatial-site equal weighted mean over 2014 and 2017 [strict annual fraction among ALL originally recorded taxa in Track − strict annual fraction in Interior Plot 1]。一年草と確認されたtaxaは1、多年草・二年草は0。残る未知taxaについては**同一speciesは全source site/zone/yearで一貫してannual=0または1**の二値と仮定し、すべての可能な配置にわたる厳密な上下限を、original incidence matrix上の線形係数で解析的に計算。未知を各plotごとに独立に選択する過度に広い欠損範囲とは異なる。未知taxaにはLactuca serriola annual/biennial、Lotus corniculatus annual/perennial、Senecio chrysanthemoides official herbarium perennial vs Kashmir field flora annualという資料不一致を保持している。

| 同一27原地点 | 未分類species=17の全ラベル組合せにおける厳密なD範囲 | Sign identified? |
|---|---:|---|
| 2014 | **−0.218468 ～ −0.110537** | sample内負方向に限定 |
| 2017 | −0.193386 ～ +0.010571 | 方向未確定 |
| 2年平均 | **−0.202726 ～ −0.053184** | sample内負方向に限定 |
| 2017−2014 2年の差 | −0.059772 ～ +0.205962 | 時間変化の方向未確定 |

重要：この範囲は **trait codingの同定不確実性** の範囲であり、フィールドからの母集団抽出不確実性・信頼区間・メタ解析prediction intervalではない。元27siteは無作為世界鉄道サンプルではなく、ほぼ同一地域の同一路線調査。station/Betweenは列車本数実測ではない。

## 源元空間地点の技術的再標本化（1999 draw）
- Source original27地点をwhole-site clusterとして復元抽選し、各回未知タクソンすべての同一speciesラベル制約の厳密な下限/上限を再計算。下限の技術 q2.5–q97.5 は **[−0.271949,−0.127971]**、上限の技術 q2.5–q97.5 は **[−0.124482,+0.022116]**。各反復の最悪上限も負だった割合は **0.911956**（=約91.2%）。したがって標本内の観測係数は未分類17種の全選択で負だが、独立field sitesをリサンプルすると95%上限は0を跨ぐ。全インドへの頻度論的CIとして扱わない。
- 2023全球traitだけの元104種段階: source 27site平均 Dのsharp bounds **[−0.202726,+0.113615]**、sign unidentified。Kew2種真正性確認の後には上限が **+0.113615→−0.053184**へ減少した。これは研究者が目標変数に合わせて生活史を恣意的に割り当てたものではなく、元ranked unknownの外部Kew corroborationに基づく追加分類である。
- LIFE CYCLEとして `Cynodon dactylon` の寄与の幅は13.686 pp、`Breea arvensis` は2.994 pp。2種を分離する以前には未分類19種全体の不確実性幅31.634 pp。この2種の原証拠により残り17種の幅は14.954 ppとなった。

## 同じ母集団での50%既知率感度と整合するか
- v28.7/v28.8の元分類済み104種だけを分母に取る対比（strict known species fraction）は、Track+Interior1の元 site-yearに trait coverage≥70%を要求すると27中3サイトに制約され、平均 −0.253582。50%閾値では24サイト、2年平均 −0.058059、サイト再標本化技術区間 [−0.159255,+0.057953]。この点推定とv28.9の全種分母によるboundsは**estimandと採用source sitesが異なる**。−0.058 vs [−0.2027,−0.0532]が一致しなくても同じ数値を比較しているわけではない。
- 欠測の原因はTrack区の形質被覆：元27site×2年の54 plot-year中、life-span known≥70%のTrackは22のみ、≥50%は51。Interior1では≥70%が48/54。123 species全体の84.6% coverageから各Track plotで70%以上が保証されるわけではない。これがstrict 3/27の主要原因。

## 他地域との対比・投稿射程
- Slovenia Šabić et al. (2026) DOI https://doi.org/10.1007/s44473-026-00155-z original 89 field sites・173 same-side BS/BT pairs・15 real line segments、strict annual fraction BS−BT **+0.054703**、source section-bootstrap technical interval [0.017619,0.086547]。対照はballast shoulder vs ballast toeであり India Track vs Interior1と深さ/幅/境界が違う。符号が違うことを同一rail-stressの反証やtrain-mediated transport効果へ直結しない。
- France Rakotoson et al. (2025) DOI https://doi.org/10.1111/wre.70025 original 394 2m² field quadrats (track 128, embankment266)。原著はtrackのannual/ruderal増加を報告するが、元完全plot×taxonとSEは著者依頼で未取得。no pooled meta-effect yet.
- Czech/Italy/Croatia-Bosnia/Slovenia four original-source functional beta differentiation results available in separate branch research/tetsu-fourstudy-trait-label-null-v285-20261118, authentic Actions #44296227526; their trait label shuffle `site CWM2D distance` is DIFFERENT estimator from annual Track-vs-Interior fraction. HR+BA are one field paper, not two replicate studies.
- 普遍性の現在地は『鉄道でどの地域も一年草化する』ではなく**環境・種プールとsource boundaryに応じて生活史選別の方向が変わり得る**という探索仮説。3大陸8独立原著で対応するeffect+SEを取得し、random-effects meta PIとholdoutを検証するという事前ゲートはまだ未達。

## 次の実装 / モデル分離
1. Strict annual gradientに同じ物理zoneのsource Francia 394 raw plotとAsia China–Laos 84を確保して、本物のwithin-study paired reference and sampling varianceを計算。Zone差を曝露変数としても扱いannual signを無理にプールしない。
2. Kashmir元source 2014と2017のsource site cluster差・未知life-span影響は既に記録済み。未分類17種のうちSenecio chrysanthemoides資料矛盾を照会し、Lactuca serriola annual/biennialを無理にannualへ固定しない。site-specific field disturbance/grass coverageを付加するまでは因果特定不可。
3. Plant-stature 4 original European field datasetsのwithin-rail taxon-fixed trait-label nullは独立の肯定的macro conditional structureであり、cross-continent taxonomic species baselineとの比較を優先。可能ならroad / non-transport matched reference同法sourceを取得する。

## 再現
- **真正ソースSOURCE SUCCESS** [GitHub Actions #44298354325](https://github.com/zuizui0223/tetsu/actions/runs/44298354325)：2021 Dryad raw SHA、2023 global Figshare source SHA、source names、Kew metadata + synthetic tests 全成功。
- Main `scripts/analyze_Kashmir_Kew_curated_strict_annual_v289.py`; 9 original synthetic math QA in `tests/test_Kashmir_global_unknown_trait_partial_id_v288.py`, plus new v28.9 provenance tests `tests/test_Kashmir_Kew_curated_lifespan_v289.py`.
- `.github/workflows/tetsu-v289-Kew-curated-Kashmir-lifespan.yml`はorig source matrixと元GPSをGitHubへ転載せず、trait evidence crosswalkと集約効果値/技術範囲のみActions artifactとして公開する。
