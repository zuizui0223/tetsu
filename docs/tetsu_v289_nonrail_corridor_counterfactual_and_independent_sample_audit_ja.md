# tetsu v28.9 — 線状人工生息地の一般性: スウェーデン道路・送電線下の公開原行列を厳格監査

更新 2026-11-19 JST。**新しい実測生態係数は、原データの真正取得・構造検証後にのみ報告**。この研究ブランチは v28.8 source-only rail vs nonrail support no-estimand の上に積む Draft で、\`main\`は変更しない。

## 研究の核心

v25.6 / v26.9 / v27.3 で、少なくとも2つの真正鉄道著者データでは地点・環境間の種の違いに「入れ子だけでない置換」が重要となり、スロベニア同一地点のBS–BT共有は、taxon-frequency/site-richness帰無を超える。MIREN mountain ROADSの16Regionでも本物のroadside–inland同一地点対応が帰無を超える。v28.0–28.4の **別地点をつなぐ同一路線効果は弱く頑健ではない**。この観察から鉄道専用効果だと判断するのでなく、人工的線状生息地全般への一般性を、**独立非鉄道道路・送電線原データ**で条件付き検証する。

### 原著に記載された外部証拠
- Dániel-Ferreira, Fourcade, Bommarco, Wissman & Öckinger (2023), **JAE**, DOI https://doi.org/10.1111/1365-2664.14378, 32の同じ4 km²スウェーデン景観、5 habitats（natural pastures, small-road verges, big-road verges, powerline corridors, field borders）、各原transect 200mに四つの1 m² quadrats、元plants recorded in 2016、元133 plant species。原著本文は植物および蝶の生息地間βで種置換が主な成分、マルハナバチでは入れ子性が相対的に重要と報告。しかし本文は**Sørensen形式β**。tetsuの既存Jaccard置換βと数値を併合しない。source original Figshare https://doi.org/10.6084/m9.figshare.22178651.v1 は公表。
- Dániel-Ferreira, Bommarco, Wissman & Öckinger (2020), **Scientific Reports**, DOI https://doi.org/10.1038/s41598-020-78090-y, **同じ32景観**、電力線有無×道路密度、pasture/small-road/big-roadの3生息地で128植物種、2016現地植物調査。元本文は送電線景観で植物種数が平均約6種多いと記載。source original Figshare https://doi.org/10.6084/m9.figshare.13286141 公開。

**重大な非独立性予測**：上の2原著は著者が重なり、調査年・32 landscape design・pasture＋2 road categoriesが一致するため、別々の「2独立研究」として効果量メタ解析に加算してはいけない。元sourceのlandscape/site identifiersを照合し、少なくとも**ひとつの原フィールドキャンペーンの重複候補**として一段上のstudy_idを共有させる。2023年は5 habitatsへの拡張であり、2020年原著で報告した効果を完全な独立検証として再利用しない。

### tetsuとは異なる元野外測定単位の注意

- Sweden 2023: 32 landscape IDs, up to 5 habitats/landscape, 4 quadrats/transect。元1m²でpresenceであり、各種群の統計解析にはtransect集計もあり得る。原1m²内部は少種なので、個plotのtrait CWMβが不識別となる可能性あり。元source species×site IDsの可否を**効果を閲覧する前に**確認。
- India Kashmir 2014/17: 31 original sites×2 years×4 zones, 100m²、比較可能27サイトのTrack taxa median4/100m²。source true structural trait CWMbeta support 0 panels、完全trait coverageでもstation0/15、between最大1/12。
- MIREN: roads only、16 Region×48 Roads×669 original site paired first-visit cohort、各原Plot100m²。こちらはKashmirと違うfield collection/country/yearsだが元原著はひとつのMIREN multinational programme。
- Slovenia: true rail BS+BT both railway microhabitats、89独立元地点/15 line sections。**非鉄道対照を持たない**のでmatched RAILvsNONRAILの効果量にはできない。

## 今回のデータ取得・QAの厳格ルール

本研究のGitHub Actionsはfigshare APIから**元source article IDs 13286141 / 22178651**、掲載DOI、元file IDs/bytes/MD5を取得。原ファイルを一時的にダウンロードして元publisher checksumがあれば照合。CSV、XLSX、ZIPの表列名、元行列サイズ、元taxon/survey/site ID列の有無を査読可能なschema metadataに縮約。元species行・元site GPS・元private property locationはCI artifactに保存しない。

まずデータ構造から次を判定する:
1. raw plant×site×habitat の真正行列があるか。original landscape ID、habitat labels、site IDがちゃんと同じ行列にあるか。**richness/summaryだけではspecies turnoverできない**。
2. 2020/2023の同一サイトが重なるか（publication数とfield study数は別）。
3. original trait harmonization可能な学名がどれだけあるか。先に原source binomial species countと元1m²/4m² field grainを固定し、可能なtrait minimum supportを判定する。
4. 実sourceのJaccard \`beta_jac=beta_jtu+beta_jne\` を再計算できるか。2023元論文のSørensen \`beta_sor=beta_sim+beta_sne\`と同じ語でも数値の直接比較をしない。min/median plot taxon richnessを先に判定し、源sourceで成立しないtrait CWM平均を作らない。
5. effect studies with genuine independent rail-nonrail AND road-pasture remain **別層として**報告。二独立鉄道研究だけでは世界的一般性に関する強いrandom-effects prediction intervalを出さない。

## 次の論文になる一般仮説

「線状人工生息地は単なる地域植生の貧弱な部分集合ではなく、地域種プールとの強い局所的な結合を残しながら、種置換を通じて独自の群集を作る。その置換が鉄道・道路・送電線とそれぞれの攪乱・基質条件でどこまで一致するか」が事前仮説。

科学的検定単位は \`publication\` ではなく \`independently surveyed field design\` を優先。source original regional structureを保つ。joint betaは各研究の現地plot内から同じ式で導く。既存のCzech species trait-label shuffled height CWM \`obs−null\` やJaccard距離 \`D_track−D_inland\` と、外部Sørensen publication coefficientsは一つのメタ平均にしない。

## 現行実装

\`scripts/audit_Figshare_Sweden_nonrail_infrastructure_source_opportunity_v289.py\`
\`tests/test_Figshare_Sweden_nonrail_source_schema_v289.py\`
\`.github/workflows/tetsu-v289-Sweden-original-nonrail-infrastructure-Figshare-source-gate.yml\`

原元著者個票は再配布せず、source checksums・source fieldsのschemaと行列適合ゲートだけを公開する。生態量は真正原行列・非独立検査後のみ報告。


## v28.9 SOURCE VERIFIEDの取得判定（2026-11-19）

**[Figshare真正原著source Actions SUCCESS #44457472437](https://github.com/zuizui0223/tetsu/actions/runs/44457472437)**。元Figshare DOI article IDs \`13286141\` と \`22178651\` と本文のpublisher file metadataを確認し、全5元public source fileの元byte取得、source file size / MD5 metadata checks と7/7 synthetic schema/private-plot exclusion QAが成功した。source元行列そのものはGitHubへ配布していない。

| 原出典 | 本当に取得したデータ構造 | 結論 |
|---|---|---|
| Dániel-Ferreira et al. 2020, figshare 13286141 | 32原景観×18列の要約CSV（元siteごと植物種名×出現行列なし） | 元field-species βを再計算できない |
| Dániel-Ferreira et al. 2023, figshare 22178651 | \`Mean_abundance.csv\` 205 rows×18 cols（Speciesラベルはあるが元site-IDなし）、\`Transect_coordinates.csv\` 133×8、\`Database_alpha_diversity.csv\` 133×18、\`Metadata.txt\` | 元真のsite×species floristic incidenceを復元できない。地点のα多様性・位置情報と種ごとの平均量を結合して真の植物種同所的共存を逆推定しない |

**新しい生態係数は1つも生成しなかった。** 元文献で報告されたSørensen replacementは引用できても、同じデータの真正site×taxon Jaccardとして再解析したと書けない。したがって本v28.9の原行列基準は **NOT_IDENTIFIABLE_ORIGINAL_SITE_TAXON_MATRIX**。2020/2023はauthor・2016景観調査が重なるので2つの独立field研究とは数えない。

**新しい次段のsource-backed外部検定**：既存tetsu v24.9/v25.0ですでに真正原行列が認証された世界MIREN山岳**道路**（8地域22道路272地点、厳格5地域13道路192地点）において、v28.5鉄道4研究と同じ真site出現固定の\`paired height/seed trait-label shuffle\`帰無を「道路際と内陸」に別々に適用する。この計画は別のDraft PRに実装し、鉄道独自の環境フィルターを道路一般の地域形質配置と識別できるかを調べる。

[元監査結果artifact #14498911726](https://github.com/zuizui0223/tetsu/actions/runs/44457472437/artifacts/14498911726)。sourceからの各植物行、個別元site座標や原元研究を再配布していない。
