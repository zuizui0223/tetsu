# tetsu v28.8 — 非欧州原著で100m²対100m²の鉄道vs周辺植生の機能的βを検定

更新 2026-11-19 JST。source-workflowが真正実行されるまでは**生物学的な新しい肯定結果を報告しない**。原著データの真正性・元地理単位・欠損処理を固定する。

## これまでのエビデンスと、新検証の独立性

- v28.5の**欧州4原著**は、元site×taxonを完全固定し、元LEDA canopy height(m)とseed mass(mg)の種への割当てだけをランダム化した帰無に対して、元地理群内の地点間 CWM 2D distance 観測−trait nullが4/4正（草丈も4/4正）。**鉄道内部**の空間分化であり、非鉄道対照の差ではない。
- v28.7のチェコ鉄道14駅対近傍CzechVeg70地点では広いnontransport文脈との距離は鉄道側が小さいが、非交通系参考元plotのHABITAT開放地明示2/70・明示的砂礫0/70。厳密な非鉄道open/gravel対照にはsupportなし。因果的均質化の証拠とはしない。
- 今回使うのは**別の1原著**Rashid et al. 2021, Kashmir (India), DOI https://doi.org/10.1111/1365-2664.13961。原Dryad https://doi.org/10.5061/dryad.0p2ngf20h、2014/2017反復。**同じ31site×2014/2017データを後年解析したHaq et al. (2024), DOI 10.1016/j.envdev.2024.101064を新しい独立原著とは数えない**。

## 新しいSOURCEにのみ基づく検証計画（解析前固定）

1. 原Dryad Matrix_plot.xlsx=122187bytes SHA256 529f3b97521933018e12623b61f531f012c67241cae2602a6675dbd534192474。真正31地点・2014/2017×4 zones 248 rows・123元植物列・欠損元species cells 4を厳密に照合。欠損plotを0 speciesに変換しない。4 zones×2年をそろえた原**駅構内15・駅間12**地点を使用する。
2. 著者原T形配置の軌道Track vs **同じ元地点・年・100m²**のInterior Plot1またはPlot2を別々に比較。Embankmentは元全zone cohortの同一site基準に使うが**非Rail対照ではない**。Interiorは元著者の周辺植生であって確認済みopen/gravel基質ではない。
3. 元Oldenburg LEDA canopy_height.txt SHA cbf5e761390d78a843512dfc1139acfc1154d46d6c1995a827354840a9f78a92、seed_mass.txt SHA 34aa12c333fc1ba6de2e560a495ac5d693aa0722a09ecf68418ad94ac1c19cef の二名法厳格一致のみ。形質欠損を0にせず、South Asiaでの欧州LEDA被覆率を記録。z(log height, log seed)座標は以前v24.5の**欧州3研究原著622一致種**で固定した log mean −0.5802365754696401 / +0.2936113186962433、SD 1.3801000545693773 / 2.1657660559995824。Indiaを見て再標準化しない。
4. 本物の各year×source zoneのknown species比を算出し、原siteを2年通して固定。あらかじめ**両LEDA形質1 plot 5種以上、元speciesの30%以上、駅構内10実地点または駅間8実地点、適格全地点内の両LEDA 20固有種以上**の閾値に届く場合のみ値を計算。未達なら\`NOT_IDENTIFIABLE_TRAIT_SUPPORT\`で停止、ゼロ生物効果とはしない。
5. 支持域を満たせば、原著の同一地理群内で元siteごとにpresence-equal CWMを作り、**異なるサイト間の平均2D / height 1D / seed 1D CWM距離**を同一元site集合のTrackとInterior間で比較。\`D_track-D_inland\`を主値、\`log(D_track/D_inland)\`を副値として出力。2原年は同siteの反復なので独立研究ではない。元siteでleave-one-out、各site×zoneでLEDA既知種数を均等化した100 technical resamplingだけ感度解析する。ペア数を独立元研究数とみなさず、technical drawの順位もfield p-valueとしない。
6. 原元species×plot/GPS/元LEDA raw種一覧はGitHubへ再配布しない。真正source読取はActions内に限定、公開aggregate provenance/被覆率/識別可能な効果のみをartifactにする。

## 広域メタ解析候補・真正元行列の追加取得ゲート

| 真正外部原著 | field units | 解析可能性を阻む事項 | DOI |
|---|---|---|---|
| Rakotoson et al. 2025, France | 2021年394植生2m² plot、Track128、Embankment266、255種 | trackとembankmentは**両方鉄道**、直接rail vs nonrailでない。データは著者request。元plant×plot行列未取得 | https://doi.org/10.1111/wre.70025 |
| Dylewski et al. 2024online/2025volume, Poland | 廃線25地点、隣接草地25地点 | **廃線**を稼働中バラストと統合しない。元raw行列は著者request | https://doi.org/10.1002/ldr.5383 |
| Li et al. 2026, China–Laos | 14駅、84植生/土壌/菌類plots、495km transect | railway-station対nearby undisturbed。本当の個別ballast・同面積plotとoriginal source matrixの取得未確認 | https://doi.org/10.1016/j.jenvman.2026.130380 |
| Winkler et al. 2026, Czech | 4地点×5鉄道/農地habitats、10m² 2021年3回 | 実plot×species matrix非公開、supp S1は学名/略号でplot occurrenceでない | https://doi.org/10.3390/land15040523 |

欧州内のtrait-label null +0.06という値と、厳格非鉄道対照との\`log(Drail/Dnonrail)\`を同じ「メタ効果量」として併合しない。4原著/5国の欧州帰無とインドの2年反復、フランスのtrack vs embankment、ポーランドのabandoned vs grassland、China–Laosのstation vs undisturbedは各々異なるestimand。

厳格な研究間メタ解析ができるのは、**真に独立した複数原著で、同じエリア・気候・調査年・field extent・物理対照・形質軸の種構成βを再計算でき、publication単位のばらつきが得られてから**。それまでは段階別のevidence mapに限定し、無理な平均CIや普遍性主張はしない。

実装: \`scripts/validate_Kashmir_out_of_Europe_rail_nonrail_height_beta_v288.py\`, \`tests/test_Kashmir_original_rail_nonrail_LEDA_beta_v288.py\`, \`.github/workflows/tetsu-v288-Kashmir-outofEurope-original-rail-nonrail-trait-coverage.yml\`.
