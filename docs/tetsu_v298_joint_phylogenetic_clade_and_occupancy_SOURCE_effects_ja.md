# v29.8 — 草丈の地理的構造は、真の植物系統群と種の出現頻度を同時に固定しても残る

解析時点：2026-11-19 JST。元publisherの鉄道植物原著4研究、独立のMIREN mountain ROADS field program 1件、Oldenburg LEDA高さ・種子重量、V.PhyloMaker2 original LCVP/WP/TPL 3樹の真正SHA/MD5を使用。

SOURCE source-runの最初の完全な原データ計算（https://github.com/zuizui0223/tetsu/actions/runs/44547195906 ）では、4鉄道出版×3treeの既存v29.6 clade75 original source係数12点と、MIREN v29.7 tree×primary/strictの元Region集合・site数・roadside/ inland とその差の基準係数を全6パネル再現した後、新しいclade75×source occupancy bin joint nullを完走した。最後のGitHub CI checkerに出力キー名「source_panels / original_panels」の不一致があり、そのCI自体はFAIL。**新数値はSOURCE計算完走・公開集計から独立検算済みだが、修正したCIの最終完走を待っている**。元真観察行・数式・source support閾値は一切変更していない。

元SOURCE集計artifact： https://github.com/zuizui0223/tetsu/actions/runs/44547195906/artifacts/14540484122 。資料は4原著出版社の元生物観察＋LEDA2原著源ハッシュ＋published元系統tree原バイナリGit blobハッシュを固定し、元個別site/taxon/GPSは公開していない。

## 数理検定の意味

元原著で観察された植物のsite×species配置はそのまま固定。系統樹 root-to-tip 元中央値75%深度で分けた既存の真単系統cladeと、元species site出現頻度に基づく5分位binを**交差**させ、2種以上が含まれる小群の内部のみで、元のheight+seed massを一体のvectorとして置換する。交換可能な種が少ないことが見かけの結果を作らないよう、元site×taxon incidenceの25%以上、異なる交換可能群5以上、元地理区分2以上（MIRENは両zoneに同25%以上かつ2実Road以上）を事前基準とした。

元RAIL4出版で系統群×種頻度群の**2種以上を持つ元交差小群**は34〜62個で、元site×plant出現の**84.3〜95.5%**が実際に交換可能。元MIREN Roadside+Inlandで適格な元Regionでは交差小群18〜58個、元site×taxon出現の**68.7〜96.8%**が交換可能。帰無が少数の種だけしか動かず機能しない状況ではない。

2 fixed seeds×199＝398 technical trait assignments /適格な元study×tree×mode。これらは種ラベルの組替え分布であり、野外field p値・独立国・全球予測区間ではない。3樹は命名体系の感度で、Rail4/ MIREN Road1の野外独立数を増やさない。

## 真正RAIL4元出版の結果：草丈CWM地点間距離（観測−帰無）

| 真のphylo tree | v29.6 系統群75%のみ制約、4独立出版等重み | v29.8 系統群75%＋出現頻度5分位を同時制約、4独立出版等重み | 新joint nullで正方向 |
|---|---:|---:|---|
| LCVP | +0.0566953 | **+0.0440228** | **4/4出版** |
| WP | +0.0554290 | **+0.0406410** | **4/4出版** |
| TPL | +0.0589976 | **+0.0464848** | **4/4出版** |

| 元TRUE独立RAIL出版 | LCVP joint | WP joint | TPL joint |
|---|---:|---:|---:|
| CZ2024, 39駅全域 | +0.0264053 | +0.0196207 | +0.0277717 |
| IT2020, 28駅全域 | +0.0203368 | +0.0119574 | +0.0187954 |
| HR+BA2017, 15駅（2国1出版） | +0.0815854 | +0.0923807 | +0.0897059 |
| SI2026, 89個100m道床地点（駅全域ではない） | +0.0477636 | +0.0386052 | +0.0496661 |

**観察**：植物系統の同cladeと観測species occupancy quintileの両方を固定すると、高さの正の分化超過は縮小するが、元4/4の符号は維持。技術順位はItalyで0.115〜0.208と強くないため、4原著とも独立な「統計的有意」とは報告しない。元野外出版の符号4/4だけで全球普遍性を証明しない。

## 独立のMIREN Mountain ROADS／Inland の真系統同時制約

元MIREN primary field supportは8 Region/22 Road/272同年・同site物理ペア。元tree直接tipに制限した実SOURCEでLCVP5 Region/ WP5 Region/TPL6 Regionがsource標本支持基準を満たした。Strictは元5 Region/13 Road/192siteから、source樹の原tip supportによりLCVP5/ WP4/ TPL5 Regionに減少。元外部研究1件に属すregionであり、5〜6独立研究ではない。

### Primary （各元site plot 5両形質種・元観察25%以上）

| tree / 実元Region | v29.8 Roadside（観測−真joint null） | 同元Site Inland（観測−真joint null） | Roadside−Inland |
|---|---:|---:|---:|
| LCVP, 5 Region | **+0.079384（5/5正）** | **+0.105060（5/5正）** | **−0.025676** |
| WP, 5 Region | **+0.077671（5/5正）** | **+0.103732（5/5正）** | **−0.026061** |
| TPL, 6 Region | **+0.061098（5/6正）** | **+0.083600（6/6正）** | **−0.022503** |

### Strict （各元site plot 10両形質種・元観察50%以上）

| tree / 実元Region | Roadside（観測−真joint null） | Inland（観測−真joint null） | Roadside−Inland |
|---|---:|---:|---:|
| LCVP, 5 Region | +0.078864（4/5正） | +0.095693（4/5正） | −0.016829 |
| WP, 4 Region | +0.073967（3/4正） | +0.096188（3/4正） | −0.022221 |
| TPL, 5 Region | +0.076786（5/5正） | +0.085363（4/5正） | −0.008577 |

Roadside−Inland joint差がPrimaryの3原樹とも平均負、元Regionを1つずつ除いた平均も負。ただし**個別regionのRoadside−Inland差の符号は混在**（Primary元LCVP 2/5負、WP3/5負、TPL4/6負）し、Strict TPLのregion除外平均は−0.02531〜+0.00949とゼロをまたぐ。従って「全道路は普遍的に機能が収斂」と結論してはいけない。

## マクロ生態学的な一文

**標準化された人工生息地での植物草丈の地域分化は、真のphylogenetic cladeとその地域の植物種出現頻度によってすべて決まるわけではない。近縁同群・似た出現頻度の種の中でも、形質と地理的な種配置には追加の非ランダムな関係が残る。** これは元欧州鉄道4原著と、同とは独立のMIREN道路・内陸植生1 field programで共通する観察学的帰結。

この共通パターン自体は、鉄道特有の適応・列車種子散布・草刈り・突風や風圧といった機構の証明ではない。元物理rail ballastと同条件のopen gravel nonrailを同調査プロトコルで比較する外部原著は未取得。日本・アフリカ・南米・豪州など新たな鉄道元field研究が必要であり、5研究を無理に一つの効果量で混合したrandom-effects prediction intervalは現時点で作らない。

## 再現性と資料
- v29.8 source helper：scripts/phylogeny_occupancy_joint_label_null_v298.py
- RAIL元study：scripts/analyze_4original_rail_joint_phylo75_occupancy_null_v298.py
- MIREN元same-site source：scripts/analyze_MIREN_roads_inland_joint_phylo75_occupancy_null_v298.py
- 10合成数理不変性テスト：tests/test_rail_road_joint_phylo_occ_source_v298.py
- 同時原著認証workflow：.github/workflows/tetsu-v298-original-RAIL-and-MIREN-road-inland-joint-clade-occupancy.yml

GitHubでの追試結果と元site×taxonの品質はSOURCE provider本体に閉じる。公開物は元source paper×treeと元MIREN region×tree/panelの集計およびSOURCE支持域だけ。
