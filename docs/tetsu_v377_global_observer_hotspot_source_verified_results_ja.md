# v377：全球鉄道網の植物記録種数予測改善は撮影者数の大きい場所だけによるか（2026-11-24）

## 真正全球観測と鉄道の出典
- 元植物原データ：iNaturalist2026-10-27 S3公開photo-linked世界全行 **288,721,938**、2025年の元research-grade/種/位置精度等条件通過**7,827,936植物レコード**、[v369 full world SOURCE SUCCESS #45375453375](https://github.com/zuizui0223/tetsu/actions/runs/45375453375)、元観測圧縮SHA256 `1c48a9646f3bf6d56ffdfdc91a1a176c173c98014830ad48fa162b2803c3b279`。元iNaturalist公開写真の種数と野外の完全植生種数は異なる。
- 元地球鉄道網：Overture2026-11-18.0全128 original transport Parquet、通常線路候補1,536,226元GERS線分、2,237,043.481kmのWGS84 measured length、[v370 source SUCCESS #45373985228](https://github.com/zuizui0223/tetsu/actions/runs/45373985228)。観測者/種の実GPS・個別GERS線分はGitHubに再掲載なし。
- 原自然地理：WorldClim2.1 BIO1 元ZIP SHA256 `00513224583665ec0f2f955a4ec252730c4deb2004cce9e793492a3f26df4dcf`、Natural Earth world original SHA256 `718ac51d68ed2bf5de4f563fff2e2e045db706c1d3222805f3e317140a64213b`。
- **[v377 genuine original world source success #45384209965](https://github.com/zuizui0223/tetsu/actions/runs/45384209965)**、合成7/7テスト成功、元v371正規化済みλ=.10・1000km地理ブロック重み・6大陸LOCOと同じM0−M1が0%除外で完全再現。元sourceの派生結果[artifact #14914249026](https://github.com/zuizui0223/tetsu/actions/runs/45384209965/artifacts/14914249026)、digest `sha256:7882437611e06027035cc4803d367d908ba5c60ca7ac46f55028cd4fa880bf7f`。

## 固定した観測努力の感度検証
v371/v375の全球植物の記録種数の増分を見た後の**事後探索的反証**。元v371で写真観測50以上・distinct撮影者5以上・植物20種以上、陸地中心、WorldClimの年平均気温有効を通過した元植物cellから、**各大陸の元distinct植物観測者数が多い上位5%、10%のcellだけを除外**。上位判定に植物種数、鉄道km、鉄道−植物の差分、残差を使っていない。同点は元グリッドIDによる決定的順序。100km・200km同じ元解析を繰り返したが、空間粒度は独立field生態研究ではない。

M0は元の植物観測数・撮影者数・WorldClim年平均気温と二乗・絶対緯度と二乗（Jaccardは元植物種数平均/差・方向も）に基づく。M1はsource格子内GERS railkmを追加。真source観測格子を落とした後、隣接景観Jaccardは両元cellが残るsource pairのみを残し、元植物相を都合よく再集計しない。6大陸を丸ごと伏せたout-of-continent block weighted ridge λ=.10のMSE差を報告する。

## 実験結果：地域記録植物種数（log1p）の予測増分

| 元植物の撮影者上位を除外 | 100km源有効植物格子／隣接Jaccard | **100km種数 ΔMSE**／改善した大陸 | 200km源有効植物格子／隣接Jaccard | **200km種数 ΔMSE**／改善した大陸 |
|---|---:|---:|---:|---:|
| **元source除外0%** | **3,455／4,736** | **+0.004463745／4/6** | 1,442／1,953 | **+0.006614244／3/6** |
| 各大陸撮影者上位5%除外 | 3,279／4,278 | **+0.004846165／4/6** | 1,368／1,771 | **+0.007035024／4/6** |
| 各大陸撮影者上位10%除外 | 3,107／3,882 | **+0.004446064／4/6** | 1,296／1,622 | **+0.007840856／4/5** |

**200km 10%除外の種数は源評価可能なwhole-continent foldが5つへ減った**。最後の+0.00784は6大陸の値と同じ条件ではないため、感度の増分を傾向の強化と解釈しない。

元地域全体の植物相Jaccardは100kmが除外0/5/10%で **+0.000147667／+0.000135489／+0.000136201** と小さいながら正（改善6/6→5/6→5/6）、200kmは **−0.000061200／−0.000018303／−0.000012033** と非正、やはりglobal植物相均質化の予測根拠は不安定であった。

## 判定
- **何が残ったか**：元公開植物の撮影者が極端に多い地域（大陸内トップ5～10%）を除いても、全球実鉄道kmを加えると地域観測植物種数への外挿予測上積みは消えない。したがって**ごく少数の写真観測ホットスポットだけが元v371の肯定的予測情報を作った**とは言いにくい。
- **何が解決していないか**：写真観測者数・撮影件数は検出率と植物分布の全ての地理的偏りを表さない。鉄道と共存する都市化・道路・人口密度・土地利用・生態系観測可能性は未解決。独立GHSL元source都市化割合で調整した[v376c全球実結果 #45384121307](https://github.com/zuizui0223/tetsu/actions/runs/45384121307)では鉄道kmによる種数予測改善は弱まりながら平均的に残るが、6/6大陸共通ではない。
- **鉄道敷の花相を直接同定したわけではない**：source100km/200km Jaccardは地域景観全体である。個別元植物・観測者GPSは源geoprivacy不明、種の在来外来・機能形質・鉄道バラストの実温度/実列車輸送をこのデータから導出していない。
- 0/5/10%除外条件は元v371肯定結果を見た後の技術的な感度解析であり、確認的に事前固定した世界実験のp値にはしない。高努力撮影地点を除外すること自体、植物相の真の代表性を保証しない。

## 再現
`docs/tetsu_v377_global_railkm_botanical_observer_hotspot_frozen_ja.md`、`scripts/analyze_global_railkm_plant_observer_hotspots_v377.py`、`tests/test_global_railkm_plant_observer_hotspots_v377.py`、`.github/workflows/tetsu-v377-global-observer-hotspot-source-robustness.yml`。元source観測行、種名前/rare GPSや鉄道原路線UUIDを本リポへ転載しない。研究Draft PR [#60](https://github.com/zuizui0223/tetsu/pull/60)、main未変更。
