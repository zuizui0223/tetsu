# v315 Thermal-affinity reference comparison — design and implementation plan

> **For agentic workers:** Use executing-plans to implement the tasks below; tests precede numerical implementation. The scientific design is frozen in this commit before community means or habitat contrasts are calculated.

2026-11-22 JST。基準コミット `7fe99dbdb5fb0167818d97f7eb74f400b7410ee3`（v314/PR44）。PR45の気候仮説A・生息地置換Bの引き継ぎを確認した。

## Goal / scientific scope

同一の外部温度指標を既存の鉄道原表へ接続し、欧州5出版研究の温度親和性と、地域の道路・都市参照に対する差を求める。カシミールは欧州指標の域外適用として分け、西オーストラリアは被覆率を示す。v314の350共通種には限定せず、各資料の確認可能な全種名を用いる。

PR45のAは「寒い地域ほど、鉄道と比較可能な近接非鉄道の温度親和性差が大きいか」である。現在、同年・同じ地点・同面積で対にできるのはカシミール1地域だけ。SIのBS/BTはいずれも鉄道内、都市参照やMIREN山地道路は気候・生息地が対応しない。Aの地域間気候交互作用は適格性不足として未推定にする。国名や緯度、駅の標高だけを気候の代用として回帰しない。今回の記述的な温度組成比較をAの検証成功へ置き換えない。

## Approaches and choice

1. 気候交互作用を既存の不一致対照で計算する案は、比較対象と気候が交絡するため採用しない。
2. SI内の温度勾配を再解析する案は原著がすでに示し、独立した地域も増やさないため採用しない。
3. 全種名・統一温度指標による複数研究の温度組成と参照差を実計算し、数値のある種・温度不関種・不明種を区別する。これを採用する。気候仮説の前提を具体的に調べ、比較可能な新しい近接対照の取得を次の識別条件として残す。

## Inputs and decisions made before effects

既存12原本とsource adapterはv314のSHA・種名・不確定マスクを継承する。スロベニアの出現表とコードブックの不一致も継承する。解析する個体数・被度ではなく、各種の記録地点頻度で重み付けする。

新しい形質原本：

- **主定義：Tichý et al. 2023、固定Zenodo v2、DOI10.5281/zenodo.7427088。** SHA256 `b3fa7c3a6e9d356e52a45472549f807967bd0c3189e58e7e8a22e0389c29cc4e`。尺度1–12。固定ファイルは8,908分類群、T数値7,402、x1,147、NA359。xは温度不関であり、欠測・0・中間値に変換しない。
- **別定義：EIVE1.0、DOI10.5281/zenodo.7534792、SM08。** SHA256 `9551b585ba13dfe3fddbd2287f65fc98ac9f4ec5f80b9ae6f6b985bed73fb4a6`。尺度0–10。14,835概念、T数値14,496。原著の構築手順では一部の地域体系のxに体系内平均を与えて統合する。その仮定を引き継ぐ別指数であり、主定義のxを埋めた同じ尺度とは扱わない。

効果計算前の種名のみの被覆確認では、主5鉄道研究のTichý数値は約40–49%、EIVE数値は約90–97%。これらの被覆率を見て両定義を採用したが、温度平均・参照差・係数はまだ計算していない。両体系は地域の原指標を共有するため独立した生態学的反復ではない。

形質側は明示された種行のみ、agg/aggr/群/節/交雑・種以下の行を除外する。異なる数値・状態が同一キーへ正規化される場合は不明にする。EIVEの種ランク欄がSpeciesでも名前にsubsp等があれば除外。新しい異名推定、種以下の平均、属・科からの補完は行わない。同じ種はすべての生息地で同じ形質値を持つ。

## Estimands

各資料jの記録地点頻度fと、数値が確定した種集合Nについて

`mu_j = sum(f_js * T_s for s in N) / sum(f_js for s in N)`。

これは**数値の温度指標がある種に限定した、地域の出現頻度加重平均**。地点ごとの群集平均を等重みにした値ではなく、種数の多い地点ほど多くの出現記録を寄与する。SK30の周辺頻度でも同じ定義が計算可能なため採用する。°C、熱耐性、微気候の実測値とはしない。

各資料で、数値・温度不関x・真の欠測・照合不能/競合の種数と頻度和の比率を出す。x比率は全確認可能種の頻度和を分母にした値と、数値+xの既知状態だけを分母にした値を併記する。xから生理的耐性や可塑性を推定しない。

## Fixed comparisons and scenarios

主5地域の参照差 `D = mu_rail - mu_reference`：CZ対MIREN道路21地点、IT対MIREN道路60地点、SI対植樹枡50区画、HRBA対植樹枡150区画、SK対Malacky都市17区画。これらは記述的な地域参照差として別々に示し、同じ処置効果として平均・メタ解析しない。

域外：IN2014の16駅Track+Embankment対同じ駅のInterior1+2、およびIN MIREN道路との別参照差。WA120は鉄道保留地の形質被覆確認のみ。各比較で鉄道・参照の分類学的不確定マスクの和集合を両側へ適用する。数字のない形質は0として入れない。

以下の5条件をすべて保存し、結果で追加・選択しない。

1. `tichy_numeric`：主定義、出現頻度加重。
2. `eive_all`：別定義、出現頻度加重。
3. `tichy_common_numeric`：両体系で数値がある同じ種集合、Tichý尺度。
4. `eive_common_numeric`：3と同じ種集合、EIVE尺度。3/4で種集合の変化と指標定義の変化を区別する。
5. `tichy_species_equal`：主定義、資料内で各記録種を1ずつ。頻度加重への感度。

共通数値集合は二つの**形質原本で数値がある種の集合**であり、鉄道と参照の両方で出現した種だけには絞らない。各生息地固有の種も、形質の条件を満たせば含める。

全候補条件で数値種10以上を表示の最低条件とする。満たさなければ平均・差はNAとし、被覆率と理由を残す。被覆率に単一の80%等の閾値を置いて、全植物相への代表性が保証されたとはしない。

カシミールの実際の対応については、上記1/2のみで、各駅の両側に数値種3以上ある対を計算する。8対以上なら、駅等重みの差の平均、中央値、正負の対数を副結果として出す。区画の反復や2017年を独立地域へ水増ししない。主たる新規結果をカシミール内の細分化に置かない。

この対の採用数・被覆率は指標体系ごとに示す。採用駅が違う場合、二つの対平均の差を尺度だけの感度とは呼ばない。3種・8対は最低表示支持条件であり、代表性の保証ではない。

## Missing-value sensitivity

EIVEの全確認可能種の平均を仮定した場合の感度幅を別に計算する。鉄道と背景の全種の重みをそれぞれ総和1にし、係数 `a_s = w_rail,s - w_reference,s` を作る。既知種は固定し、未対応・未採点の同じ種へ両側で同じ未知値 `T_s in [0,10]` を割り当てる。`sum(a_s*T_s)` の上下端を係数の符号から計算する。同じ種に生息地別の値を与えない。

これは、すべての未採点種へEIVE尺度値を割り当てられると仮定した数学的な感度包絡であり、信頼区間・実測温度範囲ではない。数値がある種の平均と推定対象が異なる。EIVE構築時のx処理や既知値の誤差までは覆わない。Tichýのxをこの計算へ便宜的に補完しない。

数値種だけのEIVE差がこの包絡内に入る必要はない。推定対象が異なるため、その点推定の誤差棒として描かない。

## Inference and reporting constraints

種を独立反復にしたp値、SK周辺頻度からの独立二項シミュレーション、種ブートストラップを行わない。本段階では全地域の数値・符号・被覆・定義感度を記述する。適応、種子輸送、気候による分布域拡大、代替生息地での個体群存続を主張しない。原行列・座標を転載せず、集計・コード・出典を保存する。

## Architecture / implementation

Python、既存NumPy/SciPy/openpyxl、PDF抽出はPyMuPDF1.24.14。旧source adapterを再利用し、形質・集計・報告を分離する。ローカルはv315_repoの独立したスナップショット、GitHubはv314から新しい研究ブランチを作る。元リポのmainへマージしない。

### Task 1: Strict trait and paired-source adapter

Files: `scripts/source_thermal_affinity_v315.py`, `tests/test_source_thermal_affinity_v315.py`, `data/manifests/v315_thermal_trait_sources.json`.

- [ ] テストを先に追加：xとNAの区別、数値0の有効性、種以下の誤ランク除外、重複値競合、同一種の単一値。
- [ ] 失敗を確認してから2原本のハッシュ・列・行数を検証する読み取りを実装。
- [ ] Interface: `load_traits(directory, normalize)` → `{'tichy': {'values':dict, 'statuses':dict, 'limits':[1,12], 'audit':dict}, 'eive': {...limits:[0,10]}, 'metadata':dict}`。
- [ ] `kashmir_pairs(data,normalize)` → in-memory list of `(rail_set,interior_set)` for the same16baseline stations; export no IDs or rawmatrix。
- [ ] ソースのみの統合確認。効果はTask2後まで実行しない。

### Task 2: Aggregation and contrasts

Files: `scripts/analyze_thermal_affinity_reference_v315.py`, `tests/test_thermal_affinity_reference_v315.py`.

- [ ] 手計算fixtureで加重/種等重み、頻度定数倍不変、同一構成で差0、不確定マスク共通除外、同一未知種の打ち消し・感度端点、尺度外値と低支持の停止、対応駅等重みを検査。
- [ ] REDを確認してから集計実装。Interface: `summarize(freq,trait,exclude=(),min_species=10,species_equal=False)`、`contrast(rail,reference,trait,exclude=(),...)`、`shared_unknown_envelope(rail,reference,values,limits,exclude=())`。
- [ ] 全条件のJSON・CSVを集計だけで出力。Aの気候交互作用は`not_identified`を明記。

### Task 3: Results and independent validation

- [ ] 固定した全条件を実行し、別担当が加重和・差・感度端点・種集合を再計算する。
- [ ] 同一尺度内の図、地域ごとの参照差、xと未対応の被覆を可視化し、日本語報告へ含める。
- [ ] 原資料マニフェストと形質由来、再現手順、限界を保存。図を目視確認。
- [ ] Draft PR、全CI、公開blob照合まで完了する。

## Review focus

- 同じ種へ生息地ごとに別の形質を与えない。
- x・0・NA・未対応を混同せず、異なる尺度を混ぜない。
- 数値のある種の条件付き平均を、全群集や実測°Cへ拡大しない。
- 同一不確定マスク・同一種集合を、直接比較する両側へ適用する。
- 気候・生息地が対応しない外部参照を、Aの地域間対照へ数えない。
