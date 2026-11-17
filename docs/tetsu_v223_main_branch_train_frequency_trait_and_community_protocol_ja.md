# tetsu v22.3 — 鉄道植生の機能的共通性と列車交通／主線・支線の検定計画・実行監査

確認日 2026-11-16。現存する研究ブランチ `research/tetsu-slovenia-identity-guild-v210-20261116` に継続、main未変更。今回追加したGIS検定コードは `scripts/audit_Slovenia_genuine_rail_main_branch_vegetation_v223.py`、合成unit testsとActions workflowは別ファイル。実原著への交通量効果はまだ計算されていない。

## 確認済み観測事実

- [スロベニア2022/23真正原著v10.3](Slovenia_real_RAIL_microhabitat_plant_functional_v102_v103.md)：89地点、15原著区間、173同じ側のBS–BTペア。一年草BS−BT平均 +0.054703、厳格therophyte +0.058871、原著cover単位−35.179775。元ゾーン幅と調査面積の同一性は未確認。
- [同原著v11](Slovenia_species_frequency_weighting_reversal_v109_v110.md)：風散布群の対BS–BT選好は等種重み+0.189388、観測サイト頻度重み−0.098125。一般的な「風適応」の肯定結果とは言えない。
- [2独立原著v18.5/v18.8](rail_flora_taxonomic_vs_functional_fixedmargin_v185_ja.md)：HR+BA 15駅×359種では機能構成重なりの観測−同気候内固定種数・出現頻度NULLは−0.025951。別研究のParma28駅×438種は−0.012497。分類学的Jaccardはそれぞれ−0.000779と+0.001574。鉄道内部帰無であり自然地域プールの対照ではない。
- [v22.0](../scripts/audit_SI_CZ_original_continuous_trait_coverage_v220.py)で本物LEDA草丈・種子重はスロベニア285観測分類群中224（78.6%）、チェコ761中525（69.0%）が両形質一致。スロベニア139の同地点・同側BS/BT対応で両形質がそれぞれ5分類群以上。**形質分散/FRicの本値は未算出**。平均の偏りと種内の遺伝的適応を混ぜない。

## 新たに確認した最重要先行研究：輸送量は既に検定されていた

Májeková J et al. 2025 *Flora* 330:152780,
DOI https://doi.org/10.1016/j.flora.2025.152780
「Climate and traffic volume control plant species richness in railway areas」。

原著はスロバキア国内**30駅×1ha固定面積**、413植物分類群（258自生、155外来）について、気候・都市環境と旅客/貨物列車輸送量を考慮。貨物列車数が多い駅で外来種（archaeophyte + neophyte）の**alpha species richness**が高いと報告した。駅間の種共有Jaccard、距離や貨物ODで調整したtaxonomic beta、生活史分散/FDisを直接評価したとの証拠は今回見つかっていない。

**論文ギャップ:** もう「貨物列車が多い駅は外来種数が多いか」という問いは新規性が低い。新規性を検証するなら、
1. 直接の貨物輸送ネットワークと実際の駅間共通種（外来/自生）のpairwise beta（供給）。
2. 同じ輸送量での道床内の生活史CWMと草丈・種子重のFDis（環境・管理選別）。
3. 同一近傍地域プールからの期待と観測された鉄道駅同士のsimilarity（収束/分化）。
3種類を別の応答として保持すること。

## 主線と支線の対応付け

スロベニア鉄道当局の同時代 `Network Statement 2023` 公開原本 https://infrastruktura.sz.si/wp-content/uploads/2023/02/NETWORK-STATEMENT-2023.pdf は路線を交通量・経済的重要性と接続役割によって**main (G)**と**regional (R)**に分類している。AKOS公開説明: https://www.akos-rs.si/en/railways/explore/public-railway-infrastructure 。OSMは `railway=rail` に `usage=main/branch`、留置/側線は `service=siding/yard/spur` を用いる https://wiki.openstreetmap.org/wiki/Key:usage 。`frequency=50` 等OSMの電化周波数は**列車運行本数でない**。

調査時期に近い地図は Geofabrik 2024-01-01 Slovenia PBF （元本280842431 bytes、SHA b572258fe000d6dfc76dc63afe5456121c5f463995e132a5dcd79957d819ae39）。元89地点の実100m WGS84原著GPSに、100m以内の確定main/branchタグが一意な場合だけ付ける。同駅ヤード等でmainとbranchが50m以内競合したらunknown。使用タグが無いまたは距離精度不足の場合もunknown。OSMの`usage`は当局の2023 G/R区分と同一保証なし、別途突合。

### v22.3実行ゲートの現状：地図が取得できない

- [初回GitHub Actions #43872122997](https://github.com/zuizui0223/tetsu/actions/runs/43872122997) は独自解析ブランチに2024 PBF元キャッシュがなく、`actions/cache/restore`段階で失敗。
- [再実行 #43872362794](https://github.com/zuizui0223/tetsu/actions/runs/43872362794) では同元ファイルの正式URLからSHA検証付きで再ダウンロードを試したが、runnerから `download.geofabrik.de:443` 接続失敗（curl exit 7）。原著植物sourceのSHAは通過している。
- **重要：main/branchの地点数・一年草率差・pairwise Jaccard効果は未算出。** 実数を捏造せず、地図原本にアクセスできない場合はfail-closed。これは生態仮説の支持・反証を意味しない。既存のv15.0原ソース成功は古いActions上の別ブランチで得られているが、現在の解析ブランチから同じPBFキャッシュを読める保証はない。

## 具体的な事前予測

- **選別（生育地効果）**：annual CWM ≈ BS−BTで正。main>branch/交通本数増加で幅（e.g. z log height, z log seed massの RaoQ/FDis）が低下する仮説。ただし路線区分は列車速度・実剪定履歴ではない。参照種プール/区域・調査面積・種数調整が必要。
- **到着（輸送効果）**：同一貨物経路/荷姿に沿う種名Jaccardが、地理・地域プール・道路等を条件付けても増すか。旅客通過数との別モデルにする。単なる同一路線/鉄道最短経路の効果を実貨物ODとみなさない。
- **駅間機能収束（機能beta）**：同じ選別強度で平均や分散が似るとしても、種名は地域種プールで違う可能性がある。逆にCWMの同方向移動だけでは、実際の駅間機能組成の重なりが増すとは限らない。各地域の独立した自然/非鉄道の種プールが必要。

主線対支線は**生態学的な類型**、路線を通過する実列車数 (trains/day) は**連続曝露**。旅客発着列車数と、貨物の通過本数・貨物重量・輸送起点終点は別変数。都市/気候/路線距離/管理強度/観測年/駅面積は必ず交絡として考慮。

### 明確な停止ルール

OSM分類できる地点が少ない、main/branchが各3以上の独立原著区間を持たない、調査年度の交通量が得られない、公開0m単位の信頼できる位置がない場合は関連係数を発表しない。3つの独立の生態的な問題をひとつの汎用「鉄道選択強度」に圧縮しない。
