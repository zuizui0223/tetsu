# tetsu v28.4 — 独立チェコ鉄道：散布戦略・生活形によって「同一路線の植物種共有」は変わるか

更新：2026-11-18 JST。真正一次資料を再取得して行った探索的・反証的な二次解析であり、2026年の現場植物採集・繁殖体捕捉を新規に実施したわけではない。直前の v28.2–v28.3 で既知の本物のチェコ35駅・元5路線区間の地理距離対応を、そのまま固定して散布形質で層化した。元生態学的測定はチェコ鉄道駅の駅全域フロラで、Sloveniaの100m BS/BT道床野外区画ではない。

## 1．原本・対象群・ソース改訂の監査

- 鉄道植物：Kutlvašr et al (2024), Preslia DOI https://doi.org/10.23855/preslia.2024.247。作者オリジナルS1植物出現Excel 79,261bytes、SHA256 df6bc0b949e7f6a5cb57755dc8f9c86187223f66ddfbef959cf0e03ef7c46b06、訂正後S2 13,533bytes SHA256 00c7cfdf5047854ba99d47dc8ec8c55b1b7df72731ce647c22b286ddef24a5a2。真正39駅・5元rail stretches、2020–22年、元RECマーク6,941 site×taxon陽性、761学名。論文本文764taxaとの差3種を捏造補填しない。
- 第三者のチェコ鉄道SR70公式台帳再生成2026年駅代表点：CSV source SHA256 fc43a1969af64f17285c7c6e3a30b00c3994d0ba1cfda6db09d8e4baaa166019。正式駅全名を原S2と一致させられた35/39本物の駅。残り4駅は不明として除外。植物の元現地plot GPSではない。
- 真正Pladias元学名・散布戦略とRaunkiaer生活形公開source: Pladias https://pladias.ibot.cas.cz/ 、Sádlo strategy feature 141とlifeform feature 136。元v14.3の原XLSX SHAはdispersion 135fff21b339b0774dd5c1ff3504240175b72be384b4813d747a2905527ec636、lifeform 6076383e6b206a74a0f23d6f8f8691c94f4c83eea8ea05a19b2d59f08e8c856cだったが、2026-11-18に再取得したURLはバイト単位で異なるXLSXを返した。最初の新source実行でdispersion 20dbb08b99e08751e57d337313bdca6cbfc0ec5dd833290363347534deea203e、lifeform 3ee134f767eb556aa08a99c1465e9616dd6df094d7e72e80fde8bd57ef209742であり、再取得ではさらに変動した。動的に生成されるXLSXのファイル梱包の違いか中身の変化かを決め打ちしない。
- 源原表の必要条件はSádlo strategy 3147行、Raunkiaer lifeform 4599行。元全761 source taxonのうち544 taxaが両trait明示で一致（元39駅全体5,081/6,941 incidences 73.2%）。元位置確認済み35駅では原S1生物記録6,142件のうち4,489件が両traitに対応（73.0869%、全元source taxa 740）。元形質ラベルの変動を検知するため、原学名→元カテゴリ値のソート後正規化をUTF8でSHA256化し、ダウンロードXLSXのZIP時刻差でなく意味的変更を検知する別source auditを追加した。意味的な対応表の正規化SHA256は 3ad00d97a88c41cb9b6c4e8775a8cc69ae0b70684a5247bcb1c4d7383e2204d1。2回目以降の別ソース取得でこのハッシュを必須にし、XLSXの生バイトは各実行時点のdigestを記録する。元v14.3の同じ学名→各trait値の完全一致までは再現データが残っていないため保証せず、v28.4で取得した版を独立の語彙固定対象として扱う。

元Pladiasの「Epilobium」または「Phragmites」と元sourceが呼ぶ散布戦略型をwind-associated proxyと分類。残りのsource戦略型をotherとする。これは**実際の種子飛距離、風による着地速度、列車風圧への実測抵抗、車両への付着能力を測定したものではない**。生活形のtherophyteは一年草のstrict寿命とは一致しない。欠損217分類群を非風散布に割り当てない。

## 2．目的と地理条件

前v28.2で凍結した元35本物の地理位置・5source stretchesの同駅i、同一区間の別駅j、異区間の別駅kの距離対応をそのまま使用。最大160km、anchorからj/kへの直線距離差<=10km、原植物名を参照駅選択のために使わない。真33station source anchor・元5source stretches（元本体のREC 761全taxaでは元同線−異線Jaccard差−0.001731）。対象33アンカーを変えず、544 taxonomy trait-known種を主対象にし、風関連と他の2群・therophyte×風・therophyte×その他・非therophyte×風・非therophyte×その他の4細群に分類して、同一j/kのgroup Jaccard(i,j)-Jaccard(i,k)を求める。

各source stretch内でanchor等重み平均→元5stretches等重みで集計。source観測分類群0同士のJaccardを都合よく0とせず未識別として除外（しかし元main群は33/33anchor全てで比較可能）。wind-minus-otherは両群が同時に実測できる**全く同じ元アンカー**上で差の差をとる。

## 3．v28.4 実源群集のJaccard差

| 元実植物・真の33source anchor／5元区間 | 同一区間−異区間のJaccard類似度 |
|---|---:|
| Pladias両traits既知全544taxa | **−0.007631** |
| Wind-associated strategy proxy | **−0.019782** |
| Other strategy proxy | **−0.004017** |
| **Wind−other：同じ本物の駅・距離組の差の差** | **−0.015765** |
| therophyte×wind | −0.013133 |
| therophyte×other | −0.008315 |
| non-therophyte×wind | −0.023604 |
| non-therophyte×other | −0.000952 |

**正の「風散布戦略の方が同一路線で共有される」という予測は、今回のプロキシでは出なかった。** 風関連proxyはむしろsource元5区間等重み平均で負。ただしこの係数の点推定を「列車が風散布種の拡散を抑制している」「風散布は距離が短い」と物理的メカニズムに翻訳しない。観測は路線系統別・地域別の種プールを調整する不確実な二次推定である。

## 4．元site×guild richnessと元species source frequencyを固定した帰無

4source guildごとに独立の元species incidence matrixを作成。2×2 checkerboardで (i)各本物の35station×guildの元種数、(ii)各元taxonの駅出現頻度を保つ。
- WORLD/whole studyに相当する「元35駅全体source frequency」ゲート：元路線内のtaxon分配を固定せず検証する。
- 元5rail stretches内source frequencyも保存する感度：これによって研究目的の「路線固有種プール」自体を保持するので**過剰調整の可能性**を明記する。どちらのoriginal sourceランダム化も現場の接続状態・実列車ODの実験的置換ではない。
- 2独立seed連鎖×79draw=158技術MCMC。独立した植物標本数、独立field p値、確証95%CIではない。

| 元群とsource frequency null | 実元同線−異線差 | null平均 | 観測−null | 実測がnullより大きい元線区 |
|---|---:|---:|---:|---:|
| Wind proxy, global35駅種頻度固定 | −0.019782 | −0.002945 | **−0.016837** | 1/5 |
| Other proxy, global35駅種頻度固定 | −0.004017 | −0.004220 | +0.000203 | 2/5 |
| **Wind−other interaction、global frequency固定** | −0.015765 | +0.001276 | **−0.017040** | 2/5 |
| Wind−other interaction、元5区間内taxon頻度も固定 | −0.015765 | +0.008070 | **−0.023835** | 1/5 |

元global frequency nullでwind−otherのTechnical upper rank 0.805（強い上側異常なし）。元区間頻度も固定では0.950（参照種pool自体の影響を過調整）。**どちらでも風戦略がその他より鉄道路線内で特に似るという証拠は見つかっていない**。ただし散布「能力」の順位としてPladias Sádlo戦略分類をそのまま読むことはできない。フランスPenone 2012の中程度のmobilityという機能クラスと一致した反証とはならない。

## 5．外部文献の照合・次のデータ取得閾値

フランスCaterina Penone et al (2012), Biological Conservation DOI https://doi.org/10.1016/j.biocon.2012.01.041：真の2路線・71地点（5×1m² quadrats/地点）を調べ、駅では連結が保持されやすく、overpassesでは中程度のmobilityを含む群で断絶効果が現れると論じている。元site×taxonの生行列・GPS・break位置の公的原本を再検証できていない。既存tetsu v14.5/v14.6ではMNHN博士論文Appendix Bの6種に関する観測数の元出典を監査し、原71×種行列を復元できないこと・元観測数が71真site別出現頻度を意味するか不明なことを証明済み。これを独立の再現した成功サンプルとして加算しない。

スロバキアMájeková et al (2025), Flora DOI https://doi.org/10.1016/j.flora.2025.152780：元30駅（1 ha/駅）の現地調査、元413種（native258/alien155）、貨物列車本数と外来種alpha richnessの正関連は**先行原著の肯定結果として尊重する**。しかし駅間beta類似度と実貨物本数の同じ元taxon×site行列＋駅別実数値の取得は未認証。原論文謝辞では貨物本数データがスロバキア鉄道会社から提供されたと記されている。したがって本人の著作外の実個別列車本数を推測で作らない。

今後の真の差別化は「列車がどの植物群に、どの散布段階で連結性を与えるか」を真正データで比較すること。散布戦略は到着確率の直接測定ではなく、局所pool・元習性・連絡線のservice, cargo OD,現地種子・発芽を個別に観測する必要がある。先行フランス研究が使用した中程度のdispersal abilityに対応する元測定軸を同じ原資料から取得できなければ、単純な風proxyとの効果量統合はしない。チェコ駅全域（3ゾーン合算、元36?ha面積不均等）とスロベニア100m道床BS/BTだけから「全球列車風圧選択」は主張できない。

## 6．再現コードと利用制約

- scripts/analyze_Czech_35stations_same_stretch_pladias_dispersal_guild_v284.py
- scripts/audit_Pladias_original_semantic_trait_fingerprint_v284.py
- tests/test_Czech_original_Pladias_dispersal_trait_matched_v284.py
- .github/workflows/tetsu-v284-Czech-35station-dispersal-lifeform-network.yml

本研究の元Pladiasのデータ利用・相当量の二次研究では、原データベースの利用規約および研究利用に関するconsortium boardの事前許諾・著者関与のガイドラインを尊重する。正式投稿前に必ず確認する。元出版社Preslia元植物taxa出現とSR70元駅代表点を新たなPRへ再配布しない。Draft PRのみ、mainは未変更。
