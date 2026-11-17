# tetsu v22.1 — チェコ鉄道文脈104分類群：ゲノム・組み立て済み転写産物・集団ゲノム候補の公開DB在庫

更新: 2026-11-16。検証：GitHub Actions [source/NBCI SUCCESS #43869305369](https://github.com/zuizui0223/tetsu/actions/runs/43869305369)。対象ブランチ `research/tetsu-slovenia-identity-guild-v210-20261116`、`main`未変更。

## 最も重要な母集団の定義

原著Presliaチェコ39駅の実植物相761分類群とCzechVeg-Openの植物×調査区原本2510340行を参照し、v17.7で「元鉄道−道路のBROADテキスト文脈」出現率標準化差 q>0、および元39駅で5駅以上・原論文2研究以上・文脈内出現5件以上の支持がある**104分類群**を抽出。対象名簿はsource-verified原著v17.4b派生データをv17.7監査コードで符号非選択的に再構成した。

**鉄道敷地内選択を実証した適応種104種、という意味ではない。** 別著者研究を留保したv19.7では種名予測の転移は再現的でない。対応する q<0 の種も元監査資料に存在し、将来の系統比較では適切に対照に含める。

## 今回の公開データベース照会結果（記載の範囲に限る）

|問い合わせ|対象104分類群中のヒット分類群数|注記|
|---|---:|---|
|NCBI Entrez Assembly, original binomial + latest|**53**|GCA/GCF等の登録候補。核ゲノム・完全性・注釈品質は未検証|
|NCBI Nucleotide TSA master|**16**|NCBI公開の組み立て済み転写産物プロジェクト。SRA RNA-Seq生リードとは異なる|
|NCBI SRA, WGS library strategy|**98**|WGSが登録された**分類群**。SRA record countsは遺伝型を持つ独立個体数ではない|
|GPS付き集団ゲノムの**個体数**|**未確定**|BioSample/SRA/VCF・原論文で個体ID×採集地×核変異を照合していない|

参照検索の対象外: Ensembl Plants、Phytozome、1KP、非NCBIリポジトリ。NCBIヒットがゼロでも実配列が「存在しない」と判断しない。種名のシノニム、自殖集団や栽培系統などのサンプリング特性も未統制。

## 保管

- [104分類群のNCBI照会CSV](../data/derived/Czech_104_broad_rail_context_NCBI_genome_TSA_WGS_inventory_v221.csv)
- [真正ソース v17.7 の監査コード](../scripts/audit_CzechVeg_rail_context_all_candidate_SIGN_balance_v177.py)
- [NCBI照会コード](../scripts/audit_rail_context_104_species_NCBI_sequence_inventory_v221.py)
- [Actions生結果のCSV/JSONアーティファクト](https://github.com/zuizui0223/tetsu/actions/runs/43869305369)

## 次の厳密な開始ゲート

1. NCBI Assemblyで候補に入った53分類群の**アセンブリアクセッション／分類学的な正確さ／核ゲノムの有無／注釈レベル／遺伝子座品質**を確認する。1種1参照配列ではその種の鉄道選択を検定できない。
2. SRA WGS実験登録98分類群について、同一種の複数BioSample・サンプル所在地の座標精度・ライブラリ方法・鉄道vs非鉄道の採集事実・時期・採集者バイアスを原著／BioSampleと1対1で突合する。原データのWGS recordは集団研究の直接証拠ではない。
3. 集団検証では地理／地域祖先・標高・気候・都市道路・土壌・除草剤履歴を調整し、事前固定した遺伝子機能候補を解析する。鉄道選択を先に仮定しない。
4. 多種比較では生命史・除草剤抵抗性を別の形質軸とし、系統樹と遺伝子系統樹を区別する。**FLC**は広く等価な単一座として扱わずBrassicaceae内の同祖性を先に確認する。**EPSPS**のコピー数は倍数性や全ゲノム重複・アセンブリ品質と切り分ける。単なる種間dN/dSを鉄道適応の証明としない。

## 論文との接続上の解釈

群集の形質分布幅（v22.0：SI/CZ真正植物原著×LEDA被覆合格）と、集団ゲノムで見つけた多型の比較は、**証拠の階層が異なる**。並列の検定として提示し、複数種で一つの選択係数を推定したとみなさない。
