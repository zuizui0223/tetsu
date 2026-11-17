# tetsu v25.3 — Kashmir 鉄道の2014→2017年「駅間類似性の変化」はどの内陸対照で変わるか

**2026-11-17に元source由来の集約JSONから追加算術監査。新しい原データ再取得・野外解析・独立検定ではない。** 元入力は `research/validation_20261117/kashmir_results.json`、Git blob SHA-1 `75dd78dd2b66a67543c248ff686c3a9dc016ea73` と完全照合。公開された元観測は Rashid et al. (2021), [JAE](https://doi.org/10.1111/1365-2664.13961) / [Dryad](https://doi.org/10.5061/dryad.0p2ngf20h)。元調査は同じ31空間地点（16駅・15駅間）を2014/2017年に再訪、4つの100 m²ゾーンで植物を記録した。ただし当該派生JSONの分析条件で除外された地点の **ID集合** は収録されず、2014と2017のNが同じでも必ず同一IDとは検証できない。

## 観測量と算術

各source年・station/between・ゾーン対ごと、`D = mean_between_sites[Jaccard_dissimilarity(rail zone)] - mean_between_sites[Jaccard_dissimilarity(interior reference)]`。負ならその組内の線路/盛土群集は内陸群集より地点間非類似度が小さい。両側種数をそろえる99反復派生値（`equalized`）も採用。ここで新規算出するのは各集約値の差 `D(2017)-D(2014)` および `[D_station(2017)-D_station(2014)] - [D_between(2017)-D_between(2014)]`。**これは集約統計からの記述的二重差であり、同じ現地点をペアにして求めたDID回帰・因果効果ではない。**

|条件|観測群数(年ごと)|非等化D変化|種数等化D変化|
|---|---:|---:|---:|
|駅 Track − Interior Plot 2|16|−0.04210|−0.07392|
|駅間 Track − Interior Plot 2|12|+0.09079|+0.03986|
|駅 Track − Interior Plot 1|15|+0.11488|+0.00288|
|駅間 Track − Interior Plot 1|15|+0.02413|−0.03683|
|駅 Embankment − Interior Plot 2|16|+0.06156|−0.13138|
|駅間 Embankment − Interior Plot 2|12|+0.11045|+0.03220|

## **駅 vs 駅間の変化の差（符号反転）**

|測定ゾーンと内陸参照|二重差（非等化）|二重差（種数等化）|
|---|---:|---:|
|Track / **Interior Plot 2**|**−0.13289**|**−0.11379**|
|Track / **Interior Plot 1**|**+0.09075**|**+0.03971**|
|Embankment / Interior Plot 2|−0.04889|−0.16358|

駅で線路側の種共有が駅間より強まった、という時間変化の方向は **遠い内陸区2を基準にしたときのみ**。内陸区1では二重差が正になり、種数等化後も逆符号が残る。同じ100 m²の現地設計であっても、内陸参照の選択が結論を左右する。ゾーン1/2の観測対象サイト集合は一致しない場合があるので、この符号差の原因を距離だけに確定しない。

## この結果から言えないこと

- 2014と2017の同一IDに限定した **paired site-level temporal change**、その統計的信頼区間、交通量の交互作用係数を推定していない。nが同数でもIDが同じとは限らない。
- 2014→2017の年次変化は、植生遷移・気候・区域管理・種数変化・土地利用・列車本数を区別しない。train frequency / passenger vs freight / main vs branch の新しい曝露データは入っていない。
- 組ごとの多くの地点ペアは同じ地点を繰り返し含む。駅と駅間はランダム割り当てではない。99希薄化drawを独立な生物標本とは数えない。
- 非等化と種数等化の両方を確認したが、二重差を確証的な列車輸送・適応・一般化した鉄道収束のp値にしない。

**次に必要な解析**：元公開Dryadのraw `Matrix_plot.xlsx` から site-ID・年・調査区・種を同時に固定して station × year × zone の**paired site-resampling（地点を丸ごとブートストラップ）** を実装する。これによって、各年に使った地点集合差と種数調整効果を切り分ける。現段階の v25.3 は、この分析で結論が参照に依存することを確定し、安易な交通量-類似性因果主張を防ぐ監査。

再現：`python -m unittest discover -s tests -p 'test_Kashmir_temporal_reference_v253.py' -v`、`python scripts/audit_Kashmir_temporal_reference_v253.py --output research/validation_20261117/kashmir_temporal_reference_audit_v253.json`。4件の合成算術/誤入力テストは開発時に実行済み。コミットした派生JSONは元Git blobを読み、同じ式を独立の計算経路で照合。