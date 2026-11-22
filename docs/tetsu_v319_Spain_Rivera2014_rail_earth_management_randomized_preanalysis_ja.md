# v319 原著数値を見る前に固定：スペインの鉄道盛土15区画無作為土壌改良試験の機構監査

2026-11-22 JST。Mell et al. (2024) review原316採用研究の元`Rivera et al 2014`行が`LTI=Waterway`・`Country=USA`、実際はSpain鉄道盛土である誤表記を原著本文から確認した後のsource再解析。これはv315の「気候による鉄道温暖種偏向」ではなく、**同一鉄道盛土内での土壌管理が植物種数・被覆の変化に与える直接介入**。予測符号は原著要旨を事前に把握しているので「新発見」「ブラインド仮説」「未公開介入」とは呼ばない。

## 真の原資料・観測単位
Rivera D et al. 2014, *PLoS ONE* 9:e101413, DOI https://doi.org/10.1371/journal.pone.0101413。原出版社はCC BY、研究全文とSupporting Table S1を公開し、原著Data Availabilityは表とsupplementaryにデータがあると記す。Europe PMC full text JATS XML `https://www.ebi.ac.uk/europepmc/webservices/rest/PMC4077791/fullTextXML` から元のDOI、元Table1/2の数値と方法を取得し、原誌サイト原文と一致することを検査する。ウェブ本文の文脈と異なる元XMLを採用しない。原sourceのSHA256、リンク、取得日時を監査。Fulltextが403/404なら強制的なスクレイピングはせずsource gateを停止。

- 真研究単位：Spain, west of Madrid、工事中の**鉄道盛土**（田園/半乾燥）、斜面14°、各6×9 m (54m²)、無作為化**5ブロック×3処置＝15野外区画**。Control=原盛土基盤、T10=上層土10cm、T30=上層土30cm。2010年春と2011年春同区画繰返し。15×2を30独立野外反復としない。
- 真比較：元原著Table2に記された2010/2011各群の**mean±SE**、植物被覆%、植物種数（richness）、群集均等度。各年 C/T10/T30 を原誌Table2の正しい処置/年コードへ割り付ける。原PDFテーブルに変形や混合単位があればブロックする。
- 主要推定量は `[(T10_2011−T10_2010)−(C_2011−C_2010)]` および同じT30の差分の差、2応答（植物被覆%ポイント、種数/区画）。これは**原群平均だけから計算できる記述的DID**、無作為介入の平均差の精密なSE/p値ではない。処置×年の相互作用F/pは元論文公表値として別に保持。元各ブロック内の共同変動と個別plot処置割当が完全再構成できなければ、新しい信頼区間、p値や階層RDAを推定しない。
- 第2の中心は**管理用土壌の追加量の飽和**。Table2の2年目T10 vs T30の植物被覆/種数の差と、不確実性を勝手に予測しないという範囲を報告。原論文の10/30cm差が有意でないという既報と照合。対照で優占の当初種や元Table S1の群別種数の変化を検証可能なら生物学的選別として分ける。
- メカニズム：土壌呼吸、β-glucosidase、phosphataseの元著者集計値（表と文中）と植物被覆の方向を区別。原著のsoil addition randomized treatment による微生物・植物の変化は土壌管理の介入に関する因果的示唆を持つが、媒介の独立実験ではない。「microbes drove cover」と書かない。列車運行本数も種子輸送も測定・操作していない。植物trait temperature / Darwinian adaptationは未評価。
- 世界の鉄道での普遍的生息地機構とも、独立したメタ効果量研究数が3になるとも扱わない。Mell Q1にもう一研究**ではなく元正しく分類すべき元スペインrandomized experiment**を追加したsource mapの1件と記す。

## source実証ゲート
- JATS XML原著DOIと主論文PMCIDの完全一致、原本のoriginal-study source縦横の3処置×2年の応答が取れる。
- 原review Supplement5真正元Excelで元`WOS_25893`国/LTIの不整合を別に記録し、元元オリジナルのSpain railwayを勝手にall meta railway poolingへ加えない。
- 元source Table2が画像ベースである場合、通常のプログラム抽出で数値を読めなければOCRや手打ちで偽の精度を作らず、原表だけを判定。
- 出力は元Table2の6個の既報処置年平均（植物被覆、richness、equitability）、元source相互作用F/pと、算術DIDの正確な再構成だけ。元JATS全文/PDF/plot個票の転載なし。
