# v14.5 — フランス71鉄道調査地点、原著博士論文Appendix Bの出現頻度を発見

2026-11-15。新規一次調査ではない。著者原表の一部を確認した二次解析。**完全な71地点×種行列は未取得。**

## 本物の出典と識別条件

Caterina Penone (2012), *Fonctionnement de la biodiversité en ville : contribution des dépendances vertes ferroviaires*, MNHN doctoral thesis, Appendix B printed pp.175–178, [original public PDF](https://www.ittecop.fr/doc/0THESES/PENONE.pdf). Related ORIGINAL same field cohort: Penone, Machon, Julliard & Le Viol (2012), *Biological Conservation* 148:126–133, [doi:10.1016/j.biocon.2012.01.041](https://doi.org/10.1016/j.biocon.2012.01.041). Authors surveyed **71 real railway verge sites** on two Paris rail lines, at five 1 m² vegetation quadrats/site, with stations and overpasses as spatial interruptions. They collected the data; tetsu only reanalyses original published margins.

The **Appendix B table says 'Number of observations*'**, not unequivocally 'number of distinct sites'. The row-level observation grain relative to the 71 sites still needs the author's protocol confirmation. Therefore calculating 71-site occupancy and pair intersection from the observations is **conditional only**, not a fully validated biological result.

## Six original Appendix B rows verified, no invented species

| Plant | Number of observations | Thesis page | Conditional 71-site fraction |
|---|---:|---:|---:|
| *Arrhenatherum elatius* | 69 | 175 | 97.18% |
| *Rubus fruticosus* | 53 | 178 | 74.65% |
| *Senecio inaequidens* | 39 | 178 | 54.93% |
| *Papaver rhoeas* | 39 | 177 | 54.93% |
| *Cirsium arvense* | 17 | 175 | 23.94% |
| *Brassica napus* | 3 | 175 | 4.23% |

The paper summarizes **186** taxa and the thesis Appendix C table for all plants says **187**, unresolved; do not silently fill the one-taxon discrepancy.

### New exactly derivable mathematical check (IF frequencies mean distinct-site incidence)

For n=71 sites, the number of unordered pairs is C(71,2)=**2,485**. For species s present at f_s distinct sites, its number of paired co-presences is exactly C(f_s,2), regardless of HOW the species are spatially allocated. Summing the six verified counts:

- Six-taxon incidence sum = **220**, contributed mean site richness **220/71=3.098592**.
- Six-taxon paired co-presence sum = **5,345**, mean over all pairs **5,345/2,485=2.150905** shared species/pair.
- *Arrhenatherum* contributes **2,346** co-present pairs, but *Brassica napus* only **3**. The ratio **782x** is a contrast in frequency-generated pair sharing, **NOT** railway dispersal efficacy.

If the original observations truly are distinct 71-site presences, 2.150905 is a lower bound on the **full-community average shared species COUNT** (not on Jaccard). Because original 'Number of observations*' has not been independently identified to that unit, that lower-bound interpretation remains conditional.

**This cannot reconstruct** pairwise Jaccard, individual site richness, connected-vs-disconnected contrast, geographic controls, wind/mowing treatments, seed arrival, establishment, or genetic adaptation. These six named rows are a verified convenience subset, not a representative species sample. No complete France field matrix or meta effect has been claimed. Do not add thesis as an eighteenth independent field study: it is the original Penone 2012 paper cohort already in the 17-source inventory.

## Next decisive data gate

Obtain original **71-site × all species** raw incidence table, geographic coordinates, interruption category, survey season and vegetation/urbanization variables from the authors or legitimate supplementary data archive. Then apply same Jaccard, geographic-matched comparison, source-frequency-fixed null as the Slovenia real original 89 field sites. Check separately the wind/mowing exposure before hypothesizing adaptation.

[Six-row source transcription](../data/literature/France_Penone2012_thesis_six_v145.csv) · [source-audited exact arithmetic](../scripts/analyze_France_Penone2012_thesis_occupancy_v145.py) · [tests](../tests/test_France_Penone2012_thesis_occupancy_v145.py) · [Actions](../.github/workflows/France-Penone2012-thesis-source-frequency-v145.yml). Original thesis PDF NOT redistributed. Keep Draft PR unmerged.
