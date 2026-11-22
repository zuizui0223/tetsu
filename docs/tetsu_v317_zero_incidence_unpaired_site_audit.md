# v317 — Zero-incidence and unpaired-site audit (2026-11-22)

The local v316 implementation passed its eight synthetic tests on rerun. A focused adversarial audit identified two **previously undetected** cases:

1. When at least one valid rail/nonrail pair existed, an additional rail-only site was silently omitted. The v317 local patch now reports `unpaired_sites` and their site keys in the run summary (rather than implying every sampled site was paired).
2. Rows with `incidence=0` were skipped before species trait/source consistency and duplicate-row checks. A zero row could therefore conceal a conflicting trait or a duplicate plot/species entry. v317 validates row identity, climate and trait provenance **before** excluding zero-incidence rows from community means.

The patched local v317 code passed **3/3 new synthetic regression tests** (unpaired site accounted, zero-row trait conflict rejected, duplicate zero-row rejected). The earlier eight v316 tests passed separately. These are synthetic QA results, not evidence of a climatic effect. The v317 code and tests are conversation artifacts; **this GitHub commit records the audit only**, and the repository's v314 script remains unchanged.

## Ecological inference gate
No original matched rail/nonrail vegetation dataset was acquired or analysed in this step. Before testing whether colder climates increase the railway thermal-affinity contrast, verify actual local plot matching, climate variation across independent regions, trait coverage, survey year and effort. Existing Slovak one-line and Slovenian within-railway data alone do not identify the proposed cross-region interaction.
