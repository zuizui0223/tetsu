"""v27.0 synthetic tests: genuine species frequency and site richness fixed null."""
import random,sys,unittest
from collections import Counter
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import audit_RAIL_ROAD_source_frequency_fixed_reference_null_v270 as v


class SourceFrequencyFixedTests(unittest.TestCase):
    def test_checkerboard_exact_all_row_column_margins(self):
        rows=[{1,2},{2,3},{1,4},{3,5}]
        before=[len(s) for s in rows]
        columns=v.source_freq(rows)
        changed=v.swaps(rows,random.Random(2026),200)
        self.assertGreater(changed,0)
        self.assertEqual([len(x) for x in rows],before)
        self.assertEqual(v.source_freq(rows),columns)

    def test_original_source_corridor_unmodified_and_local_identity_excess(self):
        # The identical same-site original vegetation has Jaccard=0,
        # but independent reference reshuffling with exact site richness
        # and species occupancy produces positive nonlocal dissimilarity.
        corridor=[{1,2},{3,4},{5,6},{7,8},{9,10}]
        reference=[set(x) for x in corridor]
        frozen=[set(x) for x in corridor]
        out=v.sampler(corridor,reference,draws=39,burnin=300,spacing=30,seed=77)
        self.assertEqual(corridor,frozen)
        self.assertEqual(out["observed"]["beta_jac"],0.)
        self.assertGreater(v.finish_null(out["observed"],out["null_draws"])["beta_jac"]["frequency_null_mean"],.1)
        self.assertTrue(out["each_original_site_reference_richness_fixed_exact"])

    def test_zero_swap_is_reported_not_fabricated(self):
        corridor=[{1},{2},{3},{4}]
        reference=[{1,2},{1,2},{1,2},{1,2}]
        out=v.sampler(corridor,reference,draws=19,burnin=20,spacing=9,seed=10)
        self.assertEqual(out["accepted_swaps_total"],0)
        for q in out["null_draws"]:
            self.assertEqual(q["beta_jac"],out["observed"]["beta_jac"])
        self.assertGreater(out["observed"]["beta_jne_nestedness"],0.)

    def test_richness_fixed_null_does_NOT_zero_out_nestedness(self):
        corridor=[{1},{2},{3},{4}]
        ref=[{1,2,5},{2,3,6},{1,4,7},{3,4,8}]
        out=v.sampler(corridor,ref,draws=19,burnin=100,spacing=10,seed=14)
        self.assertGreater(out["observed"]["beta_jne_nestedness"],0.)
        self.assertTrue(all(x["beta_jne_nestedness"]>=0 for x in out["null_draws"]))
        self.assertTrue(all(abs(x["beta_jac"]-x["beta_jtu_replacement"]-x["beta_jne_nestedness"])<1e-11 for x in out["null_draws"]))

    def test_beta_conserved_in_sample_and_invalid_draw_fails(self):
        m=v.site_mean([{1,2},{3}], [{1,2,3},{4,5}])
        self.assertAlmostEqual(m["beta_jac"],m["beta_jtu_replacement"]+m["beta_jne_nestedness"])
        with self.assertRaises(ValueError):
            v.sampler([{1}],[{1}],draws=9)
        with self.assertRaises(ValueError):
            v.site_mean([{1}],[])

    def test_null_inference_is_only_technical_rank(self):
        actual={"beta_jac":.3,"beta_jtu_replacement":.2,"beta_jne_nestedness":.1,
                "shared_taxa":2.,"corridor_only_taxa":1.,"corridor_unshared_fraction":.33}
        nul=[{**actual,"beta_jac":.4,"beta_jtu_replacement":.3,"beta_jne_nestedness":.1} for _ in range(19)]
        vout=v.finish_null(actual,nul)
        self.assertAlmostEqual(vout["beta_jac"]["observed_minus_frequency_null"],-.1)
        self.assertAlmostEqual(vout["beta_jac"]["frequency_null_mean"],.4)
        self.assertEqual(vout["beta_jac"]["upper_tail_rank_correlation_sensitive_NOT_field_p"],1.)


if __name__=="__main__":unittest.main()
