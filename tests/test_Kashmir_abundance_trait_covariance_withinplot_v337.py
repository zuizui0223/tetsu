"""v337 synthetic tests for genuine species abundance-temperature covariance null."""
import unittest,sys,random,math
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import analyze_Kashmir_abundance_trait_covariance_withinplot_v337 as m

def plot(hot,other):
    return {"poa annua":(8.,hot),"galium aparine":(4.,other)}

PAIR={
 "BT":plot(10,1),
 "BI":plot(1,10),
 "ST":plot(1,10),
 "SI":plot(10,1)
}
class SourceAbundanceNull(unittest.TestCase):
    def test_covariance_identity_preserves_community_true_source_values(self):
        q=m.cwm_from_pairs(plot(10,1),1.0)
        self.assertAlmostEqual(q["presence_CWM"],6)
        self.assertAlmostEqual(q["abundance_CWM"],84/11)
        self.assertAlmostEqual(q["abundance_alignment"],
                               q["abundance_trait_covariance_over_mean"])
        self.assertGreater(q["abundance_alignment"],0)
    def test_presence_ignores_abundance_strength(self):
        r=m.cwm_from_pairs(plot(100,1),0)
        self.assertEqual(r["presence_CWM"],r["abundance_CWM"])
        self.assertAlmostEqual(r["abundance_alignment"],0)
    def test_four_original_source_plots_additivity(self):
        res=m.means_of_pairs(m.calc_paired([PAIR],1))
        self.assertAlmostEqual(res["relative_presence"],0)
        self.assertGreater(res["relative_abundance"],0)
        self.assertAlmostEqual(res["relative_abundance_alignment"],
                               res["relative_abundance"]-res["relative_presence"])
        self.assertAlmostEqual(res["plot_alignments"]["BT"],
                               -res["plot_alignments"]["BI"])
    def test_source_abundance_shuffle_null_is_centered_near_zero(self):
        p=m.conditional_null([PAIR]*3,exponent=1.0,draws=999,seed=42)
        self.assertEqual(p["observed"]["n_matched_station_between_pairs"],3)
        self.assertLess(abs(p["technical_null"]["relative_abundance_alignment"]["null_mean"]),.2)
        self.assertGreater(p["observed"]["relative_abundance_alignment"],0)
    def test_fraction_is_technical_not_causal(self):
        r=m.conditional_null([PAIR]*3,1.0,draws=99,seed=8)
        self.assertTrue(0<r["technical_null"]["relative_abundance_alignment"]["two_sided_fraction_plus_one"]<=1)
        self.assertEqual(r["n_draws"],99)
    def test_all_community_abundance_spectra_fixed_by_shuffle(self):
        d=plot(17,2)
        before=sorted(v[1] for v in d.values())
        q=m.cwm_from_pairs(d,1.0,shuffle=True,rng=random.Random(4))
        self.assertEqual(sorted(v[1] for v in d.values()),before)
        self.assertGreater(q["effective_number_matched_species"],1)
    def test_zero_or_empty_abundance_rejected(self):
        with self.assertRaises(ValueError):m.cwm_from_pairs({},1)
        with self.assertRaises(ValueError):m.cwm_from_pairs({"A":(7,0)},1)
    def test_no_matched_field_pair_not_a_zero_result(self):
        self.assertIsNone(m.conditional_null([],1,draws=9))
    def test_fixed_panels_and_alphas(self):
        self.assertEqual(len(m.PANELS),4)
        self.assertEqual(m.SEED,337)

if __name__=="__main__":unittest.main()
