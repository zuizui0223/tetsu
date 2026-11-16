"""v11.0 source-invariant estimand and conditional label-permutation tests."""
import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from analyze_Slovenia_rail_species_occupancy_weighting_v110 import (
    group_estimator,frequency_bin,rank_logfreq_ols,
    permute_source_traits_within_lifespan_x_occurrence_band,run
)


def make_taxon(code,annual,wind,n,score):
    return {
        "taxon_code":code,
        "stratum":"strict_annual" if annual else "other_known_lifespan",
        "wind":wind,
        "n_original_locations_present_either_zone":n,
        "occupancy_preference_BS_minus_BT":score,
    }


def demo():
    a=[]
    for annual in (True,False):
        for i in range(4):
            a.append(make_taxon(str(annual)+"W"+str(i),
                                annual,True,1 if i<3 else 20,
                                1 if i<3 else -1))
            a.append(make_taxon(str(annual)+"N"+str(i),
                                annual,False,1 if i<3 else 20,0))
    return a


class SourceSpeciesWeightingTests(unittest.TestCase):
    def test_exact_frequency_bin_boundaries(self):
        self.assertEqual([frequency_bin(i) for i in (1,2,3,5,6,9,10,50)],
                         ["1-2","1-2","3-5","3-5","6-9","6-9","10+","10+"])
        with self.assertRaises(ValueError):frequency_bin(0)

    def test_species_equal_vs_occupancy_can_reverse_on_same_source_taxa(self):
        rows=demo()
        equal=group_estimator(rows,1,"species")
        weight=group_estimator(rows,1,"opportunities")
        self.assertTrue(equal["identified"])
        self.assertTrue(weight["identified"])
        self.assertEqual(equal["n_source_taxa"],weight["n_source_taxa"])
        self.assertGreater(equal["wind_minus_other_preference_lifespan_adjusted"],0)
        self.assertLess(weight["wind_minus_other_preference_lifespan_adjusted"],0)
        self.assertAlmostEqual(equal["wind_minus_other_preference_lifespan_adjusted"],.5)

    def test_model_fits_rare_common_trait_interaction(self):
        rows=demo()
        mod=rank_logfreq_ols(rows)
        self.assertEqual(mod["n_original_source_plant_taxa"],16)
        self.assertLess(mod["additional_logfreq_slope_for_wind"],0)

    def test_frequency_stratified_permutation_is_seed_deterministic(self):
        rows=demo()
        p1=permute_source_traits_within_lifespan_x_occurrence_band(rows,draws=40,seed=4)
        p2=permute_source_traits_within_lifespan_x_occurrence_band(list(reversed(rows)),draws=40,seed=4)
        self.assertEqual(p1["restricted_permutation_two_sided_centered_p"],
                         p2["restricted_permutation_two_sided_centered_p"])
        self.assertEqual(p1["n_source_lifespan_x_frequency_band_cells"],4)

    def test_wrong_weight_mode_rejected(self):
        with self.assertRaises(ValueError):
            group_estimator(demo(),1,"invented")
        self.assertFalse(group_estimator(demo(),30,"species")["identified"])

    def test_original_field_source_receipt_is_required(self):
        with self.assertRaisesRegex(ValueError,"v10.9 verified receipt required"):
            run({},{"complete":False})


if __name__=="__main__":unittest.main()
