"""v10.6 deterministic algebra and original trait missingness method checks."""
import sys
import unittest
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from analyze_Slovenia_rail_v106_annual_wind_decomposition import (
    joint_trait_counts, wind_decomposition, composition_sensitivity
)


def groups(annual_n, annual_w, other_n, other_w):
    return {"bins":{
        "strict_annual":{"n":annual_n,"wind":annual_w},
        "other_known_lifespan_including_mixed":{"n":other_n,"wind":other_w}},
        "n_original_taxa_with_both_source_traits":annual_n+other_n}


class SourceAnnualVsWindDecompositionTests(unittest.TestCase):
    def test_exact_kitagawa_identity(self):
        bs=groups(4,2,6,1)
        bt=groups(2,1,8,2)
        x=wind_decomposition(bs,bt)
        self.assertAlmostEqual(x["joint_known_wind_BS_minus_BT"],
                               x["annual_mix_component"]+x["within_lifespan_wind_component"])
        self.assertAlmostEqual(x["joint_known_wind_BS_minus_BT"], .30-.30)
        self.assertTrue(x["original_exact_algebra_identity"])

    def test_annual_composition_alone_can_produce_wind_deficit(self):
        # Both source lifespan guilds have the SAME conditional wind
        # propensity at BS and BT: 0 for annual, 1 for other.
        bs=groups(8,0,2,2)
        bt=groups(2,0,8,8)
        x=wind_decomposition(bs,bt)
        self.assertAlmostEqual(x["joint_known_wind_BS_minus_BT"], -.6)
        self.assertAlmostEqual(x["annual_mix_component"], -.6)
        self.assertAlmostEqual(x["within_lifespan_wind_component"], 0)

    def test_within_guild_change_not_confused_with_lifespan_mix(self):
        bs=groups(5,4,5,3)
        bt=groups(5,2,5,1)
        x=wind_decomposition(bs,bt)
        self.assertAlmostEqual(x["annual_mix_component"], 0)
        self.assertAlmostEqual(x["within_lifespan_wind_component"], .4)
        self.assertAlmostEqual(x["joint_known_wind_BS_minus_BT"], .4)

    def test_missing_lifespan_or_diaspore_never_imputed(self):
        tr={
            "A":{"Lifespan.categories":"annual","Diaspore.strategy":"Epilobium"},
            "B":{"Lifespan.categories":"annual","Diaspore.strategy":"Cornus"},
            "C":{"Lifespan.categories":"perennial","Diaspore.strategy":"Phragmites"},
            "D":{"Lifespan.categories":"annual/biennial","Diaspore.strategy":"Allium"},
            "E":{"Lifespan.categories":None,"Diaspore.strategy":"Epilobium"},
            "F":{"Lifespan.categories":"annual","Diaspore.strategy":None},
        }
        x=joint_trait_counts(set(tr),tr)
        self.assertEqual(x["n_original_taxa_with_both_source_traits"],4)
        self.assertEqual(x["n_source_missing_strategy"],1)
        self.assertEqual(x["n_source_missing_lifespan"],1)
        self.assertEqual(x["bins"]["strict_annual"],{"n":2,"wind":1})
        self.assertEqual(x["bins"]["other_known_lifespan_including_mixed"],
                         {"n":2,"wind":1})

    def test_unidentifiable_sparse_stratum_not_assumed_zero(self):
        self.assertIsNone(wind_decomposition(groups(0,0,5,1),groups(3,1,4,2)))
        self.assertIsNone(wind_decomposition(groups(3,1,4,2),groups(5,0,0,0)))

    def test_locked_original_v105_success_receipt_required(self):
        with self.assertRaisesRegex(ValueError,"receipt missing"):
            composition_sensitivity({}, {}, {}, {"complete":False})


if __name__=="__main__":
    unittest.main()
