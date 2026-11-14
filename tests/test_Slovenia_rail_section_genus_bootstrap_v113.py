"""v11.3 original specimen source units, positive section and genus weighting checks."""
import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from analyze_Slovenia_rail_section_genus_bootstrap_v113 import (
    weighted_taxon_records,taxon_name_map,guild_contrast,
    positive_dirichlet_weights,draws_summary,run,
)


def make_traits():
    d={}
    for annual in (True,False):
        for wind in (True,False):
            for i in range(4):
                k=f"T{int(annual)}{int(wind)}{i}"
                d[k]={"Diaspore.strategy":"Epilobium" if wind else "Allium",
                      "Lifespan.categories":"annual" if annual else "perennial"}
    return d


class StrictOriginalSectionGenusTests(unittest.TestCase):
    def test_positive_source_rail_section_mass_never_loses_rare_species(self):
        d=make_traits()
        codes=sorted(d)
        sites={
            "a":{"BS":set(codes[:7]),"BT":set(codes[4:9]),"section_number":"01","n_matched_sides":1},
            "b":{"BS":set(codes[9:]),"BT":set(codes[7:]),"section_number":"02","n_matched_sides":2},
        }
        rec=weighted_taxon_records(sites,d,{"01":.001,"02":10.0})
        self.assertEqual(len(rec),len(d))
        self.assertTrue(all(r["n_original_locations_present_either_zone"]>0 for r in rec))
        with self.assertRaises(ValueError):
            weighted_taxon_records(sites,d,{"01":0.,"02":1.})

    def test_base_species_guild_contrast_exact_equal_and_opportunities_math(self):
        d=make_traits()
        original=[]
        for code,tr in sorted(d.items()):
            w=tr["Diaspore.strategy"]=="Epilobium"
            a=tr["Lifespan.categories"]=="annual"
            original.append({
                "taxon_code":code,"wind":w,
                "stratum":"strict_annual" if a else "other_known_lifespan",
                "occupancy_preference_BS_minus_BT":(.75 if w else -.25),
                "n_original_locations_present_either_zone":1 if w else 10
            })
        r=guild_contrast(original)
        self.assertAlmostEqual(r["species"],1.)
        self.assertAlmostEqual(r["opportunities"],1.)
        self.assertFalse(r["sign_reversal"])

    def test_report_same_taxa_opposite_signs_from_frequency_only(self):
        d=make_traits(); original=[]
        for code,tr in sorted(d.items()):
            w=tr["Diaspore.strategy"]=="Epilobium"
            a=tr["Lifespan.categories"]=="annual"
            i=int(code[-1])
            original.append({
                "taxon_code":code,"wind":w,
                "stratum":"strict_annual" if a else "other_known_lifespan",
                "occupancy_preference_BS_minus_BT":1. if w and i<3 else (-1. if w else 0.),
                "n_original_locations_present_either_zone":1 if i<3 else 50,
            })
        r=guild_contrast(original)
        self.assertTrue(r["sign_reversal"])
        self.assertGreater(r["species"],0)
        self.assertLess(r["opportunities"],0)
        diag=draws_summary([r]*7)
        self.assertEqual(diag["fraction_draws_retain_sign_reversal"],1.)

    def test_genus_must_come_from_published_latin_binomial_not_Eppo_chars(self):
        original={
            "traits":[{"EPPO_code":"T1","Species":"Plantus alpha","Lifespan.categories":"annual"},
                      {"EPPO_code":"T2","Species":"Planta beta"}],
            "species_codebook":[{"EPPO_code":"T1","Latin_name":"Plantus alpha"},
                                {"EPPO_code":"T2","Latin_name":"Planta beta"}],
        }
        m,q=taxon_name_map(original,{"T1","T2"})
        self.assertEqual(m["T1"],"Plantus")
        self.assertEqual(m["T2"],"Planta")
        self.assertEqual(q["n_trait_coded_taxa_with_genus"],2)
        self.assertNotIn("T1",q["missing_taxon_codes"])

    def test_contradictory_original_source_genus_fails_closed(self):
        original={
            "traits":[{"EPPO_code":"T1","Scientific name":"Botanus alpha"}],
            "species_codebook":[{"EPPO_code":"T1","Latin_name":"Plantus alpha"}],
        }
        m,q=taxon_name_map(original,{"T1"})
        self.assertNotIn("T1",m)
        self.assertEqual(q["n_source_trait_codes_disagree_genus"],1)

    def test_positive_cluster_weights_sum_and_reproduce_fixed_seed(self):
        import random
        a=positive_dirichlet_weights(["a","b","c"],random.Random(23))
        b=positive_dirichlet_weights(["a","b","c"],random.Random(23))
        self.assertEqual(a,b)
        self.assertAlmostEqual(sum(a.values()),3.)
        self.assertTrue(all(x>0 for x in a.values()))

    def test_real_original_v11_receipt_required(self):
        with self.assertRaisesRegex(ValueError,"genuine source verified"):
            run({},{"complete":False})


if __name__=="__main__":unittest.main()
