"""v11.5 exact side/site aggregation, exponent interpolation and section swaps."""
import unittest,sys,math
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from analyze_Slovenia_rail_side_aggregation_section_flip_v115 import (
    preference_site_matrix, continuum, crossing, exact_section_label_swaps,
    section_contributions, run
)


def species(code,wind,annual,values):
    return {
        "taxon_code":code, "wind":wind,
        "stratum":"strict_annual" if annual else "other_known_lifespan",
        "site_delta":values,
        "n_original_locations_present_either_zone":len(values),
        "occupancy_preference_BS_minus_BT":sum(values.values())/len(values),
    }


class TrueSourceSidePermutationTests(unittest.TestCase):
    def test_site_side_average_differs_from_LR_union_when_sources_disagree(self):
        sides={
            "A":[
                {"section":"1","side":"L","bs":{"W","N"},"bt":{"N"}},
                {"section":"1","side":"R","bs":{"N"},"bt":{"W","N"}}
            ],
            "B":[
                {"section":"2","side":"L","bs":{"W"},"bt":set()}
            ],
        }
        tr={
            "W":{"Diaspore.strategy":"Epilobium","Lifespan.categories":"annual"},
            "N":{"Diaspore.strategy":"Allium","Lifespan.categories":"annual"},
        }
        u,_=preference_site_matrix(sides,tr,mode="union")
        s,_=preference_site_matrix(sides,tr,mode="side_average")
        w1=next(x for x in u if x["taxon_code"]=="W")
        w2=next(x for x in s if x["taxon_code"]=="W")
        self.assertAlmostEqual(w1["occupancy_preference_BS_minus_BT"],.5)
        self.assertAlmostEqual(w2["occupancy_preference_BS_minus_BT"],.5)
        self.assertEqual(w1["site_delta"]["A"],0.)
        self.assertEqual(w2["site_delta"]["A"],0.)
        # A different within-side pattern can change preference:
        sides["A"][1]["bs"].add("W")
        s2,_=preference_site_matrix(sides,tr,mode="side_average")
        u2,_=preference_site_matrix(sides,tr,mode="union")
        self.assertAlmostEqual(next(x for x in s2 if x["taxon_code"]=="W")[
            "site_delta"]["A"],.5)
        self.assertEqual(next(x for x in u2 if x["taxon_code"]=="W")[
            "site_delta"]["A"],0.)

    def test_frequency_root_derived_from_actual_same_taxa_weights(self):
        d={}
        for annual in (True,False):
            for i in range(4):
                for wind in (True,False):
                    key=str(annual)+"_"+str(wind)+"_"+str(i)
                    v={str(j):1. if wind and i<3 else
                        -1. if wind else 0. for j in range(1 if i<3 else 20)}
                    d[key]=species(key,wind,annual,v)
        rows=list(d.values())
        self.assertGreater(continuum(rows,0),0)
        self.assertLess(continuum(rows,1),0)
        z=crossing(rows)
        self.assertTrue(0<z["source_frequency_exponent_zero_crossing_0_to_1"]<1)
        self.assertLess(abs(continuum(rows,z["source_frequency_exponent_zero_crossing_0_to_1"])),1e-10)

    def test_section_linear_source_decomposition_recreates_continuum(self):
        rows=[]
        section_by_site={}
        for i in range(15):
            section_by_site[str(i)]=str(i)
        for annual in (True,False):
            for wind in (True,False):
                for t in range(4):
                    key=f"T{annual}{wind}{t}"
                    vals={str(i):float(int((i+t)%3==0)-int((i+t)%3==1))
                          for i in range(15)}
                    rows.append(species(key,wind,annual,vals))
        for power in (0,0.25,1):
            contrib=section_contributions(rows,section_by_site,power)
            self.assertEqual(len(contrib),15)
            self.assertAlmostEqual(sum(contrib.values()),continuum(rows,power))
        exact=exact_section_label_swaps(rows,section_by_site)
        self.assertEqual(exact["n_possible_section_exchangeability_flips"],32768)
        for name in ("two_sided_full_15section_flip_p_species_equal",
                     "two_sided_full_15section_flip_p_occ_weighted",
                     "two_sided_full_15section_flip_p_weighting_gap"):
            self.assertTrue(0<=exact[name]<=1)

    def test_negative_weights_and_bad_modes_not_allowed(self):
        with self.assertRaises(ValueError):
            preference_site_matrix({}, {}, mode="fake")
        with self.assertRaises(ValueError):
            continuum([], -1)

    def test_real_v110_original_source_receipt_required(self):
        with self.assertRaisesRegex(ValueError,"Locked original v11.0"):
            run({},{"complete":False})


if __name__=="__main__":unittest.main()
