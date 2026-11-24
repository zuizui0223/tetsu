"""v324 source-hierarchy synthetic tests, no real MIREN rows redistributed."""
import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import analyze_MIREN_EIVE_genus_position_width_v324 as m

TRAITS={
 "genera alpha":{"position":8,"width":2,"systems":3},
 "genera beta":{"position":7,"width":3,"systems":3},
 "genera gamma":{"position":6,"width":4,"systems":3},
 "genusb delta":{"position":5,"width":5,"systems":3},
 "genusb epsilon":{"position":4,"width":6,"systems":3},
 "genusb zeta":{"position":3,"width":7,"systems":3},
 "genusc eta":{"position":2,"width":8,"systems":3},
 "genusc theta":{"position":1,"width":9,"systems":3}
}

def fake():
    pairs=[]
    road_names=["genera alpha","genera beta","genera gamma","genusb delta","genusc eta"]
    inland=["genera gamma","genusb delta","genusb epsilon","genusb zeta","genusc theta"]
    for reg in ("AUT","CHE"):
        for rd in ("a","b"):
            for site in range(4):
                pairs.append({"region":reg,"road":rd,"site":str(site),"year":2018,
                    "roadside":set(road_names),"inland":set(inland)})
    return pairs

class GenusRoadMethods(unittest.TestCase):
    def test_full_road_region_support(self):
        pairs=fake()
        self.assertEqual(len(m.eligible_roads(pairs,TRAITS)),4)
    def test_genus_additivity_and_direct_original_road_effect(self):
        p=m.partition(m.eligible_roads(fake(),TRAITS),TRAITS,"position")
        self.assertEqual(p["regions"],2)
        self.assertEqual(p["road_ids"],4)
        self.assertEqual(p["source_first_visit_sites"],16)
        self.assertAlmostEqual(p["total_region_equal_roadside_minus_inland"],
          p["between_genus_contribution"]+p["within_genus_species_residual_contribution"])
        self.assertGreater(p["total_region_equal_roadside_minus_inland"],0)
    def test_independent_width_scale(self):
        sel=m.eligible_roads(fake(),TRAITS)
        t=m.partition(sel,TRAITS,"position")
        w=m.partition(sel,TRAITS,"width")
        self.assertLess(w["total_region_equal_roadside_minus_inland"],0)
        self.assertGreater(t["total_region_equal_roadside_minus_inland"],0)
    def test_within_genus_shuffle_uses_region_specific_species_pool(self):
        p=m.partition(m.eligible_roads(fake(),TRAITS),TRAITS,"position")
        got=m.genus_shuffle(p,draws=29,seed=2)
        self.assertEqual(got["draws"],29)
        self.assertTrue(0<got["two_sided_fraction_plus_one"]<=1)
    def test_nonreplicated_road_does_not_pass(self):
        p=fake()
        subset=[x for x in p if x["road"]=="a"]
        self.assertEqual(m.eligible_roads(subset,TRAITS),{})
if __name__=="__main__":unittest.main()
