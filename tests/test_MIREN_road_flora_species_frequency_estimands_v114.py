"""v11.4 original MIREN region species site-rate and effort safeguards."""
import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from analyze_MIREN_road_flora_species_frequency_estimands_v114 import (
    species_preference_from_site_rows,region_contrast,macro_summary,run
)


def record(site,taxon,status,plot,region="X"):
    return {"site":site,"taxon":taxon,"status":status,"plot":plot,"region":region}


class OriginalGlobalRoadFloraControlTests(unittest.TestCase):
    def test_original_pair_observation_not_false_yearly_absence(self):
        rows=[record("S1","a","Alien","1"),
              record("S2","a","Alien","1"),record("S2","a","Alien","3"),
              record("S2","b","Native","3")]
        by=species_preference_from_site_rows(rows)
        a=next(x for x in by["X"] if x["taxon"]=="a")
        b=next(x for x in by["X"] if x["taxon"]=="b")
        self.assertEqual(a["site_occupancy"],2)
        self.assertAlmostEqual(a["roadside_minus_interior_preference"],.5)
        self.assertAlmostEqual(b["roadside_minus_interior_preference"],-1.)

    def test_region_status_is_not_global_species_status(self):
        source=[record("s1","taxonX","Alien","1","A"),
                record("s2","taxonX","Native","3","B")]
        r=species_preference_from_site_rows(source)
        self.assertEqual(r["A"][0]["guild"],"Alien")
        self.assertEqual(r["B"][0]["guild"],"Native")
        with self.assertRaisesRegex(ValueError,"conflicted"):
            species_preference_from_site_rows([
                record("s1","X","Alien","1"),record("s2","X","Native","3")])

    def test_species_equal_vs_occupancy_can_flip_in_same_guilds(self):
        source=[]
        for guild in ("Alien","Native"):
            for i in range(4):
                value=(1 if i<3 else -1) if guild=="Alien" else 0
                source.append({"taxon":guild+str(i),"guild":guild,
                               "site_occupancy":1 if i<3 else 50,
                               "roadside_minus_interior_preference":value})
        res=region_contrast(source)
        self.assertTrue(res["identified"])
        self.assertGreater(res["alien_minus_native_species_equal_preference"],0)
        self.assertLess(res["alien_minus_native_occurrence_weighted_preference"],0)

    def test_regions_source_less_than_10_sites_not_macro_country(self):
        row=lambda eq,oc,n:{"identified":True,
             "alien_minus_native_species_equal_preference":eq,
             "alien_minus_native_occurrence_weighted_preference":oc,
             "n_original_matched_source_road_sites":n}
        report=macro_summary([row(1,-1,11),row(-2,-1,3),row(.5,1,20)])
        self.assertEqual(report["n_original_MIREN_road_region_codes_eligible_10plus_sites"],2)
        self.assertEqual(report["n_region_codes_with_local_weighting_sign_reversal"],1)

    def test_cannot_claim_true_MIREN_without_v99_source_receipt(self):
        with self.assertRaisesRegex(ValueError,"source receipt missing"):
            run([],{"complete":False})


if __name__=="__main__":unittest.main()
