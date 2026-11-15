"""v16.4 synthetic-only selection tests, not real railway species identities."""
import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from select_two_original_rail_20_species_for_inat_v164 import (
    rank_observed_source_floras,top_candidates,SPECIES_RE
)
class OriginalFloraTop20SelectionTests(unittest.TestCase):
    def test_binomial_species_only_and_site_occupancy_not_observation_abundance(self):
        sites=[{"Poa annua","Poa annua","Poa sp.","Setaria"},
               {"Poa annua","Bromus erectus","Poa annua subsp. annua"},
               {"Bromus erectus","Poa annua"}]
        r=rank_observed_source_floras(sites)
        self.assertEqual(r,[("Poa annua",3),("Bromus erectus",2)])
        self.assertFalse(SPECIES_RE.fullmatch("Carex"))
        self.assertFalse(SPECIES_RE.fullmatch("Poa annua var. annua"))

    def test_two_original_sources_alternate_and_skip_duplicate_species_names(self):
        flora={
            "SI":[{"Poa annua","Bromus erectus","Lolium perenne","Poa trivialis"},
                  {"Poa annua","Lolium perenne","Poa trivialis"},
                  {"Poa annua","Bromus erectus"},
                  {"Poa annua"}],
            "CZ":[{"Poa annua","Taraxacum officinale","Plantago major","Bromus inermis"},
                  {"Poa annua","Taraxacum officinale"},
                  {"Poa annua","Plantago major"},
                  {"Poa annua"}],
        }
        choose,meta=top_candidates(flora,n=4)
        self.assertEqual(len(choose),4)
        self.assertEqual(len({r["scientific_name_original_exact_binomial"] for r in choose}),4)
        self.assertEqual([x["original_source_selected_from"] for x in choose],
                         ["SI","CZ","SI","CZ"])
        self.assertEqual(choose[0]["scientific_name_original_exact_binomial"],"Poa annua")
        self.assertEqual(choose[1]["scientific_name_original_exact_binomial"],"Taraxacum officinale")
        self.assertEqual(choose[1]["other_original_country_same_source_site_count"],0)
        self.assertTrue(all(x["NOT_railway_affinity_vs_nonnear_railway"] for x in choose))

    def test_not_identifiable_species_rank_if_insufficient_true_source_names(self):
        flora={"SI":[{"Poa annua"}],"CZ":[{"Poa annua"}]}
        with self.assertRaisesRegex(ValueError,"Not enough original"):
            top_candidates(flora,n=20)

if __name__=="__main__":unittest.main()
