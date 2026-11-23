"""v338 synthetic source conditional species-abundance label-shuffle tests."""
import unittest,sys,random,math
from collections import Counter,defaultdict
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import analyze_Kashmir_species_conditioned_abundance_null_v338 as m

EIVE={
  "poa annua":{"position":8.,"width":3.,"systems":3},
  "galium aparine":{"position":4.,"width":4.,"systems":3},
  "lotus corniculatus":{"position":5.,"width":2.,"systems":1},
}
def fixture():
    values={
       ("Between A",2014,"Track"):{"poa annua":7.,"galium aparine":1.},
       ("Between A",2014,"Interior Plot 2"):{"poa annua":4.,"galium aparine":8.},
       ("Station A",2014,"Track"):{"poa annua":1.,"galium aparine":10.},
       ("Station A",2014,"Interior Plot 2"):{"poa annua":2.,"galium aparine":9.},
       ("Between A",2014,"Embankment"):{"poa annua":3.,"galium aparine":5.},
       ("Station A",2014,"Embankment"):{"poa annua":9.,"galium aparine":3.},
    }
    keys=[{"BT":("Between A",2014,"Track"),
           "BI":("Between A",2014,"Interior Plot 2"),
           "ST":("Station A",2014,"Track"),
           "SI":("Station A",2014,"Interior Plot 2")}]
    return values,keys

class SourceSpeciesAbundance(unittest.TestCase):
    def test_source_trait_crosswalk_source_system_gate(self):
        raw={("A",2014,"Track"):{0:3.,1:2.,2:4.}}
        names=["poa annua","galium aparine","lotus corniculatus"]
        self.assertEqual(len(m.site_trait_abundances(raw,names,EIVE,2014,3)["A",2014,"Track"]),2)
        self.assertEqual(len(m.site_trait_abundances(raw,names,EIVE,2014,1)["A",2014,"Track"]),3)
    def test_source_species_year_group_only(self):
        world,keys=fixture()
        groups,movable,qa=m.shuffle_setup(world,"species_all_zones")
        self.assertEqual(qa["original_source_species_groups"],2)
        self.assertEqual(qa["source_abundance_records_in_nontrivially_permutable_groups"],12)
    def test_species_zone_groups_do_not_mix_habitat(self):
        world,keys=fixture()
        groups,movable,qa=m.shuffle_setup(world,"species_within_original_zone")
        self.assertEqual(qa["original_source_species_groups"],6)
        self.assertEqual(qa["source_species_in_nontrivially_permutable_groups"],2)
        self.assertEqual(qa["source_abundance_records_in_nontrivially_permutable_groups"],12)
    def test_species_conditioned_shuffle_preserves_abundance_multisets(self):
        world,keys=fixture()
        _,movable,_=m.shuffle_setup(world,"species_all_zones")
        rng=random.Random(8)
        perm=m.reshuffle_world(world,movable,rng)
        for species in ("poa annua","galium aparine"):
            a=sorted(v[species] for v in world.values())
            b=sorted(v[species] for v in perm.values())
            self.assertEqual(a,b)
        self.assertEqual(set(perm),set(world))
    def test_species_zone_conditioned_shuffle_preserves_zone_distribution(self):
        world,keys=fixture()
        _,movable,_=m.shuffle_setup(world,"species_within_original_zone")
        perm=m.reshuffle_world(world,movable,random.Random(42))
        for species in ("poa annua","galium aparine"):
            for zone in ("Track","Interior Plot 2","Embankment"):
                a=sorted(v[species] for k,v in world.items() if k[2]==zone)
                b=sorted(v[species] for k,v in perm.items() if k[2]==zone)
                self.assertEqual(a,b)
    def test_matched_source_abundance_minus_presence_identity(self):
        world,keys=fixture()
        got=m.matched_effect(world,keys,EIVE)
        self.assertAlmostEqual(
           got["abundance_four_zone_CWM"]-got["presence_four_zone_CWM"],
           got["abundance_minus_presence_increment"])
        self.assertEqual(got["matched_pairs"],1)
    def test_permutation_null_technical_centre_nonzero_permitted(self):
        world,keys=fixture()
        got=m.technical_null(world,keys,EIVE,"species_within_original_zone",99,seed=5)
        self.assertEqual(got["draws"],99)
        self.assertTrue(0<got["technical_two_sided_fraction_centered_on_null_mean_plus_one"]<=1)
        self.assertTrue(0<got["technical_lower_tail_fraction_plus_one"]<=1)
    def test_unmatched_field_plot_not_imputed(self):
        world,keys=fixture()
        world.pop(("Station A",2014,"Track"))
        with self.assertRaises(ValueError):m.matched_effect(world,keys,EIVE)
    def test_source_frozen_panels_and_modes(self):
        self.assertEqual(len(m.PANEL_DESIGN),4)
        self.assertEqual(m.MODES,("species_all_zones","species_within_original_zone"))

if __name__=="__main__":unittest.main()
