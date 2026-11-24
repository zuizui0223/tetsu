#!/usr/bin/env python3
"""Synthetic source integrity and warm-affinity arithmetic, no field effect."""
import csv,math,sys,tempfile,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from analyze_global_plant_temperature_railway_v347 import (
    load_global,simple_binomial,summarize,pairs_kashmir,fieldmap,
    GLOBAL_K_MAX,GLOBAL_K_MIN,C4_CANDIDATES
)

class V347GlobalThermalTests(unittest.TestCase):
    def setUp(self):
        self.names=["genus "+letters(i) for i in range(24)]
        self.traits={name:(270.0+i,100) for i,name in enumerate(self.names)}
        self.source={name:1 for name in self.names}

    def test_K_to_C_and_weighting(self):
        res=summarize(self.source,self.traits,20,equal=True)
        self.assertTrue(res["estimable"])
        self.assertAlmostEqual(res["climate_tmean_C"],281.5)
        self.assertEqual(res["source_global_temperature_covered_species_n"],24)

    def test_original_observation_support_not_filled(self):
        x=dict(self.traits)
        for name in self.names[:12]:x[name]=(240.0,10)
        a=summarize(self.source,x,20)
        self.assertFalse(a["estimable"]) # 12 source species below frozen 15 threshold
        # changing original n cutoff never invents absent trait values

    def test_species_names_nonambiguous(self):
        self.assertEqual(simple_binomial(" Cynodon   dactylon "),"cynodon dactylon")
        self.assertIsNone(simple_binomial("Poa annua subsp. annua"))
        self.assertIsNone(simple_binomial("Poa cf. annua"))
        self.assertIn("cynodon dactylon",C4_CANDIDATES)

    def test_original_source_columns_required(self):
        fieldmap(["species","numrecords","tmean","genus","family"])
        with self.assertRaises(ValueError):fieldmap(["species","tmean"])

    def test_paired_stations_use_explicit_original_site_count(self):
        a=set(self.names[:12]);b=set(self.names[12:])
        pairs=[(a,b) for _ in range(16)]
        out=pairs_kashmir(pairs,self.traits,20)
        self.assertTrue(out["estimable"])
        self.assertEqual(out["original_paired_source_stations_eligible"],16)
        self.assertAlmostEqual(out["same_station_equal_species_temperature_difference_C"],-12)

def letters(i):
    a=""
    i+=1
    while i:
        i,z=divmod(i-1,26)
        a=chr(97+z)+a
    return a

if __name__=="__main__":unittest.main()
