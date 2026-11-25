#!/usr/bin/env python3
import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))

from analyze_Kashmir_C4_mechanism_v339 import (
    C4, C4_GRASS, C4_AMBIGUOUS, mean_traits, source_fraction, attenuation
)


class C4FrozenAnalysisTest(unittest.TestCase):
    def setUp(self):
        self.names=["cynodon dactylon","poa annua","euphorbia hispida"]
        self.eive={
            "cynodon dactylon":{"position":8.0,"systems":3},
            "poa annua":{"position":3.0,"systems":3},
            "euphorbia hispida":{"position":6.0,"systems":1},
        }
        self.plot={0:9.0,1:1.0,2:1.0}

    def test_strict_C4_not_infer_every_plant_C3(self):
        self.assertIn("cynodon dactylon",C4)
        self.assertNotIn("poa annua",C4)
        self.assertNotIn("euphorbia hispida",C4)
        self.assertTrue(C4_GRASS<=C4)
        self.assertFalse(C4_AMBIGUOUS & C4)

    def test_mean_trait_presence_and_source_abundance(self):
        self.assertAlmostEqual(mean_traits(self.plot,self.names,self.eive,3),5.5)
        self.assertAlmostEqual(mean_traits(self.plot,self.names,self.eive,3,"original_species_abundance"),7.5)
        self.assertAlmostEqual(mean_traits(self.plot,self.names,self.eive,3,exclude=C4),3.0)
        self.assertAlmostEqual(mean_traits(self.plot,self.names,self.eive,3,
                         "original_species_abundance",C4),3.0)

    def test_ambiguous_taxa_is_sensitivity_not_strict(self):
        self.assertAlmostEqual(mean_traits(self.plot,self.names,self.eive,1,
                         exclude=C4),4.5)
        self.assertAlmostEqual(mean_traits(self.plot,self.names,self.eive,1,
                         exclude=C4|C4_AMBIGUOUS),3.0)

    def test_no_remaining_species_is_missing_not_zero(self):
        self.assertIsNone(mean_traits({0:5.0},self.names,self.eive,3,exclude=C4))
        self.assertIsNone(mean_traits(None,self.names,self.eive,3))

    def test_fraction_source_positive_weight(self):
        self.assertAlmostEqual(source_fraction(self.plot,self.names,C4),1/3)
        self.assertAlmostEqual(source_fraction(self.plot,self.names,C4,True),9/11)

    def test_mechanism_attenuation(self):
        x=attenuation(0.2,0.05)
        self.assertAlmostEqual(x["fraction_thermal_excess_removed"],0.75)
        self.assertAlmostEqual(x["fraction_thermal_excess_remaining"],0.25)
        self.assertIsNone(attenuation(0,0.2))
        self.assertIsNone(attenuation(-1,-0.1))


if __name__=="__main__":
    unittest.main()
