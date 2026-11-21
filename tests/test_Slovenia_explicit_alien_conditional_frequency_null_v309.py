"""v309 synthetic numerical and source-frequency-null conservation tests."""
from __future__ import annotations
import unittest,sys
from pathlib import Path
from collections import Counter
import random
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from analyze_Slovenia_explicit_alien_traffic_conditional_frequency_null_v309 import (
    GROUPS, make_guild_sets,dyad_contribution_vectors,exposure_design,fit_slope,swap_communities
)
class TestAlienSharedSpeciesFrequencyNull(unittest.TestCase):
    def setUp(self):
        self.sites=[
            {"A","B"},{"A","C"},{"A","B","C"},{"B","C"}
        ]
        self.traits={
            "A":{"Nativeness.status":None, "Lifespan.categories":"annual","Diaspore.strategy":"Bidens"},
            "B":{"Nativeness.status":"invasive alien","Lifespan.categories":"perennial","Diaspore.strategy":"Epilobium"},
            "C":{"Nativeness.status":"native","Lifespan.categories":"annual/biennial","Diaspore.strategy":None},
        }
        self.dyads=[
            {"left":0,"right":1,"freight":1,"log_distance":0.7,"members":[(0,1),(0,2)]},
            {"left":0,"right":2,"freight":2,"log_distance":1.4,"members":[(0,3)]},
            {"left":1,"right":2,"freight":3,"log_distance":0.8,"members":[(1,3),(2,3)]}
        ]
    def test_explicit_only_and_unknown_never_assigned_native(self):
        g=make_guild_sets(self.sites,self.traits)
        self.assertEqual(g["explicit_alien"],{"B"})
        self.assertEqual(g["missing_nativeness"],{"A"})
        self.assertEqual(g["not_explicit_alien"],{"A","C"})
        self.assertEqual(g["strict_annual"],{"A"})
        self.assertEqual(g["wind_Epilobium_Phragmites"],{"B"})
        self.assertEqual(g["Bidens_archetype"],{"A"})
    def test_source_full_jaccard_species_contributions_and_traffic_slope_sum(self):
        guild=make_guild_sets(self.sites,self.traits)
        C=dyad_contribution_vectors(self.sites,self.dyads,guild)
        observed=np.array([
            np.mean([len(self.sites[a]&self.sites[b])/len(self.sites[a]|self.sites[b]) for a,b in d['members']])
            for d in self.dyads
        ])
        np.testing.assert_allclose(C[:,0],observed,rtol=0,atol=1e-12)
        np.testing.assert_allclose(C[:,1]+C[:,2],observed,rtol=0,atol=1e-12)
        w=np.linalg.pinv(exposure_design(self.dyads))[1]
        self.assertAlmostEqual(float(w@C[:,0]),fit_slope(self.dyads,observed),places=10)
        self.assertAlmostEqual(float(w@C[:,1]+w@C[:,2]),float(w@C[:,0]),places=10)
    def test_stratified_occupancy_and_rows_are_fixed(self):
        observed=[{"A","C"},{"B"},{"A","B"},{"C"},{"A"},{"C"}]
        orig_size=list(map(len,observed))
        orig_total=Counter(t for r in observed for t in r)
        strata=[[0,1,2],[3,4,5]]
        by=[Counter(t for i in group for t in observed[i]) for group in strata]
        work=list(map(set,observed))
        accepted=swap_communities(work,random.Random(309),700,strata)
        self.assertGreater(accepted,0)
        self.assertEqual(list(map(len,work)),orig_size)
        self.assertEqual(Counter(t for r in work for t in r),orig_total)
        self.assertEqual([Counter(t for i in group for t in work[i]) for group in strata],by)
    def test_missing_community_cannot_be_silently_imputed(self):
        bad=[set(),set(),set(),set()]
        with self.assertRaises(ValueError):
            dyad_contribution_vectors(bad,self.dyads,make_guild_sets(self.sites,self.traits))
    def test_cannot_treat_all_taxon_and_alien_only_jaccard_as_same(self):
        guild=make_guild_sets(self.sites,self.traits)
        C=dyad_contribution_vectors(self.sites,self.dyads,guild)
        # If someone recomputed Jaccard after retaining only alien B,
        # they would change the denominator. That is NOT source decomposition.
        alien_only=[s&{"B"} for s in self.sites]
        self.assertNotEqual(float(C[0,1]),np.mean([
            len(alien_only[a]&alien_only[b])/len(alien_only[a]|alien_only[b])
            for a,b in self.dyads[0]["members"] if alien_only[a]|alien_only[b]]))
if __name__=="__main__":unittest.main()
