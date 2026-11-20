"""Mathematical tests; not evidence for biological v306 outcomes."""
import sys, unittest, random
from pathlib import Path
from collections import Counter
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from analyze_Slovenia_traffic_shared_species_v306 import (
    zoneset,metric_and_contributions,fit_slope,exposure_design,
    swap_communities,contribution_analysis
)

class TestSpeciesDecomposition(unittest.TestCase):
    def setUp(self):
        self.sites=[
            {"BS":{"A","B"},"BT":{"A","C"}},
            {"BS":{"A","C"},"BT":{"B","C"}},
            {"BS":{"A","B","C"},"BT":{"A"}},
            {"BS":{"B","C"},"BT":{"B"}},
        ]
        self.dyads=[
            {"left":0,"right":1,"freight":1.,"log_distance":.7,"members":[(0,1),(0,2)]},
            {"left":0,"right":2,"freight":2.,"log_distance":1.4,"members":[(0,3)]},
            {"left":1,"right":2,"freight":3.,"log_distance":.8,"members":[(1,3),(2,3)]},
        ]

    def test_exact_jaccard_sum(self):
        arr=[zoneset(x,"BS") for x in self.sites]
        y,C=metric_and_contributions(arr,self.dyads,["A","B","C"])
        self.assertEqual(C.shape,(3,3))
        for i,d in enumerate(self.dyads):
            exp=np.mean([len(arr[a]&arr[b])/len(arr[a]|arr[b]) for a,b in d["members"]])
            self.assertAlmostEqual(y[i],exp,places=12)
        np.testing.assert_allclose(C.sum(axis=1),y,atol=1e-12)

    def test_exact_freight_coefficient_additivity(self):
        arr=[zoneset(x,"BS") for x in self.sites]
        y,C=metric_and_contributions(arr,self.dyads,["A","B","C"])
        X=exposure_design(self.dyads)
        self.assertEqual(np.linalg.matrix_rank(X),3)
        vec=np.linalg.lstsq(X,C,rcond=None)[0][1]
        self.assertAlmostEqual(float(vec.sum()),fit_slope(self.dyads,y),places=10)

    def test_unmeasured_status_not_falsely_native(self):
        src=[dict(x,source_id=str(i),section="1") for i,x in enumerate(self.sites)]
        traits={"A":{"latin_name":"Specimen A","Nativeness.status":None,
                      "Lifespan.categories":"annual","Diaspore.strategy":"Allium"},
                "B":{"latin_name":"Specimen B","Nativeness.status":"casual alien species",
                      "Lifespan.categories":"perennial","Diaspore.strategy":"Bidens"},
                "C":{"latin_name":"Specimen C","Nativeness.status":None,
                      "Lifespan.categories":"annual/biennial","Diaspore.strategy":None}}
        obs=fit_slope(self.dyads,
              metric_and_contributions([x["BS"] for x in self.sites],
                                       self.dyads,["A","B","C"])[0])
        out=contribution_analysis(src,self.dyads,traits,{},"BS",obs)
        self.assertEqual(out["source_defined_guilds"]["explicit_alien"]["n_source_species"],1)
        self.assertEqual(out["source_defined_guilds"]["strict_annual"]["n_source_species"],1)
        self.assertEqual(out["source_defined_guilds"]["source_missing_nativeness"]["n_source_species"],2)
        self.assertEqual(out["source_defined_guilds"]["author_Bidens_not_measured_rail_attachment"]["n_source_species"],1)
        self.assertAlmostEqual(out["no_urban_observed_freight_slope"],
                               sum(out["source_defined_guilds"][k]["signed_contribution_sum"]
                                   for k in ("explicit_alien","not_explicitly_alien_INCLUDES_native_unknown_crop")),places=11)

    def test_fixed_species_prevalence_and_site_richness_after_chaindraws(self):
        a=[set(x["BS"]) for x in self.sites]
        b=[set(x) for x in a]
        n=swap_communities(b,random.Random(20261120),500,[list(range(len(b)))])
        self.assertGreater(n,0)
        self.assertEqual(Counter(t for group in a for t in group),
                         Counter(t for group in b for t in group))
        self.assertEqual([len(x) for x in a],[len(x) for x in b])

    def test_zone_union_is_not_same_as_two_independent_communities(self):
        for s in self.sites:
            self.assertEqual(zoneset(s,"union"),s["BS"]|s["BT"])

if __name__=="__main__":unittest.main()
