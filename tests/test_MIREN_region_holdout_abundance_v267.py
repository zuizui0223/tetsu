"""v26.7 synthetic QA for true REGION holdout and unseen taxon transfer."""
from __future__ import annotations
import math, sys, unittest
import numpy as np
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import analyze_MIREN_source_region_holdout_abundance_v267 as m


def synthetic_rows():
    output=[]
    for r in range(4):
        for road in range(2):
            for site in range(4):
                for plot in m.PLOTS:
                    taxa=["shared0","shared1","shared2","shared3",
                          "rare"+str(r),"rare"+str(r)+"B"]
                    for t,taxon in enumerate(taxa):
                        # No source REGION effect. Only within-plot rank
                        # explains persistence; species can be unseen.
                        position=(t+site+road+(plot=="3")+r)%len(taxa)
                        rank=2*position/(len(taxa)-1)-1
                        output.append({
                            "region":"REG"+str(r),
                            "road":"road"+str(road),"site":"site"+str(site),
                            "plot":plot,"initial_year":2017,"year_gap":5,
                            "taxon":taxon,"xrank":rank,
                            "re_recorded":int(position>=3)
                        })
    return output


class OriginalRoadRegionHoldoutTests(unittest.TestCase):
    def test_no_region_feature_and_known_source_dimensions(self):
        rows=synthetic_rows()
        models,y,regions,metadata,keys=m.features_without_region(rows)
        self.assertEqual(len(regions),4)
        self.assertTrue(metadata["no_region_or_road_one_hot_features"])
        self.assertEqual(metadata["source_distinct_species"],12)
        self.assertEqual(metadata["source_species_by_plot_categories"],24)
        self.assertEqual(set(models),set(m.LABELS[:4]))
        self.assertEqual(len(y),len(rows))
        for reg,indices in regions.items():
            rest=set(range(len(rows)))-set(indices)
            self.assertTrue(all(rows[i]["region"]!=reg for i in rest))

    def test_loss_and_brier_binary(self):
        p=np.array([.1,.8,.5])
        y=np.array([0,1,1])
        ll,br=m.score_one_model(y,p)
        self.assertAlmostEqual(ll[0],-math.log(.9))
        self.assertAlmostEqual(ll[1],-math.log(.8))
        self.assertAlmostEqual(br[2],.25)
        with self.assertRaises(ValueError):
            m.score_one_model(y,np.array([0,.4,1.]))

    def test_full_region_transfer_and_unseen_taxa_falsification(self):
        rows=synthetic_rows()
        results=m.assess(rows)
        self.assertEqual(results["source_region_code_heldout_count"],4)
        self.assertTrue(results["source_heldout_folds_all_roads_and_sites_in_one_region_excluded"])
        for region,value in results["source_region_specific_results"].items():
            self.assertEqual(value["original_roads"],2)
            self.assertEqual(value["original_spatial_sites"],8)
            self.assertEqual(value["source_baseline_taxon_plot_events"],8*2*6)
            self.assertGreater(value["heldout_taxon_species_unseen_elsewhere_events"],0)
            self.assertGreater(value["heldout_taxon_species_seen_on_other_regions_events"],0)
            self.assertEqual(value["by_training_support"]["all"]["1"]["source_roads"],2)
            self.assertTrue(value["by_training_support"]["unseen_taxon"]["1"]["identified"])
        all_res=results["equal_region_predictive_gain"]["all"]["logloss"]
        self.assertEqual(all_res["M2_to_M3"]["source_region_codes_equal_weight_n"],4)
        self.assertGreater(all_res["M2_to_M3"]["mean_region_equal_roadside_gain"],.01)
        self.assertGreater(all_res["M2Z_to_M3Z"]["mean_region_equal_interior_gain"],.01)

    def test_missing_regions_fails_closed(self):
        rows=[r for r in synthetic_rows() if r["region"] in ("REG0","REG1")]
        with self.assertRaises(ValueError):m.assess(rows)

    def test_equal_site_road_weights_and_exclusion_of_unseen_taxa(self):
        # A direct deterministic nested scoring test, not MIREN observations.
        rows=[
            {"region":"R","road":"a","site":"s1","plot":"1","re_recorded":1},
            {"region":"R","road":"a","site":"s1","plot":"1","re_recorded":0},
            {"region":"R","road":"a","site":"s2","plot":"1","re_recorded":1},
            {"region":"R","road":"b","site":"s3","plot":"1","re_recorded":0},
        ]
        test=list(range(4))
        losses={k:np.array([0.,0.,1.,2.]) for k in m.LABELS}
        briers={k:np.array([0.,0.,.1,.2]) for k in m.LABELS}
        out=m.averages_from_rows(rows,test,losses,briers,
                                  [True,True,True,True],"R","1")
        # Road a site means: 0 and 1 => 0.5; road b =2; road-equal=1.25.
        self.assertAlmostEqual(out["models"]["logloss"][m.LABELS[0]],1.25)
        self.assertEqual(out["source_roads"],2)
        self.assertEqual(out["source_sites"],3)
        self.assertEqual(out["source_taxon_site_events"],4)


if __name__=="__main__":
    unittest.main()
