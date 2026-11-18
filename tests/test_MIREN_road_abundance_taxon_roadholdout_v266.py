"""v26.6 mathematical and leakage tests on synthetic MIREN-like ROAD data."""
from __future__ import annotations
import sys,unittest
import numpy as np
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import analyze_MIREN_road_abundance_after_taxon_out_of_road_v266 as m
from scipy.sparse import hstack,csr_matrix


class RoadHeldoutAbundanceTests(unittest.TestCase):
    def test_midrank_ties_and_zero_label_leakage(self):
        result=m.midrank({"Plant A":1.,"Plant B":1.,"Plant C":2.,"Plant D":3.})
        self.assertAlmostEqual(result["Plant A"],-2/3)
        self.assertAlmostEqual(result["Plant B"],-2/3)
        self.assertAlmostEqual(result["Plant C"],1/3)
        self.assertAlmostEqual(result["Plant D"],1.)
        with self.assertRaises(ValueError):m.midrank({"Taxon":2.})

    def test_training_only_taxon_zone_centering_no_test_coupling(self):
        keys=[("A","1"),("A","1"),("A","3"),("B","1")]
        rank=np.array([-1.,1.,.5,.8])
        train,test=m.centered_rank_for_fold(rank,keys,np.array([0,1,2]),np.array([3]))
        self.assertEqual(train.tolist(),[-1.,1.,0.])
        self.assertEqual(test.tolist(),[0.])
        tr,te=m.centered_rank_for_fold(rank,keys,np.array([0,2,3]),np.array([1]))
        self.assertAlmostEqual(te[0],2.)
        self.assertAlmostEqual(np.mean(tr),0.)

    def test_heldout_road_holds_all_source_sites_both_plots(self):
        rows=[]
        for reg in ("A","B"):
            for road in ("a","b","c"):
                for site in range(5):
                    for plot in m.PLOTS:
                        for taxon in range(10):
                            rank=2*((taxon+site+len(road))%10)/9-1
                            rows.append({
                                "region":reg,"road":road,"site":str(site),"plot":plot,
                                "initial_year":2017,"year_gap":5,"taxon":"Taxon"+str(taxon),
                                "xrank":rank,"re_recorded":int(rank>0)
                            })
        X,y,byroad,qa,keys=m.make_design(rows)
        self.assertEqual(len(byroad),6)
        self.assertEqual(qa["distinct_source_taxa"],10)
        self.assertEqual(qa["distinct_source_taxon_zone"],20)
        r=('A','a')
        test=np.array(byroad[r],dtype=int)
        tr=np.setdiff1d(np.arange(len(rows)),test)
        self.assertTrue(all((rows[i]["region"],rows[i]["road"])!=r for i in tr))
        self.assertEqual(len(test),5*2*10)
        p2=m.fit_prediction(X[m.LABELS[2]][tr],X[m.LABELS[2]][test],y[tr])
        p3=m.fit_prediction(X[m.LABELS[3]][tr],X[m.LABELS[3]][test],y[tr])
        b2=-np.mean(y[test]*np.log(p2)+(1-y[test])*np.log1p(-p2))
        b3=-np.mean(y[test]*np.log(p3)+(1-y[test])*np.log1p(-p3))
        self.assertLess(b3,b2)
        self.assertEqual(len(p2),len(test))


if __name__=="__main__":unittest.main()
