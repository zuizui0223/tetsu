import unittest
import numpy as np
from aggregate_world_ESA_CCI_biomass import overlap_matrix,aggregate,edges
from analyze_world_railstatus_plant_bird_biomass import fit,ACTIVE,INACTIVE,BIO,IA,II,OUT
class TestBiomass(unittest.TestCase):
    def test_overlap(self):
        x=overlap_matrix(np.array([0.,1.5,2.]),0.,2,1.)
        y=overlap_matrix(np.array([0.,1.,2.]),0.,2,1.)
        a,s,q,qa=aggregate(np.ma.array([[10.,20.],[30.,40.]],mask=[[False,True],[False,False]]),y,x)
        self.assertAlmostEqual(a.sum(),3.5);self.assertAlmostEqual(s.sum(),80.)
        self.assertAlmostEqual(a[0,1],.5)
    def test_axis(self):
        i,e=edges([.15,.05],None,0.,1.)
        self.assertEqual(list(i),[1,0]);np.testing.assert_allclose(e,[0.,.1,.2],atol=1e-8)
    def test_paired_interaction(self):
        rng=np.random.default_rng(42);r=[]
        for b in range(30):
            for i in range(8):
                a,z,v=rng.uniform(.1,2,3);p=a*v+.4*z*v+b;bird=.2*a*v+.1*z*v+b
                r.append({"continent":"X","block":b,ACTIVE:a,INACTIVE:z,BIO:v,IA:a*v,II:z*v,OUT[0]:p,OUT[1]:bird,OUT[2]:p-bird})
        q=fit(r,(ACTIVE,INACTIVE,BIO,IA,II),OUT)
        self.assertEqual(q["status"],"estimated")
        self.assertAlmostEqual(q["effects"]["active_minus_inactive_x_AGB"][OUT[2]]["estimate"],.5,places=8)
if __name__=="__main__":unittest.main()
