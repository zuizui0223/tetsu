"""Mathematical tests: toy trees/communities NEVER biological empirical results."""
import copy,sys,unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import analyze_four_rail_actual_3phylo_clade_null_python_v296 as v

class TrueTreeOriginalFieldNullMath(unittest.TestCase):
    def setUp(self):
        self.edge=np.array([[5,6],[5,7],[6,1],[6,2],[7,3],[7,4]],dtype=int)
        self.parent,self.depth,self.tree=v.original_root_depths(self.edge,np.ones(6),4,3)

    def test_original_tree_root_to_tip_distance_topology(self):
        self.assertEqual(self.tree['original_root_node'],5)
        self.assertAlmostEqual(self.tree['original_total_root_to_tip_distance_median'],2.)
        self.assertTrue(np.allclose(self.depth[1:5],2))

    def test_actual_cut_returns_monophyletic_nonoverlapping_source_species_sets(self):
        groups=v.original_true_clade_groups([1,2,3,4],self.parent,self.depth,1.)
        self.assertEqual(sorted(map(len,groups)),[2,2])
        self.assertEqual(sorted(sum([g.tolist() for g in groups],[])),[0,1,2,3])
        self.assertEqual(sorted(map(len,v.original_true_clade_groups(
            [1,2,3,4],self.parent,self.depth,0.1))),[4])

    def test_joint_trait_vectors_permuted_together_and_originals_frozen(self):
        h=np.array([[1.,10.],[2.,50.],[9.,13.],[7.,3.]])
        groups=v.original_true_clade_groups([1,2,3,4],self.parent,self.depth,1.)
        out=v.permute_paired_traits(h,groups,np.random.default_rng(23))
        self.assertEqual(sorted(map(tuple,out)),sorted(map(tuple,h)))
        np.testing.assert_array_equal(h,np.array([[1.,10.],[2.,50.],[9.,13.],[7.,3.]]))
        for g in groups:
            self.assertEqual(sorted(map(tuple,h[g])),sorted(map(tuple,out[g])))

    def test_true_author_geo_ids_not_truncated_to_one_character(self):
        P=np.array([[1,1,1,0],[1,1,1,1],[1,0,1,1]],dtype=np.int8)
        geo=np.asarray(['10','20','10'],dtype=object)
        self.assertEqual(sum(geo[0]==geo[1:]),1)
        self.assertEqual(len(set(geo)),2)
        W=P/P.sum(axis=1)[:,None]
        traits=np.array([[0.,0.],[1.,1.],[2.,2.],[3.,3.]])
        y=v.source_CWM_pairwise_beta(W,traits,np.array([0]),np.array([2]))
        self.assertTrue(y['height']>0)

    def test_no_identifiability_when_original_clades_all_immobile(self):
        P=np.eye(8,dtype=np.int8)
        groups=[np.array([i]) for i in range(8)]
        geo=np.asarray(['A','A','A','A','B','B','B','B'],dtype=object)
        qa=v.phylo_exchangeability(P,groups,geo)
        self.assertFalse(qa['predeclared_real_phylo_exchangeability_gate'])
        self.assertEqual(qa['source_site_species_incidence_exchangeable_fraction'],0)

    def test_true_occupancy_5_original_strata(self):
        P=np.zeros((10,31),dtype=np.int8)
        for i in range(10):
            P[i,(i*3)%31: min((i*3)%31+11,31)]=1
        g=v.source_occupancy_quintiles(P)
        self.assertEqual(len(g),5)
        self.assertEqual(sorted(sum([x.tolist() for x in g],[])),list(range(31)))
        self.assertTrue(all(len(x)>=6 for x in g))

    def test_frozen_original_site_matrix_between_site_CWM(self):
        P=np.array([[1,1,0,0],[0,0,1,1]],dtype=np.int8)
        H=np.array([[1.,2.],[1.,2.],[4.,6.],[4.,6.]])
        d=v.source_CWM_pairwise_beta(P/2.,H,np.array([0]),np.array([1]))
        self.assertAlmostEqual(d['height'],3.)
        self.assertAlmostEqual(d['seed'],4.)
        self.assertAlmostEqual(d['joint'],5.)

    def test_technical_draws_not_field_replicates(self):
        P=np.ones((10,20),dtype=np.int8)
        geo=np.asarray(['G0']*5+['G1']*5,dtype=object)
        h=np.array([[i*.1,(i*3%11)*.1] for i in range(20)])
        groups=[np.arange(i*2,(i+1)*2) for i in range(10)]
        a=v.source_null(P,h,groups,geo,np.array([0]),np.array([1]),0,0,1,draws=5)
        self.assertEqual(a['status'],'SOURCE_REAL_PHYLOGENETIC_CLADE_DEPTH_ESTIMABLE')
        self.assertEqual(a['n_technical_label_assignments_not_field_studies'],10)
        self.assertTrue(a['original_species_site_matrices_FROZEN'])
        self.assertEqual(v.source_null(P,h,groups,geo,np.array([0]),np.array([1]),0,0,1,draws=5),a)

    def test_tree_root_source_cycle_fails_closed(self):
        with self.assertRaises(ValueError):
            v.original_root_depths(np.array([[5,6],[6,5],[5,1],[5,2]]),np.ones(4),2,4)

if __name__=="__main__":
    unittest.main()
