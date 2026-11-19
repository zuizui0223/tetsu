"""Source-invariance synthetic tests only; never biological field results."""
import copy,math,sys,unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import phylogeny_occupancy_joint_label_null_v298 as m
from analyze_four_rail_actual_3phylo_clade_null_python_v296 import permute_paired_traits
from analyze_MIREN_roads_inland_joint_phylo75_occupancy_null_v298 import source_study_incidence

def fixture():
    # 40 original mathematical taxa arranged into 5 real tree clades with
    # eight taxon positions each. True source incidence stays unchanged.
    P=np.ones((10,40),dtype=np.uint8)
    P[0,0]=0
    g=np.asarray(["G1"]*5+["G2"]*5,dtype=object)
    clades=[np.arange(i,i+8) for i in range(0,40,8)]
    T=np.asarray([[float(i)*.2,float((i*11)%47)*.13] for i in range(40)])
    return P,g,clades,T

class AuthenticSourceJointNullMath(unittest.TestCase):
    def test_joint_source_clade_and_occupancy_partition_all_taxa_once(self):
        P,g,clades,T=fixture()
        occ=m.original_global_taxon_occupancy_quintiles(P)
        groups=m.joint_true_clade_and_occupancy_groups(clades,occ,40)
        self.assertEqual(sorted(sum((a.tolist() for a in groups),[])),list(range(40)))
        self.assertEqual(len(occ),5)
        self.assertTrue(all(len(x)>0 for x in groups))
        self.assertEqual(np.sum(P),399)

    def test_original_bivariate_trait_multiset_and_each_joint_stratum_fixed(self):
        P,g,clades,T=fixture()
        occ=m.original_global_taxon_occupancy_quintiles(P)
        bins=m.joint_true_clade_and_occupancy_groups(clades,occ,40)
        before=T.copy()
        R=permute_paired_traits(T,bins,np.random.default_rng(917))
        np.testing.assert_array_equal(T,before)
        self.assertEqual(sorted(map(tuple,R)),sorted(map(tuple,T)))
        for pair in bins:
            self.assertEqual(sorted(map(tuple,T[pair])),sorted(map(tuple,R[pair])))

    def test_true_original_species_site_incidence_stays_frozen(self):
        P,g,clades,T=fixture()
        before=P.copy()
        q=m.joint_true_clade_and_occupancy_groups(
            clades,m.original_global_taxon_occupancy_quintiles(P),40)
        _=permute_paired_traits(T,q,np.random.default_rng(123))
        np.testing.assert_array_equal(P,before)

    def test_joint_source_support_reports_true_movability_not_label_draw_p(self):
        P,g,clades,T=fixture()
        q=m.joint_true_clade_and_occupancy_groups(clades,m.original_global_taxon_occupancy_quintiles(P),40)
        support=m.original_source_joint_exchangeability(P,q,g)
        self.assertTrue(support["joint_phylogenetic_clade_AND_source_occupancy_null_identifiable"])
        self.assertGreaterEqual(support["n_source_phylo_clade_x_occupancy_strata_size_ge2"],5)
        self.assertGreaterEqual(support["fraction_original_source_site_species_incidence_exchangeable"],.25)
        self.assertEqual(support["n_original_source_geographical_field_groups_with_ge25pct_mobility"],2)

    def test_singleton_joint_intersections_fail_closed_as_not_identifiable(self):
        P,g,clades,T=fixture()
        singles=[np.asarray([i]) for i in range(40)]
        q=m.joint_true_clade_and_occupancy_groups(singles,m.original_global_taxon_occupancy_quintiles(P),40)
        qa=m.original_source_joint_exchangeability(P,q,g)
        self.assertEqual(qa["support_status"],"NOT_IDENTIFIABLE_PHYLO_OCCUPANCY_JOINT_LABEL_NULL")
        self.assertEqual(qa["n_source_phylo_clade_x_occupancy_strata_size_ge2"],0)

    def test_missing_or_double_original_phylo_tip_taxa_rejected(self):
        P,g,clades,T=fixture()
        bins=m.original_global_taxon_occupancy_quintiles(P)
        with self.assertRaises(ValueError):
            m.joint_true_clade_and_occupancy_groups(clades[:-1],bins,40)
        with self.assertRaises(ValueError):
            m.joint_true_clade_and_occupancy_groups(clades+[clades[0]],bins,40)

    def test_original_same_site_road_and_interior_site_weighting_not_split(self):
        n=40
        T=np.asarray([[i*.2,(i*3)%21*.1] for i in range(n)])
        site_rows=[]
        for i in range(2):
            a=np.ones((4,n),dtype=float)/n
            b=np.ones((4,n),dtype=float)/n
            b[0,0]=0
            b=b/b.sum(axis=1)[:,None]
            site_rows.append({"roadside":a,"interior":b})
        source={"species":[f"taxon{i}" for i in range(n)],
                "roads":site_rows,
                "original_occupancy":np.asarray([16 if i else 14 for i in range(n)])}
        P,g,zones,roads=source_study_incidence(source)
        self.assertEqual(P.shape,(16,40))
        self.assertEqual(set(zones),{"roadside","interior"})
        self.assertEqual(len(roads),2)
        self.assertFalse(any("Site" in k for k in roads))
        for key in roads:
            self.assertEqual(roads[key][0].shape,roads[key][1].shape)

    def test_one_road_cannot_be_treated_as_two_independent_field_studies(self):
        P,g,clades,T=fixture()
        groups=m.joint_true_clade_and_occupancy_groups(clades,m.original_global_taxon_occupancy_quintiles(P),40)
        q=m.original_source_joint_exchangeability(
            P,groups,g,habitat_Ps={"roadside":P[:5],"interior":P[5:]},
            original_road_strata={"road_1":[P[:5],P[5:]]})
        self.assertFalse(q["joint_phylogenetic_clade_AND_source_occupancy_null_identifiable"])
        self.assertEqual(q["n_original_roadside_and_inland_field_roads_with_ge25pct_mobility"],1)

    def test_road_and_inland_true_paired_original_source_both_check(self):
        P,g,clades,T=fixture()
        groups=m.joint_true_clade_and_occupancy_groups(clades,m.original_global_taxon_occupancy_quintiles(P),40)
        q=m.original_source_joint_exchangeability(
            P,groups,g,habitat_Ps={"roadside":P[:5],"interior":P[5:]},
            original_road_strata={"road_1":[P[:2],P[5:7]],
                                  "road_2":[P[2:5],P[7:10]]})
        self.assertTrue(q["joint_phylogenetic_clade_AND_source_occupancy_null_identifiable"])

    def test_cwm_unsupported_phylogeny_effect_not_an_assigned_zero(self):
        P,g,clades,T=fixture()
        q=m.original_source_joint_exchangeability(
            P,[np.asarray([i]) for i in range(40)],g)
        self.assertFalse(q["joint_phylogenetic_clade_AND_source_occupancy_null_identifiable"])
        self.assertNotIn("biological_effect_zero",q)

if __name__=="__main__":unittest.main()
