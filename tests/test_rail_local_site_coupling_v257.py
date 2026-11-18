"""v25.7 synthetic QA for railway same-site vs re-paired-site null."""
import random
import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import analyze_rail_local_site_coupling_v257 as m


class RailwayLocalPairingTests(unittest.TestCase):
    def test_uniform_derangement_is_bijective_and_no_fixed_sites(self):
        ids=list(range(15))
        rng=random.Random(5)
        for _ in range(100):
            d=m.derangement(ids,rng)
            self.assertEqual(set(d),set(ids))
            self.assertEqual(set(d.values()),set(ids))
            self.assertTrue(all(k!=v for k,v in d.items()))

    def test_pair_swap_for_two_sites_is_forced(self):
        self.assertEqual(m.derangement(["A","B"],random.Random(42)),{"A":"B","B":"A"})
        with self.assertRaises(ValueError):
            m.derangement(["A"],random.Random(1))

    def test_taxonomic_Jaccard_similarity_has_correct_range(self):
        self.assertEqual(m.jac_similarity({1,2},{1,2}),1.0)
        self.assertEqual(m.jac_similarity({1},{2}),0.0)
        self.assertAlmostEqual(m.jac_similarity({1,2},{2,3}),1/3)
        with self.assertRaises(ValueError):
            m.jac_similarity(set(),{3})

    def test_source_identity_coupling_detects_synthetic_local_species_pools(self):
        # Exactly same sets at each genuine site, disjoint across sites.
        # Pairing reference to another site reduces the similarity to zero.
        n=15
        sites={i:({2*i+1,2*i+2},{2*i+1,2*i+2}) for i in range(n)}
        obs=sum(m.jac_similarity(a,b) for a,b in sites.values())/n
        rng=random.Random(201)
        null=[]
        ids=list(sites)
        for _ in range(199):
            perm=m.derangement(ids,rng)
            null.append(sum(m.jac_similarity(sites[i][0],sites[perm[i]][1]) for i in ids)/n)
        out=m.diagnostic(obs,null)
        self.assertEqual(obs,1.)
        self.assertEqual(out["different_site_label_null_mean"],0.)
        self.assertEqual(out["same_site_minus_reassigned_mean"],1.)
        self.assertEqual(out["n_reassignments_equal_or_more_similar"],0)
        self.assertEqual(out["upper_tail_mapping_rank_plusone_NOT_a_field_p_value"],.005)

    def test_null_invariant_to_identical_community_lists(self):
        sites={i:({1,2},{1,2}) for i in range(5)}
        rng=random.Random(8)
        obs=1.
        null=[]
        for _ in range(99):
            perm=m.derangement(list(sites),rng)
            null.append(sum(m.jac_similarity(sites[i][0],sites[perm[i]][1]) for i in sites)/len(sites))
        out=m.diagnostic(obs,null)
        self.assertEqual(out["same_site_minus_reassigned_mean"],0.)
        self.assertEqual(out["upper_tail_mapping_rank_plusone_NOT_a_field_p_value"],1.)


if __name__=="__main__":
    unittest.main()
