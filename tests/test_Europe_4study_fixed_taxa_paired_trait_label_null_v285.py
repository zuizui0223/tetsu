"""v28.5 synthetic tests: paired original 2D trait-label shuffling with
original plant site×taxon, geography and site richness all frozen."""
import sys
import unittest
import numpy as np
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import audit_Europe_4independent_railstudies_fixed_community_trait_label_null_v285 as m


def toy_original_study():
    species=[f"taxon_{i}" for i in range(40)]
    traits={s:((i-20)/10.,((i*11)%41-20)/10.) for i,s in enumerate(species)}
    sites=[]
    for block in range(3):
        for site in range(4):
            start=block*10
            known=set(species[start:start+10])
            known.update(species[(site+11)%40:(site+11)%40+5])
            sites.append({"station":f"site{block}_{site}","group":f"region{block}",
                          "taxa":known})
    return sites,traits


class OriginalRailFixedSiteTraitShuffleTests(unittest.TestCase):
    def test_known_2D_pairwise_CWM_matches_by_hand_and_zero_alpha(self):
        W=np.eye(3)
        T=np.array([[0.,0.],[3.,4.],[6.,8.]])
        i,j=np.triu_indices(3,1)
        within=np.array([True,False,False])
        x=m.metric_values(W,T,i,j,within)
        self.assertAlmostEqual(x["cwm2D_mean_within_original_geography"],5.)
        self.assertAlmostEqual(x["cwm_height_absolute_mean_within_original_geography"],3.)
        self.assertAlmostEqual(x["cwm_seed_mass_absolute_mean_within_original_geography"],4.)
        self.assertAlmostEqual(x["source_site_within_2D_trait_variance_trace"],0.)

    def test_original_taxon_site_incidence_and_2D_trait_correlation_not_perturbed(self):
        sites,trait=toy_original_study()
        W,T,group,qa=m.build_matrix(sites,trait)
        self.assertEqual(qa["source_original_sites"],12)
        self.assertEqual(qa["genuine_original_geography_group_count"],3)
        self.assertGreaterEqual(qa["source_actual_plant_species_with_both_traits"],25)
        originalW=W.copy()
        originalT=T.copy()
        idx=np.arange(len(T))
        result=m.permutation(W,T,qa["i"],qa["j"],qa["within"],[idx],
                             seed=87,draws=59)
        self.assertEqual(result["source_permutation_draws"],59)
        self.assertTrue(result["source_trait_values_multiset_and_covariance_exactly_preserved"])
        np.testing.assert_array_equal(W,originalW)
        np.testing.assert_array_equal(T,originalT)
        self.assertTrue(all(np.isfinite(v) for values in result["draws"].values() for v in values))

    def test_source_occurrence_quintiles_keep_all_actual_taxa_in_exact_one_band(self):
        sites,trait=toy_original_study()
        W,T,geo,qa=m.build_matrix(sites,trait)
        bands=m.occupancy_quintile_strata(qa["species_site_occurrence"])
        self.assertEqual(len(bands),5)
        indices=np.concatenate(bands)
        self.assertEqual(len(indices),len(T))
        self.assertEqual(len(set(indices.tolist())),len(T))
        self.assertEqual(set(indices.tolist()),set(range(len(T))))
        counts=qa["species_site_occurrence"]
        self.assertLessEqual(max(np.max(counts[b]) for b in bands)-min(np.min(counts[b]) for b in bands),len(sites))

    def test_two_seed_sources_and_complete_study_null_metrics(self):
        sites,trait=toy_original_study()
        r=m.source_study_null(sites,trait,"CZ_Czechia",draws=51)
        self.assertEqual(r["source_study"],"CZ_Czechia")
        self.assertEqual(len(r["source_null_results"]),2)
        for mode,v in r["source_null_results"].items():
            self.assertEqual(v["n_technical_non_field_trait_shuffles"],102)
            self.assertTrue(v["source_site_species_community_matrix_is_exact_original"])
            self.assertTrue(v["source_taxon_observed_prevalence_not_permutes"])
            self.assertTrue(v["source_2D_trait_vector_correlation_preserved"])
            self.assertEqual(len(v["chains"]),2)
            for name in m.METRICS:
                outcome=v["metrics"][name]
                self.assertTrue(np.isfinite(outcome["original_observed_minus_trait_randomized_null"]))
        self.assertEqual(r,m.source_study_null(sites,trait,"CZ_Czechia",draws=51))

    def test_shuffled_bivariate_pairs_not_break_original_study_trait_distribution(self):
        sites,trait=toy_original_study()
        W,T,geo,qa=m.build_matrix(sites,trait)
        source_mean=np.mean(T,axis=0)
        source_cov=np.cov(T.T)
        # A manually drawn joint-label permutation, independent of any
        # original station species or geography assignment.
        perm=np.random.default_rng(100).permutation(len(T))
        np.testing.assert_allclose(np.mean(T[perm],axis=0),source_mean)
        np.testing.assert_allclose(np.cov(T[perm].T),source_cov)
        with self.assertRaises(ValueError):
            m.permutation(W,T,qa["i"],qa["j"],qa["within"],
                          [np.arange(len(T))],seed=2,draws=49)


if __name__=="__main__":
    unittest.main()
