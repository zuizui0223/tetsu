"""v29.3 synthetic validators; no fabricated field effect in publication results."""
import math
import sys
import unittest
import numpy as np
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import validate_Europe_original_4rail_studies_genus_trait_label_null_v293 as m

def synthetic_study():
    names=[f'Genus{chr(65+g)} species{chr(97+j)}'
           for g in range(20) for j in range(3)]
    traits={taxon:((i*7%31)*.12,(i*13%41)*.11)
            for i,taxon in enumerate(names)}
    sites=[]
    for g in range(3):
        for i in range(4):
            pool={names[(j+i*5+g*20)%60] for j in range(23)}
            sites.append({'station':f'site_{g}_{i}','group':f'field_region{g}',
                          'taxa':pool})
    return sites,traits

class GenusConstrainedRailOriginalFieldValidation(unittest.TestCase):
    def test_predeclared_source_original_study_and_four_labels(self):
        self.assertEqual(len(m.STUDIES),4)
        self.assertEqual(sum(m.N.values()),171)
        self.assertEqual(len(m.FROZEN_V285_QUINTILE_HEIGHT),4)

    def test_real_species_genus_index_matches_CWM_matrix_taxa(self):
        sites,t=synthetic_study()
        sp,groups=m.original_species_and_genus_index(sites,t)
        self.assertEqual(len(sp),len(t))
        self.assertEqual(sum(len(i) for i in groups.values()),len(sp))
        self.assertEqual(len(groups),20)
        self.assertTrue(all(len(x)==3 for x in groups.values()))

    def test_paired_bivariate_traits_are_permuted_within_original_genus(self):
        sites,t=synthetic_study()
        species,groups=m.original_species_and_genus_index(sites,t)
        T=np.array([t[k] for k in species],dtype=float)
        P=m.genus_numpy_labels(T,groups,np.random.default_rng(3))
        np.testing.assert_array_equal(T, np.array([t[k] for k in species],float))
        for g,idx in groups.items():
            self.assertEqual(sorted(map(tuple,P[idx,:])),sorted(map(tuple,T[idx,:])))
        self.assertEqual(sorted(map(tuple,P)),sorted(map(tuple,T)))

    def test_one_species_genus_immobile(self):
        T=np.array([[1.,2.],[3.,4.],[5.,6.]])
        groups={'GenusA':[0,1],'GenusB':[2]}
        P=m.genus_numpy_labels(T,groups,np.random.default_rng(5))
        np.testing.assert_array_equal(P[2],T[2])

    def test_source_support_counts_site_taxon_incidence_not_site_pairs(self):
        sites,t=synthetic_study()
        q=m.genus_source_support(sites,t)
        self.assertEqual(q['n_original_true_field_sites'],12)
        self.assertEqual(q['n_original_joint_trait_taxa'],60)
        self.assertEqual(q['n_multispecies_trait_known_genus'],20)
        self.assertEqual(q['original_source_trait_known_plot_taxon_incidence_total'],
                         sum(len(x['taxa']) for x in sites))
        self.assertTrue(q['source_genus_trait_shuffle_eligible'])

    def test_original_sites_and_species_identity_frozen_through_placebo(self):
        sites,t=synthetic_study()
        original=[{'station':x['station'],'group':x['group'],'taxa':x['taxa'].copy()}
                  for x in sites]
        r=m.source_genus_null(sites,t,m.STUDIES[0],draws=7)
        self.assertEqual(sites,original)
        self.assertEqual(r['status'],'SOURCE_GENUS_TRAIT_NULL_ESTIMABLE')
        self.assertEqual(r['n_technical_label_draws'],14)

    def test_all_true_height_seed_joint_distances_finite(self):
        sites,t=synthetic_study()
        r=m.source_genus_null(sites,t,m.STUDIES[0],draws=7)
        self.assertEqual(set(r['source_effects']),set(m.METRIC_KEYS))
        for key,x in r['source_effects'].items():
            self.assertTrue(math.isfinite(x['observed_minus_genus_label_null']))
            self.assertTrue(0<=x['technical_upper_rank_not_field_p']<=1)

    def test_frozen_source_math_repeatable(self):
        sites,t=synthetic_study()
        a=m.source_genus_null(sites,t,m.STUDIES[0],draws=7)
        b=m.source_genus_null(sites,t,m.STUDIES[0],draws=7)
        self.assertEqual(a,b)

    def test_genus_support_fails_closed_does_not_create_zero_ecology(self):
        # Single-genera source pool: no true genus trait permutation.
        site=[{'station':f's{i}','group':f'group{i//3}',
               'taxa':{f'Genus{chr(65+j)} speciesa' for j in range(10)}}
              for i in range(6)]
        t={f'Genus{chr(65+j)} speciesa':(float(j),float(j)) for j in range(10)}
        out=m.source_genus_null(site,t,m.STUDIES[0],draws=7)
        self.assertIsNone(out['source_effects'])
        self.assertEqual(out['status'],'NOT_IDENTIFIABLE_GENUS_STRATIFIED_TRAIT_SHUFFLE')

if __name__=='__main__':unittest.main()
