"""v29.4 synthetic taxonomy QA and family-constrained original field tests.

All synthetic species are mathematics fixtures, never field observations.
"""
import copy, math, sys, unittest
from pathlib import Path
from collections import Counter
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import validate_Europe_original_4rail_family_genus_taxon_null_v294 as m

def test_fixture():
    names=[f'Genus{chr(65+g)} species{chr(97+j)}' for g in range(24) for j in range(3)]
    traits={t:((i*11%37)*.13,(i*13%31)*.09) for i,t in enumerate(names)}
    taxonomy={}
    for g in range(24):
        genus=f'Genus{chr(65+g)}'
        taxonomy[genus]={'status':'EXACT_ACCEPTED_GENUS_FAMILY',
                         'family':f'Family{g//3}aceae','confidence':99}
    sites=[]
    for geo in range(3):
        for site in range(4):
            idx=geo*4+site
            taxa={names[(i+idx*6)%len(names)] for i in range(30)}
            sites.append({'station':f'site{idx}','group':f'region{geo}','taxa':taxa})
    return sites,traits,taxonomy

class GenuineOriginalFamilyNullMath(unittest.TestCase):
    def test_api_genus_accepts_only_exact_plants_with_high_confidence(self):
        ok={'matchType':'EXACT','rank':'GENUS','kingdom':'Plantae',
            'genus':'Quercus','family':'Fagaceae','confidence':99,
            'usageKey':1128,'familyKey':42}
        self.assertEqual(m.validate_gbif_genus_taxonomy('Quercus',ok)['family'],'Fagaceae')
        for k,v in [('matchType','FUZZY'),('rank','SPECIES'),
                    ('kingdom','Animalia'),('genus','Quercusx'),
                    ('family',''),('confidence',60)]:
            bad=ok.copy();bad[k]=v
            self.assertNotEqual(m.validate_gbif_genus_taxonomy('Quercus',bad)['status'],
                                'EXACT_ACCEPTED_GENUS_FAMILY')
        self.assertEqual(m.validate_gbif_genus_taxonomy('Quercus',None)['status'],
                         'NOT_IDENTIFIABLE_MALFORMED')

    def test_taxonomic_family_groups_partition_entire_original_study(self):
        sites,traits,taxonomy=test_fixture()
        names=sorted(traits)
        groups=m.family_index(names,taxonomy)
        self.assertEqual(len(groups),8)
        self.assertEqual(sorted(i for group in groups.values() for i in group),
                         list(range(len(names))))
        self.assertEqual(sorted(map(len,groups.values())),[9]*8)

    def test_ambiguous_original_genus_does_not_get_fake_family(self):
        sites,traits,taxonomy=test_fixture()
        missing={'GenusA':{'status':'NOT_IDENTIFIABLE_TAXONOMY'}}
        t=taxonomy|missing
        names=sorted(traits)
        groups=m.family_index(names,t)
        unknown=[k for k in groups if k.startswith('UNRESOLVED:')]
        self.assertEqual(len(unknown),3)
        self.assertTrue(all(len(groups[k])==1 for k in unknown))

    def test_source_family_support_not_field_pairs(self):
        sites,traits,taxonomy=test_fixture()
        gate=m.original_family_source_gate(sites,traits,taxonomy)
        self.assertTrue(gate['source_family_conditioned_null_eligible'])
        self.assertEqual(gate['n_original_true_field_sites'] if 'n_original_true_field_sites' in gate
                         else gate['n_original_sites'],12)
        self.assertGreaterEqual(gate['n_multiple_source_genera_families'],8)
        self.assertGreater(gate['source_exchangeable_multi_genus_family_incidence_fraction'],.9)

    def test_family_bivariate_vector_multiset_preserved_within_each_family(self):
        sites,traits,taxonomy=test_fixture()
        names=sorted(traits)
        group=m.family_index(names,taxonomy)
        T=np.array([traits[t] for t in names])
        arr=m.shuffle_within_family(T,group,np.random.default_rng(23))
        np.testing.assert_array_equal(T,np.array([traits[t] for t in names]))
        self.assertEqual(sorted(map(tuple,T)),sorted(map(tuple,arr)))
        for key,idx in group.items():
            self.assertEqual(sorted(map(tuple,T[idx,:])),sorted(map(tuple,arr[idx,:])))

    def test_one_singleton_unknown_is_immobile(self):
        sites,traits,taxonomy=test_fixture()
        taxonomy['GenusA']={'status':'NOT_IDENTIFIABLE_TAXONOMY'}
        names=sorted(traits)
        groups=m.family_index(names,taxonomy)
        T=np.array([traits[t] for t in names])
        R=m.shuffle_within_family(T,groups,np.random.default_rng(888))
        for i,n in enumerate(names):
            if n.startswith("GenusA "):
                np.testing.assert_array_equal(T[i,:],R[i,:])

    def test_family_source_effect_and_original_matrix_unchanged(self):
        sites,traits,taxonomy=test_fixture()
        before=copy.deepcopy(sites)
        out=m.family_null_one_study(sites,traits,taxonomy,draws=7)
        self.assertEqual(out['status'],'ORIGINAL_SOURCE_FAMILY_NULL_ESTIMATED')
        self.assertEqual(sites,before)
        self.assertEqual(set(out['effects']),set(m.METRICS))
        for v in out['effects'].values():
            self.assertTrue(math.isfinite(v['observed_minus_family_trait_null']))
            self.assertEqual(out['technical_taxon_label_shuffles_per_null'],14)

    def test_insufficient_family_knowledge_reports_not_identifiable(self):
        sites,traits,taxonomy=test_fixture()
        t={g:{'status':'NOT_IDENTIFIABLE_TAXONOMY'} for g in taxonomy}
        result=m.family_null_one_study(sites,traits,t,draws=7)
        self.assertEqual(result['status'],'NOT_IDENTIFIABLE_ORIGINAL_FAMILY_TAXONOMY_SUPPORT')
        self.assertIsNone(result['effects'])

    def test_offline_mock_preserves_original_genus_keys(self):
        sites,traits,taxonomy=test_fixture()
        observed=m.lookup_study_genera(set(taxonomy),offline=taxonomy)
        self.assertEqual(observed,taxonomy)
        with self.assertRaises(ValueError):
            m.lookup_study_genera(set(taxonomy)|{'Unrepresented'},offline=taxonomy)

    def test_repeatable_family_null_original_study_not_new_field_data(self):
        sites,traits,taxonomy=test_fixture()
        a=m.family_null_one_study(sites,traits,taxonomy,draws=7)
        b=m.family_null_one_study(sites,traits,taxonomy,draws=7)
        self.assertEqual(a,b)

if __name__=='__main__':
    unittest.main()
