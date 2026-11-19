"""Synthetic mathematical tests only. Never count invented species as empirical data."""
import copy
import math
import random
import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import validate_MIREN_genus_conditioned_trait_label_null_v292 as v

def two_roads():
    taxa=[f'Genus{chr(65+g)} species{chr(97+j)}'
          for g in range(20) for j in range(3)]
    traits={taxon:(float((i*7)%31)*.12,float((i*13)%37)*.16)
            for i,taxon in enumerate(taxa)}
    roads=[]
    for road in range(2):
        sites=[]
        for site in range(5):
            # All original sites have 30 observed species and >25% in
            # multispecies genera, with independent road/site design.
            one={taxa[(x+site*2+road)%60] for x in range(31)}
            three={taxa[(x+site+road*3+11)%60] for x in range(32)}
            sites.append({'roadside':one,'interior':three,
                          'road':'road'+str(road),'site':'site'+str(site)})
        roads.append(('road'+str(road),sites))
    return roads,traits

class TestGenusSourceNull(unittest.TestCase):
    def test_real_binomial_genus_parser_fail_closed(self):
        self.assertEqual(v.genus('Festuca rubra'),'Festuca')
        self.assertEqual(v.genus('Carex nigra'),'Carex')
        for raw in ('', 'unknown', 'Species sp.', 'Genus species subspecies'):
            with self.assertRaises(ValueError):
                v.genus(raw)

    def test_genus_groups_partition_taxa_exactly(self):
        sites,traits=two_roads()
        groups=v.genus_groups(set(traits))
        self.assertEqual(len(groups),20)
        self.assertTrue(all(len(x)==3 for x in groups))
        self.assertEqual(set().union(*map(set,groups)),set(traits))

    def test_paired_vectors_exactly_stay_in_genus(self):
        sites,traits=two_roads()
        groups=v.genus_groups(set(traits))
        original=copy.deepcopy(traits)
        shuffled=v.permute_within_genus(traits,groups,random.Random(8))
        self.assertEqual(traits,original)
        self.assertEqual(set(shuffled),set(traits))
        self.assertEqual(sorted(shuffled.values()),sorted(traits.values()))
        for members in groups:
            self.assertEqual(sorted(shuffled[t] for t in members),
                             sorted(traits[t] for t in members))

    def test_singleton_genus_species_is_immobile(self):
        traits={'GenusA alpha':(1.,2.),'GenusA beta':(5.,6.),
                'GenusB gamma':(3.,4.)}
        groups=v.genus_groups(set(traits))
        sh=v.permute_within_genus(traits,groups,random.Random(4))
        self.assertEqual(sh['GenusB gamma'],traits['GenusB gamma'])

    def test_source_taxon_incidence_unchanged(self):
        roads,traits=two_roads()
        before=copy.deepcopy(roads)
        a=v.genus_placebo(roads,traits,draws=7)
        self.assertEqual(roads,before)
        self.assertEqual(a['status'],'SOURCE_GENUS_CONDITIONED_EMPIRICAL_ESTIMABLE')
        self.assertTrue(a['original_site_taxon_presence_frequency_and_community_fixed'])

    def test_predeclared_source_support_genus_gate(self):
        roads,traits=two_roads()
        qa=v.genus_support(roads,traits)
        self.assertGreaterEqual(qa['n_multispecies_genera'],10)
        self.assertTrue(qa['genus_conditioned_null_identifiable'])
        self.assertEqual(qa['n_original_joint_trait_species_in_region'],
                         qa['n_exchangeable_joint_trait_species'])

    def test_no_effect_when_original_source_genus_support_fails(self):
        # No 10 multi-species original genera; a negative result must not be
        # invented from a mostly immovable genus trait null.
        r=[('only',[{'roadside':{'Gena alpha'},'interior':{'Genb beta'}},
                    {'roadside':{'Gena alpha'},'interior':{'Genb beta'}}])]
        t={'Gena alpha':(1.,2.),'Genb beta':(3.,4.)}
        out=v.genus_placebo(r,t,draws=7)
        self.assertEqual(out['status'],'NOT_IDENTIFIABLE_GENUS_CONDITIONED_TRAIT_SHUFFLE')
        self.assertIsNone(out['trait_beta_effects'])

    def test_actual_effects_and_null_deterministic(self):
        roads,traits=two_roads()
        a=v.genus_placebo(roads,traits,draws=7)
        b=v.genus_placebo(roads,traits,draws=7)
        self.assertEqual(a,b)
        for zone in ('roadside','interior','road_minus_interior'):
            for metric in ('joint','height','seed'):
                x=a['trait_beta_effects'][zone][metric]
                self.assertEqual(x['n_technical_label_assignments'],14)
                self.assertTrue(math.isfinite(x['observed_minus_within_genus_null']))

    def test_multi_site_weighting_not_pairs_as_independent_studies(self):
        roads,traits=two_roads()
        support=v.genus_support(roads,traits)
        self.assertEqual(support['original_site_habitat_trait_known_taxon_incidence']
                         ['roadside']['n_original_trait_known_site_taxon_incidence'],
                         sum(len(site['roadside']) for _,ss in roads for site in ss))
        self.assertEqual(len(roads),2)

if __name__=='__main__':
    unittest.main()
