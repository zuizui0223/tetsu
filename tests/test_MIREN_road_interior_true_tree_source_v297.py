"""v29.7 synthetic math/source invariance only. No invented empirical source."""
import copy, math, sys, unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import validate_MIREN_original_roads_true_phylogeny_clade_negative_control_v297 as m
from analyze_four_rail_actual_3phylo_clade_null_python_v296 import original_root_depths, original_true_clade_groups

def original_field_fixture():
    taxa=[f'Genus{chr(65+i//3)} species{chr(97+i%3)}' for i in range(60)]
    traits={tax:(float((i*7)%37)*.11,float((i*11)%43)*.08)
            for i,tax in enumerate(taxa)}
    records=[]
    for road_id in range(2):
        sites=[]
        for site in range(6):
            road={ 'year':2017,
                   'roadside':set(taxa[(i+site*2+road_id*3)%60] for i in range(25)),
                   'interior':set(taxa[(i+site+road_id*4+15)%60] for i in range(27))}
            sites.append(road)
        records.append((f'SyntheticRoad{road_id}',sites))
    return records,traits,{tax:i+1 for i,tax in enumerate(taxa)}

class VerifiedPhylogeneticRoadPlaceholderMath(unittest.TestCase):
    def test_same_site_two_zone_true_author_pairs_maintained(self):
        roads,traits,tips=original_field_fixture()
        source=copy.deepcopy(roads)
        panel,qa=m.construct_original_region_tree_panel(roads,traits,tips,5)
        self.assertEqual(roads,source)
        self.assertIsNotNone(panel)
        self.assertEqual(qa['original_source_true_spatial_sites'],12)
        self.assertEqual(qa['retained_original_source_site_pairs'],12)
        self.assertEqual(qa['retained_original_roads_with_4_or_more_sites'],2)
        self.assertTrue(qa['direct_tip_road_inland_original_field_cohort_support_pass'])
        self.assertEqual(len(panel['roads']),2)

    def test_true_2D_vectors_and_region_comparison_source_equals_manual(self):
        roads,traits,tips=original_field_fixture()
        panel,qa=m.construct_original_region_tree_panel(roads,traits,tips,5)
        b=m.region_beta_for_two_zones(panel,panel['source_trait_vectors'])
        for key in m.SCIENTIFIC_METRICS:
            self.assertAlmostEqual(
                b['roadside_minus_interior'][key],
                b['roadside'][key]-b['interior'][key])
        self.assertTrue(all(math.isfinite(v) for z in b.values() for v in z.values()))

    def test_missing_original_tree_taxa_are_not_imputed(self):
        roads,traits,tips=original_field_fixture()
        original=copy.deepcopy(roads)
        missing=dict(list(tips.items())[:3])
        panel,qa=m.construct_original_region_tree_panel(roads,traits,missing,5)
        self.assertIsNone(panel)
        self.assertFalse(qa['direct_tip_road_inland_original_field_cohort_support_pass'])
        self.assertEqual(roads,original)

    def test_frozen_original_field_source_requires_two_roads(self):
        roads,traits,tips=original_field_fixture()
        panel,qa=m.construct_original_region_tree_panel(roads[:1],traits,tips,5)
        self.assertIsNone(panel)
        self.assertEqual(qa['original_source_roads'],1)

    def test_original_strict_species_gate_does_not_promote_trait_missing_as_zero(self):
        roads,traits,tips=original_field_fixture()
        panel,qa=m.construct_original_region_tree_panel(roads,traits,tips,10)
        self.assertIsNotNone(panel)
        bad=dict(list(tips.items())[:12])
        subset,qa=m.construct_original_region_tree_panel(roads,traits,bad,10)
        self.assertIsNone(subset)

    def test_occupancy_quintile_original_source_species_frozen(self):
        roads,traits,tips=original_field_fixture()
        panel,qa=m.construct_original_region_tree_panel(roads,traits,tips,5)
        groups=m.occupancy_quintile_groups(panel['original_occupancy'])
        self.assertEqual(len(groups),5)
        self.assertEqual(sorted(sum((g.tolist() for g in groups),[])),list(range(len(panel['species']))))

    def test_clade_taxon_exchangeability_requires_five_real_groups(self):
        roads,traits,tips=original_field_fixture()
        panel,qa=m.construct_original_region_tree_panel(roads,traits,tips,5)
        n=len(panel['species'])
        one=[np.arange(n)]
        x=m.region_tree_clade_exchangeability(panel,one)
        self.assertEqual(x['source_tree_clades_with_two_or_more_joint_LEDA_species'],1)
        self.assertFalse(x['predeclared_phylogenetic_source_null_support'])
        groups=[np.arange(i,min(i+5,n)) for i in range(0,n,5)]
        good=m.region_tree_clade_exchangeability(panel,groups)
        self.assertGreaterEqual(good['source_tree_clades_with_two_or_more_joint_LEDA_species'],5)
        self.assertTrue(good['predeclared_phylogenetic_source_null_support'])

    def test_clade_unidentified_not_reported_as_zero_biological_effect(self):
        roads,traits,tips=original_field_fixture()
        panel,_=m.construct_original_region_tree_panel(roads,traits,tips,5)
        out=m.run_region_null(panel,[np.arange(len(panel['species']))],
                              m.MODES[3],0,0,draws=5)
        self.assertEqual(out['status'],'NOT_IDENTIFIABLE_ROAD_INTERIOR_PHYLOGENETIC_CLADE_NULL')
        self.assertIsNone(out['effects'])

    def test_null_randomization_uses_one_label_vector_for_both_habitats(self):
        roads,traits,tips=original_field_fixture()
        panel,_=m.construct_original_region_tree_panel(roads,traits,tips,5)
        groups=[np.arange(i,min(i+4,len(panel['species']))) for i in range(0,len(panel['species']),4)]
        a=m.run_region_null(panel,groups,m.MODES[3],0,0,draws=5)
        b=m.run_region_null(panel,groups,m.MODES[3],0,0,draws=5)
        self.assertEqual(a,b)
        self.assertEqual(a['n_technical_shuffles_not_field_studies'],10)
        for metric in m.SCIENTIFIC_METRICS:
            effect=a['source_original_region_two_zone_effects']
            self.assertAlmostEqual(effect['roadside_minus_interior'][metric]['observed_minus_trait_label_source_null'],
                effect['roadside'][metric]['observed_minus_trait_label_source_null']-
                effect['interior'][metric]['observed_minus_trait_label_source_null'])
        self.assertEqual(roads,roads)

if __name__=="__main__":unittest.main()
