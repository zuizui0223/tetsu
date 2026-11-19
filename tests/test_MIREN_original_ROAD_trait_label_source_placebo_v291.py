"""v29.1 math tests: original incidence/paired traits must remain invariant."""
import math
import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import validate_MIREN_ROAD_original_site_frozen_trait_label_counterfactual_v291 as m

def fixture():
    names=[f'Species{i}' for i in range(40)]
    traits={n:(float(i)*.14,float((i*17)%43)*.11) for i,n in enumerate(names)}
    roads=[]
    for rid in range(2):
        sites=[]
        for j in range(5):
            sites.append({'roadside':set(names[(3*rid+2*j):(3*rid+2*j+11)]),
                          'interior':set(names[(3*rid+3*j+6):(3*rid+3*j+18)]),
                          'region':'SYN','road':str(rid),'site':str(j),'year':2016})
        roads.append((f'road{rid}',sites))
    return roads,traits

class OriginalFieldSourceNullMath(unittest.TestCase):
    def test_two_sites_same_road_hand_calculation(self):
        roads=[('genuineRoadExample',[
            {'roadside':{'A'},'interior':{'A'}},
            {'roadside':{'B'},'interior':{'A'}}])]
        traits={'A':(1.,0.),'B':(4.,4.)}
        x=m.source_mean_distances(roads,traits)
        self.assertEqual(x['roadside']['joint'],5.)
        self.assertEqual(x['interior']['joint'],0.)
        self.assertEqual(x['road_minus_interior']['height'],3.)
        self.assertEqual(x['road_minus_interior']['seed'],4.)
        self.assertEqual(x['road_minus_interior']['joint'],5.)

    def test_exact_original_incidence_never_mutated(self):
        roads,traits=fixture()
        previous=[(road,[{k:(v.copy() if isinstance(v,set) else v)
                             for k,v in row.items()} for row in rows]) for road,rows in roads]
        _=m.region_null(roads,traits,draws=5)
        self.assertEqual(roads,previous)

    def test_species_frequency_quintiles_unique(self):
        roads,_=fixture()
        occ=m.selected_species_with_frequency(roads)
        groups=m.occupancy_quintile_groups(occ)
        self.assertEqual(len(groups),5)
        self.assertEqual(sum(map(len,groups)),len(occ))
        self.assertEqual(set().union(*(set(g) for g in groups)),set(occ))

    def test_2D_pair_vector_correlation_preserved(self):
        roads,traits=fixture()
        import random
        occ=m.selected_species_with_frequency(roads)
        local={k:traits[k] for k in occ}
        for groups in ([list(local)],m.occupancy_quintile_groups(occ)):
            new=m.permute_region_labels(local,groups,random.Random(777))
            self.assertEqual(sorted(new.values()),sorted(local.values()))
            self.assertEqual(set(new),set(local))
            self.assertEqual(len(new),len(local))
            self.assertAlmostEqual(sum(x[0]*x[1] for x in new.values()),
                                   sum(x[0]*x[1] for x in local.values()),places=11)

    def test_original_region_and_road_weighting(self):
        sites=[
            ('roadA',[{'roadside':{'A'},'interior':{'A'}},
                      {'roadside':{'B'},'interior':{'A'}}]),
            ('roadB',[{'roadside':{'C'},'interior':{'A'}},
                      {'roadside':{'C'},'interior':{'A'}},
                      {'roadside':{'C'},'interior':{'A'}}])]
        t={'A':(0.,0.),'B':(3.,4.),'C':(5.,12.)}
        x=m.source_mean_distances(sites,t)
        self.assertEqual(x['roadside']['joint'],2.5)
        self.assertEqual(x['interior']['joint'],0.)
        self.assertEqual(x['road_minus_interior']['height'],1.5)

    def test_both_real_habitats_and_two_source_modes(self):
        roads,traits=fixture()
        x=m.region_null(roads,traits,draws=5)
        self.assertEqual(set(x['source_null_by_mode']),set(m.MODES))
        self.assertEqual(x['n_source_original_roads'],2)
        self.assertEqual(x['n_source_field_sites'],10)
        for mode in m.MODES:
            for z in ('roadside','interior','road_minus_interior'):
                for met in ('joint','height','seed'):
                    y=x['source_null_by_mode'][mode]['effects'][z][met]
                    self.assertTrue(math.isfinite(y['observed_minus_trait_label_null']))
                    self.assertEqual(y['n_technical_assignments_not_field_studies'],10)

    def test_deterministic_null_no_fabricated_replication(self):
        roads,traits=fixture()
        self.assertEqual(m.region_null(roads,traits,draws=5),
                         m.region_null(roads,traits,draws=5))

    def test_missing_original_taxon_trait_fails(self):
        roads,traits=fixture()
        with self.assertRaises(ValueError):
            m.region_null(roads,{},draws=5)

if __name__=='__main__':unittest.main()
