"""v28.8 synthetic checks: correct beta math, original-site pairing and source gates."""
import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import validate_Kashmir_out_of_Europe_rail_nonrail_height_beta_v288 as m

def source_fixture(n):
    """Synthetic 2-year four-zone data, NOT empirical plant source."""
    names=tuple('Genus spec'+chr(97+i) for i in range(25))
    traits={name:(float(i),float(i*i%19)) for i,name in enumerate(names)}
    plots={}
    for j in range(n):
        site=f'Test {j}'
        for y in (2014,2017):
            plots[site,y,'Track']=set(range(0,12))
            plots[site,y,'Embankment']=set(range(0,12))
            plots[site,y,'Interior Plot 1']=set(range(7,25))
            plots[site,y,'Interior Plot 2']=set(range(10,22))
    return plots,names,traits

class TestKashmirIndependentSiteGates(unittest.TestCase):
    def test_2d_CWM_pair_distance_hand_verified(self):
        v=[{'joint':(0.,0.)},{'joint':(3.,4.)}]
        self.assertEqual(m.metric(v,'joint'),5.)
        self.assertEqual(m.metric(v,'height'),3.)
        self.assertEqual(m.metric(v,'seed'),4.)
        self.assertRaises(ValueError,m.metric,v,'xyz')

    def test_actual_site_reference_pairing_both_field_years(self):
        plots,names,traits=source_fixture(12)
        result=m.eligible_pairs(plots,names,traits,'station','Interior Plot 1')
        self.assertEqual(result['original_sites_with_both_field_years_and_four_zones'],12)
        self.assertEqual(result['n_original_sites_meeting_predeclared_trait_gate_both_years'],12)
        self.assertEqual(result['status'],'ESTIMABLE')
        self.assertIn('observed_rail_minus_inland_CWM_between_site_beta',result)
        self.assertEqual(len(result['equal_trait_known_species_per_paired_site_zone_technical_draws']),2)

    def test_no_estimate_when_study_support_insufficient(self):
        plots,names,traits=source_fixture(9)
        result=m.eligible_pairs(plots,names,traits,'station','Interior Plot 1')
        self.assertEqual(result['status'],'NOT_IDENTIFIABLE_TRAIT_SUPPORT')
        self.assertNotIn('observed_rail_minus_inland_CWM_between_site_beta',result)

    def test_missing_source_species_not_treated_as_absence(self):
        plots,names,traits=source_fixture(12)
        plots['Test 0',2014,'Track']=None
        result=m.eligible_pairs(plots,names,traits,'station','Interior Plot 1')
        self.assertEqual(result['original_sites_with_both_field_years_and_four_zones'],11)
        self.assertEqual(result['n_original_sites_meeting_predeclared_trait_gate_both_years'],11)

    def test_unobserved_trait_species_are_not_assigned_zero_values(self):
        plots,names,traits=source_fixture(11)
        traits={k:v for k,v in traits.items() if k in names[:3]}
        result=m.eligible_pairs(plots,names,traits,'station','Interior Plot 1')
        self.assertEqual(result['status'],'NOT_IDENTIFIABLE_TRAIT_SUPPORT')
        self.assertEqual(result['n_original_sites_meeting_predeclared_trait_gate_both_years'],0)

    def test_richness_rarefaction_same_target_n_per_original_site(self):
        a={f'one{i}' for i in range(6)};b={f'one{i}' for i in range(9)}
        t={k:(float(i),float(i*i)) for i,k in enumerate(sorted(a|b))}
        pairs=[{'traits':t,2014:{'Track':a,'reference':b}},
               {'traits':t,2014:{'Track':b,'reference':a}}]
        self.assertEqual(m.contrasts(pairs,2014,True,seed=99),m.contrasts(pairs,2014,True,seed=99))
        self.assertEqual(len(a),6)
        self.assertEqual(len(b),9)

    def test_source_sha_fails_closed(self):
        from tempfile import TemporaryDirectory
        with TemporaryDirectory() as tmp:
            path=Path(tmp)/'fake.xlsx'
            path.write_bytes(b'not genuine original source')
            with self.assertRaises(ValueError):m.strict_sha(path,m.DRYAD_SHA256)

if __name__=='__main__':unittest.main()
