import math,sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent.parent/'scripts'))
from audit_Czech_local_vegetation_exact_MILP_nonoverlap_capacity_v242 import (
    solve_full_source_global_exclusivity,evaluate_matched
)
class TestExactIntegerControlSourceFeasibility(unittest.TestCase):
    def make(self,n=4,shared=False):
        out=[]
        for i in range(n):
            choices=[]
            for k in range(6):
                pid=f'common_plot_{k}' if shared else f'plot_station_{i}_item_{k}'
                choices.append({'id':pid,'study':'independent_paper_A' if k<3 else 'independent_paper_B',
                  'd':.1+k*.1,
                  'taxa':{f'Genus taxon{chr(97+x)}' for x in range(k,k+9)}})
            out.append({'id':f'rail_station_{i}',
                'stretch':f'rail_stretch_{i}',
                'control_options':choices,
                'rail_source_trait_taxa':{f'Genus taxon{chr(97+x)}' for x in range(9)}})
        return out
    def test_disjoint_four_station_full_max_and_study_diversity(self):
        source=self.make(4,False)
        refs,qa=solve_full_source_global_exclusivity(source,min_stretches=3)
        self.assertEqual(qa['status'],'EXACT_SOURCE_BINARY_MILP_OPTIMAL')
        self.assertEqual(qa['max_exact_unique_reference_original_rail_stations'],4)
        self.assertEqual(qa['source_total_distinct_reference_plot_IDs'],20)
        self.assertEqual(len(refs),4)
        self.assertTrue(all(len({p['study'] for p in chosen})>=2 for chosen in refs.values()))
    def test_shared_six_plot_pool_cannot_serve_three_original_stations(self):
        source=self.make(4,True)
        refs,qa=solve_full_source_global_exclusivity(source,min_stretches=1)
        self.assertEqual(qa['status'],'EXACT_SOURCE_BINARY_MILP_OPTIMAL')
        self.assertEqual(qa['max_exact_unique_reference_original_rail_stations'],1)
        self.assertEqual(qa['source_total_distinct_reference_plot_IDs'],5)
    def test_three_original_stretches_constraint_can_fail_closed(self):
        source=self.make(4,True)
        refs,qa=solve_full_source_global_exclusivity(source,min_stretches=3)
        self.assertIsNone(refs)
        self.assertEqual(qa['status'],'EXACT_SOLVER_NOT_OPTIMAL_NO_MAX_SUPPORT_CLAIM')
    def test_site_beta_not_computed_if_lt10_true_disjoint_rail_stations(self):
        source=self.make(4,False)
        refs,qa=solve_full_source_global_exclusivity(source,min_stretches=3)
        a=evaluate_matched(source,refs,{f'Genus taxon{chr(97+x)}':(x*.1,x*.2) for x in range(15)},reps=5)
        self.assertFalse(a['identified'])
    def test_source_study_is_not_automatically_verified_by_plot_ID(self):
        source=self.make(4,False)
        for s in source:
            for p in s['control_options']:p['study']='single_paper_only'
        refs,qa=solve_full_source_global_exclusivity(source,min_stretches=3)
        self.assertIsNone(refs)
        self.assertEqual(qa['status'],'EXACT_SOLVER_NOT_OPTIMAL_NO_MAX_SUPPORT_CLAIM')
if __name__=='__main__':unittest.main()
