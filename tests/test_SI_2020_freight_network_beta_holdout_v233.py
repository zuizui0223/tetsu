import math,sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent.parent/'scripts'))
from validate_SI_actual_2020_freight_vegetation_beta_holdout_v233 import (
    exact_within_line_labels,build_segment_dyads,observed_adjustment,
    leave_one_segment_cv,design
)
class TestSourcePhysicalSegmentCV(unittest.TestCase):
    def make_meta(self):
        lines=['10','10','20','34','50','50','50','50','80']
        return {i:{'national_line':lines[i], 'freight_2020':(i+1)*100,
                   'passenger_2020':1200+137*i*i} for i in range(9)}
    def test_all_freight_label_permutations_are_within_true_national_lines(self):
        meta=self.make_meta()
        seq=list(exact_within_line_labels(meta))
        self.assertEqual(seq[0][0],48)
        self.assertEqual(len(seq)-1,48)
        self.assertTrue(all(len(x)==9 for _,x in seq[1:]))
        self.assertEqual(len({tuple(sorted(x.items())) for _,x in seq[1:]}),48)
        self.assertEqual({x[2] for _,x in seq[1:]},{math.log1p(300)})
    def synthetic_sites(self):
        guild={str(i):(i%2,i%2) for i in range(13)}
        out=[]
        lines=['10','10','20','34','50','50','50','50','80']
        for i in range(9):
            out.append({'source_segment':f'physical_{i}',
                'section':lines[i], 'freight':(i+1)*300,
                'passenger':5000+i*275,'year':2022+(i%2),
                'lat':45.75+i*.055,'lon':14.25+i*.11,
                'source_species_all':{str(j) for j in range(i%5,i%5+7)},
                'known_original_lifespan_diaspore_species':{str(j) for j in range(i%5,i%5+7)}})
        return out,guild
    def test_nine_true_segments_generate_36_dependent_physical_pairs(self):
        source,guild=self.synthetic_sites()
        rows,meta=build_segment_dyads(source,guild)
        self.assertEqual(len(rows),36)
        self.assertEqual(len(meta),9)
        self.assertEqual(sum(q['original_source_site_pair_n'] for q in rows),36)
    def test_held_out_segment_folds_well_formed(self):
        source,guild=self.synthetic_sites()
        rows,meta=build_segment_dyads(source,guild)
        result=observed_adjustment(rows,'taxonomic_Jaccard')
        self.assertEqual(result['status'],'DESCRIPTIVE_ONLY')
        cv=leave_one_segment_cv(rows,'taxonomic_Jaccard')
        self.assertEqual(cv['source_physical_segment_holdout_folds'],9)
        self.assertEqual(cv['source_test_pairs_per_fold'],8)
        self.assertTrue(math.isfinite(cv['gain_baseline_minus_freight_MSE']))
    def test_proper_design_includes_passenger_and_source_line_not_only_freight(self):
        rows,_=build_segment_dyads(*self.synthetic_sites())
        self.assertEqual(design(rows,False).shape,(36,5))
        self.assertEqual(design(rows,True).shape,(36,6))
if __name__=='__main__':unittest.main()
