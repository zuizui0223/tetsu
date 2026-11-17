import math,random,sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent.parent/'scripts'))
from validate_SI_89_original_track_sites_Europe82_trait_axis_v246 import (
    source_european_trait_axis,source_slovenian_stratified_null
)
from validate_Europe_3studies_82stations_original_LEDA_2D_trait_beta_v245 import (
    get_metrics,regional_margin
)
class TestFourthIndependentRealSourceGrain(unittest.TestCase):
    def test_european_622_species_trait_axes_frozen_before_slovenia(self):
        qa={'status':'SOURCE_VERIFIED_CZ_IT_HR_BA_82_INDEPENDENT_RAIL_STATIONS_2D_LEDA_INTERNAL_TRAIT_NULL',
            'source_all_three_study_taxa_shared_trait_axis_QA':{
                'source_original_LEDA_both_canopy_height_seed_mass_species_n':622,
                'species_equal_global_log_mean':[-.58,.293],
                'species_equal_global_log_SD':[1.38,2.16]},
            'original_three_study_descriptive_sign_coherence':{
                'number_original_independent_publisher_studies':3,
                'original_station_number_across_three_studies':82}}
        x,info=source_european_trait_axis(qa,
              {'Campanula punctata':10}, {'Campanula punctata':4})
        self.assertEqual(len(x),1)
        self.assertAlmostEqual(x['Campanula punctata'][0],
                 (math.log(10)+.58)/1.38)
        self.assertTrue(info['no_slovenian_species_used_to_refit_reference_trait_axis'])
        qa['source_all_three_study_taxa_shared_trait_axis_QA']['source_original_LEDA_both_canopy_height_seed_mass_species_n']=621
        with self.assertRaises(ValueError):
            source_european_trait_axis(qa,{'Campanula punctata':10},
                                             {'Campanula punctata':4})
    def test_original_route_stratified_margin_with_single_source_location(self):
        rng=random.Random(20261117)
        taxon=[f'genus spec{i}' for i in range(65)]
        vals={t:(rng.gauss(0,1),rng.gauss(0,1)) for t in taxon}
        cohort=[]
        for route in range(8):
            for station in range(2):
                cohort.append({'station':f'{route}_{station}',
                    'group':str(route),'taxa':set(rng.sample(taxon,20))})
        cohort.append({'station':'lonely','group':'8','taxa':set(rng.sample(taxon,20))})
        before=regional_margin(cohort)
        out,qa=source_slovenian_stratified_null(cohort,vals,
                                                 draws=3,burnin=6000,gap=1600)
        self.assertEqual(qa['source_n_route_sections_with_single_site_unchanged'],1)
        self.assertEqual(len(out),3)
        self.assertEqual(before,regional_margin(cohort))
        self.assertGreater(qa['total_accepted_source_swaps'],1000)
        self.assertGreater(get_metrics(cohort,vals)['functional_CWM_centroid_distance_all_stationpairs'],0)
    def test_six_route_groups_not_enough_to_support_real_section_null(self):
        rng=random.Random(4)
        species=[f'genus spec{i}' for i in range(30)]
        traits={s:(rng.random(),rng.random()) for s in species}
        c=[{'station':f'station{i}',
            'group':str(i//2),'taxa':set(rng.sample(species,10))}
              for i in range(12)]
        with self.assertRaises(ValueError):
            source_slovenian_stratified_null(c,traits,draws=1,burnin=1,gap=1)
if __name__=='__main__':unittest.main()
