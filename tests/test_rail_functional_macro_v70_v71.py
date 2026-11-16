"""All fixtures synthetic; successful tests never imply new biological results."""
from __future__ import annotations
import sys, copy, math, unittest
from pathlib import Path
from itertools import combinations
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from harvest_gift_growth_form_v70 import (
    BLOCKS_SIX,PERIODS,GIFT_sources,GBIF_exact_names,analyze as source_analyze,
    growth_group,make_source_keyset,species_binomial
)
from analyze_sixcountry_functional_rail_v71 import (
    run,validate_traits,make_pair_rows,MIN_TRAIT_COUNTRY_SOURCE_COVERAGE
)

class FakeResponse:
    def __init__(self,content):self.content=content
    def raise_for_status(self):pass
    def json(self):return self.content

class FakeClient:
    def get(self,url,params=None,timeout=None):
        q=(params or {}).get('query')
        if q=='traits_meta':return FakeResponse([{'Lvl3':'1.2.1','Trait2':'Growth_form','count':4}])
        if q=='species':
            return FakeResponse([{'work_ID':1,'work_species':'Poa annua'},
                {'work_ID':2,'work_species':'Solidago altissima'},
                {'work_ID':3,'work_species':'Pinus sylvestris'}] if params.get('startat')==0 else [])
        if q=='traits':return FakeResponse([{'work_ID':1,'trait_value':'herb','agreement':.99},
            {'work_ID':2,'trait_value':'herbaceous','agreement':.71},
            {'work_ID':3,'trait_value':'tree','agreement':.98},
            {'work_ID':4,'trait_value':'climber','agreement':.99}])
        if url.endswith('/1'):return FakeResponse({'rank':'SPECIES','kingdom':'Plantae','status':'ACCEPTED','acceptedKey':1,'canonicalName':'Poa annua'})
        if url.endswith('/2'):return FakeResponse({'rank':'SPECIES','kingdom':'Plantae','status':'ACCEPTED','acceptedKey':2,'canonicalName':'Solidago altissima'})
        if url.endswith('/3'):return FakeResponse({'rank':'SPECIES','kingdom':'Plantae','status':'ACCEPTED','acceptedKey':3,'canonicalName':'Pinus sylvestris'})
        raise AssertionError((url,params))


def fixture():
    source={'complete':True,'regions':{}}
    ra={'complete':True,'regions':{}}
    ro={'complete':True,'regions':{}}
    for k,block in enumerate(BLOCKS_SIX):
        cc=['US','FR','CA','GB','DE','AU'][k]
        cells={}
        for i in range(16):
            # four candidates in each growth class, three present in each cell;
            # missing species stay unobserved, never claimed truly absent.
            ids=[str(v) for v in (1,2,3,4) if v!=i%4+1]+[
                str(v) for v in (5,6,7,8) if v!=5+((i+k)%4)]
            periods={}
            for p in PERIODS:
                counts={z:7+((int(z)*i+k)%9) for z in ids}
                periods[p]={
                    'full_facet_source_pagination_completed':True,
                    'alien_observed_minimum_5_species_and_20_records':True,
                    'observed_alien_specieskey_counts':counts,
                    'n_indexed_GRIIS_introduced_plant_records':sum(counts.values()),
                    'indexed_all_plantae_source_records':1000+17*i+i*i,
                    'prior_2026_GBIF_indexed_source_record_count':1000+17*i+i*i}
            cells[f'{block}_x{i//4}_y{i%4}']={
                'source_pre_outcome_both_period_recording_eligible':True,'periods':periods}
        source['regions'][block]={'country':cc,'cells':cells}
        route_recs=[];road_recs=[]
        for a,b in combinations(sorted(cells),2):
            ia=int(a.split('_x')[-1].split('_')[0])*4+int(a.split('_y')[-1])
            ib=int(b.split('_x')[-1].split('_')[0])*4+int(b.split('_y')[-1])
            geo=20+abs(ia-ib)*7 + abs(ia%4-ib%4)*3
            road=geo*(1.12+.015*((ia+ib)%5))
            rail=geo*(1.25+.09*((ia*7+ib*3+k)%9))
            route_recs.append({'a':a,'b':b,'rail_to_geodesic_ratio':rail/geo,
                'cell_centroid_snap_augmented_rail_path_km':rail,
                'source_cell_centroid_geodesic_km':geo})
            road_recs.append({'a':a,'b':b,'road_to_geodesic_ratio':road/geo,
                'centroid_snap_augmented_road_path_km':road})
        ra['regions'][block]={'country':cc,'routes':{s:{'pair_routes':route_recs} for s in ('nearest_in_cell','nearest_largest_component_in_cell')}}
        ro['regions'][block]={'country':cc,'routes':{s:{'pair_paths':road_recs} for s in ('nearest_in_cell','nearest_largest_component_in_cell')}}
    traits={'complete':True,'source':{'GIFT_meta_ID':'1.2.1','GIFT_biasref':False,'GIFT_biasderiv':False,'GIFT_categorical_min_agreement':.66},
         'growth_form_by_source_gbif_species_key':{str(i):'herbaceous' if i<=4 else 'woody' for i in range(1,9)},
         'by_block':{b:{'source_GIFT_coverage_species_fraction':1.0} for b in BLOCKS_SIX}}
    return source,ra,ro,traits

class RealSourceContracts(unittest.TestCase):
    def test_categorical_ambiguity_and_taxonomy(self):
        self.assertEqual(growth_group('herb'), 'herbaceous')
        self.assertEqual(growth_group('Tree'),'woody')
        self.assertIsNone(growth_group('climber'))
        self.assertIsNone(growth_group('some wind resistance'))
        self.assertEqual(species_binomial('Poa annua'),'poa annua')
        self.assertIsNone(species_binomial('Poa annua var. minor'))
    def test_real_API_shape_mock_not_real_biology(self):
        c=FakeClient()
        groups,meta=GIFT_sources(c)
        self.assertEqual(groups['poa annua'],'herbaceous')
        self.assertEqual(groups['pinus sylvestris'],'woody')
        self.assertEqual(meta['GIFT_meta_ID'],'1.2.1')
        names,errors=GBIF_exact_names(c,['1','2','3'])
        self.assertEqual(set(names),{'1','2','3'})
    def test_frozen_source_keys_and_no_false_species_absence(self):
        f,_,_,_=fixture()
        keys,per=make_source_keyset(f)
        self.assertEqual(keys,[str(i) for i in range(1,9)])
        self.assertEqual(len(per),6)
        f['complete']=False
        with self.assertRaises(ValueError): make_source_keyset(f)

class FunctionalMacro(unittest.TestCase):
    def test_full_country_pair_group_source_study_is_descriptive(self):
        f,ra,ro,t=fixture()
        r=run(f,ra,ro,t)
        self.assertTrue(r['complete'])
        self.assertEqual(len(r['regions']),6)
        self.assertTrue(r['NO_plant_local_adaptation_or_plasticity_established'])
        self.assertTrue(r['NO_train_mediated_seed_disperal_mechanism_established'])
        z=r['regions']['us_midwest']['periods']['2018_2025']['nearest_in_cell']
        self.assertEqual(z['n_all_original_source_plant_railroad_pairs'],120)
        self.assertEqual(z['same_cell_pair_cohort_trait_sensitivity']['n_same_cell_pair_rows'],120)
        self.assertEqual(z['growth_group_eligible_cell_counts']['herbaceous'],16)
        self.assertEqual(z['growth_group_eligible_cell_counts']['woody'],16)
    def test_no_unmeasured_woody_claim_when_trait_coverage_low(self):
        f,ra,ro,t=fixture();t['by_block']['us_midwest']['source_GIFT_coverage_species_fraction']=.12
        r=run(f,ra,ro,t)
        z=r['regions']['us_midwest']['periods']['2018_2025']['nearest_in_cell']
        self.assertIsNone(z['same_cell_pair_cohort_trait_sensitivity']['descriptive_delta_rho_herbaceous_minus_woody'])
        self.assertIsNone(z['trait_group_rarefied20_source_similarity']['herbaceous']['partial_rho'])
    def test_refuse_unvetted_traits_and_mismatched_rail_road(self):
        f,ra,ro,t=fixture();t['source']['GIFT_biasref']=True
        with self.assertRaises(ValueError): run(f,ra,ro,t)
        f,ra,ro,t=fixture();ro['regions'].pop('fr_burgundy')
        with self.assertRaises(ValueError):run(f,ra,ro,t)
    def test_no_same_pairs_no_fake_comparison(self):
        f,ra,ro,t=fixture()
        for k in t['growth_form_by_source_gbif_species_key']:
            if int(k)>4:t['growth_form_by_source_gbif_species_key'][k]='herbaceous'
        r=run(f,ra,ro,t)
        z=r['regions']['de_saxony']['periods']['2018_2025']['nearest_in_cell']
        self.assertEqual(z['same_cell_pair_cohort_trait_sensitivity']['n_same_cell_pair_rows'],0)
        self.assertIsNone(z['same_cell_pair_cohort_trait_sensitivity']['descriptive_delta_rho_herbaceous_minus_woody'])

if __name__=='__main__':unittest.main()