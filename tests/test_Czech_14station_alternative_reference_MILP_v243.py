import math,sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent.parent/"scripts"))
from validate_Czech_rail_local_reference_alternative_14station_assignments_v243 import (
    solve_alternate_maximum,station_reference_profiles,stratum_loo,evaluate
)
class TestDifferentFeasibleOriginalReferences(unittest.TestCase):
    def make(self,n=16,shared=False):
        stations=[];originalgeo={}
        species=[f"Genus specimen{chr(97+k//26)}{chr(97+k%26)}" for k in range(130)]
        traits={x:(math.sin(k/10),math.cos(k/11)) for k,x in enumerate(species)}
        for i in range(n):
            reference=[]
            for j in range(9):
                id=f"shared_{j}" if shared else f"study{i}_plot{j}"
                reference.append({'id':id,'study':('Paper_A' if j<5 else 'Paper_B'),
                    'd':0.2+j*0.33,
                    'taxa':set(species[(i*9+j*4)%90:(i*9+j*4)%90+16])})
                originalgeo[id]={'lat':49+i*.025+j*.0002,'lon':13+i*.035+j*.0002}
            stations.append({'id':f"station_{i}",'stretch':f"original_rail_{i%4}",
              'rail_source_trait_taxa':set(species[(i*4)%90:(i*4)%90+35]),
              'control_options':reference})
        return stations,traits,originalgeo
    def test_exact_fourteen_stations_and_no_reused_plot(self):
        src,val,geo=self.make()
        selection,qa=solve_alternate_maximum(src,seed=45)
        self.assertIsNotNone(selection)
        self.assertEqual(qa['original_source_stations'],14)
        self.assertEqual(qa['true_source_distinct_releve_ids'],70)
        self.assertEqual(len({src[i]['stretch'] for i in selection}),4)
        allids=[p['id'] for arr in selection.values() for p in arr]
        self.assertEqual(len(allids),len(set(allids)))
        self.assertTrue(all(len({p['study'] for p in arr})>=2 for arr in selection.values()))
    def test_alternative_same_max_but_can_have_other_original_plots(self):
        src,val,geo=self.make()
        first,a=solve_alternate_maximum(src,seed=46)
        second,b=solve_alternate_maximum(src,seed=87)
        self.assertEqual(len(first),14)
        self.assertEqual(len(second),14)
        self.assertTrue(a['botany_blind_allocation_sha256'])
        self.assertTrue(b['botany_blind_allocation_sha256'])
    def test_same_species_rarefaction_and_leave_original_rail_stretch(self):
        src,val,geo=self.make()
        selection,qa=solve_alternate_maximum(src,seed=45)
        rail,ref,metrics,sourceQA=station_reference_profiles(src,selection,geo,val,rarefy_seed=22)
        self.assertEqual(metrics['n_real_rail_stations'],14)
        self.assertEqual(metrics['source_n_dependent_station_pairs'],91)
        self.assertEqual(sourceQA['n_source_original_plot_ids_same_exact_4decimal_coordinate_cell'],0)
        leave=stratum_loo(rail,ref,src,selection)
        self.assertEqual(len(leave),4)
        self.assertTrue(all(math.isfinite(x) for x in leave.values()))
    def test_infeasible_shared_nine_plots_must_not_manufacture_14_stations(self):
        src,val,geo=self.make(shared=True)
        selection,qa=solve_alternate_maximum(src,seed=45,time_limit_s=10)
        self.assertIsNone(selection)
        self.assertEqual(qa['status'],'NO_OPTIMAL_ALTERNATIVE_IDENTIFIED')
if __name__=="__main__":unittest.main()
