"""v8.9: synthetic mathematical source checks for state-specific origin
and destination fixed effects and Census-stratified simultaneous freight
node relabeling. No independently exchangeable randomized seed treatment.
"""
import sys,unittest
from pathlib import Path
import math
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from validate_US_FAF6_GRIIS_dyadic_FE_spatial_QAP_v89 import (
    MIN_STATES,MIN_ROWS,SEED,PERMUTATIONS,MODELS,PRIMARY_PAIR_CONTROLS,
    prepare_dyadic,permutation_map,predictor_vector,residual_rho,
    one_mode,RESPONSE)


def mock_source(n=44):
    states=[f"{i:02d}" for i in range(1,n+1)]
    geo={s:{"region":str((i//11)+1),"division":str((i//5)+1)}
         for i,s in enumerate(states)}
    od={}
    rows=[]
    for i,o in enumerate(states):
        for j,d in enumerate(states):
            if o==d:continue
            if (i*13+j)%17==0:od[(o,d)]=(i+j+1)*.5
            y=((i*19+j*11+(i*j)%17)%93)/1000.
            rows.append({
              "o":o,"d":d,
              "x_rail":math.log1p(od.get((o,d),0)),
              "x_truck":math.log1p((i*7+j*13)%31),
              "covars":{
                 "log_geographical_distance":math.log(2+abs(i-j)),
                 "same_Census_region":int(geo[o]["region"]==geo[d]["region"]),
                 "same_Census_division":int(geo[o]["division"]==geo[d]["division"]),
                 "log_origin_candidate_taxon_categories":math.log(15+(i*3+j*7)%49),
              },
              "responses":{
                 RESPONSE:y,
                 "raw_source_linked_first_record_fraction":y*.5+((i+j)%11)/1000.}
            })
    return rows,geo,od


class DyadicNetworkSourcelabelValidation(unittest.TestCase):
    def test_frozen_minimum_state_support_and_source_rng(self):
        self.assertEqual(MIN_STATES,40)
        self.assertEqual(MIN_ROWS,1000)
        self.assertEqual(SEED,20261112)
        self.assertEqual(PERMUTATIONS,{"Census_region_preserved":499,
            "Census_division_preserved":299,
            "unrestricted_state_codes_descriptive_only":199})
        self.assertEqual(MODELS,("RAIL","TRUCK"))
        self.assertIn("log_origin_candidate_taxon_categories",PRIMARY_PAIR_CONTROLS)
    def test_permutation_is_bijection_preserves_broad_region_and_original_source_matrix(self):
        states=[f"{i:02d}" for i in range(1,13)]
        blocks={s:str((int(s)-1)//4) for s in states}
        rng=np.random.default_rng(SEED)
        mapping=permutation_map(states,blocks,rng)
        self.assertEqual(set(mapping),set(states))
        self.assertEqual(set(mapping.values()),set(states))
        self.assertTrue(all(blocks[s]==blocks[mapping[s]] for s in states))
        self.assertTrue(any(mapping[s]!=s for s in states))
        mapping2=permutation_map(states,blocks,np.random.default_rng(SEED))
        self.assertEqual(mapping,mapping2)
    def test_same_source_OD_observed_and_permuted_direction_in_vector(self):
        a=[{"o":"01","d":"02"},{"o":"02","d":"01"}]
        od={("01","02"):100,("02","01"):2}
        x=predictor_vector(a,od,{"01":"01","02":"02"})
        r=predictor_vector(a,od,{"01":"02","02":"01"})
        self.assertAlmostEqual(x[0],math.log1p(100))
        self.assertAlmostEqual(r[0],math.log1p(2))
        self.assertEqual(sum(x),sum(r))
    def test_source_dummies_do_not_turn_1900_dyads_into_1900_independent_states(self):
        rows,geo,od=mock_source()
        self.assertGreater(len(rows),MIN_ROWS)
        ref=prepare_dyadic(rows,geo,other="TRUCK")
        self.assertEqual(len(ref["states"]),44)
        self.assertEqual(ref["n_observed_direction_pairs"],44*43)
        self.assertEqual(ref["Q"].shape[1],1+2*(44-1)+1+len(PRIMARY_PAIR_CONTROLS))
        x=predictor_vector(rows,od,{s:s for s in ref["states"]})
        r=residual_rho(x,ref["Y"][RESPONSE],ref["Q"])
        self.assertTrue(r is None or -1<=r<=1)
    def test_no_unphysical_synthetic_inference_with_too_few_states(self):
        rows,geo,od=mock_source(15)
        with self.assertRaisesRegex(ValueError,"Insufficient"):
            prepare_dyadic(rows,geo)
    def test_rho_residual_does_not_create_signal_after_perfect_projection(self):
        q=np.eye(5,2)
        x=[1,2,3,4,5]
        y=np.zeros(5)
        self.assertIsNone(residual_rho(x,y,q))
    def test_only_original_mode_names_accepted(self):
        rows,geo,od=mock_source()
        src=prepare_dyadic(rows,geo)
        with self.assertRaisesRegex(ValueError,"Only original"):
            one_mode("BIRDS",src,od)
if __name__=="__main__":unittest.main()
