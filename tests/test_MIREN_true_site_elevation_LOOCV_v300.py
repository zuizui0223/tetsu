"""v30.0 synthetic mathematics tests only, NEVER substitute for field data."""
from __future__ import annotations
import copy,sys,unittest
from pathlib import Path
from statistics import mean
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import evaluate_MIREN_native_elevation_heldout_stature_predictability_v300 as m

def fake_roads(n=9,road_n=2):
    names=[f"Genus{chr(65+i//3)} species{chr(97+i%3)}" for i in range(24)]
    T=np.asarray([[i*.15,(i*13)%25*.1] for i in range(24)])
    roadlist=[]
    for rid in range(road_n):
        a=np.zeros((n,24))
        b=np.zeros((n,24))
        for i in range(n):
            for j in range(10):
                a[i,(i*2+j+rid)%24]=1
                b[i,(i*3+j+rid)%24]=1
        a=a/a.sum(axis=1)[:,None]
        b=b/b.sum(axis=1)[:,None]
        roadlist.append({"elev_m":np.linspace(450+rid*50,2500+rid*50,n),
                         "source_original_years":[2016]*n,
                         "matrices":{"roadside":a,"interior":b},
                         "n_phylo_paired_original_sites_in_road":n})
    return roadlist,T

class TrueOriginalSiteHeldoutPredictability(unittest.TestCase):
    def test_analytic_loo_matches_explicit_n_fold_training(self):
        n=8
        e=np.asarray([300,490,760,900,1130,1490,1790,2100],dtype=float)
        x=e/1000.
        y=4.3-1.7*x+np.asarray([.04,-.02,.06,-.03,.01,.04,-.01,.01])
        W=np.eye(n)
        T=np.stack((y,np.arange(n,dtype=float)),axis=1)
        result=m.true_site_loo_cv_score(W,e,T)
        meanerr=[];linearerr=[]
        for hold in range(n):
            train=np.arange(n)!=hold
            matrix=np.column_stack((np.ones(train.sum()),x[train]))
            coef=np.linalg.lstsq(matrix,y[train],rcond=None)[0]
            pred=coef[0]+coef[1]*x[hold]
            linearerr.append((y[hold]-pred)**2)
            meanerr.append((y[hold]-y[train].mean())**2)
        reference=np.sqrt(np.mean(meanerr))-np.sqrt(np.mean(linearerr))
        self.assertAlmostEqual(result[m.METRIC],reference,places=12)
        self.assertGreater(result[m.METRIC],0)

    def test_no_self_training_and_linear_slope_in_original_site_units(self):
        e=np.linspace(100,2100,9)
        W=np.eye(9)
        y=1.3+2*(e/1000.)
        T=np.stack((y,np.linspace(0,1,9)),axis=1)
        q=m.true_site_loo_cv_score(W,e,T)
        self.assertAlmostEqual(q["held_out_elevation_RMSE"],0,places=12)
        self.assertAlmostEqual(q["original_road_slope_signed_per_km"],2.,places=12)
        self.assertTrue(q["all_original_sites_held_out_exactly_once"])

    def test_zero_elevational_data_or_too_few_sites_refused(self):
        rows,T=fake_roads(5)
        with self.assertRaises(ValueError):
            m.true_site_loo_cv_score(rows[0]["matrices"]["roadside"],
                                     rows[0]["elev_m"],T)
        rows,T=fake_roads(8)
        with self.assertRaises(ValueError):
            m.true_site_loo_cv_score(rows[0]["matrices"]["roadside"],
                                     np.ones(8)*500.,T)

    def test_unnatural_source_leverage_heldout_training_fails_closed(self):
        rows,T=fake_roads(6)
        x=np.array([200.,200.,200.,200.,200.,1800.])
        with self.assertRaises(ValueError):
            m.true_site_loo_cv_score(rows[0]["matrices"]["roadside"],x,T)

    def test_two_paired_habitats_and_road_unit_weighting(self):
        roads,T=fake_roads()
        q=m.original_region_road_equal_heldout_scores(roads,T)
        self.assertEqual(q["n_true_original_roads"],2)
        self.assertAlmostEqual(q["gain"]["roadside_minus_interior"],
            q["gain"]["roadside"]-q["gain"]["interior"],places=12)
        for z in ("roadside","interior"):
            self.assertAlmostEqual(q["mean_intercept_LOO_RMSE_by_zone"][z]-
                                  q["mean_Elevation_LOO_RMSE_by_zone"][z],
                                  q["gain"][z],places=12)

    def test_new_predictive_gate_is_predeclared_and_cannot_change_traits(self):
        roads,T=fake_roads()
        original=copy.deepcopy(roads)
        S,qa=m.original_environment_roads_with_cv_support(
            roads,{"source_true_field_sites_on_usable_elevation_roads":18})
        self.assertTrue(qa["original_region_out_of_sample_heldout_elevation_support"])
        self.assertEqual(qa["original_roads_with_ge6_true_sites_and_identifiable_heldout_elev"],2)
        for a,b in zip(roads,original):
            np.testing.assert_array_equal(a["matrices"]["roadside"],b["matrices"]["roadside"])
            np.testing.assert_array_equal(a["elev_m"],b["elev_m"])

    def test_original_source_retention_failure_reports_not_identifiable(self):
        roads,T=fake_roads()
        s,qa=m.original_environment_roads_with_cv_support(
            roads,{"source_true_field_sites_on_usable_elevation_roads":100})
        self.assertFalse(qa["original_region_out_of_sample_heldout_elevation_support"])

    def test_repeated_trait_label_technicals_do_not_invent_field_replication(self):
        roads,T=fake_roads()
        groups=[np.arange(i,i+4) for i in range(0,len(T),4)]
        s=m.source_trait_relabel_predictability_null(roads,T,groups,0,0,
               "SYNTHETIC_NOT_FIELD",draws=m.DRAWS_PER_CHAIN)
        self.assertEqual(s["original_source_original_site_heldout_predictability"]
                         ["roadside"]["n_technical_label_draws_NOT_field_studies"],398)
        self.assertTrue(s["true_field_heldout_sites_NEVER_used_in_own_road_elevation_fit"])
        self.assertAlmostEqual(s["original_source_original_site_heldout_predictability"]
                         ["roadside_minus_interior"][
                           "original_site_LOOCV_observed_minus_source_conditional_null"],
                         s["original_source_original_site_heldout_predictability"]
                         ["roadside"]["original_site_LOOCV_observed_minus_source_conditional_null"]-
                         s["original_source_original_site_heldout_predictability"]
                         ["interior"]["original_site_LOOCV_observed_minus_source_conditional_null"],places=12)

    def test_original_region_summary_does_not_pool_site_pairs_as_replicates(self):
        roads,T=fake_roads()
        groups=[np.arange(i,i+4) for i in range(0,len(T),4)]
        fake=m.source_trait_relabel_predictability_null(
            roads,T,groups,0,0,"SYNTHETIC_NOT_FIELD",draws=m.DRAWS_PER_CHAIN)
        self.assertIsNone(m.equal_region_summary({"R1":fake,"R2":fake})[
                          "source_region_equal_effects"])
        data=m.equal_region_summary({"R1":fake,"R2":fake,"R3":fake})
        self.assertEqual(data["source_supported_original_region_codes_n"],3)
        self.assertEqual(data["source_region_equal_effects"]["roadside"][
            "source_original_region_codes_evaluated"],3)

    def test_source_year_road_filter_does_not_shift_true_site_elevation(self):
        roads,T=fake_roads()
        roads[0]["source_original_years"][1]=2017
        from test_MIREN_source_elevation_phylo_occ_stature_v299 import original_year_consistent_roads
        valid=original_year_consistent_roads(roads)
        self.assertEqual(len(valid),1)
        np.testing.assert_array_equal(valid[0]["elev_m"],roads[1]["elev_m"])

if __name__=="__main__":unittest.main()
