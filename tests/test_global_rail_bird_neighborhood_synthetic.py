"""Algebraic checks on the frozen multi-source HWI estimator.

Synthetic fixture ONLY. No real occurrence data, no source access, no tuning.
"""
import json
import unittest
import numpy as np
import analyze_global_rail_bird_neighborhood_source as m

def fixture(active_interaction=.12,inactive_interaction=-.04):
    rng=np.random.default_rng(20261128)
    # All species have fixed HWI and habitat classes. Unequal case support
    # by target is intentional. Each target has four documented neighbors.
    n_species=80
    sf=rng.uniform(-.025,.025,n_species)
    n=len(m.NAMES);d=len(m.TARGET_FIELDS)
    beta=np.zeros(len(m.FIELDS))
    beta[0]=.025;beta[1]=-.015
    beta[n]=active_interaction;beta[n+1]=inactive_interaction
    beta[4:8]=[.006,-.003,.007,.002]
    beta[n+6]=-.004
    beta[2*n:2*n+d]=rng.uniform(-.004,.004,d)
    rows=[]
    for j in range(120):
        target=(j,0)
        block=("synthetic",j//5)
        target_cov=rng.normal(0,.4,d)
        # Crossings 0/1; other covariates normalized for numerical testing.
        V=rng.normal(0,.4,(4,n))
        V[:,:4]=rng.integers(0,2,(4,4))
        # Each target has both active and inactive exposure variability.
        V[:,:2]=[[0,0],[1,0],[0,1],[1,1]]
        for hab in ("h1_m1","h2_m1"):
            offset=rng.uniform(-.035,.035)
            ids=[i for i in range(n_species) if (i//20)%2==(hab=="h2_m1")]
            for i in ids:
                low=int(i>=40)
                # Keep independent species/target FE crossing, while
                # creating incomplete & imbalanced source opportunities.
                if rng.random()<.15:continue
                q=V[rng.integers(0,4)]-V.mean(axis=0)
                X=np.r_[q,low*q,low*target_cov]
                y=float(.45+sf[i]+offset+X@beta)
                assert 0<y<1
                rows.append({"tid":i,"low":low,"target":target,
                    "stratum":(target,hab,1),"hab":hab,"block":block,
                    "x":X,"y":y})
    return rows

def check_effects(test,rr,active,inactive):
    v=m.fit(rr)
    test.assertEqual(v["status"],"estimated",v)
    e=v["effects"]
    for name,expected in {
        "high_HWI_active":.025,
        "high_HWI_inactive":-.015,
        "low_HWI_active":.025+active,
        "low_HWI_inactive":-.015+inactive,
        "low_minus_high_active":active,
        "low_minus_high_inactive":inactive,
        "low_minus_high_active_minus_inactive":active-inactive,
    }.items():
        test.assertAlmostEqual(e[name]["estimate"],expected,places=7,msg=name)
    print(json.dumps({"synthetic_only":True,"injected_active":active,
        "injected_inactive":inactive,"recovered_difference":e["low_minus_high_active_minus_inactive"]["estimate"],
        "n_cases":v["n_cases"],"n_targets":v["n_targets"],"n_blocks":v["n_blocks"]}),flush=True)

class TestFrozen(unittest.TestCase):
    def test_known_positive(self):
        check_effects(self,fixture(),.12,-.04)
    def test_no_interaction(self):
        check_effects(self,fixture(0.,0.),0.,0.)
    def test_negative_interaction(self):
        check_effects(self,fixture(-.08,.03),-.08,.03)
    def test_no_status_variation(self):
        rr=fixture()
        n=len(m.NAMES)
        for r in rr:r["x"][[0,n]]=0.
        self.assertEqual(m.fit(rr)["status"],"insufficient status source variation")
    def test_alias(self):
        rr=fixture()
        n=len(m.NAMES)
        for r in rr:
            r["x"][1]=r["x"][0]
            r["x"][n+1]=r["x"][n]
        self.assertEqual(m.fit(rr)["status"],"rank deficient")

if __name__=="__main__":unittest.main()
