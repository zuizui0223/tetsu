"""v329 synthetic tests of additivity of original observed taxa loss and gain."""
import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import analyze_Kashmir_EIVE_loss_gain_identity_v329 as m

NAMES=["poa annua","galium aparine","lactuca serriola","unknown sp"]
EIVE={"poa annua":{"position":7.0,"width":4.0,"systems":3},
      "galium aparine":{"position":5.0,"width":3.0,"systems":3},
      "lactuca serriola":{"position":9.0,"width":5.0,"systems":3}}
def fixture():
    d={}
    for i in range(1,17):
        s=f"Station {i}"
        for zone in m.ZONES:
            d[s,2014,zone]={0,1}
            d[s,2017,zone]={1,2}
    for i in range(1,16):
        s=f"Between {i}"
        for zone in m.ZONES:
            d[s,2014,zone]={0,1}
            d[s,2017,zone]={1,2}
    return d
class OriginalTaxonPartition(unittest.TestCase):
    def test_exact_loss_gain_mean_identity(self):
        out=m.components({"a":7.0,"b":5.0},{"b":5.0,"c":9.0})
        self.assertAlmostEqual(out["loss_component"],-1)
        self.assertAlmostEqual(out["gain_component"],2)
        self.assertAlmostEqual(out["total_cwm_change"],1)
        self.assertEqual((out["n_lost"],out["n_gained"],out["n_persistent"]),(1,1,1))
    def test_only_loss_warm_species_reduces_cwm(self):
        out=m.components({"a":7.0,"b":5.0},{"b":5.0})
        self.assertAlmostEqual(out["gain_component"],0)
        self.assertAlmostEqual(out["loss_component"],-1)
        self.assertAlmostEqual(out["total_cwm_change"],-1)
    def test_all_new_species_not_persistence_anchor(self):
        self.assertIsNone(m.components({"a":5},{"b":7}))
    def test_missing_source_temporal_plot_not_zero(self):
        d=fixture()
        d["Station 1",2017,"Track"]=None
        out,why=m.site_components(d,NAMES,EIVE,"Station 1","Track","position",1)
        self.assertIsNone(out)
        self.assertEqual(why,"original_source_missing_taxon_cell")
    def test_one_source_site_has_one_decomposition(self):
        d=fixture()
        zones,by=m.zonal(d,NAMES,EIVE,"position",1)
        self.assertEqual(zones["Track"]["by_stratum"]["all"]["total_cwm_change"]["n"],31)
        self.assertEqual(zones["Track"]["by_stratum"]["station"]["total_cwm_change"]["n"],16)
        self.assertEqual(zones["Track"]["by_stratum"]["between"]["total_cwm_change"]["n"],15)
        self.assertAlmostEqual(zones["Track"]["by_stratum"]["all"]["total_cwm_change"]["mean"],1)
    def test_same_source_site_track_and_interior_component_difference(self):
        d=fixture()
        # Make interior no change: taxon A and B persist.
        for s in m.cohort(d):
            d[s,2017,"Interior Plot 2"]={0,1}
        zones,details=m.zonal(d,NAMES,EIVE,"position",1)
        p=m.paired_zone_comparison(details,"Interior Plot 2")
        self.assertEqual(p["n_same_source_sites"],31)
        self.assertAlmostEqual(p["all"]["track_minus_interior_2017_minus_2014_cwm"]["mean"],1)
    def test_sparse_2017_species_gate_explicit(self):
        d=fixture()
        q,why=m.site_components(d,NAMES,EIVE,"Station 1","Track","position",3)
        self.assertIsNone(q)
        self.assertEqual(why,"minimum_EIVE_matched_species")
    def test_species_system_n_gates_not_treated_as_observed_absence(self):
        eive=dict(EIVE)
        eive["lactuca serriola"]={"position":9.0,"width":5.0,"systems":1}
        r,why=m.site_components(fixture(),NAMES,eive,"Station 1","Track","position",1)
        self.assertIsNotNone(r)
        self.assertEqual(r["n_gained"],0)
if __name__=="__main__":unittest.main()
