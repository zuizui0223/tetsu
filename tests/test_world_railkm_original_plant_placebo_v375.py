"""v375 synthetic test battery for noncausal climate-matched rail map placement negatives."""
import copy,sys,unittest,math,random
from collections import Counter,defaultdict
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import validate_world_railkm_original_plant_placebo_v375 as m

def item(cont,tc,lat,x,i):
    return {"continent":cont,"bio1":tc,"abs_lat":lat,
            "key":(i,3),"rail_km":x,
            "log_rail_km":math.log1p(x),
            "observations":150,"species":55,
            "block":(i//10,0)}

class OriginalRailPlacebo(unittest.TestCase):
    def test_same_continent_temperature_and_latitude_bin_fixed(self):
        source=[item("Asia",15,35,x,i) for i,x in enumerate([0,100,20,80])]
        source+= [item("Asia",25,35,x,i+10) for i,x in enumerate([0,1,2,3])]
        source+= [item("Europe",15,35,x,i+20) for i,x in enumerate([0,7,8,9])]
        orig=copy.deepcopy(source)
        qa=m.permute_railkm_in_eligible_observed_cells(source,random.Random(3))
        self.assertEqual(qa["original_source_climate_latitude_strata"],3)
        self.assertEqual(qa["strata_with_n_at_least_three_and_nontrivial_rail_distribution"],3)
        for key in {m.shuffle_stratum(v) for v in orig}:
            a=sorted(v["rail_km"] for v in orig if m.shuffle_stratum(v)==key)
            b=sorted(v["rail_km"] for v in source if m.shuffle_stratum(v)==key)
            self.assertEqual(a,b)
    def test_rail_source_zero_stays_in_stratum_distribution(self):
        rows=[item("Oceania",21,40,x,i) for i,x in enumerate([0,0,0,55,34])]
        m.permute_railkm_in_eligible_observed_cells(rows,random.Random(7))
        self.assertEqual(sum(r["rail_km"]==0 for r in rows),3)
    def test_source_species_richness_effort_unchanged(self):
        rows=[item("Asia",15,20,x,i) for i,x in enumerate([0,5,10,40])]
        before=[(r["species"],r["observations"],r["continent"],r["key"]) for r in rows]
        m.permute_railkm_in_eligible_observed_cells(rows,random.Random(5))
        after=[(r["species"],r["observations"],r["continent"],r["key"]) for r in rows]
        self.assertEqual(before,after)
    def test_rail_log_feature_updates_in_all_cell_rows(self):
        rows=[item("Africa",15,15,x,i) for i,x in enumerate([0,10,20,40])]
        m.permute_railkm_in_eligible_observed_cells(rows,random.Random(3))
        for row in rows:
            self.assertAlmostEqual(row["log_rail_km"],math.log1p(row["rail_km"]))
    def test_original_source_small_groups_unchanged(self):
        rows=[item("Asia",7,30,3,0),item("Asia",7,30,100,1)]
        qa=m.permute_railkm_in_eligible_observed_cells(rows,random.Random(2))
        self.assertEqual(qa["source_valid_botanical_cells_railkm_movable"],0)
        self.assertEqual(rows[0]["rail_km"],3)
    def test_invalid_source_negative_source_rail_length_rejected(self):
        rows=[item("Asia",25,20,0,0)]
        rows[0]["rail_km"]=-1
        with self.assertRaises(ValueError):
            m.permute_railkm_in_eligible_observed_cells(rows,random.Random(2))
    def test_placebo_performance_fraction_is_not_pvalue(self):
        out=m.tech_summary(.7,[.1,.2,.3,.4,.5,.6])
        self.assertAlmostEqual(out[
            "source_placement_fraction_fake_at_least_as_predictive_as_real_plus_one"],1/7)
        self.assertTrue(out["exploratory_technical_negative_control_not_CI_or_causal_randomization"])
    def test_climate_band_breakpoints_stable(self):
        self.assertEqual(m.shuffle_stratum(item("Europe",7.9,35,2,1)),("Europe","cool",3))
        self.assertEqual(m.shuffle_stratum(item("Europe",8.0,35,2,1)),("Europe","temperate",3))
        self.assertEqual(m.shuffle_stratum(item("Europe",20.0,35,2,1)),("Europe","warm",3))
    def test_source_v371_frozen_actual_global_lift_checksum(self):
        self.assertAlmostEqual(m.EXPECTED[100000]["richness"],.004463745492835951)
        self.assertAlmostEqual(m.EXPECTED[200000]["jaccard"],-.00006120015031440586)
    def test_predeclared_world_placebo_count(self):
        self.assertEqual(m.DRAWS,39)
        self.assertEqual(len(m.GRID_WIDTHS),2)
if __name__=="__main__":unittest.main()
