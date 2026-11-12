"""Synthetic botanical field survey, NOT actual train traffic. Check
MIREN strict original ROAD plot first/last pair alien taxon set gains,
losses, union-standardized gains, no imputed colonization.
"""
import json,sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import pandas as pd
from analyze_MIREN_standardized_road_edge_vegetation_v99 import native_status_safe
from analyze_MIREN_roadside_alien_resurvey_turnover_v100 import (
    STRICT_PLOTS,MIN_REGION_REPEAT_SITES,per_site_turnover,summarize_change)


def fixture():
    rows=[]
    def row(site,plot,year,status,name):
        rows.append({
            "Plot_id":f"TEST.RD.{site:02d}.{plot}",
            "Region":"TEST","Road":"RD","Site":site,
            "Plot":plot,"Year":year,"Status":status,
            "Accepted.Name.MIREN":name,
            "Abundance":1,"Cover":2,
        })
    for site in range(1,12):
        # Actual 2017+2022 same three original plots; Plot1 and 3
        # have deliberately different alien taxa sets.
        for yr in (2017,2022):
            road_alien=("a","b") if yr==2017 else ("a","c","d")
            interior_alien=("p",) if yr==2017 else ("p","q")
            for sp in road_alien:row(site,"1",yr,"Alien",sp)
            for sp in interior_alien:row(site,"3",yr,"Alien",sp)
            row(site,"1",yr,"Native","native1")
            row(site,"3",yr,"Native","native3")
        # Explicit replacement (1A) should not be treated as stable plot1.
        row(site,"1A",2022,"Alien","replacement_only")
    rows.append(rows[0].copy()) # duplicated taxon row
    return pd.DataFrame(rows)


class RealPlotTurnoverLogic(unittest.TestCase):
    def test_strict_original_not_replacement_alien_sets(self):
        clean,q=native_status_safe(fixture())
        self.assertEqual(STRICT_PLOTS,("1","3"))
        self.assertEqual(MIN_REGION_REPEAT_SITES,10)
        self.assertEqual(q["n_original_nonstandard_A_replacement_plot_rows_excluded"],11)
        self.assertEqual(q["n_removed_repeated_species_same_plot_year"],1)
        x=per_site_turnover(clean)
        self.assertEqual(len(x),11)
        self.assertEqual(x["gain1"].tolist(),[2]*11)
        self.assertEqual(x["gain3"].tolist(),[1]*11)
        self.assertEqual(x["loss1"].tolist(),[1]*11)
        self.assertEqual(x["loss3"].tolist(),[0]*11)
        self.assertTrue((x["years_between_first_and_last"]==5).all())
    def test_relative_gain_is_union_based_and_not_direct_colonization(self):
        clean,q=native_status_safe(fixture())
        x=per_site_turnover(clean)
        self.assertAlmostEqual(x["gain_fraction_union1"].iloc[0],.5)
        self.assertAlmostEqual(x["gain_fraction_union3"].iloc[0],.5)
        self.assertAlmostEqual(x[
            "difference_in_gain_fraction_roadside_minus_interior"].iloc[0],0)
        self.assertEqual(x["difference_in_earliest_vs_latest_roadside_alien_advantage"].iloc[0],0)
        z=summarize_change(x)
        self.assertEqual(z["n_original_matched_first_to_last_road_sites"],11)
        self.assertEqual(z["region_code_equal_signed_mean_absolute_gain_source"][
            "n_region_codes_mean_positive"],1)
        self.assertEqual(z["region_code_equal_signed_mean_relative_gain_source"][
            "n_region_codes_mean_exact_zero"],1)
        self.assertTrue(z["no_statepair_pseudoreplication_or_local_seed_colonization_rate"])
    def test_no_plants_in_both_original_survey_plots_is_not_imputed_fake_gain_rate(self):
        d=fixture()
        # Interior original 2017/2022 all alien species removed for site 1;
        # native rows preserve plot recorded source, but denominator UNION=0.
        mask=((d.Site==1)&(d.Plot=="3")&(d.Status=="Alien"))
        d=d[~mask]
        clean,q=native_status_safe(d)
        p=per_site_turnover(clean)
        i=p[p.Site==1].iloc[0]
        self.assertTrue(pd.isna(i["gain_fraction_union3"]))
        self.assertTrue(pd.isna(i["difference_in_gain_fraction_roadside_minus_interior"]))
        self.assertEqual(i["gain3"],0)
    def test_requires_true_2year_revisits_and_matched_original_1_and_3(self):
        d=fixture()
        # Replace last original interior plot with replacement 3A:
        # region/road/site 1 lacks a 2022 original plot3.
        d.loc[(d.Site==1)&(d.Plot=="3")&(d.Year==2022),"Plot"]="3A"
        clean,q=native_status_safe(d)
        a=per_site_turnover(clean)
        self.assertEqual(len(a),10)
        self.assertNotIn(1,a.Site.to_list())
    def test_biological_change_conservation_gain_minus_loss(self):
        clean,q=native_status_safe(fixture())
        a=per_site_turnover(clean)
        actual=(a["late1"]-a["late3"])-(a["early1"]-a["early3"])
        expected=a["gain1"]-a["loss1"]-a["gain3"]+a["loss3"]
        self.assertTrue(actual.equals(expected))
    def test_no_genuine_followup_site_fails_closed(self):
        d=fixture()
        d=d[d.Year==2017]
        clean,q=native_status_safe(d)
        t=per_site_turnover(clean)
        self.assertTrue(t.empty)
        with self.assertRaisesRegex(ValueError,"No genuine"):
            summarize_change(t)


if __name__=="__main__":
    unittest.main()
