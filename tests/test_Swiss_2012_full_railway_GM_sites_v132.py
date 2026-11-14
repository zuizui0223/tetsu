"""v13.2 source original Schoenenberger and D'Andrea 2012 rail CP4 marker site parser tests."""
from __future__ import annotations
import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from analyze_Swiss_2012_full_railway_GM_sites_v132 import (
    isolate_original_site_table,extract_original_2012_table,
    describe_surveyed_original_railway,
    source_homogeneous_plant_label_randomization,
)


def test_published_table_shape():
    """Synthetic fixtures test SOURCE parsing/integrity, NOT empirical
    results. The biological source run fetches original author paper.
    """
    head="Table 1 Distribution of the sample sites and detection of transgenic and non-transgenic oilseed rape in Switzerland\n"
    rows=[
      ("Lugano",23,21),
      ("Basel St. Johann",13,12),
      ("Kleinhüningen, unloading structure (port area)",18,16),
      ("Muttenz, Auhafen (port area)",47,1),
    ]
    # 54 other Swiss original-site placeholders with true original
    # study-level 2398 test denominator but not actual labels/data.
    # Never claim these are real original other individual site values.
    rows +=[(f"synthetic Swiss rail area {i:02d}",71 if i==0 else 42,0)
            for i in range(54)]
    assert len(rows)==58
    # original four positive 101; remaining should sum2297:
    # 71 + 53*42=2297
    rows.append(("Rümlang (field), NC",5,0))
    table=head
    for i,(name,n,pos) in enumerate(rows):
        # site label x/y and Swiss canton are fixture only
        if i==2:
            table+="02.05.2012 Kleinhüningen, unloading\n"
            table+="       structure (port area) 611308 270495 BS 18 16\n"
            continue
        tabname=name
        table+=f"02.05.2012 {tabname} {601001+i} 260001 BE {n} {pos}\n"
    table+="Total 2403 50\n"
    return table


class OriginalSwiss2012RailwayFeralsMarkersTests(unittest.TestCase):
    def test_full_59_row_pdf_table_footers_and_Rumlang_control_exclusion(self):
        tab=isolate_original_site_table(test_published_table_shape())
        rows=extract_original_2012_table(tab)
        self.assertEqual(len(rows),59)
        stats=describe_surveyed_original_railway(rows)
        self.assertEqual(stats["n_original_railway_areas_surveyed_from_published_text"],79)
        self.assertEqual(stats["n_source_railway_areas_WITH_feral_Bnapus_and_test_samples"],58)
        self.assertEqual(stats["original_RAILWAY_plants_tested_excluding_field_negative_control"],2398)
        self.assertEqual(stats["original_field_negative_control_not_railway_plants_tested"],5)
        self.assertEqual(stats["original_CP4_positive_feral_OSR_plants"],50)
        self.assertEqual(len(stats["source_GM_hotspot_localities_and_actual_GM_sample_counts"]),4)
        self.assertAlmostEqual(stats["plant_individual_weighted_positive_fraction_excluding_field_control"],50/2398)
        self.assertAlmostEqual(stats["original_publication_positive_fraction_including_NC_5"],50/2403)
        self.assertAlmostEqual(stats["original_railway_area_fraction_with_any_GM_positive"],4/79)

    def test_original_railway_2012_hotspot_concentration_and_size_only_null(self):
        rows=extract_original_2012_table(isolate_original_site_table(test_published_table_shape()))
        x=source_homogeneous_plant_label_randomization(rows,draws=500,seed=20261114)
        self.assertEqual(x["source_original_58_sites_GM_positive_clusters_observed"],4)
        self.assertEqual(x["n_montecarlo_plant_label_randomizations"],500)
        self.assertAlmostEqual(x["actual_fraction_all_50_GM_plant_tests_in_top_3_positive_sites"],49/50)
        self.assertGreater(x["source_null_mean_number_sites_containing_at_least_one_GM_positive"],10)

    def test_original_negative_control_included_as_rail_produces_explicit_error(self):
        tab=test_published_table_shape()
        tab=tab.replace("Rümlang (field), NC","Rümlang railway station")
        with self.assertRaisesRegex(ValueError,"Rümlang field 5-plant NC"):
            extract_original_2012_table(isolate_original_site_table(tab))

    def test_cannot_invent_missing_source_railway_sites(self):
        tab=test_published_table_shape()
        # Do not remove source numeric summary footer; remove one actual
        # corresponding original-date-plot site record from Table1 body.
        line=next(x for x in tab.splitlines() if "synthetic Swiss rail area 25" in x)
        tab=tab.replace(line+"\n","")
        with self.assertRaisesRegex(ValueError,"Expected original 58 feral plants"):
            extract_original_2012_table(isolate_original_site_table(tab))

    def test_four_original_hotspot_phenotypes_cannot_be_changed_silently(self):
        tab=test_published_table_shape().replace("Basel St. Johann 601002 260001 BE 13 12",
                       "Basel St. Johann 601002 260001 BE 13 11")
        with self.assertRaisesRegex(ValueError,"2403/50"):
            extract_original_2012_table(isolate_original_site_table(tab))


if __name__=="__main__":unittest.main()
