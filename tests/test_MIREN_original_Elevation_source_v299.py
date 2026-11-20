"""v29.9 original field elevation metadata QA unit tests, synthetic rows only."""
from __future__ import annotations
import csv,sys,tempfile,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import audit_MIREN_original_native_elevation_source_v299 as m

class OriginalMIRENElevationQA(unittest.TestCase):
    def make_csv(self,site_altitudes=(500,700,900,1100,1300)):
        rows=[]
        for r in ("R1","R2"):
            for i,a in enumerate(site_altitudes,1):
                for z in (1,3):
                    for species in ("Plant a","Plant b"):
                        rows.append({"Region":"AAA","Road":r,"Site":str(i),"Year":"2017",
                          "Plot":str(z),"Plot_id":f"AAA.{r}.{i}.{z}",
                          "Elevation":str(a+(z==3)*10),
                          "Accepted.Name.MIREN":species})
        fieldnames=list(rows[0])
        tmp=tempfile.TemporaryDirectory()
        src=Path(tmp.name)/"original.csv"
        with src.open("w",newline="",encoding="utf-8") as handle:
            writer=csv.DictWriter(handle,fieldnames=fieldnames)
            writer.writeheader();writer.writerows(rows)
        keys={(r["Region"],r["Road"],r["Site"],int(r["Year"]),int(r["Plot"])):{"Plant a"}
              for r in rows}
        roads=[(r,[{"region":"AAA","road":r,"site":str(i),"year":2017,
                    "roadside":{"Plant a"},"interior":{"Plant b"}} for i in range(1,len(site_altitudes)+1)])
               for r in ("R1","R2")]
        return tmp,src,keys,roads

    def test_original_author_true_plot_metadata_unique_not_species_replicates(self):
        tmp,src,keys,roads=self.make_csv()
        with tmp:
            alt,qa=m.read_original_elevation_csv(src,keys)
            self.assertEqual(qa["n_source_original_plots_with_elevation_metadata"],20)
            self.assertEqual(qa["n_source_plots_with_qualified_original_Elevation"],20)
            self.assertAlmostEqual(m.original_paired_site_elevation(
                roads[0][1][0],"AAA","R1",alt),505.)

    def test_predeclared_road_gradient_gate_original_site_unit(self):
        tmp,src,keys,roads=self.make_csv()
        with tmp:
            meta,qa=m.read_original_elevation_csv(src,keys)
            passed=m.evaluate_original_region_support(roads,"AAA",meta)
            self.assertTrue(passed["source_environmental_gradient_support_PASS"])
            self.assertEqual(passed["source_original_roads"],2)
            self.assertEqual(passed["source_original_paired_sites"],10)
            self.assertEqual(passed["source_n_paired_sites_retained_with_valid_elevation"],10)
            self.assertEqual(passed["source_road_n_with_ge4_sites_and_ge150m_elevation"],2)

    def test_no_slope_when_elevation_range_too_small(self):
        tmp,src,keys,roads=self.make_csv((500,520,540,560,570))
        with tmp:
            meta,qa=m.read_original_elevation_csv(src,keys)
            info=m.evaluate_original_region_support(roads,"AAA",meta)
            self.assertFalse(info["source_environmental_gradient_support_PASS"])
            self.assertEqual(info["source_road_n_with_ge4_sites_and_ge150m_elevation"],0)

    def test_two_real_original_roads_required(self):
        tmp,src,keys,roads=self.make_csv()
        with tmp:
            meta,_=m.read_original_elevation_csv(src,keys)
            info=m.evaluate_original_region_support(roads[:1],"AAA",meta)
            self.assertFalse(info["source_environmental_gradient_support_PASS"])

    def test_elevation_plot_pair_mismatch_returns_missing_not_false_zero(self):
        tmp,src,keys,roads=self.make_csv()
        with tmp:
            meta,_=m.read_original_elevation_csv(src,keys)
            key=("AAA","R1","1",2017,3)
            meta[key]=800.
            self.assertIsNone(m.original_paired_site_elevation(roads[0][1][0],"AAA","R1",meta))
            self.assertIsNotNone(m.original_paired_site_elevation(roads[0][1][1],"AAA","R1",meta))

    def test_missing_original_plot_metadata_fails_source_cohort_80pct_gate(self):
        tmp,src,keys,roads=self.make_csv()
        with tmp:
            meta,_=m.read_original_elevation_csv(src,keys)
            for site in ("1","2","3"):
                meta[("AAA","R2",site,2017,1)]=None
            support=m.evaluate_original_region_support(roads,"AAA",meta)
            self.assertFalse(support["source_environmental_gradient_support_PASS"])

    def test_original_author_field_years_not_implicitly_one(self):
        tmp,src,keys,roads=self.make_csv()
        with tmp:
            meta,_=m.read_original_elevation_csv(src,keys)
            roads[0][1][-1]["year"]=2018
            info=m.evaluate_original_region_support(roads,"AAA",meta)
            self.assertEqual(info["source_road_n_with_mixed_original_first_eligible_survey_years"],1)
            self.assertFalse(info["source_environmental_gradient_support_PASS"])

if __name__=="__main__":
    unittest.main()
