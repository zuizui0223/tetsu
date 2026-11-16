"""v11.9 checks: no invented CzechVeg rail plots or fuzzy Pladias trait matches."""
import csv,hashlib,sys,tempfile,unittest
from pathlib import Path
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import analyze_CzechVeg_Pladias_rail_road_diaspore_v119 as q


class OriginalCzechRailSourceJoinTests(unittest.TestCase):
    def test_exact_accepted_name_normalization_no_genus_imputation(self):
        self.assertEqual(q.normalize_taxon("  Abies   alba "), "Abies alba")
        self.assertEqual(q.normalize_taxon(" ×Festulolium  braunii "), "x Festulolium braunii")
        self.assertNotEqual(q.normalize_taxon("Salix alba"),q.normalize_taxon("Salix caprea"))

    def test_source_plot_labeling_not_paper_title(self):
        original=[
            ["1","16","železniční násep","pobřeží","rail survey","2001","18a"],
            ["2","16","silnice","okraj lesa","rail survey","2005","18a"],
            ["3","16","travnatý porost","poblíž kolejiště","botanical","2010","18b"],
            ["4","NA","silnice","Brno","road","2010","18b"],
        ]
        fields=["RELEVE_NR","SURF_AREA","HABITAT","LOCALITY",
                "REFERENCE_TITLE","YEAR","SKALICKY1988"]
        with tempfile.TemporaryDirectory() as d:
            file=Path(d)/"h.csv"
            with file.open("w",encoding="utf-8",newline="") as f:
                w=csv.writer(f);w.writerow(fields);w.writerows(original)
            with patch.object(q,"EXPECTED_PLOTS",4):
                plots,audit=q.plot_labels(file)
        self.assertEqual(audit["counts"]["n_strict_rail"],1)
        self.assertEqual(audit["counts"]["n_strict_road"],2)
        self.assertEqual(audit["counts"]["n_broad_rail"],2)
        self.assertEqual(plots["1"]["strict"],"rail")
        self.assertEqual(plots["2"]["strict"],"road")
        self.assertIsNone(plots["3"]["strict"])
        self.assertNotIn("4",plots)

    def test_species_guild_uses_source_site_opportunity_not_fake_railway_seed_records(self):
        candidate={}
        observed={}
        for i in range(24):
            plot="p"+str(i)
            candidate[plot]={"strict":"rail" if i<12 else "road",
                             "broad":"rail" if i<12 else "road",
                             "area":16,"region":"1","year":"2020"}
            values=[]
            for is_annual in (True,False):
                for is_wind in (True,False):
                    for k in range(4):
                        # All source taxa observed in a real original plot
                        if (i+int(is_wind)+k)%3!=0:
                            vals="Epilobium" if is_wind else "Allium"
                            life="therophyte" if is_annual else "hemicryptophyte"
                            values.append((f"s{int(is_annual)}{int(is_wind)}{k}",vals,life))
            observed[plot]=values
        v=q.guild_panel(candidate,observed,"strict")
        self.assertTrue(v["identified"])
        self.assertEqual(v["n_source_unique_joint_trait_taxa_in_original_rail_or_road_vegetation"],16)
        self.assertEqual(v["n_original_rail_candidate_plots"],12)
        self.assertEqual(v["n_original_road_control_plots"],12)
        self.assertIn("strict_wind_minus_other_species_equal_rail_vs_road_prevalence_gap",v)

    def test_missing_original_source_file_never_bypasses_zenodo_checksum(self):
        meta={"files":[{"key":q.SPECIES_NAME,
            "checksum":"md5:"+hashlib.md5(b"demo").hexdigest(),
            "size":400000000,
            "links":{"self":"https://zenodo.org/api/records/19231152/files/x/content"}}]}
        with tempfile.TemporaryDirectory() as d:
            f=Path(d)/"fake.csv";f.write_bytes(b"demo")
            with self.assertRaisesRegex(ValueError,"bytes differs"):
                q.verify_file(meta,q.SPECIES_NAME,f)

    def test_genuine_v117_and_v118_receipts_required_before_biology(self):
        with self.assertRaisesRegex(ValueError,"source receipts"):
            q.run({}, "a","b","c","d", {"complete":False},{"complete":False})


if __name__=="__main__":unittest.main()
