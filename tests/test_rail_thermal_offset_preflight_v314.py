"""Synthetic preflight tests; no ecological source data or effect estimates."""
import csv
import importlib.util
import tempfile
import unittest
from pathlib import Path

SCRIPT=Path(__file__).resolve().parents[1]/"scripts/rail_thermal_offset_preflight_v314.py"
spec=importlib.util.spec_from_file_location("thermal_preflight",SCRIPT)
mod=importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)

class PreflightTests(unittest.TestCase):
    def test_x_is_missing_not_zero(self):
        with tempfile.TemporaryDirectory() as d:
            src=Path(d)/"in.csv";dst=Path(d)/"out.csv"
            with src.open("w",newline="") as f:
                w=csv.writer(f)
                w.writerow(["study","region","site_id","habitat","species","incidence","temperature_indicator","regional_temperature","trait_source"])
                w.writerow(["a","r","1","rail","A",1,7,12,"original"])
                w.writerow(["a","r","1","rail","B",1,"x",12,"original"])
                w.writerow(["a","r","1","nonrail","A",1,5,12,"original"])
            self.assertEqual(mod.run(src,dst),1)
            with dst.open() as f: row=next(csv.DictReader(f))
            self.assertEqual(float(row["delta_cwm"]),2)
            self.assertEqual(float(row["rail_coverage"]),0.5)
    def test_no_matched_sites_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            src=Path(d)/"in.csv";dst=Path(d)/"out.csv"
            src.write_text("study,region,site_id,habitat,species,incidence,temperature_indicator,regional_temperature,trait_source\na,r,1,rail,A,1,7,12,original\n")
            with self.assertRaises(ValueError):mod.run(src,dst)

if __name__=="__main__":unittest.main()
