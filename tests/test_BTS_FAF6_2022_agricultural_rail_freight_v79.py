"""Synthetic FAF6 CSV ZIP source-mode/commodity and noncausal provenance gates."""
import io,sys,zipfile,csv,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_BTS_FAF6_2022_agricultural_rail_freight_v79 import (
    normalize_code,source_float,records_from_faf6_zip,summarize_records,
    MODES,CROPLINKED,RESOURCE_CONTROLS)

def toy_zip():
    f=io.StringIO()
    cols=("dms_orig","dms_dest","sctg2","dms_mode","trade_type","tons_2022","value_2022")
    writer=csv.DictWriter(f,fieldnames=cols)
    writer.writeheader()
    for i in range(600):
        writer.writerow({
            "dms_orig":str(i%12+1),"dms_dest":str((i+2)%12+1),
            "sctg2":str([2,3,10,11,12,20][i%6]),
            "dms_mode":str((i//6)%2+1),
            "trade_type":"1",
            "tons_2022":str((i%8+1)*2.5),
            "value_2022":str(i%17+1)})
    buf=io.BytesIO()
    with zipfile.ZipFile(buf,"w",compression=zipfile.ZIP_DEFLATED) as z:
        z.writestr("FAF6.0_State.csv",f.getvalue())
    return buf.getvalue()

class FAF6SourceIntegrity(unittest.TestCase):
    def test_commodity_codes_and_rail_truck_are_precommitted(self):
        self.assertEqual(MODES,{"1":"TRUCK","2":"RAIL"})
        self.assertEqual(CROPLINKED,{"02","03"})
        self.assertEqual(RESOURCE_CONTROLS,{"10","11","12"})
    def test_faf6_zip_exact_predeclared_fields_and_aggregate_rail_truck(self):
        result=summarize_records(records_from_faf6_zip(toy_zip()))
        self.assertEqual(result["n_original_official_FAF6_state_OD_commodity_mode_rows"],600)
        x=result["national_trade_type_x_commodity_group_2022"]["domestic"]["crop_codes_02_03"]
        self.assertTrue(x["rail_thousand_short_tons"]>0)
        self.assertTrue(x["truck_thousand_short_tons"]>0)
        self.assertTrue(0<x["rail_fraction_among_only_rail_and_truck"]<1)
        self.assertTrue(result["NOT_actual_verified_seed_propagules_or_track_section_freight_volume"])
    def test_invalid_freight_tons_do_not_create_flow(self):
        with self.assertRaisesRegex(ValueError,"Invalid negative"):
            source_float(-1)
        with self.assertRaisesRegex(ValueError,"Malformed"):
            source_float("missing")
    def test_bad_csv_source_rejected(self):
        b=io.BytesIO()
        with zipfile.ZipFile(b,"w") as z:z.writestr("FAF6.0_State.csv","a,b,c\n1,2,3\n")
        with self.assertRaisesRegex(ValueError,"size outside|missing"):
            list(records_from_faf6_zip(b.getvalue()))
    def test_faf_codes_no_imputed_native_or_seed_content(self):
        self.assertEqual(normalize_code("2"),"02")
        self.assertEqual(normalize_code("03"),"03")
        self.assertIsNone(normalize_code("invalid"))
        self.assertNotIn("seed_species",MODES)
if __name__=="__main__":unittest.main()
