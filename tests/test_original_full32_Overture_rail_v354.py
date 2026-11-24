#!/usr/bin/env python3
"""v354: six-area original source partition and first six verified baseline guard."""
import json,sys,tempfile,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_original_full32_Overture_rail_v354 import (
    CONT_CODES,source_partition,consolidate,RAIL_RELEASE,
    WAVE1_REFERENCE,SOURCE_FRAME_SHA
)
from audit_original_iNat_world_probability_cells_v353 import selected_cells

ROOT=Path(__file__).resolve().parents[1]
FRAME=ROOT/"data/derived/v351_global_land_climate_100km_probabilistic_source_verified.json"
BASE=ROOT/WAVE1_REFERENCE

class V354GlobalSourceAuditTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.frame=json.loads(FRAME.read_text())
        cls.baseline=json.loads(BASE.read_text())

    def test_exact_32_cells_partition_six_continents(self):
        counts={c:len(source_partition(self.frame,c)) for c in CONT_CODES.values()}
        self.assertEqual(counts,{"Africa":4,"Asia":6,"Europe":4,
             "North America":6,"South America":6,"Oceania":6})
        ids=[r["cell"]["source_original_fixed_grid_index"] for c in CONT_CODES.values()
              for r in source_partition(self.frame,c)]
        self.assertEqual(len(ids),32)
        self.assertEqual(len(set(ids)),32)

    def test_first_six_previous_original_source_ids_in_frozen_full32(self):
        sourceids={r["cell"]["source_original_fixed_grid_index"] for r in selected_cells(self.frame)}
        baselineids={r["source_grid"] for r in self.baseline["source_rail_candidate_100km_cells"]}
        self.assertEqual(len(baselineids),6)
        self.assertTrue(baselineids.issubset(sourceids))
        self.assertEqual(self.baseline["original_land_grid_frame_sha256"],SOURCE_FRAME_SHA)

    def test_missing_continent_aggregate_fails_not_zero(self):
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaises(FileNotFoundError):
                consolidate(FRAME,tmp,BASE)

    def test_complete_synthetic_fallback_cannot_modify_original_six_values(self):
        baseline={r["source_grid"]:r["source_mapped_rail"] for r in self.baseline["source_rail_candidate_100km_cells"]}
        with tempfile.TemporaryDirectory() as tmp:
            for code,continent in CONT_CODES.items():
                rows=[]
                for r in source_partition(self.frame,continent):
                    idx=r["cell"]["source_original_fixed_grid_index"]
                    rail=baseline.get(idx,{
                        "source_complete":True,
                        "source_rail_candidate_in_original_land_cell":False,
                        "source_rail_geodesic_km_in_original_cell":0.0
                    })
                    rows.append({"continent":continent,
                      "source_climate_band":r["macrotemp_stratum"],
                      "source_original_grid_ix_iy":idx,
                      "source_original_BIO1_C":r["cell"]["original_macroclimate_BIO1_C"],
                      "source_Overture_rail":rail})
                obj={"source_frame_SHA":SOURCE_FRAME_SHA,"status":"COMPLETE-continent-source-partition",
                     "source_continent":continent,"source_Overture_release":RAIL_RELEASE,
                     "source_original_cells":rows}
                Path(tmp,f"v354_{code}.json").write_text(json.dumps(obj))
            result=consolidate(FRAME,tmp,BASE)
            self.assertEqual(result["total_original_source_landcells_n"],32)
            self.assertTrue(result["original_v352_six_source_baseline_reproduced"])
            self.assertEqual(result["source_rail_present_n"],3)
            self.assertTrue(result["not_global_network_length_weighted"])
            # Tamper with the prior published original-source baseline: must fail.
            path=Path(tmp,"v354_na.json")
            obj=json.loads(path.read_text())
            source_first=next(r for r in obj["source_original_cells"]
                              if r["source_original_grid_ix_iy"] in baseline)
            source_first["source_Overture_rail"]["source_rail_geodesic_km_in_original_cell"]+=1
            path.write_text(json.dumps(obj))
            with self.assertRaises(AssertionError):
                consolidate(FRAME,tmp,BASE)

if __name__=="__main__":unittest.main()
