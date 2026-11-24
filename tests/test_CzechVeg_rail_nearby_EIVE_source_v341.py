"""v341 original CzechVeg railway plot to nearby same-source vascular flora analysis tests."""
import unittest,sys,math
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import analyze_CzechVeg_rail_nearby_EIVE_source_v341 as m

TRAITS={"poa annua":{"position":7.,"width":3.,"systems":3},
        "galium aparine":{"position":5.,"width":4.,"systems":3},
        "lotus corniculatus":{"position":6.,"width":5.,"systems":1}}
def metadata(id,ref,group,lat,lon,year=2001,area=10,precision=20):
    return {"id":id,"plot_id":id,"reference":ref,"lat":lat,"lon":lon,
        "area":area,"year":year,"precision":precision,"altitude":300,
        "cov_trees":0,"cov_shrubs":0,
        "role":"rail" if group!="control" else "potential_reference",
        "habitat_type":group}

class CzechVegEIVETests(unittest.TestCase):
    def test_exact_species_name_no_hallucinated_synonym(self):
        self.assertEqual(m.source_exact_species("Poa annua"),"poa annua")
        self.assertIsNone(m.source_exact_species("Taraxacum sect. Taraxacum"))
        self.assertIsNone(m.source_exact_species("Poa spp."))
    def test_conservative_geographic_matching_outcome_blind(self):
        rows=[
          metadata("r1","SRC","railway_track_yard_or_structure",49.0,15.0),
          metadata("r2","SRC","railway_embankment_ditch_or_cut",49.0008,15.0),
          metadata("c1","SRC","control",49.0010,15.0),
          metadata("c2","SRC","control",49.0013,15.0008),
        ]
        pairs,qa=m.match_source_outcome_blind(rows,"G1_same_source_high_precision")
        self.assertEqual(len(pairs),2)
        self.assertEqual(qa["same_study_matched_pairs"],2)
        self.assertEqual(len({p["reference"]["id"] for p in pairs}),2)
    def test_different_study_not_same_study_railway_reference(self):
        rows=[
          metadata("r1","SRC","railway_track_yard_or_structure",49.0,15.0),
          metadata("c1","OTHER","control",49.001,15.0)]
        same,qa=m.match_source_outcome_blind(rows,"G2_same_source_extended")
        cross,cq=m.match_source_outcome_blind(rows,"G4_cross_source_extended")
        self.assertEqual(len(same),0)
        self.assertEqual(len(cross),1)
        self.assertEqual(cq["same_study_matched_pairs"],0)
    def test_known_trait_coverage_denominator_includes_unknown_taxon(self):
        flora={"poa annua":10.0,
               "galium aparine":20.0,
               "unresolved_source_taxon:rubus sect. rubus":None}
        q=m.compute_community(flora,TRAITS,"position",3)
        self.assertEqual(q["source_species_n"],3)
        self.assertEqual(q["trait_n"],2)
        self.assertAlmostEqual(q["trait_coverage"],2/3)
    def test_cover_weighting_uses_source_estimated_percentage(self):
        flora={"poa annua":10.0,"galium aparine":30.0}
        q=m.compute_community(flora,TRAITS,"position",3)
        self.assertEqual(q["equal_species_cwm"],6)
        self.assertAlmostEqual(q["cover_weighted_cwm"],5.5)
        self.assertEqual(q["trait_source_cover_nonmissing_n"],2)
    def test_missing_cover_not_zero(self):
        flora={"poa annua":None,"galium aparine":2.0}
        q=m.compute_community(flora,TRAITS,"position",3)
        self.assertEqual(q["trait_n"],2)
        self.assertEqual(q["trait_source_cover_nonmissing_n"],1)
        self.assertEqual(q["cover_weighted_cwm"],5)
    def test_strict_panel_does_not_impute_absent_temperature(self):
        rail=metadata("r1","SRC","railway_track_yard_or_structure",49.0,15.0)
        ctrl=metadata("c1","SRC","control",49.001,15.0)
        p={"rail":rail,"reference":ctrl,"rail_habitat":rail["habitat_type"]}
        plants={"r1":{"poa annua":1,"galium aparine":2},
                "c1":{"poa annua":2,"galium aparine":1}}
        result=m.panel_stats([p],plants,TRAITS,3,5,.5,"position")
        self.assertEqual(result["n_original_field_plot_pairs_with_EIVE"],0)
        self.assertIsNone(result["study_equal_mean_rail_minus_nearby_unlabelled"])
        self.assertEqual(result["source_qa"]["pairs_below_EIVE_numeric_trait_gate"],1)
    def test_two_original_sources_report_study_equal_not_pair_equal(self):
        rails=[metadata("r1","SRC1","railway_track_yard_or_structure",49.0,15.0),
               metadata("r2","SRC2","railway_embankment_ditch_or_cut",49.1,15.0)]
        ctrls=[metadata("c1","SRC1","control",49.001,15),
               metadata("c2","SRC2","control",49.101,15)]
        pairs=[{"rail":a,"reference":b,"rail_habitat":a["habitat_type"]}
               for a,b in zip(rails,ctrls)]
        plants={"r1":{"poa annua":1},"c1":{"galium aparine":1},
                "r2":{"galium aparine":1},"c2":{"poa annua":1}}
        res=m.panel_stats(pairs,plants,TRAITS,3,1,0,"position")
        self.assertEqual(res["n_original_studies_with_eligible_pairs"],2)
        self.assertAlmostEqual(res["study_equal_mean_rail_minus_nearby_unlabelled"],0)
    def test_source_qc_frozen_outcomes(self):
        self.assertEqual(m.EXPECTED_PLOTS,117739)
        self.assertEqual(len(m.MAIN_GATES),3)
        self.assertEqual(m.PANELS,((3,5,.5),(1,3,.3)))

if __name__=="__main__":unittest.main()
