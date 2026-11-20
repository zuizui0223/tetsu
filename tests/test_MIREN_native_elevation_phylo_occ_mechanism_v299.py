"""v29.9 synthetic math tests only, never treated as observed real MIREN data."""
from __future__ import annotations
import copy,math,sys,unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import test_MIREN_source_elevation_phylo_occ_stature_v299 as m

def real_field_structure_synthetic_fixture():
    species=[f"Genus{chr(65+i//3)} species{chr(97+i%3)}" for i in range(24)]
    names=sorted(species)
    trait=np.asarray([(float(i)*.13,float((i*13)%31)*.09)
                       for i in range(len(names))])
    source=[]
    field=[]
    altitudes={}
    for rid in ("R1","R2"):
        cells={"roadside":[],"interior":[]}
        records=[]
        for sid in range(5):
            altitude=400+sid*220
            ra={names[(i+sid*2+(rid=="R2"))%24] for i in range(13)}
            ib={names[(i+sid*3+8+(rid=="R2"))%24] for i in range(12)}
            records.append({"region":"AAA","road":rid,
                            "site":str(sid+1),"year":2017,
                            "roadside":ra,"interior":ib})
            for zone,plants in (("roadside",ra),("interior",ib)):
                arr=np.zeros(len(names))
                for plant in plants: arr[names.index(plant)]=1.0/len(plants)
                cells[zone].append(arr)
            altitudes[("AAA",rid,str(sid+1),2017,1)]=float(altitude)
            altitudes[("AAA",rid,str(sid+1),2017,3)]=float(altitude+4)
        field.append((rid,records))
        source.append({k:np.asarray(v) for k,v in cells.items()})
    return field,{"species":names,"roads":source,"source_trait_vectors":trait},altitudes,set(names)

class OriginalMIRENNativeElevationMechanismMath(unittest.TestCase):
    def test_exact_native_elevation_signed_height_slope_known_value(self):
        W=np.eye(4,dtype=float)
        T=np.asarray([[0,7],[-1,6],[-2,5],[-3,4]],dtype=float)
        a=np.asarray([100,300,500,700],dtype=float)
        z=m.road_environment_slopes(W,a,T)
        self.assertAlmostEqual(z["height_signed_slope_per_km"],-5.,places=10)
        self.assertAlmostEqual(z["height_distance_vs_elevation_gap_slope_per_km"],5.,places=10)

    def test_observed_site_taxon_matrix_equals_real_original_road_zone_rows(self):
        field,source,alt,tipset=real_field_structure_synthetic_fixture()
        original=copy.deepcopy(field)
        roads,qa=m.construct_native_environment_cohort(
            field,"AAA",source,alt,5,tipset)
        self.assertEqual(field,original)
        self.assertTrue(qa["original_environmental_altitude_support"])
        self.assertEqual(qa["source_true_field_sites_on_usable_elevation_roads"],10)
        self.assertEqual(len(roads),2)
        for entry in roads:
            self.assertEqual(len(entry["elev_m"]),5)
            self.assertTrue(np.ptp(entry["elev_m"])>150.)
            self.assertEqual(entry["matrices"]["roadside"].shape,(5,24))
            self.assertEqual(entry["matrices"]["interior"].shape,(5,24))

    def test_original_site_species_relabel_must_not_silently_reorder(self):
        field,source,alt,tipset=real_field_structure_synthetic_fixture()
        fields=copy.deepcopy(field)
        fields[0][1][0]["roadside"],fields[0][1][1]["roadside"]=(
            fields[0][1][1]["roadside"],fields[0][1][0]["roadside"])
        with self.assertRaises(ValueError):
            m.construct_native_environment_cohort(fields,"AAA",source,alt,5,tipset)

    def test_no_altitude_when_same_site_zone_elevations_disagree(self):
        field,source,alt,tipset=real_field_structure_synthetic_fixture()
        alt[("AAA","R1","1",2017,3)]=1000.
        roads,qa=m.construct_native_environment_cohort(field,"AAA",source,alt,5,tipset)
        self.assertEqual(qa["source_plant_site_pair_metadata_missing_elevation_n"],1)
        self.assertEqual(qa["source_true_field_sites_on_usable_elevation_roads"],9)

    def test_no_altitude_cohort_effect_is_not_zero(self):
        field,source,alt,tipset=real_field_structure_synthetic_fixture()
        for road in ("R1","R2"):
            for sid in (1,2,3):
                alt[("AAA",road,str(sid),2017,1)]=None
        roads,qa=m.construct_native_environment_cohort(field,"AAA",source,alt,5,tipset)
        self.assertFalse(qa["original_environmental_altitude_support"])
        self.assertLess(qa["fraction_original_phylotip_matched_sites_retained_for_elevation"],.8)

    def test_original_tree_species_universe_and_true_trait_known_pool_not_conflated(self):
        field,source,alt,tipset=real_field_structure_synthetic_fixture()
        # Existing source tree lists all 24 species, but one species is
        # absent from the final source retained leaf pool. The gate must
        # compare original phylogenetic source membership before year joining.
        missing=set(tipset)
        missing.remove(next(iter(missing)))
        try:
            roads,qa=m.construct_native_environment_cohort(field,"AAA",source,alt,5,missing)
        except ValueError:
            return
        self.assertLessEqual(qa["source_true_field_sites_on_usable_elevation_roads"],10)

    def test_mixed_year_true_road_is_flagged_and_year_consistent_sensitivity_discards(self):
        field,source,alt,tipset=real_field_structure_synthetic_fixture()
        field[0][1][-1]["year"]=2018
        # Keep original site/altitude metadata valid in its actual year.
        alt[("AAA","R1","5",2018,1)]=1280.
        alt[("AAA","R1","5",2018,3)]=1284.
        roads,qa=m.construct_native_environment_cohort(field,"AAA",source,alt,5,tipset)
        self.assertEqual(qa["source_original_phylo_roads_with_mixed_first_eligible_survey_year"],1)
        self.assertEqual(len(m.original_year_consistent_roads(roads)),1)

    def test_two_habitat_original_environment_beta_weighted_by_roads(self):
        field,source,alt,tipset=real_field_structure_synthetic_fixture()
        roads,qa=m.construct_native_environment_cohort(field,"AAA",source,alt,5,tipset)
        out=m.region_environment_source_scores(roads,source["source_trait_vectors"])
        for metric in m.METRICS:
            self.assertAlmostEqual(out["roadside_minus_interior"][metric],
                                   out["roadside"][metric]-out["interior"][metric],places=12)

    def test_trait_null_reproducible_and_exact_site_elevation_unchanged(self):
        field,source,alt,tipset=real_field_structure_synthetic_fixture()
        roads,qa=m.construct_native_environment_cohort(field,"AAA",source,alt,5,tipset)
        original=[v["elev_m"].copy() for v in roads]
        groups=[np.asarray(g,dtype=int) for g in np.array_split(np.arange(24),6)]
        a=m.label_null_environment_slopes(
            roads,source["source_trait_vectors"],groups,0,0,"original_synthetic_test",draws=5)
        b=m.label_null_environment_slopes(
            roads,source["source_trait_vectors"],groups,0,0,"original_synthetic_test",draws=5)
        self.assertEqual(a,b)
        for i,r in enumerate(roads):
            np.testing.assert_array_equal(r["elev_m"],original[i])
        for metric in m.METRICS:
            self.assertTrue(math.isfinite(a["roadside"][metric][
                "original_observed_minus_conditional_null_environmental_slope"]))
            self.assertAlmostEqual(
                a["roadside_minus_interior"][metric][
                    "original_observed_minus_conditional_null_environmental_slope"],
                a["roadside"][metric]["original_observed_minus_conditional_null_environmental_slope"]-
                a["interior"][metric]["original_observed_minus_conditional_null_environmental_slope"],places=12)

    def test_unsupported_region_not_given_fake_slope(self):
        self.assertEqual(m.source_environment_summary({})["source_status"],
                         "INSUFFICIENT_SOURCE_REGIONAL_SUPPORT")
        self.assertIsNone(m.source_environment_summary({})["region_equal_estimates"])

if __name__=="__main__":
    unittest.main()
