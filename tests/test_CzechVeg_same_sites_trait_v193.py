import importlib.util
import unittest
from pathlib import Path
p=Path(__file__).resolve().parents[1]/"scripts/analyze_CzechVeg_same_sites_regional_trait_convergence_v193.py"
spec=importlib.util.spec_from_file_location("trait193",p)
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)

class TestSameSourceTaxonomicTraitAnalysis(unittest.TestCase):
    def source_fixture(self):
        taxa=["A a","B b","C c","D d","E e","F f","G g","H h"]
        life={s:("therophyte" if i%2 else "hemicryptophyte")
              for i,s in enumerate(taxa)}
        disp={s:("Allium" if i%2 else "Epilobium") for i,s in enumerate(taxa)}
        regions=["R1","R2","R3"]
        frame={}
        plants={}
        for i,r in enumerate(regions):
            ra=f"rail_{r}"
            plants[ra]=set(taxa[:6])
            options={}
            for typ in m.CATEGORIES:
                keys=[f"{typ}_{r}_{j}" for j in range(5)]
                options[typ]=keys
                for j,k in enumerate(keys):
                    plants[k]=set(taxa[1:7]) if (j+i)%2 else set(taxa[:6])
            frame[r]=[(ra,options)]
        return frame,plants,life,disp

    def test_distinct_functional_and_taxonomic_measure(self):
        taxa={"A a","B b"}
        life={"A a":"T","B b":"T"}
        self.assertEqual(m.profile(taxa,life),{"T":1})
        self.assertAlmostEqual(m.profile_overlap({"T":.6,"H":.4},{"T":.2,"H":.8}),.6,places=12)
        self.assertEqual(m.profile_overlap({"T":1},{"H":1}),0)
        self.assertEqual(len(m.FOCUS_KEYS),16)

    def test_nothing_imputed_if_missing_trait(self):
        life={"A a":"T","B b":"H","C c":"H"}
        disp={"A a":"Allium","B b":"Epilobium"}
        plant={"A a","B b","C c","D d","E e","F f"}
        self.assertEqual(len(m.known_taxa(plant,"lifeform",life,disp)),3)
        self.assertEqual(len(m.known_taxa(plant,"both_known",life,disp)),2)
        self.assertFalse(m.site_eligible("x",{"x":plant},life,disp,"lifeform"))

    def test_identical_rail_site_for_all_reference_profiles(self):
        frame,plant,life,disp=self.source_fixture()
        a=m.result_for_frame(frame,plant,life,disp,seed=12,draws=9)
        self.assertTrue(a["identified"])
        self.assertEqual(a["n_original_phyto_regions"],3)
        self.assertEqual(a["n_region_pairs_nonindependent"],3)
        self.assertTrue(a["source_same_selected_rail_and_reference_site_ids_for_all_response_metrics"])
        for name in m.FOCUS_KEYS:
            self.assertIn(name,a["summary_statistics"])
            self.assertTrue(-1<=a["summary_statistics"][name]["mean"]<=1)

    def test_common_frame_wont_invent_same_rail_with_mismatched_reference(self):
        frame,plant,life,disp=self.source_fixture()
        all_sites=set(plant)
        meta={pid:{"region":pid.split("_")[1]} for pid in all_sites}
        arms={name:{} for name in m.CATEGORIES}
        for reg,items in frame.items():
            rid,cand=items[0]
            for group in m.CATEGORIES:arms[group][rid]=cand[group]
        selected,qa=m.eligible_frame(arms,meta,plant,life,disp,"lifeform",richness=True)
        self.assertEqual(qa["n_phyto_regions"],3)
        self.assertEqual(qa["n_study_source_rail_sites_with_all_three_reference_types_and_trait_quality"],3)
        del arms["nonanalog"]["rail_R2"]
        selected,qa=m.eligible_frame(arms,meta,plant,life,disp,"lifeform")
        self.assertEqual(qa["n_phyto_regions"],2)
        self.assertTrue(qa["all_taxonomy_and_traits_evaluated_on_exact_same_selected_site_ids"])

if __name__=="__main__":
    unittest.main()
