"""Synthetic QA: matched same taxon within a real source Site/Year in both
roadside and interior, preserving original-region folds (not actual MIREN)."""
import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import analyze_MIREN_source_taxon_matched_road_interior_region_holdout_v268 as m


def fake_rows(regions=4,roads=2,sites=5):
    rows=[]
    for region in range(regions):
        for road in range(roads):
            for site in range(sites):
                for plot in ("1","3"):
                    for taxon in ("commonA","commonB","commonC","commonD",
                                  "onlyroad" if plot=="1" else "onlyinterior"):
                        number=("commonA","commonB","commonC","commonD","onlyroad","onlyinterior").index(taxon)
                        rank=2*((number+site+road+region+(plot=="3"))%5)/4-1
                        rows.append({
                            "region":"REG"+str(region),
                            "road":"road"+str(road),
                            "site":"site"+str(site),
                            "plot":plot,
                            "taxon":taxon,
                            "initial_year":2017,
                            "year_gap":5,
                            "xrank":rank,
                            "re_recorded":int(rank>=0)
                        })
    return rows


class MatchedRoadInteriorTaxonTests(unittest.TestCase):
    def test_exact_source_taxon_and_location_match(self):
        rows=fake_rows()
        matched,qc=m.matched_taxa_at_same_original_source_site(rows)
        self.assertEqual(qc["original_unique_early_source_site_count"],4*2*5)
        self.assertEqual(qc["original_sites_with_at_least_one_2012plus_matched_taxon"],40)
        self.assertEqual(qc["source_same_site_taxon_recorded_in_both_habitats_pairs"],40*4)
        self.assertEqual(qc["original_2012plus_spatial_sites_with_zero_shared_taxa"],0)
        for site,pairs in matched.items():
            self.assertEqual(len(pairs),4)
            for i,j in pairs:
                self.assertEqual(rows[i]["taxon"],rows[j]["taxon"])
                self.assertEqual(rows[i]["plot"],"1")
                self.assertEqual(rows[j]["plot"],"3")
                self.assertEqual(rows[i]["region"],rows[j]["region"])

    def test_no_same_taxon_different_year_or_site_match(self):
        rows=fake_rows()
        for r in rows:
            if r["taxon"]=="commonA" and r["plot"]=="3" and r["site"]=="site0":
                r["year_gap"]=6
                break
        with self.assertRaises(ValueError):
            m.matched_taxa_at_same_original_source_site(rows)

    def test_original_region_out_training_prediction_and_region_weighting(self):
        rows=fake_rows()
        pairs,qa=m.matched_taxa_at_same_original_source_site(rows)
        results=m.evaluate(rows,pairs)
        self.assertEqual(set(results),{"REG0","REG1","REG2","REG3"})
        for result in results.values():
            self.assertEqual(result["source_roads_with_at_least_one_matched_site"],2)
            self.assertEqual(result["original_spatial_sites_with_shared_taxa"],10)
            self.assertEqual(result["original_2012plus_source_same_taxon_both_zones_pairs"],40)
        for threshold in m.THRESHOLDS:
            out=m.summarize(results,threshold)
            self.assertEqual(out["source_region_codes_identified"],4)
            self.assertEqual(out["source_roads_in_support"],8)
            self.assertEqual(out["source_spatial_sites_in_support"],40)
            self.assertEqual(out["source_shared_taxon_site_pairs_in_support"],160)
            self.assertTrue(out["equal_region_results"]["logloss"]["M2Z_to_M3Z"]["identified"])

    def test_original_site_first_year_and_taxon_identity_gate(self):
        rows=fake_rows(regions=3,roads=2,sites=4)
        raw=[x.copy() for x in rows]
        pair,q=m.matched_taxa_at_same_original_source_site(raw)
        self.assertEqual(q["source_same_site_taxon_recorded_in_both_habitats_pairs"],3*2*4*4)
        # A taxon only in Plot1 cannot become a paired site-taxon observation
        self.assertTrue(all(raw[i]["taxon"]!="onlyroad" and raw[j]["taxon"]!="onlyinterior"
                            for pairs in pair.values() for i,j in pairs))

    def test_no_raw_site_ids_in_macro_output_schema(self):
        rows=fake_rows()
        pairs,_=m.matched_taxa_at_same_original_source_site(rows)
        result=m.summarize(m.evaluate(rows,pairs),2)
        # Macro output has region and support counts; source site keys remain
        # internal to evaluate() and are deliberately not exported by run().
        self.assertTrue(all(set(v).issubset({
            "original_roads_in_support","original_sites_in_support",
            "original_matched_species_site_pairs_in_support",
            "original_macro_equal_road_scores"}) for v in result["region_code_results"].values()))


if __name__=="__main__":unittest.main()
