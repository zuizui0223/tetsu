import importlib.util
import pathlib
import unittest
ROOT=pathlib.Path(__file__).resolve().parents[1]
SCRIPT=ROOT/'scripts'/'audit_Swiss_route_supply_v135.py'

class RouteEvidenceTests(unittest.TestCase):
    def setUp(self):
        self.assertTrue(SCRIPT.exists(), 'route evidence parser not implemented')
        spec=importlib.util.spec_from_file_location('v135',SCRIPT)
        self.m=importlib.util.module_from_spec(spec);spec.loader.exec_module(self.m)
        self.text='''2.5.6 Monitoring von GV-Raps auf dem SBB-Gleisnetz
Risikobasierte Probenahmen zwischen dem Hafen Kleinhüningen und Jowa und Swissmill,
wobei 1314 Rapsproben erhoben wurden. Alle Proben wurden GV-negativ getestet.
Nicht-risikobasiertes Monitoring auf 30 km: 391 Rapsproben erhoben.
Dabei wurden sieben Rapsproben an einem Standort in Giornico (TI) als GV-Raps (GT73) identifiziert.
2.5.7 Untersuchungen von Importweizen auf Verunreinigung mit GV-Raps
Jowa und Swissmill: jeweils 10 Proben; von jeder Getreidemühle 150 Rapssamen.
Die Keimraten betrugen durchschnittlich 25%; GT73, MS8xRF3, MS8 und RF3.
Die Verunreinigung betrug 0.005%. Dies entspricht 16 Rapssamen pro kg.
Etwa zwei keimfähige transgene Rapssamen pro kg.
2.5.8 Badewasser
'''
    def test_matched_route_counts_and_unknown_site_denominator(self):
        out=self.m.analyze_text(self.text)
        self.assertEqual(out['route_monitoring']['targeted_route']['positive_samples'],0)
        self.assertEqual(out['route_monitoring']['targeted_route']['tested_samples'],1314)
        self.assertAlmostEqual(out['route_monitoring']['other_survey']['positive_fraction'],7/391)
        self.assertEqual(out['route_monitoring']['other_survey']['positive_locations'],1)
        self.assertIsNone(out['route_monitoring']['targeted_route']['total_locations'])
        self.assertFalse(out['field_establishment_rate_identified'])
    def test_wrong_report_number_fails_closed(self):
        with self.assertRaises(ValueError):
            self.m.analyze_text(self.text.replace('1314 Rapsproben','1315 Rapsproben'))
    def test_germination_aggregate_is_corroboration_not_independent_study(self):
        out=self.m.analyze_text(self.text)
        self.assertTrue(out['same_cargo_study_republished_not_independent_replication'])
        self.assertFalse(out['journal_Table2_cell_audit_performed'])
if __name__=='__main__':unittest.main()
