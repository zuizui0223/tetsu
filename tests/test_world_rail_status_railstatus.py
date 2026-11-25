import sys, unittest, random
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from aggregate_world_rail_by_status_railstatus import category, rollup
from analyze_world_rail_status_richness_railstatus import inject_status, shuffled, band


def row(cls, flags=()):
    return {"subtype": "rail", "class": cls,
            "rail_flags": [{"values": list(flags)}] if flags else None}


class CategoryTests(unittest.TestCase):
    def test_active_matches_v370_rule(self):
        self.assertEqual(category(row("standard_gauge")), "active_normal")
        self.assertEqual(category(row("unknown")), "active_normal")
        self.assertEqual(category(row("narrow_gauge", ("is_bridge",))), "active_normal")

    def test_inactive(self):
        self.assertEqual(category(row("standard_gauge", ("is_abandoned",))), "inactive_normal")
        self.assertEqual(category(row("narrow_gauge", ("is_disused",))), "inactive_normal")
        self.assertIsNone(category(row("standard_gauge", ("is_under_construction",))))
        self.assertIsNone(category(row("standard_gauge", ("is_abandoned", "is_under_construction"))))

    def test_urban_classes(self):
        self.assertEqual(category(row("subway")), "subway")
        self.assertEqual(category(row("tram")), "tram_lightrail")
        self.assertEqual(category(row("light_rail")), "tram_lightrail")
        self.assertIsNone(category(row("subway", ("is_abandoned",))))
        self.assertIsNone(category(row("funicular")))
        self.assertIsNone(category(row("monorail")))

    def test_rollup_200km(self):
        g = rollup({(0, 0): 1000.0, (1, 1): 2000.0, (2, 0): 500.0}, 200000)
        d = {(x["source_equal_area_grid_ix"], x["source_equal_area_grid_iy"]): x["km"] for x in g}
        self.assertEqual(d, {(0, 0): 3.0, (1, 0): 0.5})


class InjectTests(unittest.TestCase):
    def status(self, active_km):
        cats = {c: {"100000": []} for c in ("active_normal", "inactive_normal", "subway", "tram_lightrail")}
        cats["active_normal"]["100000"] = [{"source_equal_area_grid_ix": 1, "source_equal_area_grid_iy": 2, "km": active_km}]
        cats["inactive_normal"]["100000"] = [{"source_equal_area_grid_ix": 1, "source_equal_area_grid_iy": 2, "km": 4.0}]
        return {"grids_km_by_category_and_width": cats}

    def test_reproduction_gate(self):
        rows = [{"key": (1, 2), "rail_km": 10.0}, {"key": (3, 3), "rail_km": 0.0}]
        out, worst = inject_status(rows, self.status(10.0), 100000)
        self.assertEqual(worst, 0.0)
        self.assertAlmostEqual(out[0]["inactive_normal_km"], 4.0)
        self.assertEqual(out[1]["log_subway"], 0.0)
        with self.assertRaises(ValueError):
            inject_status(rows, self.status(11.0), 100000)

    def test_shuffle_keeps_band_multiset(self):
        rng = random.Random(1)
        rows = [{"continent": "Europe", "bio1": t, "abs_lat": 45 + (i % 3), "x": float(i)}
                for i, t in enumerate([5, 5, 5, 12, 12, 25])]
        sh = shuffled(rows, "x", rng)
        for b in {band(r) for r in rows}:
            orig = sorted(r["x"] for r in rows if band(r) == b)
            new = sorted(r["x"] for r in sh if band(r) == b)
            self.assertEqual(orig, new)


if __name__ == "__main__":
    unittest.main()
