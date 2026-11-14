"""Daymet weather parser tests use generated weather, NOT field observations."""
import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
from tetsu.daymet import parse_daymet,summarize_annual_climate


def synthetic_daymet(years=(2023,2024,2025),days=365):
    rows=[
        "# Synthetic Daymet-looking text; NOT real ORNL data",
        "# temperature units C",
        "year,yday,tmax (deg c),tmin (deg c),prcp (mm/day)"
    ]
    for year in years:
        for day in range(1,days+1):
            rows.append(f"{year},{day},15,5,2")
    return "\n".join(rows)+"\n"


class DaymetWeatherTests(unittest.TestCase):
    def test_complete_years_and_fixed_gdd(self):
        by_year=parse_daymet(synthetic_daymet())
        self.assertEqual(sorted(by_year),[2023,2024,2025])
        x=summarize_annual_climate(by_year[2024])
        self.assertAlmostEqual(x["gdd5_to_doy90"],450.)
        self.assertAlmostEqual(x["gdd5_to_doy120"],600.)
        self.assertAlmostEqual(x["precip_mm_to_doy120"],240.)
        self.assertAlmostEqual(x["mean_air_c_to_doy120"],10.)

    def test_missing_days_rejected(self):
        with self.assertRaises(ValueError):
            parse_daymet(synthetic_daymet(days=364))

    def test_error_or_auth_html_is_not_parsed_as_weather(self):
        with self.assertRaises(ValueError):
            parse_daymet("<html>Login required</html>")

    def test_missing_weather_column_rejected(self):
        text=synthetic_daymet().replace(
            "prcp (mm/day)","snow (mm/day)"
        )
        with self.assertRaises(ValueError):
            parse_daymet(text)

    def test_unexpected_year_rejected(self):
        with self.assertRaises(ValueError):
            parse_daymet(synthetic_daymet(years=(2022,2023,2024)))

    def test_inverted_temperature_rejected(self):
        with self.assertRaises(ValueError):
            parse_daymet(synthetic_daymet().replace(
                "2023,1,15,5,2","2023,1,-20,5,2"
            ))


if __name__=="__main__":
    unittest.main()
