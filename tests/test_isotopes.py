"""Self-tests runnable with ``python -m unittest`` or directly."""

import json
import os
import subprocess
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from isotopes import ElementTable, IsotopeError, compute_cluster, parse_formula
from isotopes.runner import run


TABLE = {
    "H": [
        {"mass_number": 1, "exact_mass": 1.0078, "abundance": 0.9},
        {"mass_number": 2, "exact_mass": 2.0141, "abundance": 0.1},
    ],
    "O": [
        {"mass_number": 16, "exact_mass": 15.995, "abundance": 0.98},
        {"mass_number": 18, "exact_mass": 17.999, "abundance": 0.02},
    ],
    "C": [
        {"mass_number": 12, "exact_mass": 12.0, "abundance": 1.0},
    ],
}


class FormulaTests(unittest.TestCase):
    def test_simple_and_omitted_subscript(self):
        self.assertEqual(parse_formula("H2O"), {"H": 2, "O": 1})

    def test_nested_parentheses(self):
        self.assertEqual(parse_formula("Ca(OH)2"), {"Ca": 1, "O": 2, "H": 2})
        self.assertEqual(parse_formula("((CH3)2CH)2"),
                         {"C": 6, "H": 14})

    def test_errors_report_position(self):
        cases = [
            ("", None),
            ("()", 2),
            ("H0", 2),
            ("H2(O", None),
            ("H)2", 2),
            ("2H", 1),
            ("X", 1),
            ("H ", 2),
        ]
        for formula, position in cases:
            with self.subTest(formula=formula):
                with self.assertRaises(IsotopeError) as ctx:
                    parse_formula(formula, {"H", "O", "C"})
                if position is not None:
                    self.assertEqual(ctx.exception.position, position)


class ElementTableTests(unittest.TestCase):
    def test_normalises_abundance(self):
        table = ElementTable({
            "X": [
                {"mass_number": 1, "exact_mass": 1.0, "abundance": 0.8999999999},
                {"mass_number": 2, "exact_mass": 2.0, "abundance": 0.1000000000},
            ]
        })
        self.assertAlmostEqual(
            sum(i.abundance for i in table.elements["X"]), 1.0)
        self.assertAlmostEqual(table.elements["X"][0].abundance, 0.9, places=9)

    def test_invalid_tables(self):
        bad = [
            {},
            {"X": []},
            {"X": [{"mass_number": 0, "exact_mass": 1.0, "abundance": 1}]},
            {"X": [{"mass_number": 1, "exact_mass": 1.0, "abundance": 0.5},
                   {"mass_number": 1, "exact_mass": 2.0, "abundance": 0.5}]},
            {"X": [{"mass_number": 1, "exact_mass": 1.0, "abundance": -0.1}]},
            {"X": [{"mass_number": 1, "exact_mass": float("nan"), "abundance": 1}]},
            {"X": [{"mass_number": 1, "exact_mass": 1.0, "abundance": 0.5}]},
        ]
        for raw in bad:
            with self.subTest(raw=raw):
                with self.assertRaises(IsotopeError):
                    ElementTable(raw)


class ClusterTests(unittest.TestCase):
    def setUp(self):
        self.table = ElementTable(TABLE)

    def test_two_hydrogen_distribution(self):
        cluster = compute_cluster("water", "H2", {"H": 2}, self.table, 1)
        peaks = {p.mass_number: p for p in cluster.peaks}
        self.assertEqual(sorted(peaks), [2, 3, 4])
        self.assertAlmostEqual(peaks[2].probability, 0.81)
        self.assertAlmostEqual(peaks[3].probability, 0.18)
        self.assertAlmostEqual(peaks[4].probability, 0.01)
        self.assertAlmostEqual(cluster.total_probability, 1.0)
        self.assertAlmostEqual(peaks[2].mean_mass, 2 * 1.0078)
        self.assertAlmostEqual(peaks[3].mean_mass, 1.0078 + 2.0141)

    def test_same_nominal_mass_is_merged_with_weighted_mean(self):
        data = {
            "X": [
                {"mass_number": 10, "exact_mass": 9.9, "abundance": 0.6},
                {"mass_number": 11, "exact_mass": 11.2, "abundance": 0.4},
            ],
            "Y": [
                {"mass_number": 10, "exact_mass": 10.1, "abundance": 0.5},
                {"mass_number": 9, "exact_mass": 8.8, "abundance": 0.5},
            ],
        }
        table = ElementTable(data)
        cluster = compute_cluster("xy", "XY", {"X": 1, "Y": 1}, table, 2)
        peak20 = next(p for p in cluster.peaks if p.mass_number == 20)
        self.assertAlmostEqual(peak20.probability, 0.3 + 0.2)
        expected_mass = (0.3 * (9.9 + 10.1) + 0.2 * (11.2 + 8.8)) / 0.5
        self.assertAlmostEqual(peak20.mean_mass, expected_mass)
        self.assertAlmostEqual(peak20.mz, expected_mass / 2)

    def test_mean_molecular_mass(self):
        cluster = compute_cluster("water", "H2O",
                                  {"H": 2, "O": 1}, self.table, 1)
        expected = 2 * (0.9 * 1.0078 + 0.1 * 2.0141) + (
            0.98 * 15.995 + 0.02 * 17.999)
        self.assertAlmostEqual(cluster.mean_molecular_mass, expected, places=10)

    def test_runner_requires_positive_charge(self):
        with self.assertRaises(IsotopeError):
            run({"elements": TABLE, "charge": 0, "formulas": "H2O"})


class CliTests(unittest.TestCase):
    def test_end_to_end(self):
        root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        with tempfile.TemporaryDirectory() as tmp:
            proc = subprocess.run(
                [sys.executable, "-m", "isotopes",
                 os.path.join(root, "examples", "input.json"),
                 "--out", tmp],
                cwd=root, capture_output=True, text=True,
            )
            self.assertEqual(proc.returncode, 0, proc.stderr)
            with open(os.path.join(tmp, "summary.json"), encoding="utf-8") as fh:
                report = json.load(fh)
            self.assertEqual(report["molecule_count"], 3)
            names = [m["formula"] for m in report["molecules"]]
            self.assertIn("Ca(OH)2", names)
            for molecule in report["molecules"]:
                path = os.path.join(tmp, molecule["csv"])
                self.assertTrue(os.path.isfile(path))


if __name__ == "__main__":
    unittest.main(verbosity=2)
