import csv
import json
import math
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from isotopes.core import compute_cluster, summarize
from isotopes.errors import IsotopeError
from isotopes.formula import parse_formula
from isotopes.table import validate_table

ELEMENTS = {
    "H": [
        {"mass_number": 1, "exact_mass": 1.007825, "abundance": 0.999885},
        {"mass_number": 2, "exact_mass": 2.014102, "abundance": 0.000115},
    ],
    "C": [
        {"mass_number": 12, "exact_mass": 12.0, "abundance": 0.9893},
        {"mass_number": 13, "exact_mass": 13.003355, "abundance": 0.0107},
    ],
    "O": [
        {"mass_number": 16, "exact_mass": 15.994915, "abundance": 0.99757},
        {"mass_number": 17, "exact_mass": 16.999132, "abundance": 0.00038},
        {"mass_number": 18, "exact_mass": 17.999160, "abundance": 0.00205},
    ],
}


def make_input(formula="CH4", charge=1, elements=None):
    return {"elements": elements or ELEMENTS, "formula": formula, "charge": charge}


class TestFormula(unittest.TestCase):
    def test_simple(self):
        self.assertEqual(parse_formula("H2O"), {"H": 2, "O": 1})

    def test_nested_parens(self):
        self.assertEqual(
            parse_formula("K4(ON(SO3)2)2"),
            {"K": 4, "O": 14, "N": 2, "S": 4},
        )

    def test_merge(self):
        self.assertEqual(parse_formula("CH3COOH"), {"C": 2, "H": 4, "O": 2})

    def test_errors(self):
        for bad in ["", "()", "H0", "H2O)", "(H2O", "2H", "H2.5O", "H2O+"]:
            with self.assertRaises(IsotopeError, msg=bad):
                parse_formula(bad)


class TestTable(unittest.TestCase):
    def test_normalize(self):
        raw = {"X": [{"mass_number": 10, "exact_mass": 10.0, "abundance": 2.0}]}
        with self.assertRaises(IsotopeError):
            validate_table(raw)

    def test_duplicate_mass_number(self):
        raw = {"X": [
            {"mass_number": 10, "exact_mass": 10.0, "abundance": 0.5},
            {"mass_number": 10, "exact_mass": 10.1, "abundance": 0.5},
        ]}
        with self.assertRaises(IsotopeError):
            validate_table(raw)

    def test_bad_values(self):
        for ab in [float("nan"), float("inf"), -0.1]:
            raw = {"X": [{"mass_number": 10, "exact_mass": 10.0, "abundance": ab}]}
            with self.assertRaises(IsotopeError, msg=str(ab)):
                validate_table(raw)


class TestCore(unittest.TestCase):
    def setUp(self):
        self.table = validate_table(ELEMENTS)

    def test_single_atom(self):
        cluster = compute_cluster({"C": 1}, self.table)
        self.assertAlmostEqual(cluster[12][0], 0.9893)
        self.assertAlmostEqual(cluster[13][0], 0.0107)
        self.assertAlmostEqual(cluster[13][1] / cluster[13][0], 13.003355)

    def test_total_probability(self):
        cluster = compute_cluster({"C": 60}, self.table)
        self.assertAlmostEqual(sum(p for p, _ in cluster.values()), 1.0, places=12)

    def test_mz_charge(self):
        cluster = compute_cluster({"H": 2, "O": 1}, self.table)
        summary = summarize(cluster, 2)
        mono = summary["peaks"][0]
        self.assertEqual(mono["mass_number"], 18)
        self.assertAlmostEqual(mono["mz"], mono["mean_exact_mass"] / 2)
        self.assertAlmostEqual(summary["average_mz"], summary["average_mass"] / 2)

    def test_unknown_element(self):
        with self.assertRaises(IsotopeError):
            compute_cluster({"Zz": 1}, self.table)


class TestCli(unittest.TestCase):
    def run_cli(self, payload, expect_ok=True):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        base = Path(tmp.name)
        inp = base / "input.json"
        inp.write_text(json.dumps(payload), encoding="utf-8")
        out = base / "out"
        proc = subprocess.run(
            [sys.executable, "-m", "isotopes", str(inp), "--out", str(out)],
            cwd=ROOT, capture_output=True, text=True,
        )
        return proc, out

    def test_success(self):
        proc, out = self.run_cli(make_input("K4(ON(SO3)2)2", 2,
                                            elements={**ELEMENTS, **{
                "K": [{"mass_number": 39, "exact_mass": 38.963707, "abundance": 0.932581},
                      {"mass_number": 41, "exact_mass": 40.961826, "abundance": 0.067419}],
                "N": [{"mass_number": 14, "exact_mass": 14.003074, "abundance": 0.99636},
                      {"mass_number": 15, "exact_mass": 15.000109, "abundance": 0.00364}],
                "S": [{"mass_number": 32, "exact_mass": 31.972071, "abundance": 0.9499},
                      {"mass_number": 33, "exact_mass": 32.971459, "abundance": 0.0075},
                      {"mass_number": 34, "exact_mass": 33.967867, "abundance": 0.0425},
                      {"mass_number": 36, "exact_mass": 35.967081, "abundance": 0.0001}],
            }}))
        self.assertEqual(proc.returncode, 0, proc.stderr)
        summary = json.loads((out / "summary.json").read_text(encoding="utf-8"))
        self.assertEqual(summary["charge"], 2)
        self.assertAlmostEqual(summary["total_probability"], 1.0, places=9)
        rows = list(csv.DictReader((out / "peaks.csv").open(encoding="utf-8")))
        self.assertEqual(len(rows), summary["peak_count"])
        mass_numbers = [int(r["mass_number"]) for r in rows]
        self.assertEqual(mass_numbers, sorted(mass_numbers))
        self.assertTrue(all(float(r["probability"]) > 0 for r in rows))

    def test_failure_no_output(self):
        proc, out = self.run_cli(make_input("Xx2"))
        self.assertNotEqual(proc.returncode, 0)
        self.assertIn("Xx", proc.stderr)
        self.assertFalse((out / "summary.json").exists())

    def test_bad_charge(self):
        proc, out = self.run_cli(make_input("H2O", charge=0))
        self.assertNotEqual(proc.returncode, 0)
        self.assertFalse((out / "summary.json").exists())


if __name__ == "__main__":
    unittest.main()
