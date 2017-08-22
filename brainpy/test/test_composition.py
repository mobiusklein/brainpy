import unittest

from brainpy.composition import parse_formula, SimpleComposition, calculate_mass
from brainpy import _has_c

if _has_c:
    from brainpy._c.composition import PyComposition


class CompositionTest(unittest.TestCase):
    def test_parse(self):
        formula = "C6H12O6"
        composition = parse_formula(formula, SimpleComposition)
        self.assertEqual(composition, {"C": 6, "H": 12, "O": 6})
        self.assertEqual(composition * 2, {"C": 6 * 2, "H": 12 * 2, "O": 6 * 2})
        if _has_c:
            composition = parse_formula(formula, PyComposition)
            self.assertEqual(composition, {"C": 6, "H": 12, "O": 6})
            self.assertEqual(composition * 2, {"C": 6 * 2, "H": 12 * 2, "O": 6 * 2})

    def test_mass(self):
        formula = "C6H12O6"
        composition = parse_formula(formula)
        self.assertAlmostEqual(composition.mass(), calculate_mass({"C": 6, "H": 12, "O": 6}))
        if _has_c:
            composition = parse_formula(formula, PyComposition)
            self.assertAlmostEqual(composition.mass(), calculate_mass({"C": 6, "H": 12, "O": 6}))
