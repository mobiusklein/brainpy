import unittest

from brainpy.composition import _parse_formula as parse_formula, calculate_mass
from brainpy import _has_c

if _has_c:
    from brainpy._c.composition import parse_formula as cparse_formula


class CompositionTest(unittest.TestCase):
    def test_parse(self):
        formula = "C6H12O6"
        composition = parse_formula(formula)
        self.assertEqual(composition, {"C": 6, "H": 12, "O": 6})
        self.assertEqual(composition * 2, {"C": 6 * 2, "H": 12 * 2, "O": 6 * 2})
        if _has_c:
            composition = cparse_formula(formula)
            self.assertEqual(composition, {"C": 6, "H": 12, "O": 6})
            self.assertEqual(composition * 2, {"C": 6 * 2, "H": 12 * 2, "O": 6 * 2})

    def test_isotope_parse(self):
        for formula in ["O1H1H1H[2]1", "O1H1H[2]1H1"]:
            composition = parse_formula(formula)
            self.assertEqual(composition, {"O": 1, "H": 2, "H[2]": 1})
            self.assertAlmostEqual(composition.mass(), 20.0246, 3)
            self.assertEqual(composition['H'], 2)
            if _has_c:
                composition = cparse_formula(formula)
                self.assertEqual(composition, {"O": 1, "H": 2, "H[2]": 1})
                self.assertAlmostEqual(composition.mass(), 20.0246, 3)
                self.assertEqual(composition['H'], 2)

    def test_mass(self):
        formula = "C6H12O6"
        composition = parse_formula(formula)
        self.assertAlmostEqual(composition.mass(), calculate_mass({"C": 6, "H": 12, "O": 6}))
        if _has_c:
            composition = cparse_formula(formula)
            self.assertAlmostEqual(composition.mass(), calculate_mass({"C": 6, "H": 12, "O": 6}))
