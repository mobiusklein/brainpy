import unittest

from brainpy.composition import _parse_formula as parse_formula, calculate_mass
from brainpy import _has_c

if _has_c:
    from brainpy._c.composition import parse_formula as cparse_formula

try:
    import faulthandler
    faulthandler.enable()
except ImportError:
    pass


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
            self.assertEqual(composition, {"O": 1, "H": 2, "H[2]": 1}, msg="%r != %r" % (
                composition, {"O": 1, "H": 2, "H[2]": 1}))
            self.assertAlmostEqual(composition.mass(), 20.0246, 3, msg="%r.mass() (%r) != %r" % (
                composition, composition.mass(), 20.0246))
            self.assertEqual(composition['H'], 2)
            if _has_c:
                composition = cparse_formula(formula)
                self.assertEqual(composition, {"O": 1, "H": 2, "H[2]": 1}, msg="%r != %r" % (
                    composition, {"O": 1, "H": 2, "H[2]": 1}))
                self.assertAlmostEqual(composition.mass(), 20.0246, 3, msg="%r.mass() (%r) != %r" % (
                    composition, composition.mass(), 20.0246))
                self.assertEqual(composition['H'], 2)

    def test_mass(self):
        formula = "C6H12O6"
        composition = parse_formula(formula)
        self.assertEqual(composition['C'], 6)
        self.assertEqual(composition['H'], 12)
        self.assertEqual(composition['O'], 6)
        self.assertAlmostEqual(composition.mass(), calculate_mass({"C": 6, "H": 12, "O": 6}),
                               msg="%s did not have the right mass (%f)" % (composition, composition.mass()))
        if _has_c:
            ccomposition = cparse_formula(formula)
            self.assertEqual(ccomposition['C'], 6)
            self.assertEqual(ccomposition['H'], 12)
            self.assertEqual(ccomposition['O'], 6)
            self.assertAlmostEqual(ccomposition.mass(), calculate_mass({"C": 6, "H": 12, "O": 6}),
                                   msg="%s did not have the right mass (%f)" % (ccomposition, ccomposition.mass()))


if __name__ == '__main__':
    unittest.main()
