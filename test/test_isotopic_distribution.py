import unittest


from brainpy import IsotopicDistribution, isotopic_variants, calculate_mass, neutral_mass, _has_c, Peak

if _has_c:
    from brainpy.brainpy import _IsotopicDistribution


class TestIsotopicDistribution(unittest.TestCase):
    def test_neutral_mass(self):
        hexnac = {'H': 13, 'C': 8, 'O': 5, 'N': 1}
        dist = isotopic_variants(hexnac)

        reference = [
            Peak(mz=203.079373, intensity=0.901656, charge=0),
            Peak(mz=204.082545, intensity=0.084376, charge=0),
            Peak(mz=205.085599, intensity=0.012998, charge=0),
            Peak(mz=206.088782, intensity=0.000969, charge=0)
        ]
        for inst, ref in zip(dist, reference):
            self.assertAlmostEqual(inst.mz, ref.mz, 3)
            self.assertAlmostEqual(inst.intensity, ref.intensity, 3)

    if _has_c:
        def test_pure_python(self):
            hexnac = {'H': 13, 'C': 8, 'O': 5, 'N': 1}
            dist = _IsotopicDistribution(hexnac, 4).aggregated_isotopic_variants()

            reference = [
                Peak(mz=203.079373, intensity=0.901656, charge=0),
                Peak(mz=204.082545, intensity=0.084376, charge=0),
                Peak(mz=205.085599, intensity=0.012998, charge=0),
                Peak(mz=206.088782, intensity=0.000969, charge=0)
            ]
            for inst, ref in zip(dist, reference):
                self.assertAlmostEqual(inst.mz, ref.mz, 3)
                self.assertAlmostEqual(inst.intensity, ref.intensity, 3)


if __name__ == '__main__':
    unittest.main()
