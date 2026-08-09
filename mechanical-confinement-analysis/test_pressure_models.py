import math
import unittest

import pressure_models as pm


class PressureModelTests(unittest.TestCase):
    def test_indentation_modulus(self):
        self.assertAlmostEqual(pm.indentation_modulus() / pm.GPA, 53.8505, places=4)

    def test_reference_geometry_2p4_nm(self):
        result = pm.evaluate_geometry(50.0, 2.4, 50.0)
        self.assertAlmostEqual(result["hamaker_empty_mpa"], 0.430324, places=6)
        self.assertAlmostEqual(result["hamaker_au_filled_mpa"], 0.356554, places=6)
        self.assertAlmostEqual(result["adhesion_scale_gpa"], 0.200833, places=6)
        self.assertAlmostEqual(result["flat_punch_two_compliant_gpa"], 1.645550, places=6)
        self.assertTrue(result["half_space_valid"])

    def test_reference_geometry_5p3_nm(self):
        result = pm.evaluate_geometry(50.0, 5.3, 50.0)
        self.assertAlmostEqual(result["hamaker_empty_mpa"], 0.039900, places=6)
        self.assertAlmostEqual(result["flat_punch_two_compliant_gpa"], 3.633923, places=6)

    def test_half_space_limit(self):
        gap = 2.5 * pm.NM
        finite = pm.hamaker_pressure(gap, 1_000 * pm.NM, pm.DEFAULTS.hamaker_empty_j)
        half_space = pm.half_space_hamaker_pressure(gap, pm.DEFAULTS.hamaker_empty_j)
        self.assertTrue(math.isclose(finite / half_space, 1.0, rel_tol=1e-7))

    def test_invalid_geometry(self):
        with self.assertRaises(ValueError):
            pm.finite_slab_kernel(0.0, 50 * pm.NM)


if __name__ == "__main__":
    unittest.main()
