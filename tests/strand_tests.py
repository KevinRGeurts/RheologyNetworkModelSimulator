"""
This module contains unit tests for the Strand, FENEStrand, and FENSStrand classes, as well as the beta function.
"""


# Standard imports
import unittest

# Local imports
from RheologyNetworkModelSimulator.strand import Strand, FENEStrand, beta, FENSStrand


class Test_beta(unittest.TestCase):
    def test_beta(self):
        z, w = 1.5, 2.5
        exp_val = 0.1963495  # Calculated by R
        act_val = beta(z, w)
        self.assertAlmostEqual(exp_val, act_val)


class Test_Strand(unittest.TestCase):
    def test_init(self):
        s = Strand(0.1, 0.2, 0.3)
        exp_val = (0.1, 0.2, 0.3)
        act_val = (s.qx, s.qy, s.qz)
        self.assertTupleEqual(exp_val, act_val)

    def test_property_getter_max_length(self):
        exp_val = None
        s = Strand(0.1, 0.2, 0.3)
        act_val = s.max_length
        self.assertEqual(exp_val, act_val)

    def test_property_getters_qx_qy_qz(self):
        exp_val = (0.1, 0.2, 0.3)
        s = Strand(0.1, 0.2, 0.3)
        act_val = (s.qx, s.qy, s.qz)
        self.assertTupleEqual(exp_val, act_val)

    def test_property_setter_qx(self):
        exp_val = 0.4
        s = Strand(0.1, 0.2, 0.3)
        s.qx = 0.4
        act_val = s.qx
        self.assertEqual(exp_val, act_val)

    def test_property_setter_qy(self):
        exp_val = 0.4
        s = Strand(0.1, 0.2, 0.3)
        s.qy = 0.4
        act_val = s.qy
        self.assertEqual(exp_val, act_val)

    def test_property_setter_qz(self):
        exp_val = 0.4
        s = Strand(0.1, 0.2, 0.3)
        s.qz = 0.4
        act_val = s.qz
        self.assertEqual(exp_val, act_val)

    def test_get_qs(self):
        exp_val = (0.1, 0.2, 0.3)
        s = Strand(0.1, 0.2, 0.3)
        act_val = s.get_qs()
        self.assertTupleEqual(exp_val, act_val)

    def test_str_len_sqr(self):
        exp_val = 0.14
        s = Strand(0.1, 0.2, 0.3)
        act_val = s.str_len_sqr()
        self.assertEqual(exp_val, act_val)

    def test_stress_XX(self):
        s = Strand(0.1, 0.2, 0.3)
        self.assertRaises(NotImplementedError, s.stress_XX)

    def test_stress_YY(self):
        s = Strand(0.1, 0.2, 0.3)
        self.assertRaises(NotImplementedError, s.stress_YY)

    def test_stress_ZZ(self):
        s = Strand(0.1, 0.2, 0.3)
        self.assertRaises(NotImplementedError, s.stress_ZZ)

    def test_stress_YX(self):
        s = Strand(0.1, 0.2, 0.3)
        self.assertRaises(NotImplementedError, s.stress_YX)

    def test_loss_rate(self):
        s = Strand(0.1, 0.2, 0.3)
        self.assertRaises(NotImplementedError, s.loss_rate)

    def test_loss_prob(self):
        s = Strand(0.1, 0.2, 0.3)
        self.assertRaises(NotImplementedError, s.loss_prob, 0.001)

    def test_generate_eq_ensemble(self):
        s = Strand(0.1, 0.2, 0.3)
        self.assertRaises(NotImplementedError, s.generate_eq_ensemble, 10)


class Test_FENEStrand(unittest.TestCase):
    def test_init(self):
        s = FENEStrand(qx=0.1, qy=0.2, qz=0.3, b=100.0, n=2.0, mm=1.0)
        exp_val = (0.1, 0.2, 0.3)
        act_val = (s.qx, s.qy, s.qz)
        self.assertTupleEqual(exp_val, act_val)
        self.assertEqual(100.0, s.b)
        self.assertEqual(2.0, s.n)
        self.assertEqual(1.0, s.mm)

    def test_strand_init_to_long(self):
        self.assertRaises(AssertionError, FENEStrand, 1.0, 1.0, 1.0)

    def test_property_getter_max_length(self):
        exp_val = 1.0
        s = FENEStrand()
        act_val = s.max_length
        self.assertEqual(exp_val, act_val)

    def test_property_setter_qx_bad(self):
        s = FENEStrand(0.1, 0.2, 0.3)
        self.assertRaises(AssertionError, setattr, s, 'qx', 0.94)

    def test_property_setter_qy_bad(self):
        s = FENEStrand(0.1, 0.2, 0.3)
        self.assertRaises(AssertionError, setattr, s, 'qy', 0.95)

    def test_property_setter_qz_bad(self):
        s = FENEStrand(0.1, 0.2, 0.3)
        self.assertRaises(AssertionError, setattr, s, 'qz', 0.98)

    def test_stress_XX(self):
        exp_val = 1.1627906976744187
        s = FENEStrand(0.1, 0.2, 0.3)
        act_val = s.stress_XX()
        self.assertEqual(exp_val, act_val)

    def test_stress_YY(self):
        exp_val = 4.651162790697675
        s = FENEStrand(0.1, 0.2, 0.3)
        act_val = s.stress_YY()
        self.assertEqual(exp_val, act_val)

    def test_stress_ZZ(self):
        exp_val = 10.465116279069768
        s = FENEStrand(0.1, 0.2, 0.3)
        act_val = s.stress_ZZ()
        self.assertEqual(exp_val, act_val)

    def test_stress_YX(self):
        exp_val = 2.3255813953488373
        s = FENEStrand(0.1, 0.2, 0.3)
        act_val = s.stress_YX()
        self.assertEqual(exp_val, act_val)

    def test_loss_rate(self):
        exp_val = 1.1627906976744187
        s = FENEStrand(0.1, 0.2, 0.3)
        act_val = s.loss_rate()
        self.assertEqual(exp_val, act_val)

    def test_loss_prob(self):
        exp_val = 0.001162114918526358
        s = FENEStrand(0.1, 0.2, 0.3)
        act_val = s.loss_prob(0.001)
        self.assertEqual(exp_val, act_val)

    def test_generate_eq_ensemble(self):
        from random import seed
        exp_val = (0.1314922734575693, 0.03212588679044972, 0.029485861493129577)
        seed(1234567890)
        e = FENEStrand().generate_eq_ensemble(10)
        self.assertEqual(10, len(e))
        act_val = (e[0].qx, e[0].qy, e[0].qz)
        self.assertTupleEqual(exp_val, act_val)


class Test_FENSStrand(unittest.TestCase):
    def test_init(self):
        s = FENSStrand(qx=0.1, qy=0.2, qz=0.3)
        exp_val = (0.1, 0.2, 0.3)
        act_val = (s.qx, s.qy, s.qz)
        self.assertTupleEqual(exp_val, act_val)

    def test_strand_init_to_long(self):
        self.assertRaises(AssertionError, FENSStrand, 100.0, 100.0, 100.0)

    def test_property_getter_max_length(self):
        exp_val = 100.0
        s = FENSStrand()
        act_val = s.max_length
        self.assertEqual(exp_val, act_val)

    def test_property_setter_qx_bad(self):
        s = FENSStrand(10.0, 20.0, 30.0)
        self.assertRaises(AssertionError, setattr, s, 'qx', 94.0)

    def test_property_setter_qy_bad(self):
        s = FENSStrand(10.0, 20.0, 30.0)
        self.assertRaises(AssertionError, setattr, s, 'qy', 95.0)

    def test_property_setter_qz_bad(self):
        s = FENSStrand(10.0, 20.0, 30.0)
        self.assertRaises(AssertionError, setattr, s, 'qz', 98.0)

    def test_stress_XX(self):
        exp_val = 0.1*0.1
        s = FENSStrand(0.1, 0.2, 0.3)
        act_val = s.stress_XX()
        self.assertEqual(exp_val, act_val)

    def test_stress_YY(self):
        exp_val = 0.2*0.2
        s = FENSStrand(0.1, 0.2, 0.3)
        act_val = s.stress_YY()
        self.assertEqual(exp_val, act_val)

    def test_stress_ZZ(self):
        exp_val = 0.3*0.3
        s = FENSStrand(0.1, 0.2, 0.3)
        act_val = s.stress_ZZ()
        self.assertEqual(exp_val, act_val)

    def test_stress_YX(self):
        exp_val = 0.1*0.2
        s = FENSStrand(0.1, 0.2, 0.3)
        act_val = s.stress_YX()
        self.assertEqual(exp_val, act_val)

    def test_loss_rate(self):
        exp_val = 1.0
        s = FENSStrand(0.1, 0.2, 0.3)
        act_val = s.loss_rate()
        self.assertEqual(exp_val, act_val)

    def test_loss_prob(self):
        exp_val = 0.001
        s = FENSStrand(0.1, 0.2, 0.3)
        act_val = s.loss_prob(0.001)
        # TODO: Investigate why this needs almost equal
        self.assertAlmostEqual(exp_val, act_val, 5)

    def test_generate_eq_ensemble(self):
        from random import seed
        exp_val = (-0.2686617612333387, 0.44611821880006447, 0.17019087101967528)
        seed(1234567890)
        e = FENSStrand().generate_eq_ensemble(10)
        self.assertEqual(10, len(e))
        act_val = (e[0].qx, e[0].qy, e[0].qz)
        self.assertTupleEqual(exp_val, act_val)


if __name__ == '__main__':
    unittest.main()
