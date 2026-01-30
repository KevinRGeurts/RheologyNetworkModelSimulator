import unittest

from RheologyNetworkModelSimulator.strand import Strand, beta

class MyTestCase(unittest.TestCase):
    def test_strand_init_to_long(self):
        self.assertRaises(AssertionError, Strand, 1.0, 1.0, 1.0)

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

    def test_property_setter_qx_bad(self):
        s = Strand(0.1, 0.2, 0.3)
        self.assertRaises(AssertionError, setattr, s, 'qx', 0.94)

    def test_property_setter_qy(self):
        exp_val = 0.4
        s = Strand(0.1, 0.2, 0.3)
        s.qy = 0.4
        act_val = s.qy
        self.assertEqual(exp_val, act_val)

    def test_property_setter_qy_bad(self):
        s = Strand(0.1, 0.2, 0.3)
        self.assertRaises(AssertionError, setattr, s, 'qy', 0.95)

    def test_property_setter_qz(self):
        exp_val = 0.4
        s = Strand(0.1, 0.2, 0.3)
        s.qz = 0.4
        act_val = s.qz
        self.assertEqual(exp_val, act_val)

    def test_property_setter_qz_bad(self):
        s = Strand(0.1, 0.2, 0.3)
        self.assertRaises(AssertionError, setattr, s, 'qz', 0.98)

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

    def test_beta(self):
        z, w = 1.5, 2.5
        exp_val = 0.1963495  # Calculated by R
        act_val = beta(z, w)
        self.assertAlmostEqual(exp_val, act_val)

    def test_stress_XX(self):
        exp_val = 1.1627906976744187
        s = Strand(0.1, 0.2, 0.3)
        act_val = s.stress_XX(100.0)
        self.assertEqual(exp_val, act_val)

    def test_stress_YY(self):
        exp_val = 4.651162790697675
        s = Strand(0.1, 0.2, 0.3)
        act_val = s.stress_YY(100.0)
        self.assertEqual(exp_val, act_val)

    def test_stress_ZZ(self):
        exp_val = 10.465116279069768
        s = Strand(0.1, 0.2, 0.3)
        act_val = s.stress_ZZ(100.0)
        self.assertEqual(exp_val, act_val)

    def test_stress_YX(self):
        exp_val = 2.3255813953488373
        s = Strand(0.1, 0.2, 0.3)
        act_val = s.stress_YX(100.0)
        self.assertEqual(exp_val, act_val)

    def test_loss_rate(self):
        exp_val = 1.1627906976744187
        s = Strand(0.1, 0.2, 0.3)
        act_val = s.loss_rate()
        self.assertEqual(exp_val, act_val)

    def test_loss_prob(self):
        exp_val = 0.001162114918526358
        s = Strand(0.1, 0.2, 0.3)
        act_val = s.loss_prob(0.001)
        self.assertEqual(exp_val, act_val)

    def test_generate_eq_ensemble(self):
        from random import seed
        exp_val = (0.1314922734575693, 0.03212588679044972, 0.029485861493129577)
        seed(1234567890)
        e = Strand().generate_eq_ensemble(10, 100.)
        self.assertEqual(10, len(e))
        act_val = (e[0].qx, e[0].qy, e[0].qz)
        self.assertTupleEqual(exp_val, act_val)

    def test_ensemble_stress(self):
        exp_val = (2.*1.1627906976744187, 2.*4.651162790697675, 2.*10.465116279069768, 2.*2.3255813953488373)
        e=(Strand(0.1, 0.2, 0.3), Strand(0.1, 0.2, 0.3))
        act_val = Strand().ensemble_stress(e, 100.0)
        self.assertTupleEqual(exp_val, act_val)


if __name__ == '__main__':
    unittest.main()
