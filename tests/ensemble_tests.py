import unittest

from RheologyNetworkModelSimulator.ensemble import Ensemble
from RheologyNetworkModelSimulator.strand import Strand

class Test_Ensemble(unittest.TestCase):
    def test_q_ave(self):
        exp_val = 0.37416573867739417
        e = Ensemble([Strand(0.1, 0.2, 0.3), Strand(0.1, 0.2, 0.3)])
        act_val = e.q_ave()
        self.assertEqual(exp_val, act_val)

    def test_len(self):
        e = Ensemble([Strand(0.1, 0.2, 0.3), Strand(0.1, 0.2, 0.3)])
        exp_val = 2
        act_val = len(e)
        self.assertEqual(exp_val, act_val)

    def test_getitem(self):
        e = Ensemble([Strand(0.1, 0.2, 0.3), Strand(0.2, 0.3, 0.1)])
        exp_val = (0.2, 0.3, 0.1)
        act_val = e[1].get_qs()
        self.assertTupleEqual(exp_val, act_val)

    def test_iter(self):
        e = Ensemble([Strand(0.1, 0.2, 0.3), Strand(0.2, 0.3, 0.1)])
        qs_list = []
        for s in e:
            qs_list.append(s.get_qs())
        exp_val = [(0.1, 0.2, 0.3), (0.2, 0.3, 0.1)]
        self.assertListEqual(exp_val, qs_list)

    def test_remove(self):
        e = Ensemble([Strand(0.1, 0.2, 0.3), Strand(0.2, 0.3, 0.1)])
        s_to_remove = e[1]
        e.remove(s_to_remove)
        exp_val = 1
        act_val = len(e)
        self.assertEqual(exp_val, act_val)
        exp_qs = (0.1, 0.2, 0.3)
        act_qs = e[0].get_qs()
        self.assertTupleEqual(exp_qs, act_qs)

    def test_append(self):
        e = Ensemble([Strand(0.1, 0.2, 0.3), Strand(0.2, 0.3, 0.1)])
        new_strand = Strand(0.3, 0.1, 0.2)
        e.append(new_strand)
        exp_val = 3
        act_val = len(e)
        self.assertEqual(exp_val, act_val)
        exp_qs = (0.3, 0.1, 0.2)
        act_qs = e[2].get_qs()
        self.assertTupleEqual(exp_qs, act_qs)

    def test_property_strands(self):
        exp_val = [Strand(0.1, 0.2, 0.3), Strand(0.2, 0.3, 0.1)]
        e = Ensemble(exp_val)
        act_val = e.strands
        self.assertListEqual(exp_val, act_val)
        # Remove a strand from the copied list returned by the property and check that the ensemble is unchanged
        act_val.remove(act_val[0])
        self.assertTrue(len(e) == 2) # Length is still 2
        self.assertListEqual(exp_val, e.strands) # Strands are unchanged


if __name__ == '__main__':
    unittest.main()
