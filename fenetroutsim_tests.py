import unittest
from fenetroutsim import *
from random import seed


class MyTestCase(unittest.TestCase):
    def test_move_strand(self):
        exp_val = (0.09570000000000001, 0.19140000000000001, 0.3258)
        sim = FeneTroutSim({})
        s = Strand(0.1, 0.2, 0.3)
        sim.move_strand(s, 0.001, 100.0, 2)
        act_val = (s.qx, s.qy, s.qz)
        self.assertTupleEqual(exp_val, act_val)

    def test_run_sim(self):
        seed(1234567890)

        exp_val = {'Qave': 0.16110272799636224, 'n/no': 0.99, 'piYX': 0.07245743211970784, 'trouton 1': 0.04737104142315748, 'trouton 2': 0.0021919416527666088}

        sim_input = {}
        sim_input['gamdot'] = 10.0
        sim_input['begstrand'] = 100
        sim_input['eps'] = 0.001
        sim_input['steps'] = int(0.01 / sim_input['eps'])
        sim_input['b'] = 100.0  # All results in J. Chem. Phys article
        sim_input['n'] = 2  # Don't change
        sim_input['mm'] = 1  # Don't change
        sim_input['outfile'] = 'testfenetrout.out'

        sim = FeneTroutSim(sim_input)
        act_val = sim.run_sim()
        print(act_val)
        self.assertAlmostEqual(exp_val['Qave'], act_val.Q_ave, 15)
        self.assertAlmostEqual(exp_val['n/no'], act_val.n_ave, 15)
        self.assertAlmostEqual(exp_val['piYX'], act_val.piYX_ave,15)
        self.assertAlmostEqual(exp_val['trouton 1'], act_val.trout1_ave, 15)
        self.assertAlmostEqual(exp_val['trouton 2'], act_val.trout2_ave, 15)


if __name__ == '__main__':
    unittest.main()
