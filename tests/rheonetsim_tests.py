"""
This module contains unit tests for:
    (1) RheoNetSim and ElongateNetSim classes
    (2) Their associated output classes, RheoNetSimOut and ElongateNetSimOut
    (3) The move_strand_fene_elongational and move_strand_fens_elongational functions
"""


# Standard imports
import unittest
from random import seed
from array import array
from math import sqrt

# Local imports
from RheologyNetworkModelSimulator.rheonetsim import ElongateNetSim, RheoNetSimOut, ElongateNetSimOut, move_strand_fene_elongational, RheoNetSim, move_strand_fens_elongational
from RheologyNetworkModelSimulator.rheonetsim import move_strand_fens_shear, move_strand_fene_shear, ShearNetSim, ShearNetSimOut
from RheologyNetworkModelSimulator.strand import FENEStrand, FENSStrand
from RheologyNetworkModelSimulator.ensemble import Ensemble
from RheologyNetworkModelSimulator.qplot import TextPlot


class Test_RheoNetSimOut(unittest.TestCase):
    def test_str(self):
        sim_out = RheoNetSimOut()
        sim_out.n_vals.append(10000)
        sim_out.Q_vals.append(0.87)
        sim_out.piYX_vals.append(10.0)
        exp_val = f"{sim_out.n_vals[0]}, {sim_out.piYX_vals[0]}, {sim_out.Q_vals[0]}"
        act_val = str(sim_out)
        self.assertEqual(exp_val, act_val)

    def test_time_step_output(self):
        sim_out = RheoNetSimOut()
        sim_out.n_vals.append(10000)
        sim_out.Q_vals.append(0.87)
        sim_out.piYX_vals.append(10.0)
        exp_val = f"n: {sim_out.n_vals[0]}, piYX: {sim_out.piYX_vals[0]}, Q: {sim_out.Q_vals[0]}"
        act_val = sim_out.time_step_output()
        self.assertEqual(exp_val, act_val)

    def test_final_output(self):
        sim_out = RheoNetSimOut()
        sim_out.n_ave = 0.95
        sim_out.Q_ave = 0.87
        sim_out.piYX_ave = 10.0
        exp_val = f"Avg n: {sim_out.n_ave}\nAvg piYX: {sim_out.piYX_ave}\nAvg Q: {sim_out.Q_ave}\n"
        act_val = sim_out.final_output()
        self.assertEqual(exp_val, act_val)

class Test_ElongateNetSimOut(unittest.TestCase):
    def test_str(self):
        sim_out = ElongateNetSimOut()
        sim_out.n_vals.append(10000)
        sim_out.Q_vals.append(0.87)
        sim_out.piYX_vals.append(10.0)
        sim_out.trout1_vals.append(8.5)
        sim_out.trout2_vals.append(0.2)
        sim_out.trout1_plt=None
        exp_val = f"{sim_out.n_vals[0]}, {sim_out.piYX_vals[0]}, {sim_out.Q_vals[0]}, {sim_out.trout1_vals[0]}, {sim_out.trout2_vals[0]}"
        act_val = str(sim_out)
        self.assertEqual(exp_val, act_val)

    def test_time_step_output(self):
        sim_out = ElongateNetSimOut()
        sim_out.n_vals.append(10000)
        sim_out.Q_vals.append(0.87)
        sim_out.piYX_vals.append(10.0)
        sim_out.trout1_vals.append(8.5)
        sim_out.trout2_vals.append(0.2)
        sim_out.trout1_plt=None
        exp_val = f"n: {sim_out.n_vals[0]}, piYX: {sim_out.piYX_vals[0]}, Q: {sim_out.Q_vals[0]}, Trout 1: {sim_out.trout1_vals[0]}, Trout 2: {sim_out.trout2_vals[0]}"
        act_val = sim_out.time_step_output()
        self.assertEqual(exp_val, act_val)

    def test_final_output(self):
        sim_out = ElongateNetSimOut()
        sim_out.trout1_vals.append(8.5)
        sim_out.trout1_vals.append(7.0)
        sim_out.time_vals.append(0.2)
        sim_out.time_vals.append(0.3)
        sim_out.n_ave = 0.95
        sim_out.Q_ave = 0.87
        sim_out.piYX_ave = 10.0
        sim_out.trout1_ave = 8.5
        sim_out.trout2_ave = 0.2

        # Make a TextPlot object for plotting elongational viscosity vs time
        _symbol = array('u')
        _symbol.append('O')
        _npts=len(sim_out.trout1_vals)
        _x=[sim_out.time_vals]
        _y=[sim_out.trout1_vals]
        _titl1='Start-up of Elongational Flow Simulation'
        _titl2='O: Elongational Viscosity vs Time'
        sim_out.trout1_plt = TextPlot(ncur=1, npts=_npts, x=_x, y=_y, symbol=_symbol, titl1=_titl1, titl2=_titl2)

        sim_out.trout1_plt=None
        exp_val = f"Avg n: {sim_out.n_ave}\nAvg piYX: {sim_out.piYX_ave}\nAvg Q: {sim_out.Q_ave}\nAvg Trout 1: {sim_out.trout1_ave}\nAvg Trout 2: {sim_out.trout2_ave}\n\n{sim_out.trout1_plt}"
        act_val = sim_out.final_output()
        self.assertEqual(exp_val, act_val)


class Test_ShearNetSimOut(unittest.TestCase):
    def test_str(self):
        sim_out = ShearNetSimOut()
        sim_out.n_vals.append(10000)
        sim_out.Q_vals.append(0.87)
        sim_out.piYX_vals.append(10.0)
        sim_out.viscosity_vals.append(8.5)
        sim_out.psi1_vals.append(0.2)
        sim_out.psi2_vals.append(0.5)
        sim_out.viscosity_plt=None
        exp_val = f"{sim_out.n_vals[0]}, {sim_out.piYX_vals[0]}, {sim_out.Q_vals[0]}, {sim_out.viscosity_vals[0]}, {sim_out.psi1_vals[0]}, {sim_out.psi2_vals[0]}"
        act_val = str(sim_out)
        self.assertEqual(exp_val, act_val)

    def test_time_step_output(self):
        sim_out = ShearNetSimOut()
        sim_out.n_vals.append(10000)
        sim_out.Q_vals.append(0.87)
        sim_out.piYX_vals.append(10.0)
        sim_out.viscosity_vals.append(8.5)
        sim_out.psi1_vals.append(0.2)
        sim_out.psi2_vals.append(0.5)
        sim_out.viscosity_plt=None
        exp_val = f"n: {sim_out.n_vals[0]}, piYX: {sim_out.piYX_vals[0]}, Q: {sim_out.Q_vals[0]}, Viscosity: {sim_out.viscosity_vals[0]}, First Normal Stress Coefficient: {sim_out.psi1_vals[0]}, Second Normal Stress Coefficient: {sim_out.psi2_vals[0]}"
        act_val = sim_out.time_step_output()
        self.assertEqual(exp_val, act_val)

    def test_final_output(self):
        sim_out = ShearNetSimOut()
        sim_out.viscosity_vals.append(8.5)
        sim_out.viscosity_vals.append(7.0)
        sim_out.time_vals.append(0.2)
        sim_out.time_vals.append(0.3)
        sim_out.n_ave = 0.95
        sim_out.Q_ave = 0.87
        sim_out.piYX_ave = 10.0
        sim_out.viscosity_ave = 8.5
        sim_out.psi1_ave = 0.2
        sim_out.psi2_ave = 0.3

        # Make a TextPlot object for plotting elongational viscosity vs time
        _symbol = array('u')
        _symbol.append('O')
        _npts=len(sim_out.viscosity_vals)
        _x=[sim_out.time_vals]
        _y=[sim_out.viscosity_vals]
        _titl1='Start-up of Shear Flow Simulation'
        _titl2='O: Shear Viscosity vs Time'
        sim_out.viscosity_plt = TextPlot(ncur=1, npts=_npts, x=_x, y=_y, symbol=_symbol, titl1=_titl1, titl2=_titl2)

        sim_out.viscosity_plt=None
        exp_val = f"Avg n: {sim_out.n_ave}\nAvg piYX: {sim_out.piYX_ave}\nAvg Q: {sim_out.Q_ave}\nAvg Viscosity: {sim_out.viscosity_ave}\nAvg First Normal Stress Coefficient: {sim_out.psi1_ave}\nAvg Second Normal Stress Coefficient: {sim_out.psi2_ave}\n\n{sim_out.viscosity_plt}"
        act_val = sim_out.final_output()
        self.assertEqual(exp_val, act_val)


class Test_RheoNetSim(unittest.TestCase):
    def test_init(self):
        # strand_mover just needs to be a callable, for this test, which is really testing that
        # _createSimOutObj raises NotImplementedError
        self.assertRaises(NotImplementedError, RheoNetSim, strand_mover=str)


class Test_ElongateNetSim(unittest.TestCase):
    def test_move_fene_strand_elongational(self):
        # def move_strand_fene_elongational(s, eps, gamdot):
        exp_val = (0.09570000000000001, 0.19140000000000001, 0.3258)
        s = FENEStrand(qx=0.1, qy=0.2, qz=0.3, b=100.0, n=2.0, mm=1.0)
        move_strand_fene_elongational(s, 0.001, 100.0)
        act_val = (s.qx, s.qy, s.qz)
        self.assertTupleEqual(exp_val, act_val)

    def test_move_fens_strand_elongational(self):
        # def move_strand_fene_elongational(s, eps, gamdot):
        exp_val = (0.09570000000000001, 0.19140000000000001, 0.3258)
        s = FENEStrand(qx=0.1, qy=0.2, qz=0.3)
        move_strand_fens_elongational(s, 0.001, 100.0)
        act_val = (s.qx, s.qy, s.qz)
        self.assertTupleEqual(exp_val, act_val)

    def test_init(self):
        s = FENEStrand()
        sim = ElongateNetSim({}, s, move_strand_fene_elongational)
        self.assertIsInstance(sim._sim_output, ElongateNetSimOut)

    def test_getOutFileHeader(self):
        s = FENEStrand()
        sim = ElongateNetSim({}, s, move_strand_fene_elongational)
        exp_val = "Time, Strands, <Q>, piYX, trouton_1, trouton_2"
        act_val = sim._getOutFileHeader()
        self.assertEqual(exp_val, act_val)

    def test_computeTimestepOutputs_before_ss(self):
        fene_b = 100.0 # All results in J. Chem. Phys article
        fene_n = 2.0 # Don't change
        fene_mm = 1.0 # Don't change
        sim_input = {}
        sim_input['gamdot'] = 10.0
        sim_input['begstrand'] = 100
        sim_input['eps'] = 0.001
        sim_input['steps'] = int(0.01 / sim_input['eps'])
        sim_input['outfile'] = 'testfenetrout.out'
        s1 = FENEStrand(qx=0.1, qy=0.2, qz=0.3, b=fene_b, n=fene_n, mm=fene_mm)
        s2 = FENEStrand(qx=0.3, qy=0.1, qz=0.2, b=fene_b, n=fene_n, mm=fene_mm)
        sim = ElongateNetSim(sim_input, s1, move_strand_fene_elongational)
        we = Ensemble([s1, s2])
        exp_val = ElongateNetSimOut()
        exp_val.n_vals.append(2)
        exp_val.piYX_vals.append(0.058139535)
        exp_val.Q_vals.append(0.374165739)
        exp_val.time_vals.append(0.007)
        exp_val.trout1_vals.append(0.003488372)
        exp_val.trout2_vals.append(-0.005813953)
        exp_val.n_ave = 0
        exp_val.piYX_ave = 0
        exp_val.Q_ave = 0
        exp_val.trout1_ave = 0
        exp_val.trout2_ave = 0
        sim._computeTimestepOutputs(we, 0.007)
        act_val = sim._sim_output
        self.assertAlmostEqual(exp_val.n_vals[0], act_val.n_vals[0], 9)
        self.assertAlmostEqual(exp_val.Q_vals[0], act_val.Q_vals[0], 9)
        self.assertAlmostEqual(exp_val.time_vals[0], act_val.time_vals[0], 9)
        self.assertAlmostEqual(exp_val.piYX_vals[0], act_val.piYX_vals[0], 9)
        self.assertAlmostEqual(exp_val.trout1_vals[0], act_val.trout1_vals[0], 9)
        self.assertAlmostEqual(exp_val.trout2_vals[0], act_val.trout2_vals[0], 9)
        self.assertAlmostEqual(exp_val.n_ave, act_val.n_ave, 9)
        self.assertAlmostEqual(exp_val.piYX_ave, act_val.piYX_ave, 9)
        self.assertAlmostEqual(exp_val.Q_ave, act_val.Q_ave, 9)
        self.assertAlmostEqual(exp_val.trout1_ave, act_val.trout1_ave, 9)
        self.assertAlmostEqual(exp_val.trout2_ave, act_val.trout2_ave, 9)

    def test_computeTimestepOutputs_after_ss(self):
        fene_b = 100.0 # All results in J. Chem. Phys article
        fene_n = 2.0 # Don't change
        fene_mm = 1.0 # Don't change
        sim_input = {}
        sim_input['gamdot'] = 10.0
        sim_input['begstrand'] = 100
        sim_input['eps'] = 0.001
        sim_input['steps'] = int(0.01 / sim_input['eps'])
        sim_input['outfile'] = 'testfenetrout.out'
        s1 = FENEStrand(qx=0.1, qy=0.2, qz=0.3, b=fene_b, n=fene_n, mm=fene_mm)
        s2 = FENEStrand(qx=0.3, qy=0.1, qz=0.2, b=fene_b, n=fene_n, mm=fene_mm)
        sim = ElongateNetSim(sim_input, s1, move_strand_fene_elongational)
        we = Ensemble([s1, s2])
        exp_val = ElongateNetSimOut()
        exp_val.n_vals.append(2)
        exp_val.piYX_vals.append(0.058139535)
        exp_val.Q_vals.append(0.374165739)
        exp_val.time_vals.append(0.008)
        exp_val.trout1_vals.append(0.003488372)
        exp_val.trout2_vals.append(-0.005813953)
        exp_val.n_ave = 2/sim_input['begstrand']
        exp_val.piYX_ave = 0.058139535
        exp_val.Q_ave = 0.374165739
        exp_val.trout1_ave = 0.003488372
        exp_val.trout2_ave = -0.005813953
        sim._computeTimestepOutputs(we, 0.008)
        act_val = sim._sim_output
        self.assertAlmostEqual(exp_val.n_vals[0], act_val.n_vals[0], 9)
        self.assertAlmostEqual(exp_val.Q_vals[0], act_val.Q_vals[0], 9)
        self.assertAlmostEqual(exp_val.time_vals[0], act_val.time_vals[0], 9)
        self.assertAlmostEqual(exp_val.piYX_vals[0], act_val.piYX_vals[0], 9)
        self.assertAlmostEqual(exp_val.n_ave, act_val.n_ave, 9)
        self.assertAlmostEqual(exp_val.piYX_ave, act_val.piYX_ave, 9)
        self.assertAlmostEqual(exp_val.Q_ave, act_val.Q_ave, 9)
        self.assertAlmostEqual(exp_val.trout1_ave, act_val.trout1_ave, 9)
        self.assertAlmostEqual(exp_val.trout2_ave, act_val.trout2_ave, 9)

    def test_computeTimestepOutputs_after_ss_finalize(self):
        fene_b = 100.0 # All results in J. Chem. Phys article
        fene_n = 2.0 # Don't change
        fene_mm = 1.0 # Don't change
        sim_input = {}
        sim_input['gamdot'] = 10.0
        sim_input['begstrand'] = 100
        sim_input['eps'] = 0.001
        sim_input['steps'] = int(0.01 / sim_input['eps'])
        sim_input['outfile'] = 'testfenetrout.out'
        s1 = FENEStrand(qx=0.1, qy=0.2, qz=0.3, b=fene_b, n=fene_n, mm=fene_mm)
        s2 = FENEStrand(qx=0.3, qy=0.1, qz=0.2, b=fene_b, n=fene_n, mm=fene_mm)
        sim = ElongateNetSim(sim_input, s1, move_strand_fene_elongational)
        we = Ensemble([s1, s2])
        exp_val = ElongateNetSimOut()
        exp_val.n_vals.append(2)
        exp_val.piYX_vals.append(0.058139535)
        exp_val.Q_vals.append(0.374165739)
        exp_val.time_vals.append(0.008)
        exp_val.trout1_vals.append(0.003488372)
        exp_val.trout2_vals.append(-0.005813953)
        exp_val.n_ave = 2/sim_input['begstrand']
        exp_val.piYX_ave = 0.058139535
        exp_val.Q_ave = 0.374165739
        exp_val.trout1_ave = 0.003488372
        exp_val.trout2_ave = -0.005813953
        # Manipulate internal state so that data is present when finalize branch is taken
        sim._sim_output.n_vals.append(2)
        sim._sim_output.piYX_vals.append(0.058139535)
        sim._sim_output.Q_vals.append(0.374165739)
        sim._sim_output.time_vals.append(0.008)
        sim._sim_output.trout1_vals.append(0.003488372)
        sim._sim_output.trout2_vals.append(-0.005813953)
        sim._sim_output.n_ave = exp_val.n_ave
        sim._sim_output.piYX_ave = exp_val.piYX_ave
        sim._sim_output.Q_ave = exp_val.Q_ave
        sim._sim_output.trout1_ave = exp_val.trout1_ave
        sim._sim_output.trout2_ave = exp_val.trout2_ave
        # Now modify exp_val to reflect finalize calculations
        exp_val.n_ave = exp_val.n_ave/round(0.25*int(0.01 / sim_input['eps']))
        exp_val.piYX_ave = exp_val.piYX_ave/round(0.25*int(0.01 / sim_input['eps']))
        exp_val.Q_ave = exp_val.Q_ave/round(0.25*int(0.01 / sim_input['eps']))
        exp_val.trout1_ave = exp_val.trout1_ave/round(0.25*int(0.01 / sim_input['eps']))
        exp_val.trout2_ave = exp_val.trout2_ave/round(0.25*int(0.01 / sim_input['eps']))
        # Now call with finalize=True
        sim._computeTimestepOutputs(we, 0.008, True)
        act_val = sim._sim_output
        self.assertAlmostEqual(exp_val.n_vals[0], act_val.n_vals[0], 9)
        self.assertAlmostEqual(exp_val.Q_vals[0], act_val.Q_vals[0], 9)
        self.assertAlmostEqual(exp_val.time_vals[0], act_val.time_vals[0], 9)
        self.assertAlmostEqual(exp_val.piYX_vals[0], act_val.piYX_vals[0], 9)
        self.assertAlmostEqual(exp_val.n_ave, act_val.n_ave, 9)
        self.assertAlmostEqual(exp_val.piYX_ave, act_val.piYX_ave, 9)
        self.assertAlmostEqual(exp_val.Q_ave, act_val.Q_ave, 9)
        self.assertAlmostEqual(exp_val.trout1_ave, act_val.trout1_ave, 9)
        self.assertAlmostEqual(exp_val.trout2_ave, act_val.trout2_ave, 9)

    def test_run_sim_elongate_fene(self):
        seed(1234567890)

        exp_val = {'Qave': 0.16110272799636224, 'n/no': 0.99, 'piYX': 0.07245743211970784, 'trouton 1': 0.04737104142315748, 'trouton 2': 0.0021919416527666088}

        fene_b = 100.0 # All results in J. Chem. Phys article
        fene_n = 2.0 # Don't change
        fene_mm = 1.0 # Don't change

        sim_input = {}
        sim_input['gamdot'] = 10.0
        sim_input['begstrand'] = 100
        sim_input['eps'] = 0.001
        sim_input['steps'] = int(0.01 / sim_input['eps'])
        sim_input['outfile'] = 'testfenetrout.out'

        _proto_strand = FENEStrand(qx=0, qy=0, qz=0, b=fene_b, n=fene_n, mm=fene_mm)
        sim = ElongateNetSim(sim_input, _proto_strand, move_strand_fene_elongational)
        act_val = sim.run_sim()
        print(act_val)
        self.assertAlmostEqual(exp_val['Qave'], act_val.Q_ave, 15)
        self.assertAlmostEqual(exp_val['n/no'], act_val.n_ave, 15)
        self.assertAlmostEqual(exp_val['piYX'], act_val.piYX_ave,15)
        self.assertAlmostEqual(exp_val['trouton 1'], act_val.trout1_ave, 15)
        self.assertAlmostEqual(exp_val['trouton 2'], act_val.trout2_ave, 15)

    def test_run_sim_elongate_fens(self):
        seed(1234567890)

        exp_val = {'Qave': 1.5689606851693894, 'n/no': 1.01, 'piYX': -0.0614371945729071, 'trouton 1': -0.0009294645015770726, 'trouton 2': -0.018659087231995483}

        sim_input = {}
        sim_input['gamdot'] = 10.0
        sim_input['begstrand'] = 100
        sim_input['eps'] = 0.001
        sim_input['steps'] = int(0.01 / sim_input['eps'])
        sim_input['outfile'] = 'testfenstrout.out'

        _proto_strand = FENSStrand(qx=0, qy=0, qz=0)
        sim = ElongateNetSim(sim_input, _proto_strand, move_strand_fens_elongational)
        act_val = sim.run_sim()
        print(act_val)
        self.assertAlmostEqual(exp_val['Qave'], act_val.Q_ave, 15)
        self.assertAlmostEqual(exp_val['n/no'], act_val.n_ave, 15)
        self.assertAlmostEqual(exp_val['piYX'], act_val.piYX_ave,15)
        self.assertAlmostEqual(exp_val['trouton 1'], act_val.trout1_ave, 15)
        self.assertAlmostEqual(exp_val['trouton 2'], act_val.trout2_ave, 15)


class Test_ShearNetSim(unittest.TestCase):
    def test_move_fene_strand_shear(self):
        s = FENEStrand(qx=0.1, qy=0.2, qz=0.3, b=100.0, n=2.0, mm=1.0)
        alpha = 1.0 - pow(sqrt(s.str_len_sqr()), s.n)
        qx_next = s.qx + alpha * 100.0 * s.qy * 0.001
        exp_val = (qx_next, 0.2, 0.3)
        move_strand_fene_shear(s, 0.001, 100.0)
        act_val = (s.qx, s.qy, s.qz)
        self.assertTupleEqual(exp_val, act_val)

    def test_move_fens_strand_shear(self):
        s = FENSStrand(qx=0.1, qy=0.2, qz=0.3)
        alpha = 1.0 - pow(sqrt(s.str_len_sqr()/s.max_length), 2.0)
        qx_next = s.qx + alpha * 100.0 * s.qy * 0.001
        exp_val = (qx_next, 0.2, 0.3)
        move_strand_fens_shear(s, 0.001, 100.0)
        act_val = (s.qx, s.qy, s.qz)
        self.assertTupleEqual(exp_val, act_val)

    def test_init(self):
        s = FENEStrand()
        sim = ShearNetSim({}, s, move_strand_fene_shear)
        self.assertIsInstance(sim._sim_output, ShearNetSimOut)

    def test_getOutFileHeader(self):
        s = FENEStrand()
        sim = ShearNetSim({}, s, move_strand_fene_shear)
        exp_val = "Time, Strands, <Q>, piYX, viscosity, first normal stress coefficient, second normal stress coefficient"
        act_val = sim._getOutFileHeader()
        self.assertEqual(exp_val, act_val)


if __name__ == '__main__':
    unittest.main()
