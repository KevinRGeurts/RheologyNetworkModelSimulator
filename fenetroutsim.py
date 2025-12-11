"""
This module defines the class FeneTroutSime which simulates a FENE (finitely extensible non-linear elastic) network model in elongational viscometric flow.

References:
(1) Geurts, K.R. and L.E. Wedgewood, "A finitely extensible network strand model with nonlinear backbone
forces and entanglement kinetics," J. Chem. Phys., 1-January-1997, 106(1), pp. 339-346.
(2) Biller and Petruccione, J. Chem. Phys., 89(1), pp.577-582, 1988.  

Exported classes:
    FeneTroutSim: Class simulates a FENE (finitely extensible non-linear elastic) network model in elongational viscometric flow.

Exported functions:
    None

Exported exceptions:
    None
"""


# Standard imports
from copy import deepcopy
from math import log, copysign
from array import array

# Local imports
from strand import *
from qplot import TextPlot


class FeneTroutSim:
    """
    Simulates a FENE (finitely extensible non-linear elastic)network model in elongational viscometric flow.

    Uses the method of Biller and Petruccione, J. Chem. Phys., 89(1), pp.577-582, 1988.  

    Nondimensionalization Scheme:
        Strand internal coordinates are divided by Qo, the maximum strand length.
        Loss rate is multiplied by the loss rate constant, Lambdao.
        Time step size is divided by the loss rate constant, Lamdao.
        Stress tensor is divided by NokT.
    """
    def __init__(self, si={}):
        """
        Create an object used to run simulation.

        :param si: Dictionary of non-dimensional input values for the simulation as follows:
            gamdot = elongation rate, float
            begstrand = number of strands in equilibrium ensemble, int
            eps = time step size, float
            steps = total time steps, int
            b = FENE Parameter, float
            n = Exponent in Prefactor, int
            mm = Loss rate exponent is m*n, int
            outfile = file for writing simulation steps output, string
        """
        assert(type(si)==dict)
        self.sim_input = si

    def run_sim(self):
        """
        Execute the simulation.

        :return: A dictionary of non-dimensional results from the simulation, as follows:

            'n/no' = Number of entangled strands in the network ensemble relative to the equilibrium number, int
            'piYX' = YX-component of the total stress tensor, float
            'trouton 1' = Trouton viscosity, float
            'trouton 2' = Trouton viscosity, float
            'Qave' = Average length of a strand in the network ensemble, float
            'trouton 1 plot' = Trouton viscosity plot, qplot.TextPlot object
        """
        #  Open a file for writing output from steps
        out_f = open(self.sim_input['outfile'], 'w')

        # Write a header line into the files with column headers
        out_f.write('Step, Strands, <Q>, piYX, trouton_1, trouton_2\n')

        # The equilibrium ensemble is used for selecting new strands to join the network
        ee = generate_eq_ensemble(self.sim_input['begstrand'], self.sim_input['b'])
        # TODO: Fix so that uses any path provided by sim_input['outfil']
        write_ensemble_to_file(ee, 'ensemblestrands.csv')

        # Copy the equilibrium ensemble as the initial working ensemble
        we = deepcopy(ee)

        # Create some variables so that have the right scope
        q_ave, piYX, trout1, trout2 = 0.0, 0.0, 0.0, 0.0

        # Arrays to hold all values of trout1 and time for plotting later
        trout1_vals = array('f')
        time_vals = array('f')
        
        # For each time step of the simulation ...
        for ts in range(1, self.sim_input['steps']):

            # Advance each strand over one time step
            for s in we:
                self.move_strand(s, self.sim_input['eps'], self.sim_input['gamdot'], self.sim_input['n'])

            # Check survival of each strand
            # Make a copy to iterate over. TODO: Try list comprehension instead? Or a list ot strands to remove afterwards.
            for s in list(we):
                p1 = uniform(0.0, 1.0)
                if p1 <= s.loss_prob(self.sim_input['eps'], self.sim_input['n'], self.sim_input['mm']):
                    we.remove(s)

            # Create new strands to enter the network, from the equilibrium ensemble
            for s in ee:
                p1 = uniform(0.0, 1.0)
                if p1 <= s.loss_prob(self.sim_input['eps'], self.sim_input['n'], self.sim_input['mm']):
                    ns = Strand(s.qx, s.qy, s.qz)
                    # Move the new strand for the portion of the part of a time step for which it existed
                    p1 = uniform(0.0, 1.0)
                    y = -log(1. - p1 * ns.loss_prob(self.sim_input['eps'], self.sim_input['n'], self.sim_input['mm'])) / ns.loss_rate(self.sim_input['n'], self.sim_input['mm'])
                    self.move_strand(ns, y, self.sim_input['gamdot'], self.sim_input['n'])

                    we.append(ns)

            # Compute some output values for this time step
            q_ave = ensemble_q_ave(we)
            stress = ensemble_stress(we, self.sim_input['b'])
            piYX = stress[3] / self.sim_input['begstrand']
            trout1 = (stress[2] - stress[0]) / self.sim_input['begstrand'] / self.sim_input['gamdot']
            trout2 = (stress[1] - stress[0]) / self.sim_input['begstrand'] / self.sim_input['gamdot']

            # Print some output to the stdout
            print('Step: %i, Strands: %i, <Q>: %f, piYX: %f, trouton_1: %f, trouton_2: %f\n' % (
            ts, len(we), q_ave, piYX, trout1, trout2))

            # Write some output to a file
            out_f.write('%i, %i, %f, %f, %f, %f\n' % (ts, len(we), q_ave, piYX, trout1, trout2))

            # Append the trout1 value to the list for plotting later
            trout1_vals.append(trout1)
            time_vals.append(ts * self.sim_input['eps'])

        out_f.close()

        # Build a dictionary of results to return
        results = {}
        results['n/no'] = len(we)/len(ee)
        results['piYX'] = piYX
        results['trouton 1'] = trout1
        results['trouton 2'] = trout2
        results['Qave'] = q_ave
        symbol=array('u')
        symbol.append('O')
        results['trouton 1 plot'] = TextPlot(ncur=1, npts=len(trout1_vals), x=[time_vals], y=[trout1_vals],
                                             symbol=symbol, titl1='FENE Network Start-up of Elongational Flow Simulation',
                                             titl2='O: Elongational Viscosity vs Time')

        return results

    # Appy Euler integration to advance a strand s, one time step eps, under elongation rate gamdot.
    # n parameters the non-affine motion
    def move_strand(self, s, eps, gamdot, n):
        """
        Apply Euler integration to advance the internal coordinates of a network strand by one time step.

        :param s: Network Strand to advance, as Strand object
        :param eps: Time step size, as float
        :param gamdot: Elongation rate, as float
        :param n: Exponent on Q/Qo in the non-affine motion prefactor. Should be equal to 2, as float

        :return: None. Strand s's internal coordinates are updated upon return.
        """

        alpha = 1.0 - pow(sqrt(s.str_len_sqr()), n)

        # Actual movement of the strand by Euler integration
        s.qx = s.qx - alpha * gamdot * .5 * s.qx * eps
        s.qy = s.qy - alpha * gamdot * .5 * s.qy * eps
        s.qz = s.qz + alpha * gamdot * s.qz * eps

        if s.str_len_sqr() > 1.0:
            s.qz = copysign(sqrt(0.9999 - s.qy * s.qy - s.qx * s.qx), s.qz)

