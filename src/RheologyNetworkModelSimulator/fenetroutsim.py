"""
This module defines the class FeneTroutSime which simulates a FENE (finitely extensible non-linear elastic) network model in elongational viscometric flow.

References:
(1) Geurts, K.R. and L.E. Wedgewood, "A finitely extensible network strand model with nonlinear backbone
forces and entanglement kinetics," J. Chem. Phys., 1-January-1997, 106(1), pp. 339-346.
(2) Biller and Petruccione, J. Chem. Phys., 89(1), pp.577-582, 1988.  

Exported classes:
    FeneTroutSim: Class simulates a FENE (finitely extensible non-linear elastic) network model in elongational viscometric flow.
    FeneTroutSimOut: Class attributes are the results of FeneTroutSim.run_sim() method.

Exported functions:
    None

Exported exceptions:
    None
"""


# Standard imports
from copy import deepcopy
from math import log, copysign
from array import array
from random import uniform
from math import sqrt

# Local imports
from RheologyNetworkModelSimulator.strand import Strand, generate_eq_ensemble, write_ensemble_to_file, ensemble_stress, ensemble_q_ave
from RheologyNetworkModelSimulator.qplot import TextPlot


class FeneTroutSimOut:
    """
    A class to accumulate (over time steps) and then return time-averaged output from FeneTroutSim.run_sim().
    """
    def __init__(self):
        """
        n_ave = n/no = Time-averaged, non-dimensional (relative to equilibrium number) number of entangled strands in the network ensemble, int
        piYX_ave = Time-averaged, non-dimensional YX-component of the total stress tensor, float
        trout1_ave = Time-averaged, non-dimensional Trouton 1 viscosity, float
        trout3_ave = Time-averaged, non-dimensional Trouton 2 viscosity, float
        Q_ave = Time-averaged, non-dimensinal ensemble-average network strand length, float
        trout1_plt = TextPlot object plot of Trouton 1 viscosity
        """
        self.n_ave=0
        self.piYX_ave=0
        self.trout1_ave=0
        self.trout2_ave=0
        self.Q_ave=0
        self.trout1_plt=None

    def __str__(self):
        """
        Convert object to str, for example, to print results.
        :return: Class attributes converted to a string representation, suitable for printing results, as str
        """
        result = f"n/no = {self.n_ave}\nPi-YX = {self.piYX_ave}\ntrout1 = {self.trout1_ave}\ntrout2 = {self.trout2_ave}\nQave = {self.Q_ave}\n\n{self.trout1_plt}\n"
        return result


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
        self._time_frac_ss = 0.25 # The fraction of simulation time the simulation is assumed to be at stead state

    def run_sim(self):
        """
        Execute the simulation.
        :return: Results from the simulation, as a FeneTroutSimOut object
        """
        results = FeneTroutSimOut()
        # Results will be averaged over this length of time
        t_res = round(self.sim_input['steps']*self._time_frac_ss)
        
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

            # If we are far enough along (that is, in the last 1/3 of the simulation time steps)
            # to likely be at steady state, store some information into results,
            # which will later be used to generate time-averaged values
            if ts >= round(self.sim_input['steps']*(1.-self._time_frac_ss)):
                results.n_ave+= len(we)/len(ee)
                results.piYX_ave+= piYX
                results.trout1_ave+= trout1
                results.trout2_ave+= trout2
                results.Q_ave+= q_ave

        out_f.close()

        # Finalize results to return
        results.n_ave=results.n_ave/t_res
        results.piYX_ave=results.piYX_ave/t_res
        results.trout1_ave=results.trout1_ave/t_res
        results.trout2_ave=results.trout2_ave/t_res
        results.Q_ave=results.Q_ave/t_res
        symbol=array('u')
        symbol.append('O')
        results.trout1_plt = TextPlot(ncur=1, npts=len(trout1_vals), x=[time_vals], y=[trout1_vals],
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

