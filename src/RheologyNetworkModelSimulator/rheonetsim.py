"""
This module defines classes for simulating network models for polymer melt rheology in viscometric flows.

References:
(1) Geurts, K.R. and L.E. Wedgewood, "A finitely extensible network strand model with nonlinear backbone
forces and entanglement kinetics," J. Chem. Phys., 1-January-1997, 106(1), pp. 339-346.
(2) Biller and Petruccione, J. Chem. Phys., 89(1), pp.577-582, 1988.  

Exported classes:
    RheoNetSim: Base class for simulators of rheology network models in viscometric flows.
    ElongateNetSim: Subclass of RheoNetSim that simulates a network strand model in elongational viscometric flow.
    RheoNetSimOut: Class to accumulate (over time steps) and then return time-averaged output from RheoNetSim.run_sim() method.
    ElongateNetSimOut: Subclass of RheoNetSimOut to accumulate (over time steps) and then return time-averaged output from RheoNetSim.run_sim() method,
                       with additional output values appropriate for elongational flow simulation.

    Note that:
    (1) A XNetSim class uses it's factory method _createSimOutObj() to create the correct type of XSimOut object.
    (2) The XSimOut classes are populated by calls to XNetSim._computeTimestepOutputs() method.

Exported functions:
    move_strand_fene_elongational: Appy Euler integration to advance a FENE strand s, one time step eps, under elongation rate gamdot.
    move_strand_fens_elongational: Appy Euler integration to advance a FENS strand s, one time step eps, under elongation rate gamdot.

Exported exceptions:
    None
"""


# Standard imports
from copy import deepcopy
from math import log, copysign
from array import array
from random import uniform
from math import sqrt
from symtable import Symbol

# Local imports
from RheologyNetworkModelSimulator.strand import Strand, FENEStrand
from RheologyNetworkModelSimulator.ensemble import Ensemble
from RheologyNetworkModelSimulator.qplot import TextPlot


class RheoNetSimOut:
    """
    This is a base class for subclasses which accumulate (over time steps) and then returns time-averaged output from RheoNetSim.run_sim().
    """
    def __init__(self):
        """
        Time-averaged output values from simulation:
            n_ave = n/no = Time-averaged, non-dimensional (relative to equilibrium number) number of entangled strands in the network ensemble, int
            piYX_ave = Time-averaged, non-dimensional YX-component of the total stress tensor, float
            Q_ave = Time-averaged, non-dimensinal ensemble-average network strand length, float
        Arrays of values over time steps in simulation output:
            n_vals = Array of strand number values in simulation output, array of int
            piYX_vals = Array of piYX values in simulation output, array of float
            Q_vals = Array of Q values in simulation output, array of float
            time_vals = Array of time values in simulation output, array of float
        """
        self.n_ave=0
        self.piYX_ave=0
        self.Q_ave=0
        self.n_vals = array('i')
        self.piYX_vals = array('f')
        self.Q_vals = array('f')
        self.time_vals = array('f')

    def __str__(self):
        """
        Convert object to str, for example, to write results to file.
        :return: Class attributes converted to a string representation, suitable for writing results to a file, as str
        """
        result = f"{self.n_vals[len(self.n_vals)-1]}, {self.piYX_vals[len(self.piYX_vals)-1]}, {self.Q_vals[len(self.Q_vals)-1]}"
        return result

    def time_step_output(self):
        """
        Generate a string representation of labeled output values for the latest time step, for example, to print to show progress of simulation.
        :return: Class attributes converted to a string representation, suitable for printing time step results to show progress of simulation, as str
        """
        result = f"n: {self.n_vals[len(self.n_vals)-1]}, piYX: {self.piYX_vals[len(self.piYX_vals)-1]}, Q: {self.Q_vals[len(self.Q_vals)-1]}"
        return result

    def final_output(self):
        """
        Generate a string representation of labeled output values for the final time-averaged results, for example, to print at the end of simulation.
        :return: Class attributes converted to a string representation, suitable for printing final results at end of simulation, as str
        """
        result = f"Avg n: {self.n_ave}\nAvg piYX: {self.piYX_ave}\nAvg Q: {self.Q_ave}\n"
        return result


class ElongateNetSimOut(RheoNetSimOut):
    """
    A class to accumulate (over time steps) and then return time-averaged output from ElongateNetSim.run_sim().
    """
    def __init__(self):
        """
        Time-averaged output values from simulation:
            trout1_ave = Time-averaged, non-dimensional Trouton 1 viscosity, float
            trout2_ave = Time-averaged, non-dimensional Trouton 2 viscosity, float
        Arrays of values over time steps in simulation output:
            trout1_vals = Array of Trouton 1 viscosity values in simulation output, array of float
            trout2_vals = Array of Trouton 2 viscosity values in simulation output, array of float
        TextPlot objects from simulation:
            trout1_plt = TextPlot object plot of Trouton 1 viscosity
        """
        super().__init__()
        self.trout1_ave=0
        self.trout2_ave=0
        self.trout1_vals= array('f')
        self.trout2_vals= array('f')
        self.trout1_plt=None

    def __str__(self):
        """
        Convert object to str, for example, to write results to file.
        :return: Class attributes converted to a string representation, suitable for writing results to a file, as str
        """
        result = super().__str__()
        result += f", {self.trout1_vals[len(self.trout1_vals)-1]}, {self.trout2_vals[len(self.trout2_vals)-1]}"
        return result

    def time_step_output(self):
        """
        Generate a string representation of labeled output values for the latest time step, for example, to print to show progress of simulation.
        :return: Class attributes converted to a string representation, suitable for printing time step results to show progress of simulation, as str
        """
        result = super().time_step_output()
        result += f", Trout 1: {self.trout1_vals[len(self.trout1_vals)-1]}, Trout 2: {self.trout2_vals[len(self.trout2_vals)-1]}"
        return result

    def final_output(self):
        """
        Generate a string representation of labeled output values for the final time-averaged results, for example, to print at the end of simulation.
        :return: Class attributes converted to a string representation, suitable for printing final results at end of simulation, as str
        """
        result = super().final_output()
        result += f"Avg Trout 1: {self.trout1_ave}\nAvg Trout 2: {self.trout2_ave}\n\n{self.trout1_plt}"
        return result


class ShearNetSimOut(RheoNetSimOut):
    """
    A class to accumulate (over time steps) and then return time-averaged output from ShearNetSim.run_sim().
    """
    def __init__(self):
        """
        Time-averaged output values from simulation:
            viscosity_ave = Time-averaged, non-dimensional viscosity, float
            psi1_ave = Time-averaged, non-dimensional first normal stress coefficient, float
            psi2_ave = Time-averaged, non-dimensional second normal stress coefficient, float
        Arrays of values over time steps in simulation output:
            viscosity_vals = Array of viscosity values in simulation output, array of float
            psi1_vals = Array of first normal stress coefficient values in simulation output, array of float
            psi2_vals = Array of second normal stress coefficient values in simulation output, array of float
        TextPlot objects from simulation:
            visocisty_plt = TextPlot object plot of viscosity
        """
        super().__init__()
        self.viscosity_ave=0
        self.psi1_ave=0
        self.psi2_ave=0
        self.viscosity_vals = array('f')
        self.psi1_vals = array('f')
        self.psi2_vals = array('f')
        self.viscosity_plt=None

    def __str__(self):
        """
        Convert object to str, for example, to write results to file.
        :return: Class attributes converted to a string representation, suitable for writing results to a file, as str
        """
        result = super().__str__()
        result += f", {self.viscosity_vals[len(self.viscosity_vals)-1]}, {self.psi1_vals[len(self.psi1_vals)-1]}, {self.psi2_vals[len(self.psi2_vals)-1]}"
        return result

    def time_step_output(self):
        """
        Generate a string representation of labeled output values for the latest time step, for example, to print to show progress of simulation.
        :return: Class attributes converted to a string representation, suitable for printing time step results to show progress of simulation, as str
        """
        result = super().time_step_output()
        result += f", Viscosity: {self.viscosity_vals[len(self.viscosity_vals)-1]}, First Normal Stress Coefficient: {self.psi1_vals[len(self.psi1_vals)-1]}, Second Normal Stress Coefficient: {self.psi2_vals[len(self.psi2_vals)-1]}"
        return result

    def final_output(self):
        """
        Generate a string representation of labeled output values for the final time-averaged results, for example, to print at the end of simulation.
        :return: Class attributes converted to a string representation, suitable for printing final results at end of simulation, as str
        """
        result = super().final_output()
        result += f"Avg Viscosity: {self.viscosity_ave}\nAvg First Normal Stress Coefficient: {self.psi1_ave}\nAvg Second Normal Stress Coefficient: {self.psi2_ave}\n\n{self.viscosity_plt}"
        return result


class RheoNetSim:
    """
    Base class for simulators of rheology network models in viscometric flows.

    Uses the method of Biller and Petruccione, J. Chem. Phys., 89(1), pp.577-582, 1988.  

    Nondimensionalization Scheme:
        Strand internal coordinates are non-dimensionalized consistent with proto_strand argument to __init__().
        Loss rate is multiplied by the loss rate constant, Lambdao.
        Time step size is divided by the loss rate constant, Lamdao.
        Stress tensor is divided by NokT.
    """
    def __init__(self, si={}, proto_strand=Strand(), strand_mover=None):
        """
        Create an object used to run simulation.

        :param si: Dictionary of non-dimensional input values for the simulation as follows:
            gamdot = shear rate or elongation rate, float
            begstrand = number of strands in equilibrium ensemble, int
            eps = time step size, float
            steps = total time steps, int
            outfile = file for writing simulation steps output, string
        :param proto_strand: A prototype Strand object used to create new strands in the network ensemble, Strand subclass object
        :param strand_mover: A function that moves a Strand object according to the rheological model and viscometric flow being simulated, function
                             Function signature: strand_mover(strand, eps, gamdot, n), where
                                strand is a Strand object to move,
                                eps is the time step size, float,
                                gamdot is the shear or elongation rate, float,
        """
        assert(type(si)==dict)
        self.sim_input = si
        self._time_frac_ss = 0.25 # The fraction of simulation time the simulation is assumed to be at steady state
        assert(isinstance(proto_strand, Strand))
        self._proto_strand = proto_strand
        assert(callable(strand_mover))
        self._strand_mover = strand_mover
        self._sim_output = self._createSimOutObj()

    def _createSimOutObj(self):
        """
        Create and return an output object appropriate for this simulation class. Must be implemented by subclasses.
        Will raise NotImplementedError if called from this base class.
        :return: An output object appropriate for this simulation class, RheoNetSimOut subclass object
        """
        raise NotImplementedError("Subclasses must implement _createSimOutObj method")
        return RheoNetSimOut()

    def _getOutFileHeader(self):
        """
        Return a string to write as a header line into the output file, that is, the column headers. Should be extended as needed by subclasses.
        This would typically be done by calling super() to get the base class functionality, and then appending additional column headers by
        += f", {}" the returned string.
        :return: String for output file header, that is, the column headers, as string
        """
        result = f"Time, Strands, <Q>, piYX"
        return result

    def _computeTimestepOutputs(self, we=None, time_value=0, finalize=False):
        """
        Compute output values for a single time step. Values should be appended to the "vals" arrays of the output object.
        Should be extended as needed by subclasses. Subclasses should call super() to get base class functionality.
        :param we: Working ensemble at current time step, Ensemble object
        :param time_value: Current time value, float
        :param finalize: If True, finalize any calculations needed at the end of the simulation, bool
                         Typically, this would be used to compute time-averaged values from accumulated sums.
        :return: None. Output values are appended to the "vals" arrays of the output object.
        """
        assert(isinstance(we, Ensemble))

        if not finalize:
            # Compute and store output values for this time step, and accumulate sums for time-averaged values later

            n_val = len(we)
            self._sim_output.n_vals.append(n_val)

            stress = we.ensemble_stress()
            piYX = stress[3] / self.sim_input['begstrand']
            self._sim_output.piYX_vals.append(piYX)

            q_val = we.q_ave()
            self._sim_output.Q_vals.append(q_val)

            self._sim_output.time_vals.append(time_value)

            # If we are far enough along to assume steady state, store some information into results,
            # which will later be used to generate time-averaged values
            if time_value >= self.sim_input['eps']*self.sim_input['steps']*(1.-self._time_frac_ss):
                self._sim_output.n_ave += n_val/self.sim_input['begstrand']
                self._sim_output.piYX_ave += piYX
                self._sim_output.Q_ave += q_val

        else:
            # Finalize any calculations needed at the end of the simulation.
            # Here we will compute time-averaged values from accumulated sums, by dividing by the number of time steps used in the averaging.
            t_res = round(self.sim_input['steps']*self._time_frac_ss)
            self._sim_output.n_ave = self._sim_output.n_ave / t_res
            self._sim_output.piYX_ave = self._sim_output.piYX_ave / t_res
            self._sim_output.Q_ave = self._sim_output.Q_ave / t_res

        return None
        
    def run_sim(self):
        """
        Execute the simulation, following a template design.
        :return: Results from the simulation, as a RheoNetSimOut object
        """
        #  Open a file for writing output from steps
        out_f = open(self.sim_input['outfile'], 'w')

        # Write a header line into the files with column headers
        header_str = self._getOutFileHeader() + '\n'
        out_f.write(f"{header_str}")

        # The equilibrium ensemble is used for selecting new strands to join the network.
        # It is generated using the prototype strand provided when the simulation object was created.
        ee = Ensemble(self._proto_strand.generate_eq_ensemble(self.sim_input['begstrand']))
        # TODO: Fix so that uses any path provided by sim_input['outfil']
        ee.write_ensemble_to_file('eq_ensemble_strands.csv')

        # Make a deep copy of the equilibrium  ensemble as the initial working ensemble.
        # It is important that this is a deep copy, so that the Strand objects themselves are not shared between ee and we.
        # Otherwise, as the simulation proceeded, the Strand objects in ee would be modified (e.g., deformed by movement) as well as those in we,
        # and it would no longer represent an equilibrium ensemble.
        we=deepcopy(ee)

        # For each time step of the simulation ...
        for ts in range(1, self.sim_input['steps']):

            # Advance each strand over one time step
            for s in we:
                self._strand_mover(s, self.sim_input['eps'], self.sim_input['gamdot'])

            # Check survival of each strand.
            # Since we.strands is a copy of the list of strands in we, but the Strand objects in the list themselves are shared,
            # the call to we.remove(s) below removes the correct Strand object from we, but the iteration over s in we.strands is unaffected.
            # TODO: Try list comprehension instead? Or a list ot strands to remove afterwards.
            for s in we.strands:
                p1 = uniform(0.0, 1.0)
                if p1 <= s.loss_prob(self.sim_input['eps']):
                    we.remove(s)

            # Create new strands to enter the network, from the equilibrium ensemble
            for s in ee:
                p1 = uniform(0.0, 1.0)
                if p1 <= s.loss_prob(self.sim_input['eps']):
                    # Here we create a new Strand object, ns, to add to the working ensemble we. It starts as a copy of
                    # strand s from the equilibrium ensemble, and thus starts with the same internal coordinates as strand s.
                    # However, from then on, it deforms independently.
                    ns = deepcopy(s)
                    # Move the new strand for the portion of the part of a time step for which it existed
                    p1 = uniform(0.0, 1.0)
                    y = -log(1. - p1 * ns.loss_prob(self.sim_input['eps'])) / ns.loss_rate()
                    self._strand_mover(ns, y, self.sim_input['gamdot'])

                    we.append(ns)

            # Compute some output values for this time step
            self._computeTimestepOutputs(we, ts * self.sim_input['eps'])

            # Print some output to the stdout
            print(f"Time: {round(ts * self.sim_input['eps'], 9)}, {self._sim_output.time_step_output()}")

            # Write some output to a file
            out_f.write(f"{round(ts * self.sim_input['eps'], 9)}, {self._sim_output}\n")

        out_f.close()

        # Finalize results to return
        self._computeTimestepOutputs(we, 0, finalize=True)

        return self._sim_output


class ElongateNetSim(RheoNetSim):
    """
    Simulates a network strand model in elongational viscometric flow.

    Uses the method of Biller and Petruccione, J. Chem. Phys., 89(1), pp.577-582, 1988.  

    Nondimensionalization Scheme:
        Strand internal coordinates are non-dimensionalized consistent with proto_strand argument to __init__().
        Loss rate is multiplied by the loss rate constant, Lambdao.
        Time step size is divided by the loss rate constant, Lamdao.
        Stress tensor is divided by NokT.
    """
    def __init__(self, si={}, proto_strand=Strand(), strand_mover=None):
        """
        Create an object used to run simulation.

        :param si: Dictionary of non-dimensional input values for the simulation as follows:
            gamdot = elongation rate, float
            begstrand = number of strands in equilibrium ensemble, int
            eps = time step size, float
            steps = total time steps, int
            outfile = file for writing simulation steps output, string
        """
        super().__init__(si, proto_strand, strand_mover)

    def _createSimOutObj(self):
        """
        Create and return an output object appropriate for this simulation class.
        :return: An output object appropriate for this simulation class, RheoNetSimOut subclass object
        """
        return ElongateNetSimOut()

    def _getOutFileHeader(self):
        """
        Return a string to write as a header line into the output file, that is, the column headers.
        :return: String for output file header, that is, the column headers, as string
        """
        result = super()._getOutFileHeader()
        result += f", trouton_1, trouton_2"
        return result

    def _computeTimestepOutputs(self, we=None, time_value=0, finalize=False):
        """
        Compute output values for a single time step.
        :param we: Working ensemble at current time step, Ensemble object
        :param time_value: Current time value, float
        :param finalize: If True, finalize any calculations needed at the end of the simulation, bool
                         Typically, this would be used to compute time-averaged values from accumulated sums.
        :return: None. Output values are appended to the "vals" arrays of the output object.
        """
        assert(isinstance(we, Ensemble))
        super()._computeTimestepOutputs(we, time_value, finalize)

        if not finalize:
            # Compute and store output values for this time step, and accumulate sums for time-averaged values later

            stress = we.ensemble_stress()
            trout1 = (stress[2] - stress[0]) / self.sim_input['begstrand'] / self.sim_input['gamdot']
            self._sim_output.trout1_vals.append(trout1)
            trout2 = (stress[1] - stress[0]) / self.sim_input['begstrand'] / self.sim_input['gamdot']
            self._sim_output.trout2_vals.append(trout2)


            # If we are far enough along to assume steady state, store some information into results,
            # which will later be used to generate time-averaged values
            if time_value >= self.sim_input['eps']*self.sim_input['steps']*(1.-self._time_frac_ss):
                self._sim_output.trout1_ave+=trout1
                self._sim_output.trout2_ave+=trout2

        else:
            # Finalize any calculations needed at the end of the simulation.
            # Here we will compute time-averaged values from accumulated sums, by dividing by the number of time steps used in the averaging.
            t_res = round(self.sim_input['steps']*self._time_frac_ss)
            self._sim_output.trout1_ave = self._sim_output.trout1_ave / t_res
            self._sim_output.trout2_ave = self._sim_output.trout2_ave / t_res

            # Make a TextPlot object for plotting elongational viscosity vs time
            _symbol = array('u')
            _symbol.append('O')
            _npts=len(self._sim_output.trout1_vals)
            _x=[self._sim_output.time_vals]
            _y=[self._sim_output.trout1_vals]
            _titl1='Start-up of Elongational Flow Simulation'
            _titl2='O: Elongational Viscosity vs Time'
            self._sim_output.trout1_plt = TextPlot(ncur=1, npts=_npts, x=_x, y=_y, symbol=_symbol, titl1=_titl1, titl2=_titl2)

        return None


class ShearNetSim(RheoNetSim):
    """
    Simulates a network strand model in shear viscometric flow.

    Uses the method of Biller and Petruccione, J. Chem. Phys., 89(1), pp.577-582, 1988.  

    Nondimensionalization Scheme:
        Strand internal coordinates are non-dimensionalized consistent with proto_strand argument to __init__().
        Loss rate is multiplied by the loss rate constant, Lambdao.
        Time step size is divided by the loss rate constant, Lamdao.
        Stress tensor is divided by NokT.
    """
    def __init__(self, si={}, proto_strand=Strand(), strand_mover=None):
        """
        Create an object used to run simulation.

        :param si: Dictionary of non-dimensional input values for the simulation as follows:
            gamdot = shear rate, float
            begstrand = number of strands in equilibrium ensemble, int
            eps = time step size, float
            steps = total time steps, int
            outfile = file for writing simulation steps output, string
        """
        super().__init__(si, proto_strand, strand_mover)

    def _createSimOutObj(self):
        """
        Create and return an output object appropriate for this simulation class.
        :return: An output object appropriate for this simulation class, RheoNetSimOut subclass object
        """
        return ShearNetSimOut()

    def _getOutFileHeader(self):
        """
        Return a string to write as a header line into the output file, that is, the column headers.
        :return: String for output file header, that is, the column headers, as string
        """
        result = super()._getOutFileHeader()
        result += f", viscosity, first normal stress coefficient, second normal stress coefficient"
        return result

    def _computeTimestepOutputs(self, we=None, time_value=0, finalize=False):
        """
        Compute output values for a single time step.
        :param we: Working ensemble at current time step, Ensemble object
        :param time_value: Current time value, float
        :param finalize: If True, finalize any calculations needed at the end of the simulation, bool
                         Typically, this would be used to compute time-averaged values from accumulated sums.
        :return: None. Output values are appended to the "vals" arrays of the output object.
        """
        assert(isinstance(we, Ensemble))
        super()._computeTimestepOutputs(we, time_value, finalize)

        if not finalize:
            # Compute and store output values for this time step, and accumulate sums for time-averaged values later

            stress = we.ensemble_stress()
            viscosity = stress[3] / self.sim_input['begstrand'] / self.sim_input['gamdot']
            self._sim_output.viscosity_vals.append(viscosity)
            psi1 = (stress[0] - stress[1]) / self.sim_input['begstrand'] / self.sim_input['gamdot'] / self.sim_input['gamdot']
            self._sim_output.psi1_vals.append(psi1)
            psi2 = (stress[1] - stress[2]) / self.sim_input['begstrand'] / self.sim_input['gamdot'] / self.sim_input['gamdot']
            self._sim_output.psi2_vals.append(psi2)

            # If we are far enough along to assume steady state, store some information into results,
            # which will later be used to generate time-averaged values
            if time_value >= self.sim_input['eps']*self.sim_input['steps']*(1.-self._time_frac_ss):
                self._sim_output.viscosity_ave+=viscosity
                self._sim_output.psi1_ave+=psi1
                self._sim_output.psi2_ave+=psi2

        else:
            # Finalize any calculations needed at the end of the simulation.
            # Here we will compute time-averaged values from accumulated sums, by dividing by the number of time steps used in the averaging.
            t_res = round(self.sim_input['steps']*self._time_frac_ss)
            self._sim_output.viscosity_ave = self._sim_output.viscosity_ave / t_res
            self._sim_output.psi1_ave = self._sim_output.psi1_ave / t_res
            self._sim_output.psi2_ave = self._sim_output.psi2_ave / t_res

            # Make a TextPlot object for plotting viscosity vs time
            _symbol = array('u')
            _symbol.append('O')
            _npts=len(self._sim_output.viscosity_vals)
            _x=[self._sim_output.time_vals]
            _y=[self._sim_output.viscosity_vals]
            _titl1='Start-up of Shear Flow Simulation'
            _titl2='O: Viscosity vs Time'
            self._sim_output.viscosity_plt = TextPlot(ncur=1, npts=_npts, x=_x, y=_y, symbol=_symbol, titl1=_titl1, titl2=_titl2)

        return None


# Appy Euler integration to advance a strand s, one time step eps, under elongation rate gamdot.
def move_strand_fene_elongational(s, eps, gamdot):
    """
    Apply Euler integration to advance the internal coordinates of a network strand by one time step.
    :param s: Network Strand to advance, as Strand object
    :param eps: Time step size, as float
    :param gamdot: Elongation rate, as float
    :return: None. Strand s's internal coordinates are updated upon return.
    """

    alpha = 1.0 - pow(sqrt(s.str_len_sqr()), s.n)

    # Actual movement of the strand by Euler integration
    s.qx = s.qx - alpha * gamdot * .5 * s.qx * eps
    s.qy = s.qy - alpha * gamdot * .5 * s.qy * eps
    s.qz = s.qz + alpha * gamdot * s.qz * eps

    if s.max_length is not None:
        if sqrt(s.str_len_sqr()) > s.max_length:
            s.qz = copysign(sqrt(0.9999*s.max_length - s.qy * s.qy - s.qx * s.qx), s.qz)


# Appy Euler integration to advance a strand s, one time step eps, under elongation rate gamdot.
def move_strand_fens_elongational(s, eps, gamdot):
    """
    Apply Euler integration to advance the internal coordinates of a network strand by one time step.
    :param s: Network Strand to advance, as Strand object
    :param eps: Time step size, as float
    :param gamdot: Elongation rate, as float
    :return: None. Strand s's internal coordinates are updated upon return.
    """

    alpha = 1.0 - pow(sqrt(s.str_len_sqr()/s.max_length), 2.0)

    # Actual movement of the strand by Euler integration
    s.qx = s.qx - alpha * gamdot * .5 * s.qx * eps
    s.qy = s.qy - alpha * gamdot * .5 * s.qy * eps
    s.qz = s.qz + alpha * gamdot * s.qz * eps

    if s.max_length is not None:
        if sqrt(s.str_len_sqr()) > s.max_length:
            s.qz = copysign(sqrt(0.9999*s.max_length - s.qy * s.qy - s.qx * s.qx), s.qz)


# Appy Euler integration to advance a strand s, one time step eps, under shear rate gamdot.
def move_strand_fens_shear(s, eps, gamdot):
    """
    Apply Euler integration to advance the internal coordinates of a network strand by one time step.
    :param s: Network Strand to advance, as Strand object
    :param eps: Time step size, as float
    :param gamdot: Shear rate, as float
    :return: None. Strand s's internal coordinates are updated upon return.
    """

    alpha = 1.0 - pow(sqrt(s.str_len_sqr()/s.max_length), 2.0)

    # Actual movement of the strand by Euler integration
    s.qx = s.qx + alpha * gamdot * s.qy * eps

    if s.max_length is not None:
        if sqrt(s.str_len_sqr()) > s.max_length:
            s.qz = copysign(sqrt(0.9999*s.max_length - s.qy * s.qy - s.qx * s.qx), s.qz)


# Appy Euler integration to advance a strand s, one time step eps, under shear rate gamdot.
def move_strand_fene_shear(s, eps, gamdot):
    """
    Apply Euler integration to advance the internal coordinates of a network strand by one time step.
    :param s: Network Strand to advance, as Strand object
    :param eps: Time step size, as float
    :param gamdot: Shear rate, as float
    :return: None. Strand s's internal coordinates are updated upon return.
    """

    alpha = 1.0 - pow(sqrt(s.str_len_sqr()), s.n)

    # Actual movement of the strand by Euler integration
    s.qx = s.qx + alpha * gamdot * s.qy * eps

    if s.max_length is not None:
        if sqrt(s.str_len_sqr()) > s.max_length:
            s.qz = copysign(sqrt(0.9999*s.max_length - s.qy * s.qy - s.qx * s.qx), s.qz)
