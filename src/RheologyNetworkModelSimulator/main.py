"""
This module's __main__ collects simulation input interactively from the user, executes a FeneTroutSim, and prints some results.

Exported classes:
    None

Exported functions:
    __main__: Collects input, runs simulation, prints results output

Exported exceptions:
    None
"""


# standard imports

# Environment package imports
from UserResponseCollector.UserQueryCommand import askForInt, askForFloat, askForPathSave, askForPathSave, askForStr

# local imports
from fenetroutsim import *


if __name__ == '__main__':
    # Get input from the user

    print('--------------------------------------------------')
    print('-----     FENE Network Model Simulator       -----')
    print('-----  (decreasingly affine strand motion)   -----')
    print('--------------------------------------------------')

    gamdot = askForFloat('Enter the elongation rate (suggest 1.0)', minimum=0.0)
    begstrand = askForInt('Enter the number of strands in the equilibrium ensemble (suggest 10000)', minimum=1)
    eps = askForFloat('Enter the time step size (suggest 0.001)', minimum=1e-5, maximum=0.1)
    endt = askForFloat('Enter the end time of the simulation (suggest 5.0). Steady-state results will be averaged over the last quarter of this time.', minimum=10*eps)
    steps = int(endt / eps)
    # Suggested value for b matches all results in J. Chem. Phys. article
    fene_b = askForFloat('Enter the FENE Parameter b (suggest 100.0)', minimum=1.0)
    outfil = askForPathSave('Path of output file')
    banner = askForStr('Type an identifier for ouput')

    # Package the required simulation input into a dictionary
    sim_input = {}
    sim_input['gamdot'] = gamdot
    sim_input['begstrand'] = begstrand
    sim_input['eps'] = eps
    sim_input['steps'] = steps
    sim_input['b'] = fene_b
    sim_input['n'] = 2.0 # Don't change
    sim_input['mm'] = 1.0  # Don't change
    sim_input['outfile'] = outfil
    # TODO: Investigate if this is actually used
    sim_input['banner'] = ''

    # Initialize the simulation
    sim = FeneTroutSim(sim_input)

    # Execute the simulation
    so=sim.run_sim()

    # Print some output from the simulation
    print(so)

