# standard imports
from array import array

# local imports
from fenetroutsim import *

if __name__ == '__main__':
    # Get input from the user

    # print('--------------------------------------------------')
    # print('-----  LEW/KRG FENE Network Model Simulator  -----')
    # print('-----  (decreasingly affine strand motion)   -----')
    # print('--------------------------------------------------')
    # gamdot = float(input('elongation rate. . .'))
    # begstrand = int(input('# strands in equilibrium ensemble. . .'))
    # eps = float(input('time step size. . .'))
    # steps = int(input('total time steps. . .'))
    # qdex = float(input('Q-Dist extension factor. . .'))
    # b = float(input('FENE Parameter. . .'))
    # n = int(input('Exponent in Prefactor. . .'))
    # mm = int(input('Loss rate exponent is m*n.  Enter m. . .'))
    # outfil = input('Name of output file . . .')
    # banner = input('Idenifier for output. . .')

    sim_input = {}
    sim_input['gamdot'] = 1
    sim_input['begstrand'] = 10000 # Should be 10000 for accurate results
    sim_input['eps'] = 0.001 # Should be 0.001 for accurate results
    sim_input['steps'] = int(5 / sim_input['eps'])
    sim_input['b'] = 100.0  # All results in J. Chem. Phys article
    sim_input['n'] = 2  # Don't change
    sim_input['mm'] = 1  # Don't change
    sim_input['outfile'] = 'fenetrout.out'

    sim = FeneTroutSim(sim_input)

    so=sim.run_sim()

    # Print some output from the simulation

    print(so['trouton 1 plot'])

