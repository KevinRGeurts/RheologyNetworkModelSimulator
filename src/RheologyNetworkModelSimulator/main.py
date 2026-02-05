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
from UserResponseCollector.UserQueryCommand import askForInt, askForFloat, askForPathSave, askForPathSave, askForStr, askForMenuSelection

# local imports
from RheologyNetworkModelSimulator.rheonetsim import ElongateNetSim, move_strand_fene_elongational, move_strand_fens_elongational
from RheologyNetworkModelSimulator.strand import FENEStrand, FENSStrand


def debug():
    """
    Run a debugging scenario.
    """
    return None


def Fene_Trout():
    """
    Run an elongatonal flow simulation using the FENE network model.
    """
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
    sim_input['outfile'] = outfil
    # TODO: Investigate if this is actually used
    sim_input['banner'] = ''

    # Initialize the simulation
    _proto_strand = FENEStrand(qx=0, qy=0, qz=0, b=fene_b, n=2.0, mm=1.0) # Don't change n or mm
    sim = ElongateNetSim(sim_input, _proto_strand, move_strand_fene_elongational)

    # Execute the simulation
    so=sim.run_sim()

    # Print some output from the simulation
    print(f"\n**Final Simulation Results**\n{so.final_output()}")

    return None


def Fens_Trout():
    """
    Run an elongatonal flow simulation using the FENS network model.
    """
    gamdot = askForFloat('Enter the elongation rate (suggest 1.0)', minimum=0.0)
    begstrand = askForInt('Enter the number of strands in the equilibrium ensemble (suggest 10000)', minimum=1)
    eps = askForFloat('Enter the time step size (suggest 0.001)', minimum=1e-5, maximum=0.1)
    endt = askForFloat('Enter the end time of the simulation (suggest 5.0). Steady-state results will be averaged over the last quarter of this time.', minimum=10*eps)
    steps = int(endt / eps)
    outfil = askForPathSave('Path of output file')
    banner = askForStr('Type an identifier for ouput')

    # Package the required simulation input into a dictionary
    sim_input = {}
    sim_input['gamdot'] = gamdot
    sim_input['begstrand'] = begstrand
    sim_input['eps'] = eps
    sim_input['steps'] = steps
    sim_input['outfile'] = outfil
    # TODO: Investigate if this is actually used
    sim_input['banner'] = ''

    # Initialize the simulation
    _proto_strand = FENSStrand(qx=0, qy=0, qz=0)
    sim = ElongateNetSim(sim_input, _proto_strand, move_strand_fens_elongational)

    # Execute the simulation
    so=sim.run_sim()

    # Print some output from the simulation
    print(f"\n**Final Simulation Results**\n{so.final_output()}")

    return None


if __name__ == '__main__':
    """
    Query the user for how they wish to use the Network Model Simulator, and then launch that usage.
    This includes a "debug" usage to set up what ever situation is needed for debugging, since I can't seem to reliably debug unit tests.
    """
    print('---------------------------------------------')
    print('-----     Network Model Simulator       -----')
    print('---------------------------------------------')
    
    # Build a query for the user to obtain their choice of how to user the simulator
    query_preface = 'How do you want to use the simulator?'
    query_dic = {'q':'Quit', '1':'FENE Elongational', '2':'FENS Elongational', 'd':'Debug Scenario'}
    response = askForMenuSelection(query_preface, query_dic)
    
    while response != 'q':
        
        match response:
            
            case '1':
                Fene_Trout()
                
            case '2':
                Fens_Trout()
                
            case 'd':
                debug()
        
        print('--------------------')
        response = askForMenuSelection(query_preface, query_dic)



