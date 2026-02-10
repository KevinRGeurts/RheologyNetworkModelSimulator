"""
This module's __main__ first asks the user how they want to use the Network Model Simulator, then
collects simulation input interactively from the user, executes a simulation, and prints some results.

Exported classes:
    None

Exported functions:
    __main__: Requests simulation type, collects input, runs simulation, prints results output
    Fene_Elongational: Run an elongatonal flow simulation using the FENE network model.
    Fens_Elongational: Run an elongatonal flow simulation using the FENS network model
    Fene_Shear: Run a shear flow simulation using the FENE network model.
    Fens_Shear: Run a shear flow simulation using the FENS network model.
    debug: Run a debugging scenario (currently does nothing).

Exported exceptions:
    None
"""


# local imports
from UserResponseCollector.UserQueryCommand import askForInt, askForFloat, askForPathSave, askForPathSave, askForStr, askForMenuSelection
from RheologyNetworkModelSimulator.rheonetsim import ElongateNetSim, move_strand_fene_elongational, move_strand_fens_elongational
from RheologyNetworkModelSimulator.rheonetsim import ShearNetSim, move_strand_fene_shear, move_strand_fens_shear
from RheologyNetworkModelSimulator.strand import FENEStrand, FENSStrand


def debug():
    """
    Run a debugging scenario.
    """
    return None


def Fene_Elongational():
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
    sim_input['banner'] = banner

    # Initialize the simulation
    _proto_strand = FENEStrand(qx=0, qy=0, qz=0, b=fene_b, n=2.0, mm=1.0) # Don't change n or mm
    sim = ElongateNetSim(sim_input, _proto_strand, move_strand_fene_elongational)

    # Execute the simulation
    so=sim.run_sim()

    # Print some output from the simulation
    print(f"\n**Final Simulation Results**\n{so.final_output()}")

    return None


def Fene_Shear():
    """
    Run a shear flow simulation using the FENE network model.
    """
    gamdot = askForFloat('Enter the shear rate (suggest 30.0)', minimum=0.0)
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
    sim_input['banner'] = banner

    # Initialize the simulation
    _proto_strand = FENEStrand(qx=0, qy=0, qz=0, b=fene_b, n=2.0, mm=1.0) # Don't change n or mm
    sim = ShearNetSim(sim_input, _proto_strand, move_strand_fene_shear)

    # Execute the simulation
    so=sim.run_sim()

    # Print some output from the simulation
    print(f"\n**Final Simulation Results**\n{so.final_output()}")

    return None

def Fens_Elongational():
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
    sim_input['banner'] = banner

    # Initialize the simulation
    _proto_strand = FENSStrand(qx=0, qy=0, qz=0)
    sim = ElongateNetSim(sim_input, _proto_strand, move_strand_fens_elongational)

    # Execute the simulation
    so=sim.run_sim()

    # Print some output from the simulation
    print(f"\n**Final Simulation Results**\n{so.final_output()}")

    return None


def Fens_Shear():
    """
    Run a shear flow simulation using the FENS network model.
    """
    gamdot = askForFloat('Enter the shear rate (suggest 30.0)', minimum=0.0)
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
    sim_input['banner'] = banner

    # Initialize the simulation
    _proto_strand = FENSStrand(qx=0, qy=0, qz=0)
    sim = ShearNetSim(sim_input, _proto_strand, move_strand_fens_shear)

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
    query_dic = {'q':'Quit', '1':'FENE Elongational', '2':'FENS Elongational', '3':'FENE Shear', '4':'FENS Shear','d':'Debug Scenario'}
    response = askForMenuSelection(query_preface, query_dic)
    
    while response != 'q':
        
        match response:
            
            case '1':
                Fene_Elongational()
                
            case '2':
                Fens_Elongational()

            case '3':
                Fene_Shear()

            case '4':
                Fens_Shear()
                
            case 'd':
                debug()
        
        print('--------------------')
        response = askForMenuSelection(query_preface, query_dic)
