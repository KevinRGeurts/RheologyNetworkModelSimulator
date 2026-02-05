# RheologyNetworkModelSimulator

Source code: [GitHub](https://github.com/KevinRGeurts/RheologyNetworkModelSimulator)
---
RheologyNetworkModelSimulator is a python implementation of a rheology simulator for polymer network models. In the current version it can
can simulate the FENE (finitely extensible non-linear elastic) network model in elongational viscometric flow,
and the FENS (finitely extensible network strand) network model in elongational viscometric flow.

The FENE network model is based on the work of Geurts and Wedgewood (1997).
The FENS network model is based on the work of Wedgewood and Geurts (1995).

## References:
1. Geurts, K.R. and L.E. Wedgewood, "A finitely extensible network strand model with nonlinear backbone forces and entanglement kinetics," J. Chem. Phys., 1-January-1997, 106(1), pp. 339-346.
2. Wedgewood, L.E. and K.R. Geurts, "A non-affine network model for polymer melts," Rheol. Acta 34, 196 (1995).
3. Biller and Petruccione, J. Chem. Phys., 89(1), pp.577-582, 1988.  

## Requirements
- UserResponseCollector>=1.1.0: [GitHub](https://github.com/KevinRGeurts/UserResponseCollector), [PyPi](https://pypi.org/project/UserResponseCollector/)

## Usage
To run a simulation interactively:
```
python -m RheologyNetworkModelSimulator.main
```
First the type of simulation will be requested from the console, and then the required simulation input will be requested from the console.
When the simulation is completed, select results will be printed to the console.

Additional results will appear in the files:
1. eq__ensemble__strands.csv: A listing for each network strand in the initial ensemble of its internal coordinates and its length.
2. _filename_: In the output file path/name requested in the simulation input, will be a listing at each simulation time step of the number of entangled strands, the average strand length, the YX-component of the stress tensor, the Trouton 1 viscosity value, and the Trouton 2 viscosity value.

## Unittests
Unit tests for RheologyNetworkModelSimulator have filenames ending with _tests.py. To run the unit tests,
type ```python -m unittest discover -s .\..\tests -p "*_tests.py" -v``` in a terminal window in the project directory.

## License
MIT License. See the LICENSE file for details