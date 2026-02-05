"""
This module defines the class Ensemble, which represents an ensemble (collection) of network strands in a network model for polymer melt rheology.

Exported classes:
    Ensemble: Class to represent an ensemble of network strands in a polymer network.

Exported functions:


Exported exceptions:
    None
"""


# Standard imports
from math import sqrt


# Local imports
from RheologyNetworkModelSimulator.strand import Strand


class Ensemble(object):
    """
    This class represents an ensemble of network strands in a polymer network. This class is iterable over its Strand objects.
    """
    def __init__(self, strands):
        """
        Construct an Ensemble with a given list of Strand objects.
        :param strands: List of Strand objects representing the ensemble.
        """
        assert(isinstance(strands, list))
        assert(isinstance(strands[0], Strand))  # Test the first element to (roughly) confirm list of Strand objects.
        self._strands = list(strands)

    @property
    def strands(self):
        """
        Get a copy of the list of Strand objects in the ensemble.
        :return: A copy of the list of Strand objects in the ensemble.
        """
        return list(self._strands)
        
    def write_ensemble_to_file(self, filename):
        """
        Write a CSV formatted text file with the internal coordinates and length of each strand  in list.
        :param filename: Name of file to write to, as string.
        :return: None.
        """
        f = open(filename, 'w')
        f.write('qx, qy, qz, length\n')
        i = 1
        for s in self._strands:
            f.write('%i, %f, %f, %f, %f\n' % (i, s.qx, s.qy, s.qz, sqrt(s.str_len_sqr())))
            i = i +1
        f.close()
        return None

    def q_ave(self):
        """
        Compute the average non-dimensional length of a strand in the network ensemble.
        :return: Average non-dimensional length of a strand in the network ensemble, as float
        Strand length is non-dimensionalized consistent self._strands.
        """
        total_q = 0.0
        for s in self._strands:
            total_q = total_q + sqrt(s.str_len_sqr())
        return total_q / len(self._strands)

    def ensemble_stress(self):
        """
        Compute the components of the total stress tensor for the network,
        by summing over the strands in the network ensemble.  Stress tensor is divided by NokT.
        :param e: List of strands (the "ensemble"), as [Strand objects].
        :return: Tuple (XX, YY, ZZ, YX) of components of total stress tensor of the network, divided by NokT, as (float, float, float, float)
        """
        piXX = 0.0
        piYY = 0.0
        piZZ = 0.0
        piYX = 0.0
        for s in self._strands:
            piXX = piXX + s.stress_XX()
            piYY = piYY + s.stress_YY()
            piZZ = piZZ + s.stress_ZZ()
            piYX = piYX + s.stress_YX()
        return (piXX, piYY, piZZ, piYX)

    # Implement methods so that an Ensemble can be treated as a container type. See Section 3.3.7. Emulating container types in
    # https://docs.python.org/3/reference/datamodel.html#classgetitem-versus-getitem
    
    def __iter__(self):
        return iter(self._strands)
    
    def __len__(self):
        return len(self._strands)
    
    def __getitem__(self, item):
        if item >= len(self):
            raise IndexError("Ensemble index out of range")
        return self._strands[item]

    def remove(self, value):
        for i in self._strands:
            if i == value:
                self._strands.remove(i)
                break
        return

    def append(self, value):
        self._strands.append(value)
        return
