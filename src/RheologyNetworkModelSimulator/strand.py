"""
This module defines the class Strand, which represents a FENE (finitely extensible non-linear elastic) strand in a network model for polymer melt rheology.

Exported classes:
    Strand: Class to represent a FENE strand in a polymer network.

Exported functions:
    beta: Function to compute the mathematical beta function.
    generate_eq_ensemble: Function to create an equilibrium ensemble of strands.
    write_ensemble_to_file: Function to write an ensemble of strands to a CSV file.
    ensemble_stress: Function to compute the total stress tensor components for an ensemble of strands.
    ensemble_q_ave: Function to compute the average strand length in an ensemble.

Exported exceptions:
    None
"""

# Standard imports
from math import pow, sqrt, exp, lgamma, pi, sin, cos
from random import uniform

class Strand:
    """
    Represent a strand in a network model for polymer melt rheology.

    Reference: Geurts, K.R. and L.E. Wedgewood, "A finitely extensible network strand model with nonlinear backbone
    forces and entanglement kinetics," J. Chem. Phys., 1-January-1997, 106(1), pp. 339-346.

    Nondimensionalization Scheme:
        Strand internal coordinates are divided by Qo, the maximum strand length ??Hookean has no max length??.
        Loss rate is multiplied by the loss rate constant, Lambdao.
        Time step size is divided by the loss rate constant, Lamdao.
        Stress tensor is divided by NokT.
    """
    def __init__(self, qx=0.0, qy=0.0, qz=0.0):
        """
        Construct a Strand with given X, Y, and Z internal coordinate lengths, all non-dimensionalized by
        dividing by the maximum strand length Qo ??Hookean has no max length??.

        :param qx: Internal X-coordinate of Strand length, non-dimensionalizec by dividing by Qo, the maximum strand length, as float
        :param qy: Internal Y-coordinate of Strand length, non-dimensionalizec by dividing by Qo, the maximum strand length, as float
        :param qz: Internal Z-coordinate of Strand length, non-dimensionalizec by dividing by Qo, the maximum strand length, as float
        """
        self._max_length = None
        self._qx=qx
        self._qy=qy
        self._qz=qz

    @property
    def max_length(self):
        return self._max_length

    @property
    def qx(self):
        return self._qx

    # TODO: The setters for the internal strand coordinates enforce the length constraint immediately.
    # But be aware that since they can be called independently, there could be a situation where the user
    # plans to change the other coordinates as well, so that the final length is valid, but when the first to second
    # setter is called, that it isn't valid yet.

    @qx.setter
    def qx(self, value):
        if self._max_length is not None:
            assert(sqrt(value*value + self.qy*self.qy + self.qz*self.qz) <= self._max_length)
        self._qx = value

    @property
    def qy(self):
        return self._qy

    @qy.setter
    def qy(self, value):
        if self._max_length is not None:
            assert(sqrt(self.qx*self.qx + value*value + self.qz*self.qz) <= self._max_length)
        self._qy = value

    @property
    def qz(self):
        return self._qz

    @qz.setter
    def qz(self, value):
        if self._max_length is not None:
            assert(sqrt(self.qx*self.qx + self.qy*self.qy + value*value) <= self._max_length)
        self._qz = value

    def get_qs(self):
        """
        Get the internal coordinates of the Strand, all non-dimensionalized by dividing by the maximum strand length Qo.
        :return: Tuple (qx, qy, qz) of internal coordinates of the Strand, as (float, float, float)
        """
        return (self.qx, self.qy, self.qz)

    def str_len_sqr(self):
        """
        Compute the length of the Strand, squared. Non-dimensionalized by dividing by maximum strand length Qo, squared.
        :return: Strand length squared, as float
        """
        return pow(self.qx, 2)+pow(self.qy, 2)+pow(self.qz, 2)

    def loss_rate(self):
        """
        Compute the loss_rate for the Strand. Loss rate is multiplied by the loss rate constant, Lambdao.
        Must be implemented by subclasses.
        :return: Loss rate for strand, multiplied by the loss rate constant, Lambdao, as float
        """
        raise NotImplementedError("Subclasses must implement loss_rate method.")
        return None

    def loss_prob(self, eps):
        """
        Compute the probability of the Strand disentangling from the network during a time step of length eps.
        :param eps: Time step size, divided by the loss rate constant, Lamdao, as float
        :return: Probability of the Strand disentangling from the network during a time step of length eps, as float
        """
        return 1.0 - exp(-1.0 * self.loss_rate() * eps)

    # XX component of total stress tensor
    def stress_XX(self):
        """
        Compute the XX-component of the strand's contrribution to the total stress tensor of the network.
        Stress tensor is divided by NokT. Must be implemented by subclasses.
        :return: XX-component of the strand's contribution to the total stress tensor of the network, as float
        """
        raise NotImplementedError("Subclasses must implement stress_XX method.")
        return None

    # YY component of total stress tensor
    def stress_YY(self):
        """
        Compute the YY-component of the strand's contribution to the total stress tensor of the network.
        Stress tensor is divided by NokT. Must be implemented by subclasses.
        :return: YY-component of the strand's contrribution to the total stress tensor of the network, as float
        """
        raise NotImplementedError("Subclasses must implement stress_YY method.")
        return None

    # ZZ component of total stress tensor
    def stress_ZZ(self):
        """
        Compute the ZZ-component of the strand's contribution to the total stress tensor of the network.
        Stress tensor is divided by NokT. Must be implemented by subclasses.
        :return: ZZ-component of the strand's contrribution to the total stress tensor of the network, as float
        """
        raise NotImplementedError("Subclasses must implement stress_ZZ method.")
        return None

    # YX component of total stress tensor
    def stress_YX(self):
        """
        Compute the YX-component of the strand's contrribution to the total stress tensor of the network.
        Stress tensor is divided by NokT. Must be implemented by subclasses.
        :return: YX-component of the strand's contrribution to the total stress tensor of the network, as float
        """
        raise NotImplementedError("Subclasses must implement stress_YX method.")
        return None

    def generate_eq_ensemble(self, n=1):
        """
        Create an ensemble of size n number of strands, where the distribution of strand lengths will
        fit the equilibrium distribution. Must be implemented by subclasses.
        :param n: Number of strands to generate for the ensemble, as int
        :return: List of strands in the equilibrium ensemble, as [Strand objects]
        """
        raise NotImplementedError("Subclasses must implement generate_eq_ensemble method.")
        return None


class FENSStrand(Strand):
    """
    Represent a FENS (finitely extensible network strand) strand in a network model for polymer melt rheology.

    References:
    
    (1) Wedgewood, L.E. and K.R. Geurts, Rheol. Acta 34, 196 (1995).
    (2) Geurts, K.R. and L.E. Wedgewood, "A finitely extensible network strand model with nonlinear backbone
        forces and entanglement kinetics," J. Chem. Phys., 1-January-1997, 106(1), pp. 339-346.

    Non-affine Motion prefactor is  1-(Q/Qo)^n. No strand may exceed the nondimensional length of one.
    The Hookean (linear) force law is used. The loss rate is constant.

    Nondimensionalization Scheme:
        Strand internal coordinates are divided by Qo, the maximum strand length.
        Loss rate is multiplied by the loss rate constant, Lambdao.
        Time step size is divided by the loss rate constant, Lamdao.
        Stress tensor is divided by NokT.
    """
    def __init__(self, qx=0.0, qy=0.0, qz=0.0):
        """
        Construct a Strand with given X, Y, and Z internal coordinate lengths, all non-dimensionalized by
        dividing by the maximum strand length Qo.
        :param qx: Internal X-coordinate of Strand length, non-dimensionalizec by dividing by Qo, the maximum strand length, as float
        :param qy: Internal Y-coordinate of Strand length, non-dimensionalizec by dividing by Qo, the maximum strand length, as float
        :param qz: Internal Z-coordinate of Strand length, non-dimensionalizec by dividing by Qo, the maximum strand length, as float
        """
        super().__init__(qx, qy, qz)
        self._max_length = 1.0 # Override max length to 1.0 for FENS Strand
        assert(sqrt(qx*qx + qy*qy + qz*qz) <= self._max_length)

    def loss_rate(self):
        """
        Compute the loss_rate for the Strand. Loss rate is multiplied by the loss rate constant, Lambdao.
        :return: Loss rate for strand, multiplied by the loss rate constant, Lambdao, as float
        """
        return 1.0

    # XX component of total stress tensor
    def stress_XX(self):
        """
        Compute the XX-component of the strand's contrribution to the total stress tensor of the network.
        Stress tensor is divided by NokT.
        :return: XX-component of the strand's contrribution to the total stress tensor of the network, as float
        """
        return self.qx * self.qx

    # YY component of total stress tensor
    def stress_YY(self):
        """
        Compute the YY-component of the strand's contrribution to the total stress tensor of the network.
        Stress tensor is divided by NokT.
        :return: YY-component of the strand's contrribution to the total stress tensor of the network, as float
        """
        return self.qy * self.qy

    # ZZ component of total stress tensor
    def stress_ZZ(self):
        """
        Compute the ZZ-component of the strand's contrribution to the total stress tensor of the network.
        Stress tensor is divided by NokT.
        :return: ZZ-component of the strand's contrribution to the total stress tensor of the network, as float
        """
        return self.qz * self.qz

    # YX component of total stress tensor
    def stress_YX(self):
        """
        Compute the YX-component of the strand's contrribution to the total stress tensor of the network.
        Stress tensor is divided by NokT.
        :return: YX-component of the strand's contrribution to the total stress tensor of the network, as float
        """
        return self.qy * self.qx

    def generate_eq_ensemble(self, n=1):
        """
        Create an ensemble of size n number of strands, where the distribution of strand lengths will
        fit the equilibrium distribution.
        :param n: Number of strands to generate for the ensemble, as int
        :return: List of strands in the equilibrium ensemble, as [Strand objects]
        """
        e = []
        # jeq: Normalization constant for equilibriium Hookean Dumbbell distribution function DPL2 Eq. (L) of Table 11.5-1
        jeq = (1. / (2. * pi))
        for s in range(n):
            p1 = 1.0
            p2 = 0.0
            # r, theta, phi: Internal coordinates of a network strand in spherical coordinates used to find equilibrium distribution.
            r, theta, phi = 0.0, 0.0, 0.0
            while p1>p2:
                r = uniform(0.0,1.0)
                theta = pi * uniform(0.0, 1.0)
                phi = 2.0 * pi * uniform(0.0, 1.0)
                p1 = uniform(0.0, 1.0)
                p2 = jeq * ((1.0 - r * r) ** (1.0 / 2.0)) * r * r * sin(theta)
            x = r * sin(theta) * cos(phi)
            y = r * sin(theta) * sin(phi)
            z = r * cos(theta)
            e.append(FENSStrand(x, y, z))
        return e


class FENEStrand(Strand):
    """
    Represent a FENE (finitely extensible non-linear elastic) strand in a network model for polymer melt rheology.

    Reference: Geurts, K.R. and L.E. Wedgewood, "A finitely extensible network strand model with nonlinear backbone
    forces and entanglement kinetics," J. Chem. Phys., 1-January-1997, 106(1), pp. 339-346.

    Non-affine Motion prefactor is  1-(Q/Qo)^n. No strand may exceed the nondimensional length of one.
    The FENE force law is used. The loss rate is (1/Lambdao)(1)/(1-(Q/Qo)^n*m).

    Nondimensionalization Scheme:
        Strand internal coordinates are divided by Qo, the maximum strand length.
        Loss rate is multiplied by the loss rate constant, Lambdao.
        Time step size is divided by the loss rate constant, Lamdao.
        Stress tensor is divided by NokT.
    """
    def __init__(self, qx=0.0, qy=0.0, qz=0.0, b=100.0, n=2.0, mm=1.0):
        """
        Construct a Strand with given X, Y, and Z internal coordinate lengths, all non-dimensionalized by
        dividing by the maximum strand length Qo.
        :param qx: Internal X-coordinate of Strand length, non-dimensionalizec by dividing by Qo, the maximum strand length, as float
        :param qy: Internal Y-coordinate of Strand length, non-dimensionalizec by dividing by Qo, the maximum strand length, as float
        :param qz: Internal Z-coordinate of Strand length, non-dimensionalizec by dividing by Qo, the maximum strand length, as float
        :param b: The FENE Parameter (HQo^2/kT), or non-dimensional maximum strand length, as float        
        :param n: Exponent on Q/Qo in the non-affine motion prefactor. Should be equal to 2, as float
        :param mm: mm*n is exponent on (Q/Qo) in the loss rate. mm should be equal to 1, as float
        """
        super().__init__(qx, qy, qz)
        self._max_length = 1.0 # Override max length to 1.0 for FENE Strand
        assert(sqrt(qx*qx + qy*qy + qz*qz) <= self._max_length)
        assert(b>0.0)
        assert(n==2.0)
        assert(mm==1.0)
        self._b = b
        self._n = n
        self._mm = mm

    @property
    def b(self):
        return self._b

    @property
    def n(self):
        return self._n

    @property
    def mm(self):
        return self._mm

    def loss_rate(self):
        """
        Compute the loss_rate for the Strand. Loss rate is multiplied by the loss rate constant, Lambdao.
        :return: Loss rate for strand, multiplied by the loss rate constant, Lambdao, as float
        """
        return 1.0 / (1.0 - pow(sqrt(self.str_len_sqr()), (self.mm * self.n)))

    # XX component of total stress tensor
    def stress_XX(self):
        """
        Compute the XX-component of the strand's contrribution to the total stress tensor of the network.
        Stress tensor is divided by NokT.
        :return: XX-component of the strand's contrribution to the total stress tensor of the network, as float
        """
        return self.b * self.qx * self.qx / (1.0 - self.str_len_sqr())

    # YY component of total stress tensor
    def stress_YY(self):
        """
        Compute the YY-component of the strand's contrribution to the total stress tensor of the network.
        Stress tensor is divided by NokT.
        :return: YY-component of the strand's contrribution to the total stress tensor of the network, as float
        """
        return self.b * self.qy * self.qy / (1.0 - self.str_len_sqr())

    # ZZ component of total stress tensor
    def stress_ZZ(self):
        """
        Compute the ZZ-component of the strand's contrribution to the total stress tensor of the network.
        Stress tensor is divided by NokT.
        :return: ZZ-component of the strand's contrribution to the total stress tensor of the network, as float
        """
        return self.b * self.qz * self.qz / (1.0 - self.str_len_sqr())

    # YX component of total stress tensor
    def stress_YX(self):
        """
        Compute the YX-component of the strand's contrribution to the total stress tensor of the network.
        Stress tensor is divided by NokT.
        :return: YX-component of the strand's contrribution to the total stress tensor of the network, as float
        """
        return self.b * self.qy * self.qx / (1.0 - self.str_len_sqr())

    def generate_eq_ensemble(self, n=1):
        """
        Create an ensemble of size n number of strands, where the distribution of strand lengths will
        fit the equilibrium distribution.
        :param n: Number of strands to generate for the ensemble, as int
        :return: List of strands in the equilibrium ensemble, as [Strand objects]
        """
        e = []
        # jeq: Normalization constant for equilibriium FENE Dumbbell distribution function DPL2 Eq. (L) of Table 11.5-1
        jeq = (1. / (2. * pi * beta(3. / 2., (self.b + 2.) / 2.)))
        for s in range(n):
            p1 = 1.0
            p2 = 0.0
            # r, theta, phi: Internal coordinates of a network strand in spherical coordinates used to find equilibrium distribution.
            r, theta, phi = 0.0, 0.0, 0.0
            while p1>p2:
                r = uniform(0.0,1.0)
                theta = pi * uniform(0.0, 1.0)
                phi = 2.0 * pi * uniform(0.0, 1.0)
                p1 = uniform(0.0, 1.0)
                p2 = jeq * ((1.0 - r * r) ** (self.b / 2.0)) * r * r * sin(theta)
            x = r * sin(theta) * cos(phi)
            y = r * sin(theta) * sin(phi)
            z = r * cos(theta)
            e.append(FENEStrand(x, y, z, self.b, self.n, self.mm))
        return e


def beta(z, w):
    """
    Compute beta function of z and w.

    :param z: Parameter z of beta function, as float
    :param w: Parameter w of beta function, as float
    :return: Beta function of z and w, as float
    """
    return exp(lgamma(z)+lgamma(w)-lgamma(z+w))