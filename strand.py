from math import pow, sqrt, exp, lgamma, pi, sin, cos
from random import uniform


class Strand:
    """
    Represent a FENE (finitely extensible non-linear elastic) strand in a network model.

    Non-affine Motion prefactor is  1-(Q/Qo)^n. No strand may exceed the nondimensional length of one.
    The FENE force law is used. The loss rate is (1/Lambdao)(1)/(1-(Q/Qo)^n*m).

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
        # TODO: Enforce that Strand length should not exceed 1
        self.qx=qx
        self.qy=qy
        self.qz=qz

    def str_len_sqr(self):
        """
        Compute the length of the Strand, squared. Non-dimensionalized by dividing by maximum strand length Qo, squared.
        :return: Strand length squared, as float
        """
        return pow(self.qx, 2)+pow(self.qy, 2)+pow(self.qz, 2)

    def loss_rate(self, n=2.0, mm=1.0):
        """
        Compute the loss_rate for the Strand. Loss rate is multiplied by the loss rate constant, Lambdao.

        :param n: Exponent on Q/Qo in the non-affine motion prefactor. Should be equal to 2, as float
        :param mm: mm*n is exponent on (Q/Qo) in the loss rate. mm should be equal to 1, as float
        :return: Loss rate for strand, multiplied by the loss rate constant, Lambdao, as float
        """
        return 1.0 / (1.0 - pow(sqrt(self.str_len_sqr()), (mm * n)))

    def loss_prob(self, eps, n=2.0, mm=1.0):
        """
        Compute the probability of the Strand disentangling from the network during a time step of length eps.

        :param eps: Time step size, divided by the loss rate constant, Lamdao, as float
        :param n: Expoent on Q/Qo in the non-affine motion prefactor. Should be equal to 2, as float
        :param mm: mm*n is exponent on (Q/Qo) in the loss rate. mm should be equal to 1, as float
        :return: Probability of the Strand disentangling from the network during a time step of length eps, as float
        """
        return 1.0 - exp(-1.0 * self.loss_rate(n, mm) * eps)

    # XX component of total stress tensor
    def stress_XX(self, b):
        """
        Compute the XX-component of the strand's contrribution to the total stress tensor of the network.
        Stress tensor is divided by NokT.

        :param b: The FENE Parameter (HQo^2/kT), or non-dimensional maximum strand length, as float
        :return: XX-component of the strand's contrribution to the total stress tensor of the network, as float
        """
        return b * self.qx * self.qx / (1.0 - self.str_len_sqr())

    # YY component of total stress tensor
    def stress_YY(self, b):
        """
        Compute the YY-component of the strand's contrribution to the total stress tensor of the network.
        Stress tensor is divided by NokT.

        :param b: The FENE Parameter (HQo^2/kT), or non-dimensional maximum strand length, as float
        :return: YY-component of the strand's contrribution to the total stress tensor of the network, as float
        """
        return b * self.qy * self.qy / (1.0 - self.str_len_sqr())

    # ZZ component of total stress tensor
    def stress_ZZ(self, b):
        """
        Compute the ZZ-component of the strand's contrribution to the total stress tensor of the network.
        Stress tensor is divided by NokT.

        :param b: The FENE Parameter (HQo^2/kT), or non-dimensional maximum strand length, as float
        :return: ZZ-component of the strand's contrribution to the total stress tensor of the network, as float
        """
        return b * self.qz * self.qz / (1.0 - self.str_len_sqr())


    # YX component of total stress tensor
    def stress_YX(self, b):
        """
        Compute the YX-component of the strand's contrribution to the total stress tensor of the network.
        Stress tensor is divided by NokT.

        :param b: The FENE Parameter (HQo^2/kT), or non-dimensional maximum strand length, as float
        :return: YX-component of the strand's contrribution to the total stress tensor of the network, as float
        """
        return b * self.qy * self.qx / (1.0 - self.str_len_sqr())


def beta(z, w):
    """
    Compute beta function of z and w.

    :param z: Parameter z of beta function, as float
    :param w: Parameter w of beta function, as float
    :return: Beta function of z and w, as float
    """
    return exp(lgamma(z)+lgamma(w)-lgamma(z+w))


# Generate and equilibrium ensemble of strands, returned as a list, using a rejection technique
#             b               The FENE Parameter (HQo^2/kT)
def generate_eq_ensemble(n=1, b=1.0):
    """
    Create an ensemble of size n number of strands, where the distribution of strand lengths will
    fit the equilibrium distribution.

    :param n: Number of strands to generate for the ensemble, as int
    :param b: The FENE Parameter (HQo^2/kT), or non-dimensional maximum strand length, as float
    :return: List of strands in the equilibrium ensemble, as [Strand objects]
    """
    e = []
    # jeq: Normalization constant for equilibriium FENE Dumbbell distribution function DPL2 Eq. (L) of Table 11.5-1
    jeq = (1. / (2. * pi * beta(3. / 2., (b + 2.) / 2.)))
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
            p2 = jeq * ((1.0 - r * r) ** (b / 2.0)) * r * r * sin(theta)
        x = r * sin(theta) * cos(phi)
        y = r * sin(theta) * sin(phi)
        z = r * cos(theta)
        e.append(Strand(x, y, z))
    return e


def write_ensemble_to_file(e, filename):
    """
    Write a CSV formatted text file with the internal coordinates and length of each strand  in list.

    :param e: List of strands (the "ensemble") to output to file, as [Strand objects].
    :param filename: Name of file to write to, as string.
    :return: None.
    """
    import os
    f = open(filename, 'w')
    f.write('qx, qy, qz, length\n')
    i = 1
    for s in e:
        f.write('%i, %f, %f, %f, %f\n' % (i, s.qx, s.qy, s.qz, sqrt(s.str_len_sqr())))
        i = i +1
    f.close()


# Compute stress component members for ensemble
def ensemble_stress(e, b):
    """
    Compute the components of the total stress tensor for the network,
    by summing over the strands in the network ensemble.  Stress tensor is divided by NokT.

    :param e: List of strands (the "ensemble"), as [Strand objects].
    :param b: The FENE Parameter (HQo^2/kT), or non-dimensional maximum strand length, as float
    :return: Tuple (XX, YY, ZZ, YX) of components of total stress tensor of the network, divided by NokT, as (float, float, float, float)
    """
    piXX = 0.0
    piYY = 0.0
    piZZ = 0.0
    piYX = 0.0
    for s in e:
        piXX = piXX + s.stress_XX(b)
        piYY = piYY + s.stress_YY(b)
        piZZ = piZZ + s.stress_ZZ(b)
        piYX = piYX + s.stress_YX(b)
    return (piXX, piYY, piZZ, piYX)


# Compute average strand length for the ensemble
def ensemble_q_ave(e):
    """
    Compute the average length of a strand in the network ensemble, divided by the maximum strand length.

    :param e: List of strands (the "ensemble"), as [Strand objects].
    :return: Average length of a strand in the network ensemble, divided by the maximum strand length, as float
    """
    total_q = 0.0
    for s in e:
        total_q = total_q + sqrt(s.str_len_sqr())
    return total_q / len(e)