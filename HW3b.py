import math
from math import pi, sqrt
from SimpsonNumInt import Simpson


# we will be using Table A9 for comparison

def FZ(m, u, Km):
    """Function FZ computes the area below the t-distribution"""

    fcn = (1 + ((u ** 2) / m)) ** (-m - 1) / 2
    Fz = Km * Simpson(fcn, 0, -100, u, npoints=21)
    return Fz


def gamma(alpha):
    """ Gamma Function computes the function for a positive m"""

    if alpha % 1 == 0:
        g = 1
        for i in range(1, int(alpha)):
            g *= i
        return g

    def fn(args):
        t, a = args
        return math.exp(-t) * math.pow(t, a - 1)

    g = Simpson(fn, alpha, 0, 50, 100000)
    return g


def km(m, g):
    """ Function to find Km constant as a function of gamma and degrees of freedom"""

    Km = g * ((0.5 * m) + 0.5) / ((sqrt(m * pi)) * g * (0.5 * m))
    return Km


def main():
    """Main function for finding the  probability function"""

    getout = False
    while getout is False:
        m = input("Degrees of freedom (integer): ")
        u = input("Upper integration limit (float):")
        m = int(m)
        u = float(u)
        Km = km(m, gamma(m))
        Fz = FZ(m, u, Km)
        print("F({:0.3f})={:0.3f}".format(u, Fz))
        getout = input("Go Again (Y/N)?") == "N"
    pass


if __name__ == '__main__':
    main()
