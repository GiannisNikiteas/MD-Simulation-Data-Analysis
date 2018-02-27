import numpy as np
from PathHandling import *


def rdf(power, par_a):
    """
    A function meant to compare multiple RDF files and find the intersects.

    :param power: Potential strength n, of the RDF
    :param par_a: Parameter A, of the potential
    :return: numpy array that contains the RDF points
    """
    self.n_str = str(int(power))
    A = '{:.2f}'.format(par_a)
    HIST = 'Hist' + self.n_str + '~' + A + '.txt'
    num_lines = sum(1 for line in open(HIST))
    Hist = np.loadtxt(HIST, skiprows=1, delimiter='\n')
    Hist = np.delete(Hist, 99)
    return Hist


def find_nearest(array, value):
    """
    Finds the closest value for an array of floats
    :param array: Array to search [INPUT]
    :param value: Value to find [INPUT]
    :return: Index of value in the array [OUTPUT]
    """
    idx = (np.abs(array-value)).argmin()
    return idx


path = OSPaths()
path.dir('Density 0.5', '5000')

rdf1= rdf(6, 0)
rdf2 = rdf(12, 0)

final = np.intersect1d(rdf1, rdf2)
print(find_nearest(rdf1, final[1]))



