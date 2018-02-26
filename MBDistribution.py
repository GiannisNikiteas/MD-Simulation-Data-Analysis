import numpy as np
from scipy.interpolate import interp1d
from scipy.special import erf


RANGE_OF_SPEEDS_MONTE_CARLO = np.arange(0, 5, 0.001)
K_BOLTZMANN = 1.
TEMPERATURE = 1.4
MASS = 1.
NUMBER_OF_PARTICLES = 10**3
NUMBER_OF_DIMENSIONS = 3


def mb_speed_distribution():
    """
    Maxwell-Boltzmann's speed distribution
    Takes in speeds range, temperature and mass.
    :return: numpy array of M-B speed distribution
    """
    term1 = (MASS / (2 * np.pi * K_BOLTZMANN * TEMPERATURE)) ** 1.5
    term2 = 4 * np.pi * RANGE_OF_SPEEDS_MONTE_CARLO ** 2
    term3 = np.exp(-MASS * RANGE_OF_SPEEDS_MONTE_CARLO ** 2 / (2 * K_BOLTZMANN * TEMPERATURE))
    return term1 * term2 * term3


def mb_cdf():
    """
    Maxwell-Boltzmann's speed distribution Cumulative Distribution Function (CDF)
    Takes in speeds range, temperature and mass.
    :return: numpy array of CDF
    """
    term0 = np.sqrt(K_BOLTZMANN * TEMPERATURE / MASS)
    term1 = erf(RANGE_OF_SPEEDS_MONTE_CARLO / (np.sqrt(2) * term0))
    term2 = np.sqrt(2 / np.pi) * RANGE_OF_SPEEDS_MONTE_CARLO
    term3 = np.exp(-RANGE_OF_SPEEDS_MONTE_CARLO ** 2 / (2 * term0 ** 2)) / term0
    return term1 - term2 * term3


class VelGen:
    def __init__(self):
        # Generates random speeds according to
        # Maxwell-Boltzmann's speed distribution
        if TEMPERATURE == 0:
            self.speeds = np.zeros(NUMBER_OF_PARTICLES)
        else:
            # Create interpolation function to CDF
            inv_cdf = interp1d(mb_cdf(), RANGE_OF_SPEEDS_MONTE_CARLO)
            # Generate random speeds
            self.speeds = inv_cdf(np.random.random(NUMBER_OF_PARTICLES))

    def get_velocities(self):
        """
        Converts random speeds into their x, y, z components.
        Takes in random speeds
        :return: numpy arrays of speeds and their x,y,z components
        """
        theta = np.arccos(np.random.uniform(-1, 1, NUMBER_OF_PARTICLES))
        phi = np.random.uniform(0, 2 * np.pi, NUMBER_OF_PARTICLES)
        vx = self.speeds * np.cos(phi)
        vy = self.speeds * np.sin(phi)
        vz = self.speeds * np.cos(theta)
        if NUMBER_OF_DIMENSIONS == 3:
            vx *= np.sin(theta)
            vy *= np.sin(theta)
            return np.array([vx, vy, vz], dtype=np.double)
        return np.array([vx, vy], dtype=np.double)


a = VelGen()
a, b, c = a.get_velocities()
np.savetxt('../MD simulation/MD simulation/vx.txt', a, delimiter='\n')
np.savetxt('../MD simulation/MD simulation/vy.txt', b, delimiter='\n')
np.savetxt('../MD simulation/MD simulation/vz.txt', c, delimiter='\n')




















