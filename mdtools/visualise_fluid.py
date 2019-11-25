from mdtools.stat_quantities import FileNaming
import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


class ParticleVisualisation(FileNaming):
    def __init__(self, steps, particles):
        super().__init__(steps, particles)
        # Do some initialisation here if needed

    def particle_plot(self, sim_name, rho, t, power=None, par_a=None):
        """
        Creates a 3D plot for the particles in the fluid.
        The colormap depicts the z-position of the particle.

        @:param sim_name: simulation name used as the prefix in the log files
        @:param rho: Density
        @:param t: Temperature
        @:param power: Pair potential strength
        @:param par_a: Softening @:parameter
        @:return: Nothing. Simply adds a plot on the corresponding canvas
        """
        file_id = self.file_searcher(rho, t, power, par_a)
        data = f"{sim_name}Positions_Velocities{file_id}.txt"

        rx, ry, rz = np.loadtxt(data, usecols=(0, 1, 2), delimiter='\t',
                                comments='#', unpack=True)
        name = f"\N{GREEK SMALL LETTER RHO}: {self.rho_str} T: {self.t_str}" \
               f" n: {self.n_str} A: {self.a_str}"

        fig = plt.figure('3D Scatter Plot')
        ax = fig.add_subplot(111, projection='3d')
        s = ax.scatter(rx, ry, rz, c=rz, cmap='gnuplot_r', label=name)
        plt.colorbar(s)
        plt.legend(loc='best', fancybox=True)

    def vector_field(self, sim_name, rho, t, power=None, par_a=None):
        """
        Creates a 2D projection of the of the loaded files of the fluid for
        position and velocities.

        @:param sim_name: simulation name used as the prefix in the log files
        @:param rho: Density
        @:param t: Temperature
        @:param power: Pair potential strength
        @:param par_a: Softening @:parameter
        @:return: Nothing. Simply adds a plot on the corresponding canvas
        """
        file_id = self.file_searcher(rho, t, power, par_a)
        data = f"{sim_name}Positions_Velocities{file_id}.txt"

        rx, ry, rz, vx, vy, vz = np.loadtxt(data,
                                            # redundant
                                            usecols=(0, 1, 2, 3, 4, 5),
                                            delimiter='\t',
                                            comments='#',
                                            unpack=True)

        name = f"\N{GREEK SMALL LETTER RHO}: {self.rho_str} T: {self.t_str}" \
               f" n: {self.n_str} A: {self.a_str}"

        plt.figure('2D Vector Field of particles')
        q = plt.quiver(rx, ry, vx, vy, rz, pivot='mid',
                       cmap=cm.get_cmap('gnuplot_r'), alpha=0.75, label=name)
        # plt.scatter(rx, ry, alpha=0.4, label=name)
        plt.colorbar(q)
        plt.legend(loc="best")

    # 3D visualisation of the fluid with vector arrows
    def vector_field_3d(self, sim_name, rho, t, power=None, par_a=None):
        """
        Creates a 3D projection based on the last iteration of the MD algorithm
        of the fluids last position and velocities on a vector map.

        @:param sim_name: simulation name used as the prefix in the log files
        @:param rho: Density
        @:param t: Temperature
        @:param power: Pair potential strength
        @:param par_a: Softening @:parameter
        @:return: Nothing. Simply adds a plot on the corresponding canvas
        """
        file_id = self.file_searcher(rho, t, power, par_a)
        data = f"{sim_name}Positions_Velocities{file_id}.txt"

        rx, ry, rz, vx, vy, vz = np.loadtxt(data,
                                            # redundant
                                            usecols=(0, 1, 2, 3, 4, 5),
                                            delimiter='\t',
                                            comments='#',
                                            unpack=True)

        name = f"\N{GREEK SMALL LETTER RHO}: {self.rho_str} T: {self.t_str}" \
               f" n: {self.n_str} A: {self.a_str}"
        fig = plt.figure('3D Vector Field of particles')
        ax = fig.gca(projection='3d')
        # v = np.sqrt(vx**2 + vy**2 + vz**2)
        # x[startAt:endBefore:skip]
        stride = 1
        ax.view_init(elev=18, azim=30)  # camera elevation and angle
        ax.dist = 8  # camera distance
        q = ax.quiver(rx[::stride], ry[::stride], rz[::stride],
                      vx[::stride], vy[::stride], vz[::stride],
                      cmap=cm.get_cmap('viridis'), label=name, alpha=0.7,
                      normalize=True, pivot='middle')
        q.set_array(rz[::stride])
        # m = cm.ScalarMappable(cmap=cm.jet)
        # m.set_array(rz)
        fig.colorbar(q, cmap=cm.get_cmap('viridis'))
        plt.legend(loc='best')
