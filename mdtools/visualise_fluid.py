from mdtools.stat_quantities import FileNaming
import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation


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
        data = f"{sim_name}Positions_Velocities{file_id}.log"

        rx, ry, rz = np.loadtxt(data, usecols=(0, 1, 2), delimiter='\t',
                                comments='#', unpack=True)
        name = self.get_label(file_id)
        fig = plt.figure('3D Scatter Plot')
        ax = fig.add_subplot(111, projection='3d')
        s = ax.scatter(rx, ry, rz, c=rz, cmap='viridis_r', label=name)
        ax.legend(loc='best', fancybox=True)
        fig.colorbar(s)

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
        data = f"{sim_name}Positions_Velocities{file_id}.log"

        rx, ry, rz, vx, vy, vz = np.loadtxt(data,
                                            # redundant
                                            usecols=(0, 1, 2, 3, 4, 5),
                                            delimiter='\t',
                                            comments='#',
                                            unpack=True)

        name = self.get_label(file_id)
        plt.figure('2D Vector Field of particles')
        q = plt.quiver(rx, ry, vx, vy, rz, pivot='mid',
                       cmap=cm.get_cmap('viridis_r'), alpha=0.75, label=name)
        plt.legend(loc="best")
        plt.colorbar(q)

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
        data = f"{sim_name}Positions_Velocities{file_id}.log"

        rx, ry, rz, vx, vy, vz = np.loadtxt(data,
                                            # redundant
                                            usecols=(0, 1, 2, 3, 4, 5),
                                            delimiter='\t',
                                            comments='#',
                                            unpack=True)

        name = self.get_label(file_id)
        fig = plt.figure('3D Vector Field of particles')
        ax = fig.gca(projection='3d')
        stride = 1
        ax.view_init(elev=18, azim=30)  # camera elevation and angle
        ax.dist = 8  # camera distance
        q = ax.quiver(rx[::stride], ry[::stride], rz[::stride],
                      vx[::stride], vy[::stride], vz[::stride],
                      cmap=cm.get_cmap('viridis'), label=name, alpha=0.7,
                      normalize=True)
        q.set_array(rz[::stride])
        fig.colorbar(q, cmap=cm.get_cmap('viridis'))
        plt.legend(loc='best')

    def animation3D(self, sim_name, rho, t, power=None, par_a=None, save=False):

        file_id = self.file_searcher(rho, t, power, par_a)

        x_all = np.loadtxt(f"{sim_name}x_data{file_id}.log")
        y_all = np.loadtxt(f"{sim_name}y_data{file_id}.log")
        z_all = np.loadtxt(f"{sim_name}z_data{file_id}.log")

        fig = plt.figure(figsize=(10, 10))
        ax = fig.add_subplot(111, projection="3d")
        graph = ax.scatter(x_all[0], y_all[0], z_all[0], color='orange')
        text = fig.text(0, 1, "TEXT", va='top')  # displays the frames

        box_length = (int(self.p_str) / rho) ** (1./3.)

        ax.set_xlim3d(0, box_length)
        ax.set_ylim3d(0, box_length)
        ax.set_zlim3d(0, box_length)

        # Mock function for animation frame generation
        def update_figure(num):
            x = x_all[num]
            y = y_all[num]
            z = z_all[num]
            # TODO: NN calculation between all particles & -> assign colormap
            text.set_text(f"Frame: {num}")
            graph._offsets3d = (x, y, z)
            return graph

        ani = FuncAnimation(fig, update_figure,
                            frames=int(self.steps_str), interval=0, blit=False)
        if save:
            ani.save(f"{sim_name}.mp4", fps=60)

        plt.show()
