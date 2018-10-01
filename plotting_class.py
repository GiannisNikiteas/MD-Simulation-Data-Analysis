import numpy as np
import matplotlib.cm as cm
import scipy.stats as stats
from scipy import interpolate
import matplotlib.pyplot as plt
import statistics as stat
from mpl_toolkits.mplot3d import Axes3D


class FilePlotting:
    """
    Class responsible for Plotting the data generated by the MD simulation
    software written in C++.
    Constructor purposely does not have any arguments
    """

    def __init__(self, __step, __particles):
        self.step_str = str(__step)
        self.particles_str = str(__particles)
        self.rho_str = None
        self.t_str = None
        self.n_str = None
        self.a_str = None
        self.dif_coef = np.array([])
        self.reduced_dif_coef = np.array([])
        self.dif_err = np.array([])
        self.reduced_dif_err = np.array([])
        self.dif_y_int = np.array([])
        self.reduced_dif_y_int = np.array([])
        self.line_style = ['solid', 'dashed', 'dotted', 'dashdot']  # TODO: Fix with itertools
        self.interpolated_data = []             # TODO: this is a butcher's approach to solving the problem
        self.dr = None
        # TODO: get rid of all these with itertools, This is python not C++
        # This is an iterator for the color array
        self.p, self.c = 0, 0
        self.j = 0  # stride for MSD
        self.v = 0  # index for MSD
        self.line_it = 0    # Index iterator for line styles

    def file_searcher(self, rho, t, n, a=None):
        """
        Constructs the file signature of the MD simulation in order for the information to be read
        :param rho: density
        :param t:   temperature
        :param n:   potential strength
        :param a:   softness parameter
        :return: string with file identifier
        """
        p_str = self.particles_str
        self.rho_str = "{:.4f}".format(rho)
        self.t_str = "{:.4f}".format(t)
        self.n_str = str(n)
        alpha = None

        if a is not None:
            self.a_str = "{:.5f}".format(a)
            alpha = '_A_' + self.a_str
        else:
            self.a_str = ""
            alpha = self.a_str

        name_id = '_step_' + self.step_str + '_particles_' + p_str + '_rho_' + \
                  self.rho_str + '_T_' + self.t_str + '_n_' + self.n_str + alpha
        return name_id

    def energy_plots(self, rho, t, power, par_a):
        """
        Creates a figure where the average kinetic, potential and total energy are displayed.
        Separately and in a combined graph.
        :param rho: Density
        :param t: Temperature
        :param power: Pair potential strength
        :param par_a: Softening parameter
        :return: Nothing. Simply adds a plot on the corresponding canvas
        """
        file_id = self.file_searcher(rho, t, power, par_a)
        data = "Data" + file_id + ".txt"

        num_lines = 0
        # Measuring the line number in a file ignoring comments
        for line in open(data):
            if line.startswith('#'):
                continue
            num_lines += 1

        pot_en, kin_en = np.loadtxt(data, usecols=(3, 4), delimiter='\t', comments='#', unpack=True)
        tot_en = pot_en + kin_en
        #  Plots the Energies
        step = 0.005
        time = num_lines * step
        x = np.linspace(0, time, num_lines)
        fig = plt.figure('Energy Plots')

        kin_f = plt.subplot2grid((3, 2), (0, 0), colspan=1)
        pot_f = plt.subplot2grid((3, 2), (1, 0), colspan=1)
        tot_f = plt.subplot2grid((3, 2), (2, 0), colspan=1)
        all_f = plt.subplot2grid((3, 2), (0, 1), rowspan=3)

        # k, u, t = kin_en[500], pot_en[500], tot_en[500]
        kin_f.plot(x, kin_en, 'r')
        kin_f.locator_params(axis='y', nbins=4), kin_f.set_ylim(ymax=4)
        pot_f.plot(x, pot_en, 'g')
        pot_f.locator_params(axis='y', nbins=3)
        pot_f.set_ylabel("Energy units", size=16)
        tot_f.plot(x, tot_en, 'b')
        tot_f.locator_params(axis='y', nbins=4)
        tot_f.set_ylim(ymax=6)

        # x_r = time / 2 - time / 4
        # Kin.set_title('Individual Plots n = %d' %power, fontsize=17)
        kin_f.set_title('Individual Plots', fontsize=17)
        all_f.set_title('Energy Contributions', fontsize=17)
        all_f.set_xlabel(r"Time $t$", fontsize=16)

        tot_f.set_xlabel(r"Time $t$", fontsize=16)
        fig.subplots_adjust(hspace=0)

        # Tick correction
        for ax in [kin_f, pot_f]:
            plt.setp(ax.get_xticklabels(), visible=False)
            # The y-ticks will overlap with "hspace=0", so we'll hide the bottom tick
            ax.set_yticks(ax.get_yticks()[1:])

        all_f.plot(x, kin_en, 'r', x, pot_en, 'g', x, tot_en, 'b')
        all_f.set_ylim(ymax=5)

    def potential_data(self, rho, t, power, par_a):
        """
        Creates plots for the visualisation of the average potential energy of the fluid.
        :param rho: Density
        :param t: Temperature
        :param power: Pair potential strength
        :param par_a: Softening parameter
        :return: Nothing. Simply adds a plot on the corresponding canvas
        """
        file_id = self.file_searcher(rho, t, power, par_a)
        data = "Data" + file_id + ".txt"

        num_lines = 0
        for line in open(data):
            if line.startswith('#'):
                continue
            num_lines += 1

        rho_list, u = np.loadtxt(data, usecols=(1, 3), delimiter='\t', comments='#', unpack=True)

        #  Plots the Energies
        name = "rho: " + self.rho_str + "T: " + \
               self.t_str + "n: " + self.n_str + "A: " + self.a_str
        step = 0.005
        time = num_lines * step
        x = np.linspace(0, time, num_lines)
        plt.figure('Potential Plots of Data')
        plt.plot(rho_list, u, label=name)
        plt.legend(loc='best', fancybox=True)

    def particle_plot(self, rho, t, power, par_a):
        """
        Creates a 3D plot for the particles in the fluid.
        The colormap depicts the z-position of the particle.
        :param rho: Density
        :param t: Temperature
        :param power: Pair potential strength
        :param par_a: Softening parameter
        :return: Nothing. Simply adds a plot on the corresponding canvas
        """
        file_id = self.file_searcher(rho, t, power, par_a)
        data = "Positions_Velocities" + file_id + ".txt"

        rx, ry, rz = np.loadtxt(data, usecols=(0, 1, 2), delimiter='\t',
                                comments='#', unpack=True)
        name = "rho: " + self.rho_str + "T: " + \
               self.t_str + "n: " + self.n_str + "A: " + self.a_str
        fig = plt.figure('3D Scatter Plot')
        ax = fig.add_subplot(111, projection='3d')
        s = ax.scatter(rx, ry, rz, c=rz, cmap='gnuplot_r', label=name)
        plt.colorbar(s)
        plt.legend(loc='best', fancybox=True)

    def vector_field(self, rho, t, power, par_a):
        """
        Creates a 2D projection of the of the loaded files of the fluid for
        position and velocities
        :param rho: Density
        :param t: Temperature
        :param power: Pair potential strength
        :param par_a: Softening parameter
        :return: Nothing. Simply adds a plot on the corresponding canvas
        """
        file_id = self.file_searcher(rho, t, power, par_a)
        data = "Positions_Velocities" + file_id + ".txt"

        rx, ry, rz, vx, vy, vz = np.loadtxt(data,
                                            usecols=(0, 1, 2, 3, 4, 5),  # redundant
                                            delimiter='\t',
                                            comments='#',
                                            unpack=True)

        plt.figure('2D Vector Field of particles')
        name = "rho: " + self.rho_str + "T: " + \
               self.t_str + "n: " + self.n_str + "A: " + self.a_str
        q = plt.quiver(rx, ry, vx, vy, rz, pivot='mid',
                       cmap=cm.gnuplot_r, alpha=0.75, label=name)
        # plt.scatter(rx, ry, alpha=0.4, label=name)
        plt.colorbar(q)
        plt.legend(loc="best")

    # 3D visualisation of the fluid with vector arrows
    def vector_field_3d(self, rho, t, power, par_a):
        """
        Creates a 3D projection based on the last iteration of the MD algorithm
        of the fluids last position and velocities on a vector map
        :param rho: Density
        :param t: Temperature
        :param power: Pair potential strength
        :param par_a: Softening parameter
        :return: Nothing. Simply adds a plot on the corresponding canvas
        """
        file_id = self.file_searcher(rho, t, power, par_a)
        data = "Positions_Velocities" + file_id + ".txt"

        rx, ry, rz, vx, vy, vz = np.loadtxt(data,
                                            usecols=(0, 1, 2, 3, 4, 5),  # redundant
                                            delimiter='\t',
                                            comments='#',
                                            unpack=True)
        fig = plt.figure('3D Vector Field of particles')
        ax = fig.gca(projection='3d')
        # v = np.sqrt(vx**2 + vy**2 + vz**2)
        name = "rho: " + self.rho_str + "T: " + \
               self.t_str + "n: " + self.n_str + "A: " + self.a_str
        # x[startAt:endBefore:skip]
        stride = 1
        ax.view_init(elev=18, azim=30)  # camera elevation and angle
        ax.dist = 8  # camera distance
        q = ax.quiver(rx[::stride], ry[::stride], rz[::stride],
                      vx[::stride], vy[::stride], vz[::stride],
                      cmap=cm.jet, label=name, alpha=0.7,
                      normalize=True, pivot='middle')
        q.set_array(rz[::stride])
        # m = cm.ScalarMappable(cmap=cm.jet)
        # m.set_array(rz)
        fig.colorbar(q, cmap=cm.jet)
        plt.legend(loc='best')

    # Radial Distribution Function
    def rdf(self, rho, t, power, par_a, iso_scale=False, show_iso=False):
        """
        Creates a plot for the Radial Distribution Function of the fluid, which depicts
        the microscopic density fluctuations of the molecules as a function of distance.
        The parameters iso_scale and show_iso are optional and should normally be set to False,
        unless the very specific isosbestic behaviour is examined.


        :param rho: Density
        :param t: Temperature
        :param power: Pair potential strength
        :param par_a: Softening parameter
        :param iso_scale: Optional, scales the values of r, the radius, based on the isosbestic model that
                          has been created
        :param show_iso: Optional, shows the location of the theoretical isosbestic point between different rdfs
        :return: The numpy.array for the RDF data

        """
        file_id = self.file_searcher(rho, t, power, par_a)
        data = "RDF" + file_id + ".txt"
        num_lines = sum(1 for line in open(data))
        rdf = np.loadtxt(data, delimiter="\t", usecols=1, comments="#")

        # Number of particles, Number of bins
        particles, bins = int(self.particles_str), 500
        # Cut off radius
        rg = 3.0   # ((particles / rho) ** (1. / 3.)) / 2.
        dr = rg / bins

        # Initial normalisation factor of

        # Normalise the raw RDF for long distances
        # for i in range(1, bins-1):
        #     norm = particles * \
        #            int(self.step_str) * \
        #            float(self.rho_str) *\
        #            2*np.pi/3.0 *\
        #            ((rg*i/bins + dr/2.0)**3 - (rg*i/bins - dr/2.0)**3)
        #     rdf[i] /= norm

        # r range, r=0 is intentionally neglected due to division by 0
        # num_lines-1 because one of them is a comment
        num_lines -= 2
        r = np.linspace(1, num_lines, num_lines)
        r = np.multiply(r, dr)
        a_tilde = par_a * rho ** (1./3.)    # Scale a

        # Isomorphic scaling of r for the isomorph plane
        if iso_scale is True:
            r = np.multiply(r, rho ** (1. / 3.))  # Scale r

        plt.figure('Radial Distribution Function')

        # Plotting isosbestic point
        max_scaling = np.max(rdf)  # Scaling the ymax
        if show_iso is True:    # Show isosbestic point
            iso = np.sqrt(1 - par_a ** 2)   # TODO: this is not correct, missing a factor probably, revise theory!
            # iso = np.sqrt(rho ** (2./3) - a_tilde ** 2)   # This is probably wrong
            plt.plot([iso, iso], [0, max_scaling + 0.1], '--', color='red')

        name = "rho: " + self.rho_str + " T: " + self.t_str + " n: " + self.n_str + " A: " + self.a_str
        plt.plot(r, rdf, '-o', markersize=4, label=name)

        # Plot labels
        plt.xlabel(r"$r$", fontsize=16)
        plt.ylabel(r"$g(r)$", fontsize=16)
        # Plotting isosbestic location of points

        # Line through y = 1
        plt.plot([0, r[-1]], [1, 1], '--', color='black', linewidth=0.5)
        # Plot limits and legends
        plt.xlim(xmin=0, xmax=3)
        plt.ylim(ymin=0, ymax=max_scaling + 0.1)
        plt.legend(loc="best", fancybox=True, prop={'size': 8})
        self.c += 1

        # return the plotting lists
        return r, rdf

    # Velocity Autocorrelation Function
    def vaf(self, rho, t, power, par_a):
        """
        Creates a figure for the Velocity Autocorrelation Function of the fluid, which illustrates
        if the fluid remains a coupled system through time (correlated) or it uncouples.
        :param rho: Density
        :param t: Temperature
        :param power: Pair potential strength
        :param par_a: Softening parameter
        :return: Nothing. Simply adds a plot on the corresponding canvas
        """
        file_id = self.file_searcher(rho, t, power, par_a)
        data = "Data" + file_id + ".txt"

        cr = np.loadtxt(data, usecols=8, delimiter='\t', unpack=True, comments='#')

        num_lines = 0
        for line in open(data):
            if line.startswith('#'):
                continue
            num_lines += 1

        time_step = 0.005 / np.sqrt(t)
        time = time_step * num_lines
        x = np.linspace(0, time, num_lines)
        t_tilde = x * rho ** (1./3.) * t ** 0.5
        name = ""
        if num_lines < 100000:
            name = "rho: " + self.rho_str + "T: " + \
                   self.t_str + "n: " + self.n_str + "A: " + self.a_str

        plt.figure('Velocity Autocorrelation Function')
        y = np.full(num_lines, 0)
        # xx = np.full(num_lines, time)
        # yy = np.linspace(5, -0.5, num_lines)
        # plt.plot(xx, yy, '--', color='black')
        plt.plot(t_tilde, y, '--', color='black')
        plt.plot(t_tilde, cr, label=name)
        plt.xlabel(r"Time $t$", fontsize=16)
        plt.ylabel(r"$C_v$", fontsize=16)
        plt.ylim(ymax=5, ymin=-0.5)
        plt.xlim(xmin=t_tilde[0], xmax=t_tilde[-1])
        plt.legend(loc="best", ncol=1, borderpad=0.1,
                   labelspacing=0.01, columnspacing=0.01, fancybox=True,
                   fontsize=12)

    # Mean Square Displacement
    def msd(self, rho, t, power, par_a):
        """
        Creates a figure which depicts the Mean Square Displacement for our fluid.
        According to diffusion theory the slope of the MSD corresponds to the inverse of the
        diffusion coefficient
        :param rho: Density
        :param t: Temperature
        :param power: Pair potential strength
        :param par_a: Softening parameter
        :return: Nothing. Simply adds a plot on the corresponding canvas
        """
        file_id = self.file_searcher(rho, t, power, par_a)
        data = "Data" + file_id + ".txt"

        rho_list, msd_data = np.loadtxt(data, usecols=(1, 7), delimiter='\t', unpack=True)

        num_lines = 0
        for line in open(data):
            if line.startswith('#'):
                continue
            num_lines += 1

        step = 0.005 / np.sqrt(t)
        limit = step * num_lines
        step = int(0.6 * num_lines)
        x = np.linspace(0, num_lines-1, num=num_lines)

        if par_a >= 0:
            x_sliced = x[step:]
            msd_sliced = msd_data[step:]
            gradient, intercept, r_value, p_value, std_err = stats.linregress(x_sliced, msd_sliced)

            self.reduced_dif_coef = np.append(self.reduced_dif_coef, gradient)
            self.reduced_dif_err = np.append(self.reduced_dif_err, std_err)
            self.reduced_dif_y_int = np.append(self.reduced_dif_y_int, intercept)

        # Regular coefs are calculated independent of the if loop
        gradient, intercept, r_value, p_value, std_err = stats.linregress(x, msd_data)
        self.dif_coef = np.append(self.dif_coef, gradient)
        self.dif_err = np.append(self.dif_err, std_err)
        self.dif_y_int = np.append(self.dif_y_int, intercept)

        # TODO: this should be output to a log file for reference
        print('Diffusion coef: ', gradient, '\n',
              'y-intercept: ', intercept, '\n',
              'R value: ', r_value, '\n',
              'Fit Error: ', std_err)

        name = "rho: " + self.rho_str + " T: " + \
               self.t_str + " n: " + self.n_str + " a: " + self.a_str
        plt.figure('Mean Square Displacement')
        plt.plot(x, msd_data, label=name)
        plt.xlabel(r"$t$", fontsize=16)
        plt.ylabel(r"$MSD$", fontsize=16)
        # plt.xlim(xmin=0, xmax=x[num_lines - 1])
        # plt.ylim(ymin=0, ymax=msd_data[num_lines - 1])
        plt.legend(loc="best", fancybox=True)
        print("@ index: ", np.argmax(msd_data), " value: ", max(msd_data))

    # Pressure C
    def pc(self, rho, t, power, par_a):
        file_id = self.file_searcher(rho, t, power, par_a)
        pc_name = "Data" + file_id + ".txt"

        pc_data = np.loadtxt(pc_name, usecols=5, delimiter='\t',
                             comments='#', unpack=True)
        num_lines = 0
        for line in open(pc_name):
            if line.startswith('#'):
                continue
            num_lines += 1

        time = num_lines * 0.005
        xxx = np.linspace(0, time, num=num_lines)

        name = "rho: " + self.rho_str + "T: " + \
               self.t_str + "n: " + self.n_str + "A: " + self.a_str
        plt.figure('Configurational Pressure')

        plt.plot(xxx, pc_data, label=name)
        plt.xlabel(r"Time $t$", size=18)
        plt.ylabel(r"Configurational Pressure $P_C$", size=18)
        plt.legend(loc="best", prop={'size': 12},
                   borderpad=0.2, labelspacing=0.2, handlelength=1)

    # Averages
    # TODO: Look how AVG files are named and fix
    # TODO: Probably these methods will not be static after fixing
    def avg_pressure(self, rho, t, power):
        """
        Plots the average Configurational (virial) pressure of the fluid throughout
        the entire simulation
        :param rho:
        :param t:
        :param power:
        :return:
        """
        file_id = self.file_searcher(rho, t, power)
        pc_name = "AVGdata" + file_id + ".txt"
        name = "rho: " + self.rho_str + "T: " + self.t_str + "n: " + self.n_str

        num_lines = sum(1 for line in open(pc_name))
        a, pc = np.loadtxt(pc_name, delimiter='\t', comments='#', usecols=(0, 5), unpack=True)

        plt.figure('Average Pressure')
        plt.plot(a, pc, '-o', label=name, markersize=3)
        plt.xlim(xmin=0, xmax=4.0)
        # plt.title("Configurational Pressure for multiple Potentials")
        plt.xlabel(r"$A$", size=16)
        plt.ylabel(r"$P_c$", size=16)
        plt.legend(loc="best")

    def avg_kin(self, rho, t, power):
        """
        Plots the average kinetic energy of the fluid throughout the entire simulation.
        The kinetic energy is supposed to be constant since the fluid is placed in an
        isothermal container, hence the deviations in kinetic energy can be used as an
        error metric
        :param rho: Density
        :param t: Temperature
        :param power: Pair potential strength
        :return: Nothing. Simply adds a plot on the corresponding canvas
        """
        file_id = self.file_searcher(rho, t, power)
        k_name = "AVGdata" + file_id + ".txt"
        name = "rho: " + self.rho_str + "T: " + self.t_str + "n: " + self.n_str

        a, k = np.loadtxt(k_name, delimiter='\t', comments='#', usecols=(0, 2), unpack=True)

        plt.figure('Average Kinetic Energy')
        plt.plot(a, k, label=name)
        plt.title("Kinetic Energy vs multiple values of A")
        plt.xlabel(r"Parameter A")
        plt.ylabel(r"Kinetic Energy $K$")
        plt.legend(loc="best")

    def avg_pot(self, rho, t, power):
        file_id = self.file_searcher(rho, t, power)
        u = "AVGdata" + file_id + ".txt"
        name = "rho: " + self.rho_str + "T: " + self.t_str + "n: " + self.n_str

        a, k = np.loadtxt(u, delimiter='\t', comments='#', usecols=(0, 3), unpack=True)

        plt.figure('Average Potential Energy')
        plt.plot(a, k, label=name)
        plt.title("Potential Energy against parameter A")
        plt.xlabel(r"Parameter A")
        plt.ylabel(r"Potential Energy $U$")
        plt.legend(loc="best")

    def avg_en(self, rho, t, power):
        """
        Plots the average total energy of the fluid throughout the entire simulation
        :param rho: Density
        :param t: Temperature
        :param power: Pair potential strength
        :return: Nothing. Simply adds a plot on the corresponding canvas
        """
        file_id = self.file_searcher(rho, t, power)
        e_name = "AVGdata" + file_id + ".txt"
        # name = "rho: " + self.rho_str + "T: " + self.t_str + "n: " + self.n_str

        num_lines = sum(1 for line in open(e_name))

        a, K, U, T = np.loadtxt(e_name, delimiter='\t', comments='#', usecols=(0, 2, 3, 4), unpack=True)

        fig = plt.figure('Average Energies')
        kin_f = plt.subplot2grid((3, 2), (0, 0), colspan=1)
        pot_f = plt.subplot2grid((3, 2), (1, 0), colspan=1)
        tot_f = plt.subplot2grid((3, 2), (2, 0), colspan=1)
        all_f = plt.subplot2grid((3, 2), (0, 1), rowspan=3)

        kin_f.plot(a, K, color='r')
        kin_f.set_ylim(ymin=2.0), kin_f.locator_params(axis='y', nbins=5, prune='lower')
        pot_f.plot(a, U, color='g')
        pot_f.locator_params(axis='y', nbins=4), pot_f.set_ylabel("Energy units")
        tot_f.plot(a, T, color='c')
        tot_f.locator_params(axis='y', nbins=3), tot_f.locator_params(axis='x', nbins=4)
        tot_f.set_xlabel(r"Parameter $A$")

        kin_f.set_title('Individual Plots for n = %d' % power, fontsize=14)
        all_f.set_title('Total Energy, Kinetic and Potential')
        all_f.set_xlabel(r"Parameter $A$")

        fig.subplots_adjust(hspace=0)

        for ax in [kin_f, pot_f]:
            plt.setp(ax.get_xticklabels(), visible=False)
            # The y-ticks will overlap with "hspace=0", so we'll hide the bottom tick
            ax.set_yticks(ax.get_yticks()[1:])

        all_f.plot(a, K, 'r', a, U, 'g', a, T, 'c')
        all_f.set_ylim(ymax=6)
        all_f.locator_params(axis='x', nbins=4)

    def diffusion_plot(self, rho, t, power, my_list):
        """
        A graph of the Diffusion coefficients D against a list of parameter A values
        :param rho: Density
        :param t: Temperature
        :param power: Potential strength parameter n
        :param my_list: List of parameter A coefficients
        :return: Figure of D vs A for a given number of iterations
        """
        # I recommend not to fuck with this list, when using the GUI
        # otherwise comment out for custom implementation
        # my_list = [0, 0.25, 0.50, 0.65, 0.7, 0.75, 0.8, 0.85, 0.90, 0.97, 1.00, 1.1, 1.2,
        #    1.25, 1.50, 1.75, 2.00, 2.25, 2.50, 2.75, 4.00]
        for i in my_list:
            self.msd(rho, t, power, i)
            print("-----------------------------")
        name = "n: " + str(power)
        # steps = ['5k', '10k', '12.5k', '15k', '20k']

        plt.figure('Diffusion coefficients D vs A')
        plt.plot(my_list, self.dif_coef, '--o', label=name, markersize=3.5)
        plt.xlabel(r"$a$", fontsize=16)
        plt.ylabel(r"$D$", fontsize=16)
        plt.legend(loc="best", fancybox=True, ncol=2)

        self.dif_coef, self.dif_err, self.dif_y_int = np.array([]), np.array([]), np.array([])
        self.reduced_dif_coef, self.reduced_dif_err, self.reduced_dif_y_int = np.array([]), np.array([]), np.array([])
        plt.ylim(ymin=0)
        plt.xlim(xmin=0)
        self.j += 15
        self.v += 1

    # No data plots
    def potential(self, power, par_a, show_inf_p=False):
        """
        A plot of the pair potential used in the simulations.
        It includes the functionality of displaying the isosbestic points
        of the potential
        :param power: Pair potential strength
        :param par_a: Softening parameter
        :param show_inf_p: Show the x-axis origin, for cases where A=0
                           the value has to be False for the graphs to display properly
        :return: Nothing. Simply adds a plot on the corresponding canvas
        """
        self.n_str = str(power)
        a = str(float(par_a))

        x = np.linspace(0, 3, 150)  # TODO:  0.5
        phi = 1 / pow((x ** 2 + par_a ** 2), power / 2)
        plt.figure('Potential')

        name = "n: " + self.n_str + " A: " + a
        # name = "A: " + A
        if par_a <= 1. and self.p == 0:
            iso = np.sqrt(1 - par_a ** 2)
            sak = 1 / pow((iso ** 2 + par_a ** 2), power / 2.0)
            plt.scatter(iso, sak, marker='o', color='magenta', label='Isosbestic point')

        plt.plot(x, phi, label=name)  # , linestyle=self.line_style[self.line_it], color='black')
        plt.xlim(xmin=x[0], xmax=2)
        plt.ylim(ymin=0)

        plt.xlabel(r"$r$", size=16)
        plt.ylabel(r"$\Phi$", size=16)
        plt.legend(loc="best", fancybox=True, ncol=1)
        f, v = [], []
        if show_inf_p is True:
            temp = np.sqrt(par_a / (1. + power))
            f.append(temp)
            temp = 1 / pow((temp ** 2 + par_a), power / 2.)
            v.append(temp)
            plt.scatter(f, v, color='red', marker='o', s=13, label='Inflection point')
            plt.locator_params(axis='x', nbins=5)
        self.p += 1
        self.line_it += 1

    # In experimental stage
    def scaled_potential(self, rho, power, par_a):
        """
        A Scaled potential for the isosbestic point theory
        :param rho:
        :param power:
        :param par_a:
        :return:
        """
        self.n_str = str(power)
        a = str(float(par_a))

        r = np.linspace(0, 3, 150)  # TODO:  0.5
        a_tilde = par_a * rho ** (1./3.)    # Scale a
        r_tilde = np.multiply(r, rho ** (1./3.))  # Scale r

        phi = rho ** (power/3.) * (1 / pow((r_tilde ** 2 + a_tilde ** 2), power/2.))
        plt.figure('Scaled Potential')
        plt.plot(r_tilde, phi, )

    def force(self, power, par_a):
        """
        Plots the force experienced by the molecules
        :param power: Potential strength
        :param par_a: Softening parameter
        :return: Nothing, Simply adds a plot to the corresponding canvas
        """
        self.n_str = str(power)
        a = str(float(par_a))

        x = np.linspace(0.0, 3, 300)  # division with 0
        phi = power * x * (pow((x ** 2 + par_a), (-power / 2 - 1)))
        plt.figure('Force')
        name = " A: " + a
        plt.plot(x, phi, label=name, linestyle=self.line_style[self.line_it], color='black')
        plt.xlabel(r"$r$", size=16)
        plt.ylabel(r"$f$", size=16)
        plt.legend(loc="best", fancybox=True)  # fontsize=16
        f, v = [], []
        temp = np.sqrt(par_a / (1. + power))
        f.append(temp)
        temp = power * temp * (pow((temp ** 2 + par_a), (-power / 2 - 1)))
        v.append(temp)
        plt.scatter(f, v, color='red', marker='o', s=13, label='Inflection point')
        plt.xlim(0, 2)
        plt.ylim(ymin=0)

        self.line_it += 1

    def rdf_2(self, power, par_a):
        """
        A theoretical plot for the RDF in the ideal gas limit approximation.
        The integral terms of the equation have not been included
        :param power: Potential strength
        :param par_a: Softening parameter
        :return: Nothing, Simply adds a figure to the corresponding canvas
        """
        self.n_str = str(power)
        a = str(float(par_a))
        name = "n: " + self.n_str + " A: " + a

        x = np.linspace(0, 5, 300)
        G = np.exp(-1. / ((x ** 2 + par_a) ** (power / 2.0)) / 1.4)  # T_0 = 1.4 = thermostat temperature
        plt.figure('RDF for Ideal gas')
        plt.plot(x, G, label=name)
        plt.xlim(xmin=0, xmax=2.5)
        plt.ylim(ymin=0)
        plt.legend(loc='best', fancybox=True)

        self.c += 1

    def vel_dist(self, rho, t, power, par_a):
        """
        Plots the velocity distributions for the X, Y, Z and the combined velocity vector
        for the last savec position of the fluid
        :param rho: Density
        :param t: Temperature
        :param power: Pair potential strength
        :param par_a: Softening parameter
        :return: Nothing. Simply adds a plot on the corresponding canvas
        """
        file_id = self.file_searcher(rho, t, power, par_a)
        data = "Positions_Velocities" + file_id + ".txt"
        name = "n: " + self.n_str + " A: " + self.a_str
        vx, vy, vz = np.loadtxt(data,
                                usecols=(3, 4, 5),
                                delimiter='\t',
                                comments='#',
                                unpack=True)
        v = np.sqrt(np.square(vx) + np.square(vy) + np.square(vz))

        # n, bins, patches = plt.hist(v, 100, normed=1, label=name)
        xmin, xmax = 0, max(v) + 1
        lnspc = np.linspace(xmin, xmax, len(v))
        # m, var, skew, kurt = stats.maxwell.stats(moments='mvsk')
        mean, std = stats.maxwell.fit(v, loc=0, scale=1)
        pdf_mb = stats.maxwell.pdf(lnspc, mean, std)

        fig = plt.figure('Velocity Dist Vx, Vy, Vz, V')

        vx_plot = plt.subplot2grid((2, 3), (0, 0), colspan=1)
        vy_plot = plt.subplot2grid((2, 3), (0, 1), colspan=1)
        vz_plot = plt.subplot2grid((2, 3), (0, 2), colspan=1)
        v_plot = plt.subplot2grid((2, 3), (1, 0), colspan=3)

        n, bins, patches = vx_plot.hist(vx, 150, density=1, label='vx')
        n, bins, patches = vy_plot.hist(vy, 150, density=1, label='vy')
        n, bins, patches = vz_plot.hist(vz, 150, density=1, label='vz')
        n, bins, patches = v_plot.hist(v, 150, density=1, label='v')

        plt.plot(lnspc, pdf_mb, label='Theory')
        plt.xlim(xmin=0)
        plt.title(self.n_str + '~' + self.a_str)
        plt.legend(loc='best', fancybox=True)

    def rdf_interpolate(self, rho, t, power, par_a):
        """
        It interpolates linearly between the data provided for the RDF which in turn
        makes possible to find the intersection point between the curves.
        :param rho: Density
        :param t: Temperature
        :param power: Pair potential strength
        :param par_a: Softening parameter
        :return: A numpy.array of the interpolated RDF data
        """
        file_id = self.file_searcher(rho, t, power, par_a)
        data = "RDF" + file_id + ".txt"
        # Measure the number of lines in the file in order to normalise
        # the distance r to the correct units
        num_lines = sum(1 for line in open(data))
        # Load the RDF data
        rdf = np.loadtxt(data, delimiter="\t", usecols=1, comments="#")

        # Number of particles, Number of bins
        particles, bins = int(self.particles_str), 500
        # Cut off radius
        rg = 3.0
        self.dr = rg / bins

        num_lines -= 2
        r = np.linspace(1, num_lines, num_lines)
        r = np.multiply(r, self.dr)
        max_scaling = np.max(rdf)  # Scaling the plot to ymax

        name = "rho: " + self.rho_str + " T: " + self.t_str + " n: " + self.n_str + " A: " + self.a_str

        # Make and interpolating function
        f = interpolate.interp1d(r, rdf, kind='linear')
        # Create a more accurate radius array
        r_interp = np.linspace(r[0], r[-1], 1000)

        # Use the interpolation function
        rdf_interp = f(r_interp)
        # begin ????
        self.dr = rg / 1000
        self.interpolated_data.append(rdf_interp)
        # end ????
        plt.figure('Interpolated RDF')
        plt.plot(r_interp, rdf_interp, '-o', label='interpolation '+name)
        plt.plot([0, r[-1]], [1, 1], '--', color='black', linewidth=0.5)
        # Plot limits and legends
        plt.xlim(xmin=0, xmax=3)
        plt.ylim(ymin=0, ymax=max_scaling + 0.1)
        # Plot labels
        plt.xlabel(r"$r$", fontsize=16)
        plt.ylabel(r"$g(r)$", fontsize=16)
        plt.legend(loc="best", fancybox=True, prop={'size': 8})
        return r_interp, rdf_interp

    def find_intersect(self, rho, t, n_list, par_a):
        rdf_interp_list = []
        rdf_list = []
        r_interp = []
        r = []
        # Get all RDFs with different n values into a single container
        for n in n_list:
            # Get normal rdfs
            r, rdf_temp = self.rdf(rho, t, n, par_a)
            # Get interpolated rdf
            r_interp, rdf_interp_temp = self.rdf_interpolate(rho, t, n, par_a)
            rdf_interp_list.append(rdf_interp_temp)
            rdf_list.append(rdf_temp)

        # RDF and RDF_interp need different for-loops since they use diff # bins
        # Normal rdf
        cwd = "/home/gn/Code/Python/MD-Simulation-Data-Analysis"
        f = open(f"{cwd}/intersect_rho_{rho}_T_{t}_A_{par_a}.log", "w")
        f.write("r\tstdev\n")
        for i in range(len(rdf_list[0])):
            # Define and clear list per iteration
            values = []
            for n in range(len(rdf_list)):
                # Pass the values of the same i to a list
                values.append(rdf_list[n][i])
            # Calculate stdev and store
            f.write(f"{r[i]}\t{stat.pstdev(values)}\n")
        f.close()

        f = open(f"{cwd}/logs/intersect_interpolated_rho_{rho}_T_{t}_A_{par_a}.log", "w")
        f.write("r\tstdev\n")
        for i in range(len(rdf_interp_list[0])):
            # Define and clear list per iteration
            values = []
            for n in range(len(rdf_interp_list)):
                # Pass the values of the same i to a list
                values.append(rdf_interp_list[n][i])
            # Calculate stdev and store
            f.write(f"{r_interp[i]}\t{stat.pstdev(values)}\n")
        f.close()




        # Loop through the container of the RDFs to find intersection points



    # def find_intersect(self, interp_list):
    #     """
    #     Finds the intersection of our interpollated data, assuming that our interpolated data
    #     have a reasonable sampling frequency so that the interpolation is not coarse.
    #     :param interp_list: a list of lists, or a 2D array
    #     :return:
    #     """
    #     # TODO: this needs serious fixing, code below just for demonstrative purposes in meeting
    #     l0 = interp_list[2]
    #     l1 = interp_list[3]
    #     l2 = interp_list[2]
    #     l3 = interp_list[3]
    #     list_r = []
    #     tolerance = 0.0005
    #     for i in range(len(l0)):
    #         if (l0[i] > 0.5 and l1[i] > 0.5) and abs(l0[i] - l1[i]) < tolerance:
    #             list_r.append(i*self.dr)
    #     print(list_r)
    #     return list_r

    @staticmethod
    def savefig(figname):
        return plt.savefig(figname, bbox_inches='tight', pad_inches=0)


