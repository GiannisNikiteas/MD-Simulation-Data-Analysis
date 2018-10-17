import numpy as np
import matplotlib.cm as cm
import scipy.stats as stats
from scipy import interpolate
import matplotlib.pyplot as plt
import statistics as stat
from mpl_toolkits.mplot3d import Axes3D


class FileNaming(object):
    def __init__(self, steps, particles):
        self.steps_str = str(steps)
        self.p_str = str(particles)
        self.rho_str = None
        self.t_str = None
        self.n_str = None
        self.a_str = None

    def file_searcher(self, rho, t, n, alpha=None):
        """
        Constructs the file signature of the MD simulation in order for the information to be read
        :param rho: density
        :param t:   temperature
        :param n:   potential strength
        :param alpha:   softness parameter
        :return: string with file identifier
        """
        self.rho_str = "{:.4f}".format(rho)
        self.t_str = "{:.4f}".format(t)
        self.n_str = str(n)
        a = None
        if alpha is not None:
            self.a_str = "{:.5f}".format(alpha)
            a = '_A_' + self.a_str
        else:
            self.a_str = ""
            a = self.a_str

        name_id = f"_step_{self.steps_str}_particles_{self.p_str}_rho_{self.rho_str}_T_{self.t_str}_n_{self.n_str}{a}"
        return name_id


class StatQ(FileNaming):

    def __init__(self, steps, particles):
        super().__init__(steps, particles)
        self.dif_coef = np.array([])
        self.dif_err = np.array([])
        self.dif_y_int = np.array([])
        self.line_style = ['solid', 'dashed', 'dotted', 'dashdot']
        self.interpolated_data = []
        self.dr = None
        # This is an iterator for the color array
        self.p, self.c = 0, 0
        self.j = 0  # stride for MSD
        self.v = 0  # index for MSD
        self.line_it = 0  # Index iterator for line styles

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
        data = f"RDF{file_id}.txt"
        num_lines = sum(1 for line in open(data))
        rdf = np.loadtxt(data, delimiter="\t", usecols=1, comments="#")

        # Number of particles, Number of bins
        particles, bins = int(self.p_str), 500
        # Cut off radius
        rg = 3.0
        dr = rg / bins

        # r range, r=0 is intentionally neglected due to division by 0
        # num_lines-1 because one of them is a comment
        num_lines -= 2
        r = np.linspace(1, num_lines, num_lines)
        r = np.multiply(r, dr)
        a_tilde = par_a * rho ** (1. / 3.)  # Scale a

        # Isomorphic scaling of r for the isomorph plane
        if iso_scale is True:
            r = np.multiply(r, rho ** (1. / 3.))  # Scale r

        plt.figure('Radial Distribution Function')

        # Plotting isosbestic point
        max_scaling = np.max(rdf)  # Scaling the ymax
        if show_iso is True:  # Show isosbestic point
            # TODO: this is not correct, missing a factor probably, revise theory!
            iso = np.sqrt(1 - par_a ** 2)
            # iso = np.sqrt(rho ** (2./3) - a_tilde ** 2)   # This is probably wrong
            plt.plot([iso, iso], [0, max_scaling + 0.1], '--', color='red')

        name = f"rho: {self.rho_str} T: {self.t_str} n: {self.n_str} A: {self.a_str}"
        plt.plot(r, rdf, '-', markersize=4, label=name)

        # Plot labels
        plt.xlabel(r"$r$", fontsize=16)
        plt.ylabel(r"$g(r)$", fontsize=16)
        # Plotting isosbestic location of points

        # Line through y = 1
        plt.plot([0, r[-1]], [1, 1], '--', color='black', linewidth=0.5)
        # Plot limits and legends
        plt.xlim(left=0, right=3)
        plt.ylim(bottom=0, top=max_scaling + 0.1)
        plt.legend(loc="best", fancybox=True, prop={'size': 8})
        print("@ index: ", np.argmax(rdf), " value: ", max(rdf))

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

        cr = np.loadtxt(data, usecols=8, delimiter='\t',
                        unpack=True, comments='#')

        # Measure the number of lines with data in file
        num_lines = 0
        for line in open(data):
            if line.startswith('#'):
                continue
            num_lines += 1

        time_step = 0.005 / np.sqrt(t)
        time = time_step * num_lines
        x = np.linspace(0, time, num_lines)
        t_tilde = x * rho ** (1. / 3.) * t ** 0.5
        name = ""
        if num_lines < 100000:
            name = f"rho: {self.rho_str} T: {self.t_str} n: {self.n_str} A: {self.a_str}"

        plt.figure('Velocity Autocorrelation Function')
        y = np.full(num_lines, 0)
        # xx = np.full(num_lines, time)
        # yy = np.linspace(5, -0.5, num_lines)
        # plt.plot(xx, yy, '--', color='black')
        plt.plot(t_tilde, y, '--', color='black')
        plt.plot(t_tilde, cr, label=name)
        plt.xlabel(r"Time $t$", fontsize=16)
        plt.ylabel(r"$C_v$", fontsize=16)
        # plt.ylim(top=5, bottom=-0.5)
        plt.xlim(left=t_tilde[0], right=t_tilde[-1])
        plt.legend(loc="best", ncol=1)

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
        data = f"Data{file_id}.txt"

        msd_data = np.loadtxt(
            data, usecols=7, delimiter='\t', unpack=True)

        # Measure lines filled with data
        num_lines = 0
        for line in open(data):
            if line.startswith('#'):
                continue
            num_lines += 1

        step = 0.005 / np.sqrt(t)
        x = np.linspace(0, num_lines - 1, num=num_lines)

        # Regular coefs are calculated independent of the if loop
        gradient, intercept, r_value, p_value, std_err = stats.linregress(
            x, msd_data)
        self.dif_coef = np.append(self.dif_coef, gradient)
        self.dif_err = np.append(self.dif_err, std_err)
        self.dif_y_int = np.append(self.dif_y_int, intercept)

        # TODO: this should be output to a log file for reference
        # print('Diffusion coef: ', gradient, '\n',
        #       'y-intercept: ', intercept, '\n',
        #       'R value: ', r_value, '\n',
        #       'Fit Error: ', std_err)

        name = f"rho: {self.rho_str} T: {self.t_str} n: {self.n_str} A: {self.a_str}"

        plt.figure('Mean Square Displacement')
        plt.plot(x, msd_data, label=name)
        plt.xlim(left=0, right=x[-1])
        plt.xlabel(r"$t$", fontsize=16)
        plt.ylabel(r"$MSD$", fontsize=16)
        plt.legend(loc="best", fancybox=True)

        return msd_data

    def diffusion_plot(self, rho, t, power, my_list):
        """
        A graph of the Diffusion coefficients D against a list of parameter A values
        :param rho: Density
        :param t: Temperature
        :param power: Potential strength parameter n
        :param my_list: List of parameter A coefficients
        :return: Figure of D vs A for a given number of iterations
        """
        # I recommend this list of values for a
        # my_list = [0, 0.25, 0.50, 0.65, 0.7, 0.75, 0.8, 0.85, 0.90, 0.97, 1.00, 1.1, 1.2,
        #    1.25, 1.50, 1.75, 2.00, 2.25, 2.50, 2.75, 4.00]
        for i in my_list:
            self.msd(rho, t, power, i)
            print("-----------------------------")
        name = f"n: {str(power)}"

        plt.figure('Diffusion coefficients D vs A')
        plt.plot(my_list, self.dif_coef, '--o', label=name, markersize=3.5)
        plt.xlabel(r"$a$", fontsize=16)
        plt.ylabel(r"$D$", fontsize=16)
        plt.legend(loc="best", fancybox=True, ncol=2)

        # Resetting the Best Fit model
        self.dif_coef, self.dif_err, self.dif_y_int = np.array(
            []), np.array([]), np.array([])
        plt.ylim(bottom=0)
        plt.xlim(left=0)
        self.j += 15
        self.v += 1


class StateProperties(FileNaming):
    def __init__(self, step, particles):
        super().__init__(step, particles)
        self.line_style = ['solid', 'dashed', 'dotted',
                           'dashdot']  # TODO: Fix with itertools
        self.p, self.c = 0, 0
        self.line_it = 0  # Index iterator for line styles

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
        data = f"Data{file_id}.txt"

        num_lines = 0
        # Measuring the line number in a file ignoring comments
        for line in open(data):
            if line.startswith('#'):
                continue
            num_lines += 1

        pot_en, kin_en = np.loadtxt(data, usecols=(
            3, 4), delimiter='\t', comments='#', unpack=True)
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

        rho_list, u = np.loadtxt(data, usecols=(
            1, 3), delimiter='\t', comments='#', unpack=True)

        #  Plots the Energies
        name = "rho: " + self.rho_str + "T: " + \
               self.t_str + "n: " + self.n_str + "A: " + self.a_str
        step = 0.005
        time = num_lines * step
        x = np.linspace(0, time, num_lines)
        plt.figure('Potential Plots of Data')
        plt.plot(rho_list, u, label=name)
        plt.legend(loc='best', fancybox=True)


if __name__ == "__main__":
    import os

    os.chdir("/home/gn/Desktop/test_data")

    obj = StatQ(15000, 1000)
    n = [6, 8, 10, 12]
    a = [0, 0.25, 0.50, 0.75, 0.8, 0.90, 1.00, 1.1,
         1.25, 1.50, 1.75, 2.00, 2.25, 2.50, 4.00]
    for i in n:
        obj.rdf(0.5, 0.5, i, 0)

        # obj.msd(0.5, 0.5, 12, i)
        # obj.vaf(0.5, 0.5, 12, i)

    plt.show()
