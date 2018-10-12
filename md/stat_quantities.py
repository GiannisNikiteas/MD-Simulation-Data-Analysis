import numpy as np
import matplotlib.cm as cm
import scipy.stats as stats
from scipy import interpolate
import matplotlib.pyplot as plt
import statistics as stat
from mpl_toolkits.mplot3d import Axes3D
#
# class FileNaming:
#     def __init__(self, step, particles):
#         self.step_str = str(step)
#         self.p_str = str(particles)
#         self.rho_str = None
#         self.t_str = None
#         self.n_str = None
#         self.alpha = None
#         self.a_str = None
#
#     def file_searcher(self, rho, t, n, alpha=None):
#         self.rho_str = "{:.4f}".format(rho)
#         self.t_str = "{:.4f}".format(t)
#         self.n_str = str(n)
#         self.alpha = alpha
#         a = None
#         if self.alpha is not None:
#             self.a_str = "{:.5f}".format(self.alpha)
#             a = '_A_' + self.a_str
#         else:
#             self.a_str = ""
#             a = self.a_str
#
#         name_id = f"_step_{self.step_str}_particles_{self.p_str}_rho_{self.rho_str}_T_{self.t_str}_n_{self.n_str}{a}"
#         return name_id


class StatQ:

    def __init__(self, __step, __particles):
        self.step_str = str(__step)
        self.p_str = str(__particles)
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
        self.line_style = ['solid', 'dashed', 'dotted', 'dashdot']
        self.interpolated_data = []
        self.dr = None
        # This is an iterator for the color array
        self.p, self.c = 0, 0
        self.j = 0  # stride for MSD
        self.v = 0  # index for MSD
        self.line_it = 0    # Index iterator for line styles

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

        name_id = f"_step_{self.step_str}_particles_{self.p_str}_rho_{self.rho_str}_T_{self.t_str}_n_{self.n_str}{a}"
        return name_id

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
        plt.xlim(left=0, right=3)
        plt.ylim(bottom=0, top=max_scaling + 0.1)
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
        plt.ylim(top=5, bottom=-0.5)
        plt.xlim(left=t_tilde[0], right=t_tilde[-1])
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
        data = f"Data{file_id}.txt"

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


if __name__ == "__main__":
    import os
    os.chdir("/home/gn/Code/MD-simulation/examples/example_data")
    obj = StatQ(10000, 1000)
    obj.rdf(0.5, 0.5, 6, 0.5)
    plt.show()
