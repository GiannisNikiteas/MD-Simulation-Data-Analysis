import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt


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
        Constructs the file signature of the MD simulation
        in order for the information to be read

        :param rho: density
        :param t:   temperature
        :param n:   potential strength
        :param alpha:   softness parameter
        :return: string with file identifier
        """
        self.rho_str = "{:.4f}".format(rho)
        self.t_str = "{:.4f}".format(t)
        self.n_str = "{:.2f}".format(n)
        a = None
        if alpha is not None:
            self.a_str = "{:.5f}".format(alpha)
            a = '_A_' + self.a_str
        else:
            self.a_str = ""
            a = self.a_str

        name_id = f"_step_{self.steps_str}_particles_{self.p_str}" \
                  f"_rho_{self.rho_str}_T_{self.t_str}_n_{self.n_str}{a}"
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
        self.step = 0.005

        # This is an iterator for the color array
        self.p, self.c = 0, 0
        self.j = 0  # stride for MSD
        self.v = 0  # index for MSD
        self.line_it = 0  # Index iterator for line styles

        # RDF shared variables
        self.r = []         # container for the rdf-x data
        self.rdf_data = []  # container for the rdf y-data
        self.rdf_bins = 0  # number of lines filled with data in the rdf file
        self.rg = 3.0   # cut-off radius
        self.dr = 0     # distance increment in the radius
        self.iso = 0    # x-location for the theoretical isosbestic points

    # Radial Distribution Function
    def rdf(self, rho, t, power, par_a, iso_scale=False):
        """
        Reads the data corresponding to the Radial Distribution Function from
        file.

        @param rho: Density
        @param t: Temperature
        @param power: Pair potential strength
        @param par_a: Softening parameter
        @param iso_scale: Optional, scales the values of r, the radius, 
                          based on the isosbestic model that has been created
        @return: The numpy.array for the RDF data

        """
        file_id = self.file_searcher(rho, t, power, par_a)
        data = f"RDF{file_id}.txt"
        # self.rdf_bins = sum(1 for line in open(data))
        self.rdf_data = np.loadtxt(data, delimiter="\t",
                                   usecols=1, comments="#")
        # Calculate the number of bins present in RDF
        self.rdf_bins = int(len(self.rdf_data))

        # Bin width in r-units
        self.dr = self.rg / self.rdf_bins

        # r=0 is intentionally neglected due to division by 0
        self.r = np.linspace(1, self.rdf_bins, self.rdf_bins)
        self.r = np.multiply(self.r, self.dr)

        # Isomorphic scaling of r for the isomorph plane
        if iso_scale is True:
            self.r = np.multiply(self.r, rho ** (1. / 3.))  # Scale r

        # return the plotting lists
        return self.r, self.rdf_data

    def rdf_plot(self, rho, t, power, par_a,
                 iso_scale=False,
                 show_label=True,
                 **kwargs):
        """
        Creates a plot for the Radial Distribution Function of the fluid, 
        which depicts the microscopic density fluctuations of the molecules 
        as a function of distance.
        The parameters iso_scale is optional and should normally 
        be set to False.

        @param rho: Density
        @param t: Temperature
        @param power: Pair potential strength
        @param par_a: Softening parameter
        @param iso_scale: Optional, scales the values of r, the radius, 
                          based on the isosbestic model that has been created
        @param show_label:
        """
        self.rdf(rho, t, power, par_a, iso_scale)
        plt.figure('Interpolated RDF')

        max_scaling = np.max(self.rdf_data)  # Scaling the ymax

        # Naming the curves
        name = ""
        if show_label is True:
            name = f"\N{GREEK SMALL LETTER RHO}: {self.rho_str}" \
                   f" T: {self.t_str} n: {self.n_str} A: {self.a_str}"

        plt.plot(self.r, self.rdf_data, '-',
                 markersize=4, label=name, **kwargs)

        # Plot labels
        plt.xlabel(r"$r$")
        plt.ylabel(r"$g(r)$")

        # Line through y = 1
        plt.plot([0, self.r[-1]], [1, 1], '--', color='black', linewidth=0.5)

        # Plot limits and legends
        plt.xlim(left=0, right=self.rg)
        plt.ylim(bottom=0, top=max_scaling + 0.1)
        if show_label is True:
            plt.legend(loc="best", fancybox=True, prop={'size': 8})

    # Velocity Autocorrelation Function
    def vaf(self, rho, t, power, par_a, iso_scale=False, **kwargs):
        """
        Creates a figure for the Velocity Autocorrelation Function of the fluid,
        which illustrates if the fluid remains a coupled system
        through time (correlated) or it uncouples.

        @param rho: Density
        @param t: Temperature
        @param power: Pair potential strength
        @param par_a: Softening parameter
        @param iso_scale:
        @type iso_scale:
        @return: Nothing. Simply adds a plot on the corresponding canvas
        """
        file_id = self.file_searcher(rho, t, power, par_a)
        data = "Data" + file_id + ".txt"

        cr = np.loadtxt(data, usecols=8, delimiter='\t',
                        unpack=True, comments='#')

        num_lines = int(len(cr))
        time_step = self.step / np.sqrt(t)
        time_max = time_step * num_lines
        time = np.linspace(0, time_max, num_lines)

        name = f"\N{GREEK SMALL LETTER RHO}: {self.rho_str} T: {self.t_str}"\
               f"n: {self.n_str} A: {self.a_str}"

        plt.figure('Velocity Autocorrelation Function')
        y = np.full(num_lines, 0)

        if iso_scale is True:
            time = time * (rho ** (1.0 / 3.0)) * (t ** 0.5)

        plt.plot(time, y, '--', color='black')
        plt.plot(time, cr, label=name, **kwargs)
        plt.xlabel(r"Time $t$", fontsize=16)
        plt.ylabel(r"$C_r$", fontsize=16)

        plt.xlim(left=time[0], right=time[-1])
        plt.legend(loc="best", ncol=1)

    # Mean Square Displacement
    def msd(self, rho, t, power, par_a, **kwargs):
        """
        Plots the Mean Square Displacement for our fluid.
        According to diffusion theory the slope of the MSD corresponds 
        to the inverse of the diffusion coefficient.

        @param rho: Density
        @param t: Temperature
        @param power: Pair potential strength
        @param par_a: Softening parameter
        @return: msd list
        """
        file_id = self.file_searcher(rho, t, power, par_a)
        data = f"Data{file_id}.txt"

        msd_data = np.loadtxt(data, usecols=7,
                              delimiter='\t', unpack=True)

        num_lines = int(len(msd_data))
        step = self.step / np.sqrt(t)
        x = np.linspace(0, num_lines - 1, num=num_lines)

        # Perform a linear fit to the MSD data and get fit parameters
        grad, intercept, rms, p_val, std = stats.linregress(x, msd_data)
        self.dif_coef = np.append(self.dif_coef, grad)
        self.dif_err = np.append(self.dif_err, std)
        self.dif_y_int = np.append(self.dif_y_int, intercept)

        name = file_id.replace("_", " ")
        plt.figure('Mean Square Displacement')
        plt.plot(x, msd_data, label=name, **kwargs)
        plt.xlim(left=0, right=x[-1])
        plt.xlabel(r"$t$")
        plt.ylabel(r"$MSD$")
        plt.legend(loc="best", fancybox=True)

        return msd_data

    def diffusion_plot(self, rho, t, power, my_list):
        """
        A graph of the Diffusion coefficients D against a list of parameter A

        @param rho: Density
        @param t: Temperature
        @param power: Potential strength parameter n
        @param my_list: List of parameter A coefficients
        @return: Figure of D vs A for a given number of iterations
        """

        for i in my_list:
            self.msd(rho, t, power, i)
            print("-----------------------------")
        name = f"n: {str(power)}"

        plt.figure('Diffusion coefficients D vs A')
        plt.plot(my_list, self.dif_coef, '--o', label=name, markersize=3.5)
        plt.xlabel(r"$a$")
        plt.ylabel(r"$D$")
        plt.legend(loc="best", fancybox=True, ncol=2)

        # Resetting the Best Fit model
        self.dif_coef, self.dif_err, self.dif_y_int = np.array(
            []), np.array([]), np.array([])
        plt.ylim(bottom=0)
        plt.xlim(left=0)
        self.j += 15
        self.v += 1

    def vel_dist(self, rho, t, power, par_a):
        """
        Plots the velocity distributions for the X, Y, Z and 
        the combined velocity vector for the last saved position of the fluid.

        @param rho: Density
        @param t: Temperature
        @param power: Pair potential strength
        @param par_a: Softening parameter
        @return: Nothing. Simply adds a plot on the corresponding canvas
        """
        file_id = self.file_searcher(rho, t, power, par_a)
        data = "Positions_Velocities" + file_id + ".txt"
        name = f"n: {self.n_str} A: {self.a_str}"
        vx, vy, vz = np.loadtxt(data,
                                usecols=(3, 4, 5),
                                delimiter='\t',
                                comments='#',
                                unpack=True)
        v = np.sqrt(np.square(vx) + np.square(vy) + np.square(vz))

        xmin, xmax = 0, max(v) + 1
        x = np.linspace(xmin, xmax, len(v))
        # m, var, skew, kurt = stats.maxwell.stats(moments='mvsk')
        mean, std = stats.maxwell.fit(v, loc=0, scale=1)
        pdf_mb = stats.maxwell.pdf(x, mean, std)

        fig = plt.figure('Velocity Dist Vx, Vy, Vz, V')

        vx_plot = plt.subplot2grid((2, 3), (0, 0), colspan=1)
        vy_plot = plt.subplot2grid((2, 3), (0, 1), colspan=1)
        vz_plot = plt.subplot2grid((2, 3), (0, 2), colspan=1)
        v_plot = plt.subplot2grid((2, 3), (1, 0), colspan=3)

        vx_plot.hist(vx, 150, density=1)
        vx_plot.set_title(r'$v_x$')
        vy_plot.hist(vy, 150, density=1)
        vy_plot.set_title(r'$v_y$')
        vz_plot.hist(vz, 150, density=1)
        vz_plot.set_title(r'$v_z$')

        v_plot.hist(v, 150, density=1, label='v')
        v_plot.plot(x, pdf_mb, label='Theory')

        plt.xlim(left=0)
        plt.title(self.n_str + ' ' + self.a_str)
        plt.legend(loc='best', fancybox=True)

    def sf(self, rho, t, power, par_a):
        file_id = self.file_searcher(rho, t, power, par_a)
        data = "Data" + file_id + ".txt"
        name = f"n: {self.n_str} A: {self.a_str}"
        
        sf = np.loadtxt(data, usecols=(9, 10, 11),
                        delimiter='\t', comments='#', unpack=True)

        x = np.arange(1, len(sf[0])+1)

        fig, ax = plt.subplots(3, 1, figsize=(12, 6))

        for (i, sf_i) in enumerate(sf):
            ax[i].scatter(x, sf[i], label=i)
            ax[i].set(xlabel='Steps', ylabel='SF')
            ax[i].legend(loc='best')
        ax[0].set_title(f'Structure Factor {file_id}')

# %%
if __name__ == "__main__":

    import os
    from mdtools.stat_quantities import *
    import matplotlib.pyplot as plt

    plt.style.use('default')

