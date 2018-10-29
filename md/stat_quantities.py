import numpy as np
import scipy.stats as stats
from scipy import interpolate
from scipy.signal import argrelextrema, butter, filtfilt
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

        # RDF shared variables
        self.r = []         # container for the rdf-x data
        self.rdf_data = []  # container for the rdf y-data
        self.num_lines_rdf = 0  # number of lines filled with data in the rdf file
        self.rg = 3.0   # cut-off radius
        self.dr = 0     # distance increment in the radius
        self.iso = 0    # x-location for the theoretical isosbestic points
        # Containers for interpolated data
        self.r_interp = []
        self.rdf_interp = []
        self.rdf_interp_smooth = []

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
        self.num_lines_rdf = sum(1 for line in open(data))
        self.rdf_data = np.loadtxt(
            data, delimiter="\t", usecols=1, comments="#")

        # Number of particles, Number of bins
        particles, bins = int(self.p_str), 500
        # Cut off radius
        self.rg = 3.0
        self.dr = self.rg / bins

        # r range, r=0 is intentionally neglected due to division by 0
        # num_lines-2 because of comments and headers
        self.num_lines_rdf -= 2
        self.r = np.linspace(1, self.num_lines_rdf, self.num_lines_rdf)
        self.r = np.multiply(self.r, self.dr)

        # Isomorphic scaling of r for the isomorph plane
        if iso_scale is True:
            self.r = np.multiply(self.r, rho ** (1. / 3.))  # Scale r

        # Plotting isosbestic point
        if show_iso is True:  # Show isosbestic point
            # TODO: this is not correct, missing a factor probably, revise theory!
            self.iso = np.sqrt(1 - par_a ** 2)

        # return the plotting lists
        return self.r, self.rdf_data

    def rdf_plot(self, rho, t, power, par_a, iso_scale=False, show_iso=False):
        """
        Generates a plot of the Radial Distribution Data after it calls the rdf method


        :param rho: Density
        :param t: Temperature
        :param power: Pair potential strength
        :param par_a: Softening parameter
        :param iso_scale: Optional, scales the values of r, the radius, based on the isosbestic model that
                          has been created
        :param show_iso: Optional, shows the location of the theoretical isosbestic point between different rdfs
        """
        self.rdf(rho, t, power, par_a, iso_scale, show_iso)
        plt.figure('Radial Distribution Function')

        max_scaling = np.max(self.rdf_data)  # Scaling the ymax

        # Plotting isosbestic location of points
        if show_iso is True:
            plt.plot([self.iso, self.iso], [
                     0, max_scaling + 0.1], '--', color='red')

        # Naming the curves
        name = f"rho: {self.rho_str} T: {self.t_str} n: {self.n_str} A: {self.a_str}"

        plt.plot(self.r, self.rdf_data, '-', markersize=4, label=name)

        # Plot labels
        plt.xlabel(r"$r$", fontsize=16)
        plt.ylabel(r"$g(r)$", fontsize=16)

        # Line through y = 1
        plt.plot([0, self.r[-1]], [1, 1], '--', color='black', linewidth=0.5)

        # Plot limits and legends
        plt.xlim(left=0, right=self.rg)
        plt.ylim(bottom=0, top=max_scaling + 0.1)
        plt.legend(loc="best", fancybox=True, prop={'size': 8})
        print("@ index: ", np.argmax(self.rdf_data),
              " value: ", max(self.rdf_data))

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

    def rdf_interpolate(self, rho, t, power, par_a, range_refinement=2000, ignore_zeroes=True):
        """
        It smooths the data using a forward-backward low-pass filter.
        Then it interpolates linearly between the data provided for the RDF which in turn
        makes possible to find the intersection point between the curves.


        :param rho: Density
        :param t: Temperature
        :param power: Pair potential strength
        :param par_a: Softening parameter
        :param range_refinement: The accuracy of the interpolated data
        :param ignore_zeroes: Ignores zeroes and close to zero values in the specified RDF array
        :return: 3 numpy.arrays of the interpolated and smoothed RDF data
        """
        self.rdf(rho, t, power, par_a)
        # Smooth the data before interpolating with a forward backward filter
        # First create a lowpass butterworth filter
        # TODO: fix the parameters for the filter
        b, a = butter(3, 0.09)

        # Make the interpolating functions
        f = interpolate.interp1d(self.r, self.rdf_data, kind='linear')

        # The loop creates a filter that ignores zeroes
        # which in turn allows a more accurate smoothing of the data
        rdf_smooth = None
        if ignore_zeroes is True:
            # Get the non-zero values and make a filter
            non_zero = self.rdf_data[np.nonzero(self.rdf_data)[0]]
            rdf_smooth = filtfilt(b, a, non_zero)
            # Locate the zero values in the rdf data
            zero_idx = np.where(self.rdf_data == 0)[0]
            # Knowing that the RDF values == 0 are always
            # at the start of the array allows us to append
            # the two arrays to each other
            rdf_smooth = np.concatenate((self.rdf_data[zero_idx], rdf_smooth))

        else:
            rdf_smooth = filtfilt(b, a, self.rdf_data)

        f_smooth = interpolate.interp1d(self.r, rdf_smooth)

        # Create a radius array with increased precision (number of bins)
        self.r_interp = np.linspace(self.r[0], self.r[-1], range_refinement)

        # Use the interpolation functions
        self.rdf_interp = f(self.r_interp)
        self.rdf_interp_smooth = f_smooth(self.r_interp)

        # Generating new separation unit distance for the x-axis
        self.dr = self.rg / range_refinement

        # Passing interpolated data to be stored later in file
        self.interpolated_data.append(self.rdf_interp)
        return self.r_interp, self.rdf_interp, self.rdf_interp_smooth

    def rdf_interpolate_plot(self, rho, t, power, par_a,
                             range_refinement=2000):
        """
        Generates a plot of the interpolated data after it calls rdf_interpolate


        :param rho: Density
        :param t: Temperature
        :param power: Pair potential strength
        :param par_a: Softening parameter
        :param range_refinement: The accuracy of the interpolated data
        """
        self.rdf_interpolate(rho, t, power, par_a, range_refinement)

        # Naming the curves
        name = f"rho: {self.rho_str} T: {self.t_str} n: {self.n_str} A: {self.a_str}"
        max_scaling = np.max(self.rdf_data)  # Scaling the plot to ymax

        plt.figure('Interpolated RDF')

        plt.plot(self.r_interp, self.rdf_interp,
                 '-', label='interpolation ' + name)
        plt.plot([0, self.r_interp[-1]], [1, 1],
                 '--', color='black', linewidth=0.5)

        # Plot limits and legends
        plt.xlim(left=0, right=3)
        plt.ylim(bottom=0, top=max_scaling + 0.1)

        # Plot labels
        plt.xlabel(r"$r$", fontsize=16)
        plt.ylabel(r"$g(r)$", fontsize=16)
        plt.legend(loc="best", fancybox=True, prop={'size': 8})

    def rdf_interpolate_smooth_plot(self, rho, t, power, par_a,
                                    range_refinement=2000):
        """
        Generates a plot of the smoothed and interpolated
        data after it calls rdf_interpolate


        :param rho: Density
        :param t: Temperature
        :param power: Pair potential strength
        :param par_a: Softening parameter
        :param range_refinement: The accuracy of the interpolated data
        """
        self.rdf_interpolate(rho, t, power, par_a, range_refinement)

        # Naming the curves
        name = f"rho: {self.rho_str} T: {self.t_str} n: {self.n_str} A: {self.a_str}"
        max_scaling = np.max(self.rdf_data)  # Scaling the plot to ymax

        plt.figure('Interpolated RDF')

        plt.plot(self.r_interp, self.rdf_interp_smooth,
                 '-.', label='smooth interp ' + name)
        plt.plot([0, self.r_interp[-1]], [1, 1],
                 '--', color='black', linewidth=0.5)

        # Plot limits and legends
        plt.xlim(left=0, right=3)
        plt.ylim(bottom=0, top=max_scaling + 0.1)

        # Plot labels
        plt.xlabel(r"$r$", fontsize=16)
        plt.ylabel(r"$g(r)$", fontsize=16)
        plt.legend(loc="best", fancybox=True, prop={'size': 8})

    def rdf_intersect(self, rho, t, power_list, par_a,
                      range_refinement=2000, r_lower=0, r_higher=-1, intersections=1):
        """
        Finds the points of intersection in RDF functions by looping through multiple ns.
        The RDF data are first smoothed with a forward-backward filter in order to reduce
        the noise and then the data are then interpolated linearly.
        An intersection between curves is quantified as the minimum in the standard
        deviation at the same index gr[i] but for multiple curves. At an intersection,
        the std is theoretically at a minimum.

        The points of inflection are identified on the last curve passed on the method,
        which in turn splits the RDF data into regions. In every region between a max-min
        an intersection point is sought (although for some extreme cases that is not entirely true).

        Worth remembering is that smoothing the data results into the introduction of unwanted
        inflection points, neat neighbourhoods where the data do not fluctuate much.

        # TODO: make the interpolation/smoothing act to all the non-zero values
                pass arg ignore_zeroes=True and handle with if-statement


        :param rho: Density
        :param t: Temperature
        :param power_list: List of pair potential strengths
        :param par_a: Softening parameter
        :param range_refinement: The accuracy of the interpolated data
        :param r_lower: Specifies the lower bound of the array that
                        the method will look for intersections.
                        It should be provided in terms of indices
        :param r_higher: Specifies the upper bound of the array that
                         the method will look for intersections.
                         It should be provided in terms of indices
        :param intersections: Number of intersections to look for between a max and a min
        :return: A plot with the interpolated and smoothed RDF
                 along with the local max, min, and isosbestic (intersection) points
        """
        # Passing the smoothed RDF data into a list of lists
        rdf_interp_list = []
        for n in power_list:
            # calls rdf internally, and initialises r_interp list
            self.rdf_interpolate(rho, t, n, par_a, range_refinement)
            # Selecting only the required range of smooth interpolated rdf values
            sliced_rdf_data = self.rdf_interp_smooth[r_lower:r_higher]
            rdf_interp_list.append(sliced_rdf_data)

        # Transpose the RDF array for easier file output
        rdf_interp_list = np.transpose(rdf_interp_list)

        # Calculate the average and standard deviation, across n runs at the same index
        mean_list = np.mean(rdf_interp_list, axis=1)
        std_list = np.std(rdf_interp_list, axis=1)

        # Get the coordinates for the local maxima and intersections in the provided range
        # The max, min plotted here correspond to the last n in the power_list
        x_max, y_max, x_min, y_min, idx_max, idx_min = self.find_local_min_max(
            self.r_interp[r_lower:r_higher], self.rdf_interp_smooth[r_lower:r_higher])

        # Plotting the intersection results into the interpolated RDF canvas
        plt.figure('Interpolated RDF')

        # Merge and sort the arrays containing the indices of the inflection points
        idx_max_min = np.concatenate((idx_max[0], (idx_min[0])))
        idx_max_min.sort()

        # Loop through all the max-min combinations and spot the
        # intersect of the curves based on the std
        for i in range(1, len(idx_max_min)):
            # Get the index of the k smallest stds, which in theory
            # should correspond to the point where the curves intersect
            idx_intersect = np.argpartition(
                std_list[idx_max_min[i-1]:idx_max_min[i]], intersections)
            # Adjust index to match with global array index
            idx_intersect += idx_max_min[i-1]

            # Get the mean for the intersection points
            mean_scatter = mean_list[idx_intersect[:intersections]]

            # Add the lower boundary index to the std,
            # to correspond to the actual r index
            idx_intersect += r_lower

            # Get the r-values for the corresponding RDF averaged values
            r_data = [self.r_interp[index]
                      for index in idx_intersect[:intersections]]

            # scatter_label = [f"x:{r_data[i]}, y:{mean_scatter[i]}" r]
            plt.scatter(r_data, mean_scatter, marker='x', color='red')

        plt.scatter(x_max, y_max, color='orange')
        plt.scatter(x_min, y_min, color='green')

        # Plotting the data curves of the interpolated data
        for n in power_list:
            self.rdf_interpolate_smooth_plot(
                rho, t, n, par_a, range_refinement)
            # self.rdf_interpolate_plot(
            #     rho, t, n, par_a, range_refinement)

    @staticmethod
    def find_local_min_max(x, y):
        """
        Finds the local minima and maxima and their indices. Assuming that x and y
        have the same dimensions


        :param x: The x variable array
        :param y: The y variable array
        :return: x_max, y_max, x_min, y_min, idx_max, idx_min.
                 The x-y coordinates for the maxima and minima.
                 Along with their indices so that the max, min values
                 can be further used
        """
        # Find the index of the local maxima and minima
        # this index can then be used to find the x-y values of the points
        idx_local_max = argrelextrema(y, np.greater)
        idx_local_min = argrelextrema(y, np.less)

        # Get the x-value for the local maxima and minima
        x_local_max = x[idx_local_max[0]]
        x_local_min = x[idx_local_min[0]]
        y_local_max = y[idx_local_max[0]]
        y_local_min = y[idx_local_min[0]]

        return x_local_max, y_local_max, x_local_min, y_local_min, idx_local_max, idx_local_min


if __name__ == "__main__":
    import os

    os.chdir("/home/gn/Desktop/test_data")

    obj = StatQ(15000, 1000)
    n = [6, 7, 8, 9, 10]
    a = [0, 0.25, 0.50, 0.75, 0.8, 0.90, 1.00, 1.1,
         1.25, 1.50, 1.75, 2.00, 2.25, 2.50, 4.00]

    obj.rdf_intersect(0.5, 0.5, n, 0, r_lower=550,
                      r_higher=-1, intersections=1)

    # for i in n:
    # obj.rdf_plot(0.5, 0.5, i, 0.25)
    # Shows the interpolated data
    # obj.rdf_interpolate_plot(0.5, 0.5, i, 0.25)
    # obj.msd(0.5, 0.5, i, 0.25)
    # obj.vaf(0.5, 0.5, i, 0.25)

    plt.show()
