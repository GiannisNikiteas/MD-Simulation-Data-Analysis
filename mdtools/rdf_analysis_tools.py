import numpy as np
from md import StatQ
from scipy import interpolate
import matplotlib.pyplot as plt
from scipy.signal import argrelextrema, butter, filtfilt
import itertools

MARKERS = itertools.cycle(("o", "v", "s", "p", "P", "*", "+", "x", "d"))


class RDFAnalysis(StatQ):

    def __init__(self, steps, particles):
        super().__init__(steps, particles)
        # Containers for interpolated data
        self.r_interp = []
        self.rdf_interp = []
        self.rdf_interp_smooth = []

    def rdf_interpolate(self, rho, t, power, par_a,
                        range_refinement=2000,
                        iso_scale=False,
                        ignore_zeroes=True):
        """
        It smooths the data using a forward-backward low-pass filter.
        Then it interpolates linearly between the data provided for the RDF which in turn
        makes possible to find the intersection point between the curves.


        :param rho: Density
        :param t: Temperature
        :param power: Pair potential strength
        :param par_a: Softening parameter
        :param range_refinement: The accuracy of the interpolated data
        :param iso_scale: Choosing whether the RDF data will be scaled onto an isomorphic surface
        :param ignore_zeroes: Ignores zeroes and close to zero values in the specified RDF array
        :return: 3 numpy.arrays of the interpolated and smoothed RDF data
        """
        self.rdf(rho, t, power, par_a, iso_scale)
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
        # TODO: we lose one index at some point could be bc of slicing the upper boundary
        # Use the interpolation functions
        self.rdf_interp = f(self.r_interp)
        self.rdf_interp_smooth = f_smooth(self.r_interp)

        # Generating new separation unit distance for the x-axis
        self.dr = self.rg / range_refinement

        # Passing interpolated data to be stored later in file
        self.interpolated_data.append(self.rdf_interp)
        return self.r_interp, self.rdf_interp, self.rdf_interp_smooth

    def rdf_interpolate_plot(self, rho, t, power, par_a,
                             range_refinement=2000,
                             iso_scale=False, **kwargs):
        """
        Generates a plot of the interpolated data after it calls rdf_interpolate


        :param rho: Density
        :param t: Temperature
        :param power: Pair potential strength
        :param par_a: Softening parameter
        :param range_refinement: The accuracy of the interpolated data
        """
        self.rdf_interpolate(rho, t, power, par_a, range_refinement, iso_scale)

        # Naming the curves
        name = f"rho: {self.rho_str} T: {self.t_str} n: {self.n_str} A: {self.a_str}"
        max_scaling = np.max(self.rdf_data)  # Scaling the plot to ymax

        plt.figure('Interpolated RDF')

        plt.plot(self.r_interp, self.rdf_interp,
                 '-', label='interpolation ' + name, **kwargs)
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
                                    range_refinement=2000,
                                    iso_scale=False,
                                    show_label=True,
                                    **kwargs):
        """
        Generates a plot of the smoothed and interpolated
        data after it calls rdf_interpolate


        :param rho: Density
        :param t: Temperature
        :param power: Pair potential strength
        :param par_a: Softening parameter
        :param range_refinement: The accuracy of the interpolated data
        """
        self.rdf_interpolate(rho, t, power, par_a, range_refinement, iso_scale)

        # Naming the curves
        name = ''
        if show_label is True:
            name = f"rho: {self.rho_str} T: {self.t_str}" \
                   f" n: {self.n_str} A: {self.a_str}"
        max_scaling = np.max(self.rdf_data)  # Scaling the plot to ymax

        plt.figure('Interpolated RDF')

        plt.plot(self.r_interp, self.rdf_interp_smooth,
                 '-.', label='smooth interp ' + name, **kwargs)
        plt.plot([0, self.r_interp[-1]], [1, 1],
                 '--', color='black', linewidth=0.5)

        # Plot limits and legends
        plt.xlim(left=0, right=3)
        plt.ylim(bottom=0, top=max_scaling + 0.1)

        # Plot labels
        plt.xlabel(r"$r$", fontsize=16)
        plt.ylabel(r"$g(r)$", fontsize=16)
        if show_label is True:
            plt.legend(loc="best", fancybox=True)

    def rdf_intersect(self, rho, t, power_list, par_a,
                      range_refinement=2000,
                      r_lower=0, r_higher=-1,
                      intersections=1):
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
        inflection points, near neighbourhoods where the data do not fluctuate much.


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

        # Get the coordinates for the local max and min in the provided range,
        # with the intent of looking between a min-max for intersection points
        # The max, min returned here correspond to the last curve for n in the power_list
        __, __, __, __, idx_max, idx_min = self.find_local_min_max(
            self.r_interp[r_lower:r_higher], self.rdf_interp_smooth[r_lower:r_higher])

        # Plotting the intersection results into the interpolated RDF canvas
        plt.figure('Interpolated RDF')

        # Merge and sort the arrays containing the indices of the inflection points
        idx_max_min = np.concatenate((idx_max, idx_min))
        idx_max_min.sort()

        # Loop through all the max-min combinations and spot the
        # intersect of the curves based on the std
        r_iso_list, mean_iso_list = [], []
        used_max_idx, used_min_idx = [], []
        for i in range(1, len(idx_max_min)):
            # Filter out values that have small fluctuations in the y-axis
            # and/or are closely located in the x-axis
            tolerance = 0.05
            if (abs(mean_list[idx_max_min[i]] - mean_list[idx_max_min[i-1]]) > tolerance) and \
                    (idx_max_min[i] - idx_max_min[i-1] >= range_refinement/20):
                # Get the index of the k smallest stds, which in theory
                # should correspond to the point where the curves intersect
                idx_intersect = np.argpartition(
                    std_list[idx_max_min[i-1]:idx_max_min[i]], intersections)
                # Adjust index to match with global(interpolated) array index
                idx_intersect += idx_max_min[i-1]

                # Get the mean for the intersection points
                mean_scatter = mean_list[idx_intersect[:intersections]]

                # Add the lower boundary index to the std,
                # to correspond to the actual r index
                idx_intersect += r_lower

                # Get the r-values for the corresponding RDF averaged values
                r_iso = [self.r_interp[i]
                         for i in idx_intersect[:intersections]]

                # Plot only the max used
                # Extract the indices of the max/min that is used
                if (idx_max_min[i] in idx_max) and (idx_max_min[i] not in used_max_idx):
                    used_max_idx.append(idx_max_min[i])
                if idx_max_min[i-1] in idx_max and (idx_max_min[i-1] not in used_max_idx):
                    used_max_idx.append(idx_max_min[i-1])
                if idx_max_min[i] in idx_min and (idx_max_min[i] not in used_min_idx):
                    used_min_idx.append(idx_max_min[i])
                if idx_max_min[i-1] in idx_min and (idx_max_min[i-1] not in used_min_idx):
                    used_min_idx.append(idx_max_min[i-1])

                # BUG: using [0] does not account if multiple intersections are chosen
                # Storing the x-values for export, r_iso is a list of a single element
                r_iso_list.append(r_iso[0])
                # Storing the y-value for the intersection
                mean_iso_list.append(mean_scatter[0])

        # Plot the intersection points
        plt.scatter(r_iso_list, mean_iso_list, marker='x', color='red', s=100)

        # Plotting the local maxima and minima of the RDF
        # Uses the extracted indices to match them to
        # the smoothed RDF data
        used_x_max = [self.r_interp[i+r_lower] for i in used_max_idx]
        used_y_max = [self.rdf_interp_smooth[i+r_lower] for i in used_max_idx]
        used_x_min = [self.r_interp[i+r_lower] for i in used_min_idx]
        used_y_min = [self.rdf_interp_smooth[i+r_lower] for i in used_min_idx]

        # Plot the local max, min that contain an intersection between them
        plt.scatter(used_x_max, used_y_max, color='orange')
        plt.scatter(used_x_min, used_y_min, color='green')

        # Plotting the data curves of the interpolated data
        for n in power_list:
            self.rdf_interpolate_smooth_plot(
                rho, t, n, par_a, range_refinement,
                show_label=False)
            self.rdf_plot(
                rho, t, n, par_a,
                show_label=False)

        return r_iso_list, mean_iso_list

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

        # Realigning for convenience
        idx_local_max = idx_local_max[0]
        idx_local_min = idx_local_min[0]

        # Fetch the last global minimum
        g_min = np.where(y == y.min())[-1]

        # Test to see if global min is already present
        # TODO: element comparisson will be deprecated
        # https://stackoverflow.com/questions/40659212/futurewarning-elementwise-comparison-failed-returning-scalar-but-in-the-futur?rq=1
        if g_min not in idx_local_min:
            idx_local_min = np.append(idx_local_min, g_min)
            idx_local_min.sort()

        # Get the x-value for the local maxima and minima
        x_local_max = x[idx_local_max]
        x_local_min = x[idx_local_min]
        y_local_max = y[idx_local_max]
        y_local_min = y[idx_local_min]

        return x_local_max, y_local_max, x_local_min, y_local_min, idx_local_max, idx_local_min

    def get_intersections_to_file(self, rho_list, t_list, n_list, a_list,
                                  filename,
                                  delimiter='\t'):
        """
        Writes the isosbestic points coordinates (r_iso and rdf_iso) to two different files


        :param rho_list: list of densities
        :param t_list: list of temperatures
        :param n_list: list of pair potential strenghts
        :param a_list: list of softening parameters
        :param filename: Output filename/ directory
        :param delimiter: string that separates the data in the files
        """

        x_iso = f"{filename}/r_iso.dat"
        y_iso = f"{filename}/rdf_iso.dat"
        with open(x_iso, 'w+') as f_x, open(y_iso, 'w+') as f_y:
            f_x.write('rho\tT\ta\tr_iso\n')
            f_y.write('rho\tT\ta\tr_iso\n')
            for rho in rho_list:
                for t in t_list:
                    for a in a_list:
                        r_iso, rdf_iso = self.rdf_intersect(
                            rho, t, n_list, a, r_lower=100)

                        line_x = f"{rho}{delimiter}" \
                                 f"{t}{delimiter}" \
                                 f"{a}{delimiter}{r_iso}\n"
                        line_x = line_x.replace('[', '')
                        line_x = line_x.replace(']', '')
                        line_x = line_x.replace(', ', delimiter)

                        line_y = f"{rho}{delimiter}{t}" \
                                 f"{delimiter}{a}" \
                                 f"{delimiter}{rdf_iso}\n"
                        line_y = line_y.replace('[', '')
                        line_y = line_y.replace(']', '')
                        line_y = line_y.replace(', ', delimiter)

                        f_x.write(line_x)
                        f_y.write(line_y)

    def plot_intersection(self, rho, t):
        # Read rho and T from file if it matches rho and T read
        data = "r_iso.dat"
        rho_list, t_list, a_list, r_iso_list = np.loadtxt(
            data, usecols=(0, 1, 2, 3), unpack=True, skiprows=1)

        indices = [i for i, x in enumerate(t_list) if x == t]
        rho_list = rho_list[indices]
        t_list = t_list[indices]
        a_list = a_list[indices]
        r_iso_list = r_iso_list[indices]

        indices = [i for i, x in enumerate(rho_list) if x == rho]
        rho_list = rho_list[indices]
        t_list = t_list[indices]
        a_list = a_list[indices]
        r_iso_list = r_iso_list[indices]
        plt.figure(f"r_iso Intersections with T: {t}")
        name = fr"$\rho$: {rho}  T: {t}"
        plt.scatter(a_list, r_iso_list, label=name)
        plt.legend(loc="best")
