import matplotlib.pyplot as plt
import numpy as np
from mdtools.stat_quantities import FileNaming


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
        kin_f.locator_params(axis='y', nbins=4), kin_f.set_ylim(top=4)
        pot_f.plot(x, pot_en, 'g')
        pot_f.locator_params(axis='y', nbins=3)
        pot_f.set_ylabel("Energy units", size=16)
        tot_f.plot(x, tot_en, 'b')
        tot_f.locator_params(axis='y', nbins=4)
        tot_f.set_ylim(top=6)

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
        all_f.set_ylim(top=5)

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
        a, pc = np.loadtxt(pc_name, delimiter='\t',
                           comments='#', usecols=(0, 5), unpack=True)

        plt.figure('Average Pressure')
        plt.plot(a, pc, '-o', label=name, markersize=3)
        plt.xlim(left=0, right=4.0)
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

        a, k = np.loadtxt(k_name, delimiter='\t', comments='#',
                          usecols=(0, 2), unpack=True)

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

        a, k = np.loadtxt(u, delimiter='\t', comments='#',
                          usecols=(0, 3), unpack=True)

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

        a, K, U, T = np.loadtxt(e_name, delimiter='\t',
                                comments='#', usecols=(0, 2, 3, 4), unpack=True)

        fig = plt.figure('Average Energies')
        kin_f = plt.subplot2grid((3, 2), (0, 0), colspan=1)
        pot_f = plt.subplot2grid((3, 2), (1, 0), colspan=1)
        tot_f = plt.subplot2grid((3, 2), (2, 0), colspan=1)
        all_f = plt.subplot2grid((3, 2), (0, 1), rowspan=3)

        kin_f.plot(a, K, color='r')
        kin_f.set_ylim(bottom=2.0), kin_f.locator_params(
            axis='y', nbins=5, prune='lower')
        pot_f.plot(a, U, color='g')
        pot_f.locator_params(
            axis='y', nbins=4), pot_f.set_ylabel("Energy units")
        tot_f.plot(a, T, color='c')
        tot_f.locator_params(axis='y', nbins=3), tot_f.locator_params(
            axis='x', nbins=4)
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
        all_f.set_ylim(top=6)
        all_f.locator_params(axis='x', nbins=4)
