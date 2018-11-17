class Isomorph:
    """
    Isomorph state generator for fluid transitioning from MD to continuum limit
    """

    def __init__(self, rho_r, t_r, a_r, t_out):
        """
        :param t_r: Reference Temperature
        :param rho_r: Reference density
        :param a_r: Reference A parameter
        :param t_out: Range of temperatures the isomorph will be developed
        """
        self.t_r = t_r  # Reference T
        self.rho_r = rho_r  # Reference density
        self.a_r = a_r  # Reference A par
        self.t_out = t_out   # LIST Isomorph T
        self.rho2_list = []  # LIST
        self.a2_list = []    # LIST

    @staticmethod
    def get_rho(rho1, t1, t2, n):
        return rho1 * (t2 / t1) ** (3.0 / n)

    @staticmethod
    def get_a(a1, rho1, rho2):
        a2 = a1 * (rho1 / rho2) ** (1.0 / 3.0)
        return a2

    def gen_line(self, n):
        """
        :param n: Potential power strength of the pair potential
        :return: Output for Density and A of the isomorph, along the given rho, A and T reference point
                 and the using T_OUT as a range of values for the isomorph
        """
        self.rho2_list = []
        self.a2_list = []
        # Extract T2, from the range of temperatures
        for t2 in self.t_out:
            # Generate the rho2, for T1, T2, rho1, a1
            rho_out = self.get_rho(self.rho_r, self.t_r, t2, n)
            # Generate a2 for T1, T2, rho1, rho2, a1
            a_out = self.get_a(self.a_r, self.rho_r, rho_out)
            # Pass rho2, a2 to lists
            self.rho2_list.append(rho_out)
            self.a2_list.append(a_out)
        return self.rho2_list, self.a2_list

    def get_iso_point(self, t2, n):
        """
        Generate a single isomorphic point by providing a new temperature t2.

        :param t2:  Temperature where the isomorphic point will be calculated at
        :param n:  Potential strength
        :return:  x, y, z, coordinates of the isomorphic point
        """
        rho2 = self.get_rho(self.rho_r, self.t_r, t2, n)
        a2 = self.get_a(self.a_r, self.rho_r, rho2)
        return rho2, t2, a2
