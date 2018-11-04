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
        self.rho_out = []  # LIST
        self.a_out = []    # LIST

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
        self.rho_out = []
        self.a_out = []
        for i in range(len(self.t_out)):
            t = self.t_out[i]
            rho_out = self.get_rho(self.rho_r, self.t_r, t, n)
            a_out = self.get_a(self.a_r, self.rho_r, rho_out)
            self.rho_out.append(rho_out)
            self.a_out.append(a_out)
        return self.rho_out, self.a_out
