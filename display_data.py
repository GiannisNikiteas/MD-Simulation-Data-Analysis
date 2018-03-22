from PathHandling import *
from plotting_class import *


class Isomorph:
    def __init__(self, t_r, rho_r, a_r, t_out):
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
        for i in range(len(self.t_out)):
            t = self.t_out[i]
            rho_out = self.get_rho(self.rho_r, self.t_r, t, n)
            a_out = self.get_a(self.a_r, self.rho_r, rho_out)
            self.rho_out.append(rho_out)
            self.a_out.append(a_out)
        return self.rho_out, self.a_out


if __name__ == '__main__':
    os.chdir('../../Archives of Data/')
    obj = FilePlotting(10000, 1000)  # steps, particles

    my_list = [0, 0.25, 0.50, 0.65, 0.7, 0.75, 0.8, 0.85, 0.90, 0.97, 1.00, 1.1, 1.2,
               1.25, 1.50, 1.75, 2.00, 2.25, 2.50, 2.75, 4.00]
    a_list = [0, 0.25, 0.50, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25]
    n_list = [6, 8, 10, 12]

    # RDF Isomorph space creation
    t_iso = np.linspace(0.5, 3, 5)
    isomorph_line = Isomorph(0.5, 0.5, 0.5, t_iso)
    rho_iso, a_iso = isomorph_line.gen_line(8)

    # for i in range(len(t_iso)):
    #     obj.vaf(rho_iso[i], t_iso[i], 8, a_iso[i])
    #     obj.scaled_potential(rho_iso[i], 8, a_iso[i])
    #     obj.rdf(rho_iso[i], t_iso[i], 8, a_iso[i])

    for i in n_list:
        obj.potential(i, 1)
    #     obj.scaled_potential(0.5, i, 0.5)
        # obj.avg_pressure(1, 1, i)
        # obj.rdf(0.5, 0.5, i, 1)
        # obj.diffusion_plot(0.5, 1.0, 12, a_list)
        # obj.avg_pressure(i)
        # obj.RDF2(6, i)

    plt.tight_layout()
    plt.show()
