from PathHandling import *
from plotting_class import *
from isomorphs import *

if __name__ == '__main__':
    os.chdir('../../Archives of Data/')
    obj = FilePlotting(10000, 1000)  # steps, particles

    my_list = [0, 0.25, 0.50, 0.65, 0.7, 0.75, 0.8, 0.85, 0.90, 0.97, 1.00, 1.1, 1.2,
               1.25, 1.50, 1.75, 2.00, 2.25, 2.50, 2.75, 4.00]
    a_list = [0.50, 0.75, 1.0, 1.25]
    n_list = [6, 8, 10, 12]

    # RDF Isomorph space creation
    t_iso = np.linspace(0.5, 3, 5)
    isomorph_line = Isomorph(0.5, 0.5, 0.5, t_iso)
    rho_iso, a_iso = isomorph_line.gen_line(8)

    # for i in range(len(t_iso)):
    #     obj.vaf(rho_iso[i], t_iso[i], 8, a_iso[i])
    #     obj.scaled_potential(rho_iso[i], 8, a_iso[i])
    #     obj.rdf(rho_iso[i], t_iso[i], 8, a_iso[i])

    # for i in n_list:
    #     obj.rdf(0.5, 0.5, i, 0., iso_scale=False)
    for i in n_list:
        obj.rdf(0.5, 0.5, i, 0, iso_scale=False, show_iso=True)
        # obj.potential(i, 1)
    #     obj.scaled_potential(0.5, i, 0.5)
        # obj.avg_pressure(1, 1, i)
        # obj.diffusion_plot(0.5, 1.0, 12, a_list)
        # obj.avg_pressure(i)
        # obj.RDF2(6, i)

    plt.tight_layout()
    plt.show()
