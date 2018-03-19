from PathHandling import *
from plotting_class import *

os.chdir('../../Archives of Data/')
obj = FilePlotting(10000, 1000)  # steps, particles

my_list = [0, 0.25, 0.50, 0.65, 0.7, 0.75, 0.8, 0.85, 0.90, 0.97, 1.00, 1.1, 1.2,
           1.25, 1.50, 1.75, 2.00, 2.25, 2.50, 2.75, 4.00]
a_list = [0, 0.25, 0.50, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25]
n_list = [6, 8, 10, 12]
temperature = np.linspace(0.6, 3, 5)
obj.radial_dist_func(0.5, 0.5, 8, 0.5)
obj.radial_dist_func(0.8409, 2.0, 8, 0.42045)
# obj.radial_dist_func(0.5, 0.5, 12, 0.5)
# obj.radial_dist_func(0.52, 0.6, 12, 0.4925)
# obj.radial_dist_func(0.62, 1.2, 12, 0.4648)
# obj.radial_dist_func(0.69, 1.8, 12, 0.4494)
# obj.radial_dist_func(0.74, 2.4, 12, 0.4387)
# obj.radial_dist_func(0.78, 3.0, 12, 0.4306)

# for i in n_list:
    # obj.avg_pressure(1, 1, i)
    # obj.radial_dist_func(1.5, 2, i, 1.5)
    #obj.diffusion_plot(0.5, 1.0, 12, a_list)
    #obj.avg_pressure(i)
    #obj.RDF2(6, i)





#                       rho  T   n   a
# obj.radial_dist_func(0.5, 0.5, 8, 0.5)
# obj.radial_dist_func(0.8409, 2, 8, 0.42045)
# obj.velocity_autocorrelation_func(0.5, 0.5, 8, 0.5)
# obj.velocity_autocorrelation_func(0.8409, 2, 8, 0.42045)

plt.tight_layout()
plt.show()
