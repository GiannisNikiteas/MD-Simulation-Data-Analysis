from src.plotting_class import *
import os

os.chdir('/home/gn/Desktop/test_data')
obj = FilePlotting(10000, 1000)  # steps, particles

my_list = [0, 0.25, 0.50, 0.7, 0.75, 0.8, 0.85, 0.90, 0.95,
           1.00, 1.1, 1.25, 1.50, 1.75, 2.00]
a_list = [0, 0.50, 0.75, 1.0, 1.25, 2.0]
n_list = [6, 8, 10, 12]
t_list = np.linspace(0.001, 0.01, 10)

for i in n_list:
    obj.rdf(0.2, 0.5, i, 0.0)
    # obj.rdf_interpolate(0.5, 1, i, 0.5)
    # obj.msd(0.5, 0.5, i, 0.5)
    # obj.rdf(1, 1, i, 0, iso_scale=False, show_iso=False)
    # obj.vaf(0.5, 1, 6, i)

# obj.rdf(0.5, 0.5, 6, 0)
# obj.particle_plot(0.5, 0.5, 6, 0)
# obj.find_intersect(obj.interpolated_data)
# obj.rdf(0.5, 0.5, 6, 0)
# obj.rdf(0.5, 1., 12, 0)


plt.show()
