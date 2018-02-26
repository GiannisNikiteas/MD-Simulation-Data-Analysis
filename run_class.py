from PathHandling import *
from plotting_class import *
import path

# TODO: Adjust RDF plotting method not start from 0. Start from next list entry since 0 is not measured
# TODO: make methods in plotting class return strings with plotnames
#       in order to save figures
# TODO: Obvious CON when plotting multiple plots on the same figure
# -> the final image will have the name of the last plot called (not much you can do about it
path = OSPaths()
path.dir('T_0.5', 'Density_0.5', '10000')
obj = FilePlotting()  # constructor is left empty

my_list = [0, 0.25, 0.50, 0.65, 0.7, 0.75, 0.8, 0.85, 0.90, 0.97, 1.00, 1.1, 1.2,
           1.25, 1.50, 1.75, 2.00, 2.25, 2.50, 2.75, 4.00]
a_list = [0, 0.25, 0.5, 0.75, 0.97, 1, 2, 3, 4]
n_list = [6, 8, 10, 12]

# obj.vel_dist(12, 4)
# for i in a_list:
#     # obj.avg_pressure(i)
#     obj.radial_dist_func(12, i)
#     # obj.avg_pressure(i)
#     # obj.RDF2(6, i)
#     # obj.savefig('rdf_a0_n6-12.pdf')


obj.potential(12,1), obj.potential(12, 0.75)
# obj.potential_data(12, 0.89)
# obj.radial_dist_func(12, 0.89)
# os.chdir('../T_1.0_step_10000')
# obj.potential_data(12, 1)
# obj.radial_dist_func(12, 1)
plt.tight_layout()
plt.show()
