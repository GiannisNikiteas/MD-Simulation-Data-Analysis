from PathHandling import *
from plotting_class import *

os.chdir('../../Archives of Data/')
obj = FilePlotting(10000, 1000)  # steps, particles

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

#                     rho  T   n   a
obj.radial_dist_func(0.5, 0.5, 8, 0.5)
obj.radial_dist_func(0.8409, 2, 8, 0.42045)
o

plt.tight_layout()
plt.show()
