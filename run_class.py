from PathHandling import *
from plotting_class import *

path = OSPaths()
path.dir('Density 0.5', '10000')
obj = FilePlotting()  # constructor is left empty

my_list = [0, 0.25, 0.50, 0.65, 0.7, 0.75, 0.8, 0.85, 0.90, 1.00, 1.1, 1.2, 1.25, 1.50, 1.75, 2.00, 2.25, 2.50,
           2.75, 4.00]

n_list = [6, 8, 10, 12]

for i in n_list:
    obj.radial_dist_func(i, 0.75)

plt.tight_layout()

plt.show()
