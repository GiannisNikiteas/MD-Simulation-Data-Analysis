from PathHandling import *
from plotting_class import *


os.chdir('C:/Code/C++/MD simulation/Archives of Data/')
obj = FilePlotting(50000, 1000)  # steps, particles

my_list = [0, 0.25, 0.50, 0.7, 0.75, 0.8, 0.85, 0.90, 0.95,
           1.00, 1.1, 1.2, 1.25, 1.50, 1.75, 2.00]
a_list = [0, 0.50, 0.75, 1.0, 1.25, 2.0]
n_list = [6, 8, 10, 12]

for i in a_list:
    obj.rdf(0.5, 1, 6, i, iso_scale=False, show_iso=False)
    # obj.vaf(0.5, 1, 6, i)


# obj.rdf(0.5, 0.5, 8, 0.5, iso_scale=True)
# obj.rdf(0.8409, 2.0, 8, 0.42045, iso_scale=True)
# obj.diffusion_plot(0.5, 0.5, 6, a_list)
# obj.potential(6, 1.5)
# obj.potential(8, 1.5)
# obj.potential(12, 1.5)

plt.show()
