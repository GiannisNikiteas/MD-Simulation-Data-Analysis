import os
import numpy as np
from md import StatQ, Isomorph
import matplotlib.pyplot as plt

os.chdir("/home/gn/Desktop/test_data/isomorph")

# Generate the isomorph rho, T and A
# using a reference point of 0.5, 0.5 0.5

t_iso = np.linspace(0.5, 2.0, 5)
iso = Isomorph(0.5, 0.5, 0.5, t_iso)
rho_iso, a_iso = iso.gen_line(8)

obj = StatQ(15000, 10**3)

for i in range(len(t_iso)):
    obj.vaf(rho_iso[i], t_iso[i], 8, a_iso[i], iso_scale=True)

# Find the isosbestic points
os.chdir("/home/gn/Desktop/test_data")
n = [6, 7, 8, 9, 10, 12]
rho = [0.3, 0.5, 1.0, 1.5]
t = [0.5, 1.0, 2.0]
a = [0, 0.25, 0.50, 0.7, 0.75, 0.8, 0.90]
obj.get_intersections_to_file(rho, t, n, a, "/home/gn/Desktop")

plt.show()
