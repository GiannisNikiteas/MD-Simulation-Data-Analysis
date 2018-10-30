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
    obj.rdf_interpolate_smooth_plot(
        rho_iso[i], t_iso[i], 8, a_iso[i], iso_scale=True)

plt.show()
