# %%
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
from mdtools import StatQ, Isomorph, RDFAnalysis

# Find the isosbestic points
os.chdir("/home/gn/Desktop/test_data/long")

n = [8, 10, 12]
rho = [0.20, 0.3, 0.5, 1.0]
t = [0.5, 1.0, 1.5, 2.0]
a = [0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]   # 0.1, 0.5

# Generating the RDF intersection data
# ! this generates the RDF Interpolate plot
rdf_stat = RDFAnalysis(35000, 1000)
r_iso_dir = "/home/gn/Desktop/test_data/long"
if not os.path.exists(f"{r_iso_dir}/r_iso.dat"):
    rdf_stat.get_intersections_to_file(rho, t, n, a, r_iso_dir)


# Plot the intersection points r_iso vs a
for temp in t:
    for density in rho:
        rdf_stat.plot_intersection(density, temp)
    az = np.linspace(0, 1, 500)
    tr = (1 - az**2)**(1./2.)
    plt.plot(az, tr, color='red', linewidth='3')
    plt.show()
