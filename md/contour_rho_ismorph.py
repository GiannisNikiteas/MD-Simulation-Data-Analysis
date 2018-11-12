import matplotlib.pyplot as plt
from isomorphs import *
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

from md import Isomorph

"""
    The program plots the isomorphic surfaces generated for a fluid.
    A reference Temperature and a parameter are given. 
    The program generates a list of isomorphic points (t_iso, a_iso) for the given references parameters.
    The program proceeds to calculate the isomorphic lines for a series of reference densities rho_ref.
    Output is a contour plot for of the isomorph plane.
"""
rho_ref = np.linspace(0.2, 1.0, 10)   # horizontal par
n = 8
# List of isomorphic temperatures:
# num, adjusts horizontal refinement of surface
t_iso_line = np.linspace(0.2, 2.0, 20)
# Lists for isomorphic density, parameter a and temperature:
rho_iso = np.empty((0, len(t_iso_line)))
a_iso = np.empty((0, len(t_iso_line)))
t_iso = np.empty((0, len(t_iso_line)))

# Plot contours
fig = plt.figure('3D contour plot')
ax = fig.gca(projection='3d')

# Generate matrix of isomorph quantities for various rhos
for i in rho_ref:
    isomorph_line = Isomorph(i, 0.2, 0.5, t_iso_line)  # This is the reference point
    temp_rho, temp_a = isomorph_line.gen_line(n)
    rho_iso = np.append(rho_iso, [temp_rho], axis=0)
    a_iso = np.append(a_iso, [temp_a], axis=0)
    t_iso = np.append(t_iso, [t_iso_line], axis=0)
    # ax.scatter(i, 1, 1, marker='+', color='red')

ax.plot_wireframe(rho_iso, t_iso, a_iso, alpha=0.7, color='cyan')
# Contour projections on 3D plot, enable if needed for visualisation
# cset = ax.contourf(rho_iso, a_iso, t_iso, zdir='z', offset=np.amin(t_iso), cmap=cm.coolwarm)  # m.coolwarm
# cset = ax.contourf(rho_iso, a_iso, t_iso, zdir='x', offset=np.amin(rho_iso), cmap=cm.coolwarm, alpha=0.5)
# cset = ax.contourf(rho_iso, a_iso, t_iso, zdir='y', offset=np.amax(a_iso), cmap=cm.coolwarm, alpha=0.5)

# Labels and Limits
ax.set_xlabel(r'$\rho$')
# ax.set_xlim(0, np.amax(rho_iso))  # todo: test why max > rho_ref max
# ax.set_xlim(0, 5)  # todo: test why max > rho_ref max
ax.set_ylabel(r'T')
# ax.set_ylim(0, np.amax(t_iso))
ax.set_zlabel(r'a')
# ax.set_zlim(np.amin(a_iso), np.amax(a_iso))

# Grid background color set to white
ax.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
ax.w_yaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
ax.w_zaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))

# Contour projections plots, enable if needed
# plt.figure('Separate contours')
# plt.subplot(131)    # X. Y, Z
# cset = plt.contourf(rho_iso, t_iso, a_iso)
# plt.xlabel(r'$\rho$')
# plt.ylabel(r'T')
# plt.subplot(132)    # Y, Z, X
# cset = plt.contourf(t_iso, a_iso, rho_iso)
# plt.xlabel(r'T')
# plt.ylabel(r'a')
# plt.subplot(133)    # X, Z, Y
# cset = plt.contourf(rho_iso, a_iso, t_iso)
# plt.xlabel(r'$\rho$')
# plt.ylabel(r'a')

plt.tight_layout(pad=0.1, h_pad=0, w_pad=0)
plt.show()
