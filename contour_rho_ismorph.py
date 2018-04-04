from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
from matplotlib import cm
from isomorphs import *
import numpy as np


"""
    The program plots the isomorphic surfaces generated for a fluid.
    A reference Temperature and a parameter are given. 
    The program generates a list of isomorphic points (t_iso, a_iso) for the given references parameters.
    The program proceeds to calculate the isomorphic lines for a series of reference densities rho_ref.
    Output is a contour plot for of the isomorph plane.
"""
rho_ref = np.linspace(0.1, 5, 20)   # horizontal par
n = 8
# List of isomorphic temperatures:
# num, adjusts horizontal refinement of surface
t_iso_line = np.linspace(0.01, 5, 5000)
# Lists for isomorphic density, parameter a and temperature:
rho_iso = np.empty((0, len(t_iso_line)))
a_iso = np.empty((0, len(t_iso_line)))
t_iso = np.empty((0, len(t_iso_line)))

# Plot contours
fig = plt.figure('3D contour plot')
ax = fig.gca(projection='3d')

# Generate matrix of isomorph quantities for various rhos
# todo: This could potentially be wrong
# todo: The wireframe here is lines through contour planes at some points (don't know why these points)
# todo: or could be simply wrong. Think through the how does T get fed into the isomorph
# todo: and what t1 and t2 are used to generate the data
# todo: examine how I end up with higher densities than max of rho_ref
for i in range(len(rho_ref)):
    isomorph_line = Isomorph(1, rho_ref[i], 1, t_iso_line)  # This is the reference point, T=1,rho?,a=1
    temp_rho, temp_a = isomorph_line.gen_line(n)
    rho_iso = np.append(rho_iso, [temp_rho], axis=0)
    a_iso = np.append(a_iso, [temp_a], axis=0)
    t_iso = np.append(t_iso, [t_iso_line], axis=0)
    ax.scatter(rho_ref[i], 1, 1, marker='+', color='red')

ax.plot_wireframe(rho_iso, t_iso, a_iso, alpha=0.7, color='cyan')
# Contour projections on 3D plot, enable if needed for visualisation
# cset = ax.contourf(rho_iso, a_iso, t_iso, zdir='z', offset=np.amin(t_iso), cmap=cm.coolwarm)  # m.coolwarm
# cset = ax.contourf(rho_iso, a_iso, t_iso, zdir='x', offset=np.amin(rho_iso), cmap=cm.coolwarm, alpha=0.5)
# cset = ax.contourf(rho_iso, a_iso, t_iso, zdir='y', offset=np.amax(a_iso), cmap=cm.coolwarm, alpha=0.5)

# Labels and Limits
ax.set_xlabel(r'$\rho$')
ax.set_xlim(0, np.amax(rho_iso))
ax.set_ylabel(r'T')
ax.set_ylim(0, np.amax(t_iso))
ax.set_zlabel(r'a')
ax.set_zlim(np.amin(a_iso), np.amax(a_iso))

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
