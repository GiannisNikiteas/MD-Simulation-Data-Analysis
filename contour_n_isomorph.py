from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
from matplotlib import cm
from isomorphs import *
import numpy as np
import itertools


"""
 Potential strength list:
 The larger the n-range, the wider the surface diagonal
 For small n e.g. 8-12, the isomorph contours are concave and monotonic for an increasing T
 As n increases e.g. 25-30, the isomorphs, become more linear, with a more prompt increase of rho, a vs T
"""
n_list = np.arange(6, 19, dtype=int)

# List of isomorphic temperatures:
# num, adjusts horizontal refinement of surface
t_iso_line = np.linspace(0.05, 5, 15)
# Lists for isomorphic density, parameter a and temperature:
rho_iso = np.empty((0, len(t_iso_line)))
a_iso = np.empty((0, len(t_iso_line)))
t_iso = np.empty((0, len(t_iso_line)))
# Reference seed with rho=0.5, T=0.5, a=0.5:
# Note: if reference T > than any T in t_iso_line, then a "flipping" of the isomorph is observed
# Explanation not yet given!
isomorph_line = Isomorph(0.5, 0.5, 0.5, t_iso_line)

# Generate matrix of isomorph quantities for various n s
for i in range(len(n_list)):
    temp_rho, temp_a = isomorph_line.gen_line(n_list[i])
    rho_iso = np.append(rho_iso, [temp_rho], axis=0)
    a_iso = np.append(a_iso, [temp_a], axis=0)
    t_iso = np.append(t_iso, [t_iso_line], axis=0)

# Plot contours
fig = plt.figure('3D contour plot')
ax = fig.gca(projection='3d')
for i in range(len(n_list)):
    # ax.scatter(rho_iso[i], a_iso[i], t_iso[i], alpha=0.9, linestyle='-',
    #            marker=marker[i], label='n: ' + str(n_list[i]))
    ax.plot_surface(rho_iso, a_iso, t_iso, alpha=0.6, label='n: ' + str(n_list[i]))

# ax.contourf(rho_iso, a_iso, t_iso, alpha=0.9)
# Projections of X, Y, Z onto corresponding planes in 3D plot
# cset = ax.contourf(rho_iso, a_iso, t_iso, zdir='z', offset=np.amin(t_iso), cmap=cm.coolwarm)  # m.coolwarm
# cset = ax.contourf(rho_iso, a_iso, t_iso, zdir='x', offset=np.amin(rho_iso), cmap=cm.coolwarm, alpha=0.5)
# cset = ax.contourf(rho_iso, a_iso, t_iso, zdir='y', offset=np.amax(a_iso), cmap=cm.coolwarm, alpha=0.5)

# Labels and Limits
ax.set_xlabel(r'$\rho$')
# ax.set_xlim(np.amin(rho_iso), np.amax(rho_iso))
ax.set_ylabel(r'a')
# ax.set_ylim(np.amin(a_iso), np.amax(a_iso))
ax.set_zlabel(r'T')
# ax.set_zlim(np.amin(t_iso), np.amax(t_iso))
# ax.legend(loc='best', fancybox=True, fontsize='small')

# plt.figure('Separate contours')
# plt.subplot(131)    # X. Y, Z
# cset = plt.contourf(rho_iso, a_iso, t_iso, cmap='coolwarm')
# plt.xlabel(r'$\rho$')
# plt.ylabel(r'a')
# plt.subplot(132)    # Y, Z, X
# cset = plt.contourf(a_iso, t_iso, rho_iso, cmap='coolwarm')
# plt.xlabel(r'a')
# plt.ylabel(r'T')
# plt.subplot(133)    # X, Z, Y
# cset = plt.contourf(rho_iso, t_iso, a_iso, cmap='coolwarm')
# plt.xlabel(r'$\rho$')
# plt.ylabel(r'T')

plt.tight_layout(pad=0.1, h_pad=0, w_pad=0)
plt.show()
