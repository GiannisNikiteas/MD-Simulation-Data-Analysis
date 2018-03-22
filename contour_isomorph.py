from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
from matplotlib import cm
from isomorphs import *
import numpy as np


# Potential strength list
n_list = range(6, 25)
# List of isomorphic temperatures
t_iso_line = np.linspace(0.1, 5, 100)
# Lists for isomorphic density and a
rho_iso = np.empty((0, len(t_iso_line)))
a_iso = np.empty((0, len(t_iso_line)))
t_iso = np.empty((0, len(t_iso_line)))
# Reference seed with rho=0.5, T=0.5, a=0.5
isomorph_line = Isomorph(0.1, 0.5, 0.5, t_iso_line)

# Generate matrix of isomorph quantities for various n s
for i in range(len(n_list)):
    temp_rho, temp_a = isomorph_line.gen_line(n_list[i])
    rho_iso = np.append(rho_iso, [temp_rho], axis=0)
    a_iso = np.append(a_iso, [temp_a], axis=0)
    t_iso = np.append(t_iso, [t_iso_line], axis=0)

# Plot contours
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot_surface(rho_iso, a_iso, t_iso, rstride=8, cstride=8, alpha=0.3)
# cset = ax.contourf(rho_iso, a_iso, t_iso, zdir='z', cmap=cm.coolwarm)
# cset = ax.contourf(rho_iso, a_iso, t_iso, zdir='x', cmap=cm.coolwarm)
# cset = ax.contourf(rho_iso, a_iso, t_iso, zdir='y', cmap=cm.coolwarm)

ax.set_xlabel(r'\rho')
# ax.set_xlim(-40, 40)
ax.set_ylabel(r'a')
# ax.set_ylim(-40, 40)
ax.set_zlabel(r'T')
# ax.set_zlim(-100, 100)

plt.show()
