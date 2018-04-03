from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np

"""
    Functional form of the isomorph generation.
    The file generates a the contour for all the possible parameter a values, for a given T0 and a0.
"""


def iso_function(rho, t, rho0=1.5, t0=1.0, a0=0.5, n=8):
    a = a0 * (rho/rho0 * (t0/t) ** (-6/n)) ** (1.0/3.0)
    return a


fig = plt.figure('3D contour plot, with function')
ax = fig.gca(projection='3d')

rho = np.arange(0.1, 5.0, 0.05)
t = np.arange(0.1, 5.0, 0.05)
RHO, T = np.meshgrid(rho, t)
a = np.array([iso_function(rho, t) for rho, t in zip(np.ravel(RHO), np.ravel(T))])
A = a.reshape(RHO.shape)
ax.plot_wireframe(RHO, T, A, color='blue', alpha=0.5, rstride=2, cstride=2)

################################################################
rho = np.arange(0.1, 5.0, 0.05)
t = np.arange(0.1, 5.0, 0.05)
RHO, T = np.meshgrid(rho, t)
a = np.array([iso_function(rho, t, rho0=1, t0=1, a0=1) for rho, t in zip(np.ravel(RHO), np.ravel(T))])
A = a.reshape(RHO.shape)
ax.plot_wireframe(RHO, T, A, color='red', alpha=0.8, rstride=2, cstride=2)

################################################################
# Labels and Axis limits
ax.set_xlabel(r'$\rho$')
# ax.set_xlim(np.amin(RHO), np.amax(RHO))
ax.set_ylabel(r'T')
# ax.set_ylim(np.amin(T), np.amax(T))
ax.set_zlabel(r'a')
# ax.set_zlim(np.amin(A), np.amax(A))

# Grid background color set to white
ax.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
ax.w_yaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
ax.w_zaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))

plt.show()
