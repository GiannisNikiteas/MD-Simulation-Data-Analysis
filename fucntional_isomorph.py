from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np

"""
    Functional form of the isomorph generation.
    The file generates a the contour for all the possible parameter a values, for a given T0 and a0.
"""


def iso_function(rho, t, rho0=1.0, t0=1.0, a0=0.5, n=8):
    a = a0 * (rho/rho0 * (t0/t) ** (-6/n)) ** (1.0/3.0)
    return a


# Arrays of X, Y planes
rho = np.arange(0.1, 5.0, 0.05)
t = np.arange(0.1, 5.0, 0.05)
RHO, T = np.meshgrid(rho, t)

fig = plt.figure('3D contour plot, with function')
ax = fig.gca(projection='3d')
# Simple contour plot, USE FUNCTION plot_isomorph_contour instead
# a = np.array([iso_function(rho, t) for rho, t in zip(np.ravel(RHO), np.ravel(T))])
# A = a.reshape(RHO.shape)
# ax.plot_wireframe(RHO, T, A, color='blue', alpha=0.5, rstride=2, cstride=2)
################################################################


def label_name(_a0, _n, _rho=None, _t0=None):
    rho_str = '{:.4f}'.format(_rho)
    t_str = '{:.4f}'.format(_t0)
    a_str = '{:.5f}'.format(_a0)
    n_str = str(_n)

    name_id = r'$\rho$: ' + rho_str + 'T: ' + t_str + 'a: ' + a_str + 'n: ' + n_str
    if (_rho is None) or (_t0 is None):  # or (_rho and _t0 is None)
        # No need for case handling when _t0=None and _rho!=None
        name_id = 'a: ' + a_str + 'n: ' + n_str

    return name_id


def plot_isomorph_contour(_rho0, _t0, _a0, _n):
    """
    Uses rho and t meshes from outer scope along with RHO, T 2D arrays to plot a wireframe contour
    :param _rho0: reference density
    :param _t0: reference temperature
    :param _a0: reference parameter A
    :param _n: reference potential strength parameter
    """
    a = np.array([iso_function(rho, t, rho0=_rho0, t0=_t0, a0=_a0, n=_n)
                  for rho, t in zip(np.ravel(RHO), np.ravel(T))])
    A = a.reshape(RHO.shape)
    ax.plot_wireframe(RHO, T, A, color='red', alpha=0.8, rstride=2, cstride=2)


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
