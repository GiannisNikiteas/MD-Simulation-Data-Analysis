from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np

"""
Functional form of the isomorph generation.
The file generates a the contour for all the possible parameter a values, for a given T0 and a0.
"""

jet = plt.get_cmap('jet')


def iso_function(rho, t, rho0=1.0, t0=1.0, a0=0.5, n=8):
    a = a0 * (rho/rho0) ** (1./3.) * ((t0/t) ** (6/n)) ** (1.0/3.0)
    return a


# Arrays of X, Y planes
rho = np.arange(0.1, 5.0, 0.05)
t = np.arange(0.1, 5.0, 0.05)
RHO, T = np.meshgrid(rho, t)
# n_list = np.arange(6, 13, step=3, dtype=int)
n_list = [6, 12]


#############################################################################
def label_name(_a0, _n, _rho=None, _t0=None):
    rho_str = None
    t_str = None
    a_str = '{:.5f}'.format(_a0)
    n_str = str(_n)
    name_id = None
    if (_rho is None) or (_t0 is None):  # or (_rho and _t0 is None)
        # No need for case handling when _t0=None and _rho!=None
        name_id = 'a: ' + a_str + ' n: ' + n_str
    else:
        rho_str = '{:.4f}'.format(_rho)
        t_str = '{:.4f}'.format(_t0)
        name_id = r'$\rho$: ' + rho_str + 'T: ' + t_str + 'a: ' + a_str + 'n: ' + n_str

    return name_id


def plot_isomorph_contour(_rho0, _t0, _a0, _n, show_projections=False, transparency=1.0):
    """
    Uses rho and t meshes from outer scope along with RHO, T 2D arrays to plot a wireframe contour.
    Can be used intuitivly in loops to generate rho vs T vs a vs n, contours
    :param _rho0: reference density
    :param _t0: reference temperature
    :param _a0: reference parameter A
    :param _n: reference potential strength parameter
    :param show_projections: Shows X, Y, Z projections of the contours
    :param transparency: alpha parameter for the plots to adjust transparency
    """
    a = np.array([iso_function(rho, t, rho0=_rho0, t0=_t0, a0=_a0, n=_n)
                  for rho, t in zip(np.ravel(RHO), np.ravel(T))])
    A = a.reshape(RHO.shape)
    name = label_name(_a0, _n)
    w = ax.plot_wireframe(RHO, T, A, alpha=transparency,
                          rstride=4, cstride=4, label=name)
    # w = ax.surface_plot(RHO, T, A, alpha=transparency, rstride=4, cstride=4, label=name)

    # Projections of X, Y, Z onto corresponding planes in 3D plot
    if show_projections is True:
        # Plots projections for every single contour. Can get confusing when plotting against n as well
        cset = ax.contourf(RHO, T, A, zdir='z', offset=0,
                           cmap=plt.get_cmap('jet'), alpha=0.5)       # offset=np.amin(A),
        cset = ax.contourf(RHO, T, A, zdir='x', offset=0,
                           cmap=plt.get_cmap('coolwarm'), alpha=0.5)  # offset=np.amin(RHO)
        cset = ax.contourf(RHO, T, A, zdir='y', offset=5,
                           cmap=plt.get_cmap('coolwarm'), alpha=0.5)  # offset=np.amax(T),
    return w


if __name__ == "__main__":
    ############################################################################
    # PLOTTING RHO vs T vs A vs N with COLORBAR
    ############################################################################
    # Custom colormap
    color_map = jet(np.linspace(0, 1, len(n_list)))
    fig = plt.figure('3D contour plot, with function')
    ax = fig.gca(projection='3d')
    canvas = None

    for i in range(len(n_list)):
        canvas = plot_isomorph_contour(
            _rho0=0.5, _t0=0.05, _a0=0.5, _n=n_list[i], transparency=0.5)
        canvas.set_color(color_map[i])

    m = cm.ScalarMappable(cmap=jet, norm=canvas.norm)
    m.set_array(n_list)
    # fig.colorbar(m, cmap=cm.jet)
    ax.legend(loc='best', fancybox=True, fontsize='small')

    #############################################################################
    # PLOTTING RHO vs T vs A vs N with VARYING RHO
    #############################################################################
    fig = plt.figure(r'3D contour plot, with function var $\rho$')
    ax = fig.gca(projection='3d')

    rho = t = None
    rho = t = np.arange(0.01, 5, 0.05)
    RHO, T = np.meshgrid(rho, t)
    rho_ref = np.linspace(0.1, 5, 20)   # list of different reference densities
    n = 8
    g = None
    color_map_rho = jet(np.linspace(0, 1, len(rho_ref)))
    for i in range(len(rho_ref)):
        r_temp = rho_ref[i]
        g = plot_isomorph_contour(
            _rho0=r_temp, _t0=1.0, _a0=1, _n=8, transparency=0.5)
        g.set_color(color_map_rho[i])

    m = cm.ScalarMappable(cmap=jet, norm=g.norm)
    m.set_array(rho)
    fig.colorbar(m, cmap=jet)
    ax.legend(loc='best', fancybox=True, fontsize='small')

    #############################################################################
    # PLOTTING RHO vs T vs A vs N with DIFFERENT A
    #############################################################################
    color_map = jet(np.linspace(0, 1, len(n_list)))
    fig = plt.figure('3D contour plot, with function var A')
    ax = fig.gca(projection='3d')

    p = plot_isomorph_contour(0.5, 1, 0.5, 8, transparency=0.6)
    p.set_color('blue')
    p = plot_isomorph_contour(0.5, 1, 1.0, 8, transparency=0.6)
    p.set_color('red')
    ax.legend(loc='best', fancybox=True, fontsize='small')

    #############################################################################
    # LABELS, AXIS AND VISUALISATION IMPROVEMENTS
    #############################################################################
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

    fig.tight_layout()
    fig.show()
    plt.show()
