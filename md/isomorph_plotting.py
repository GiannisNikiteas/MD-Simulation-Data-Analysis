import os
import numpy as np
import pickle as pl
import matplotlib.pyplot as plt
from md import Isomorph
from mpl_toolkits.mplot3d import Axes3D


"""
PARAMETER A:

1. A0 affects the curve of the isomorphic surface, with a0 = 0, being a flat surface and
a0 = 1.0 being more curved as T becomes smaller.

2. It also affects, the position on the z-axis (a), with higher values of A0 shifting the
surface higher up. When A0 = 0, then the figure is collapsed on the rho-T axis

3. The surfaces have the exact same footprint on the X-Y axis


    rho2 increases as we go down the isomorphic line
    T2 also increases
    A2 decreases

    if a_r == 0; then all the a2s are forced to be zero (0) hence causing the surface to be flat
"""


def isomorphic_surface_array(rho_list, t_list, n_list, a_list, t2, figname, save_fig=False, save_pickle=False):
    """

    @param rho_list:
    @type rho_list: ndarray
    @param t_list:
    @type t_list: ndarray
    @param n_list:
    @type n_list: ndarray
    @param a_list:
    @type a_list: ndarray
    @param t2:
    @type t2: ndarray
    @param figname:
    @type figname: basestring
    @param save_fig:
    @type save_fig: bool
    @param save_pickle:
    @type save_pickle: bool
    @return:
    @rtype: ndarray[ndarray]
    """

    # Generate the 3D plot canvas here
    fig = plt.figure(figname)
    ax = Axes3D(fig)
    canvas = None

    edg_col = ["none", "black", "cyan", "white"]
    edg_it = 0
    for a_r in a_list:
        for t_r in t_list:
            # Varying the reference density in the line
            # Merging multiple rho lines creates a surface
            for n in n_list:
                rho_iso = np.empty((0, len(t2)))
                t_iso = np.empty((0, len(t2)))
                a_iso = np.empty((0, len(t2)))

                label_title = fr"n: {n}, $T_0$: {t_r:.1f}, $A_0$: {a_r:.1f}"
                for i, rho_r in enumerate(rho_list):
                    iso = Isomorph(rho_r, t_r, a_r, t2)
                    # Generate isomorphic line
                    rho2, a2 = iso.gen_line(n)

                    # Creating the 3D mesh grid
                    rho_iso = np.append(rho_iso, [rho2], axis=0)
                    t_iso = np.append(t_iso, [t2], axis=0)
                    a_iso = np.append(a_iso, [a2], axis=0)

                surf = ax.plot_surface(rho_iso, t_iso, a_iso, label=label_title,
                                       alpha=0.9, edgecolor=edg_col[edg_it])

                # This is a fix for the bug that does not allow legends to 3D surfaces
                surf._facecolors2d = surf._facecolors3d
                surf._edgecolors2d = surf._edgecolors3d
                ax.legend(loc="best")

        edg_it += 1

    ax.set_xlabel(r'$\rho$')
    ax.set_ylabel(r'T')
    ax.set_zlabel(r'a')
    # ax.view_init(elev=44, azim=-128.5)

    if save_fig is True:
        plt.savefig(f"{figname}.pdf")

    # Save figures for loading at a later time
    if save_pickle is True:
        pl.dump(fig, open(f"{figname}.pickle", "wb"))


def plot_all_surfaces():
    """
    Generates all the plots of the ismorphic surfaces for the MD fluid.
    It also pickles the figures, so that they can be loaded at a later time.
    """
    # Leave these unchanged
    # Density reference
    rho_list = np.linspace(0.2, 1, 10)
    # Input temperature over which the isomorphic points are generated
    t2 = np.linspace(0.2, 2, 20)

    # PLOT VARIOUS TEMPERATURES T0
    # Reference parameters
    t_list = np.linspace(0.2, 1, 5)
    a_list = [0.5]
    n_list = [8]
    isomorphic_surface_array(rho_list, t_list, n_list,
                             a_list, t2, "Isomorphs_with_varying_T0")

    # PLOT VARIOUS PARAMETERS A0
    a_list = np.linspace(0.5, 1, 3)
    t_list = [0.5]
    n_list = [8]
    isomorphic_surface_array(rho_list, t_list, n_list,
                             a_list, t2, "Isomorphs_with_varying_A0")

    # PLOT VARIOUS PAIR POTENTIAL STRENGTHS n
    # TODO: increase color intensity cyan, light green, DISABLE SHADING, wireframe color and decrease thickness
    a_list = [0.5]
    t_list = [0.5]
    n_list = list(range(8, 15, 2))
    isomorphic_surface_array(rho_list, t_list, n_list,
                             a_list, t2, "Isomorphs_with_varying_n")

    # PLOT VARIOUS A0 AND T0
    t_list = np.linspace(0.2, 1, 3)
    a_list = np.linspace(0, 1, 3)
    n_list = [8]
    isomorphic_surface_array(rho_list, t_list, n_list,
                             a_list, t2, "Isomorphs_with_varying_A0_and_T0")

    # PLOT VARYING A0 and n
    t_list = [0.5]
    a_list = np.linspace(0, 1, 3)
    n_list = list(range(8, 15, 2))
    isomorphic_surface_array(rho_list, t_list, n_list,
                             a_list, t2, "Isomorphs_with_varying_A0_and_n")

    # PLOT VARYING T0 and n
    t_list = np.linspace(0.2, 1, 3)
    a_list = [0.5]
    n_list = list(range(8, 13, 2))
    isomorphic_surface_array(rho_list, t_list, n_list,
                             a_list, t2, "Isomorphs_with_varying_T0_and_n")

    # PLOT WITH TEMPERATURE TO INFINITY
    t_list = np.linspace(0.2, 10, 15)
    a_list = [0.5]
    n_list = [8]
    isomorphic_surface_array(rho_list, t_list, n_list,
                             a_list, t2, "Isomorphs_T0_to_infinity")


# Not working with 3DAxes
def load_figures(fig_names):
    """
    Uses a list of the figure names to load them into a list
    @param fig_names:
    @type fig_names:
    @return: A list containing all the figures
    @rtype: list
    """
    fig_list = []
    for i, name in enumerate(fig_names):
        fig_list.append(pl.load(open(f"{name}.pickle", "rb")))
    return fig_list


if __name__ == "__main__":
    import os
    import matplotlib.pyplot as plt

    # Change directory to where files will be saved and loaded from
    os.chdir("/home/gn/Desktop/surface_figures")
    # List that contains a copy of the figure names, 
    # that will be loaded into the load_pickle function
    # to resume the figures to their normal state without having to re-plot them
    fig_names = ["Isomorphs_with_varying_T0",
                 "Isomorphs_with_varying_A0",
                 "Isomorphs_with_varying_n",
                 "Isomorphs_with_varying_A0_and_T0",
                 "Isomorphs_with_varying_A0_and_n",
                 "Isomorphs_with_varying_T0_and_n",
                 "Isomorphs_T0_to_infinity"]

    plot_all_surfaces()
    plt.show()
