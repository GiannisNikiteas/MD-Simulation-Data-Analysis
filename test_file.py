import matplotlib.pyplot as plt
from isomorphs import *
import numpy as np


T_LINE = np.linspace(0.05, 5, 100)
rho_ref = np.arange(0.25, 2, 0.25)
n = 8
fig = plt.figure('3D line plot')
ax = fig.gca(projection='3d')


def contour_gen(rho_r, t_r, a_r, rho_list=None, t_list=None, a_list=None):
    arg_list = [rho_list, t_list, a_list]
    no_none_entry = [i for i, entry in enumerate(arg_list) if entry is not None]
    ref_list = arg_list[no_none_entry[0]]
    rho_iso = np.empty((0, len(T_LINE)))
    a_iso = np.empty((0, len(T_LINE)))
    t_iso = np.empty((0, len(T_LINE)))

    for i in ref_list:
        isomorph_line = Isomorph(i, t_r, a_r, T_LINE)
        temp_rho, temp_a = isomorph_line.gen_line(n)
        rho = np.append(rho_iso, [temp_rho], axis=0)
        a = np.append(a_iso, [temp_a], axis=0)
        t = np.append(t_iso, [T_LINE], axis=0)
        return rho, t, a


rho_iso, t_iso, a_iso = contour_gen(0.5, 0.5, 0.5, rho_list=rho_ref)
print(np.shape(a_iso))


ax.set_ylabel(r'$\rho$')
ax.set_xlabel(r'T')
ax.set_zlabel(r'a')

# Grid background color set to white
ax.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
ax.w_yaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
ax.w_zaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))

plt.tight_layout()
plt.show()

