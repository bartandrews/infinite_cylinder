import numpy as np
import matplotlib
import csv
import matplotlib.pyplot as plt
import sys
from itertools import product
import matplotlib.gridspec as gridspec
from matplotlib.patches import Polygon
import matplotlib.patches as patches
from numpy.polynomial.polynomial import polyfit
import os
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.ticker as ticker
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import ConnectionPatch

plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}')
# matplotlib.verbose.level = 'debug-annoying'

if __name__ == '__main__':

    fig = plt.figure(figsize=(6, 4))
    gs = gridspec.GridSpec(2, 2, hspace=0, wspace=0.5)

    name = "pi_by_4"
    model = "FerBBH2"
    nu = (1, 2)

    ####################################################################################################################

    ax1 = plt.subplot(gs[0])

    plot_dir = f'/home/bart/PycharmProjects/infinite_cylinder/data/O_I_n_flow/{model}/diag_initial_state'
    f90_1 = f'O_I_n_flow_{model}_chi_100_n_{nu[0]}_{nu[1]}_LxMUC_6_Ly_6_custom.dat.90.{name}_part_1'
    f90_2 = f'O_I_n_flow_{model}_chi_100_n_{nu[0]}_{nu[1]}_LxMUC_6_Ly_6_custom.dat.90.{name}_part_2'
    p90_1 = os.path.join(plot_dir, f90_1)
    p90_2 = os.path.join(plot_dir, f90_2)

    t = []
    mag = []
    ang = []

    for i in [p90_1, p90_2]:
        with open(i, 'r') as csvfile:
            plots = csv.reader(csvfile, delimiter='\t')
            for row in plots:
                t.append(float(row[0]))
                mag.append(float(row[2]))
                ang.append(float(row[3]))

    ax1.plot(t, mag, 's', marker='x', color='k', markersize=3)
    # ax1.set_xlabel("$t_{3,4}/t_{1,2}$", fontsize=11)
    ax1.set_ylabel("$|z_1|$", fontsize=11)
    ax1.set_xlim([0, 4])
    ax1.axvspan(0, 1, alpha=0.5, color='slateblue')
    ax1.axvspan(1, 4, alpha=0.5, color='gold')

    ax2 = plt.subplot(gs[2], sharex=ax1)
    ax1.tick_params('x', direction='in', bottom=True)
    ax2.tick_params('x', direction='in', top=True)
    plt.setp(ax1.get_xticklabels(), visible=False)
    ax2.plot(t, ang, 's', marker='x', color='k', markersize=3)
    ax2.set_ylabel("$\\angle z_1 / (\pi/4)$", fontsize=11)
    ax2.set_xlabel("$t_{3,4}/t_{1,2}$", fontsize=11)
    ax2.set_xlim([0, 4])
    ax2.axvspan(0, 1, alpha=0.5, color='slateblue')
    ax2.axvspan(1, 4, alpha=0.5, color='gold')

    ####################################################################################################################

    ax3 = plt.subplot(gs[1])

    plot_dir = f'/home/bart/PycharmProjects/infinite_cylinder/data/O_I_n_flow/{model}/diag_initial_state'
    f180_1 = f'O_I_n_flow_{model}_chi_100_n_{nu[0]}_{nu[1]}_LxMUC_6_Ly_6_custom.dat.180.{name}_part_1'
    f180_2 = f'O_I_n_flow_{model}_chi_100_n_{nu[0]}_{nu[1]}_LxMUC_6_Ly_6_custom.dat.180.{name}_part_2'
    p180_1 = os.path.join(plot_dir, f180_1)
    p180_2 = os.path.join(plot_dir, f180_2)

    t = []
    mag = []
    ang = []

    for i in [p180_1, p180_2]:
        with open(i, 'r') as csvfile:
            plots = csv.reader(csvfile, delimiter='\t')
            for row in plots:
                t.append(float(row[0]))
                mag.append(float(row[2]))
                ang.append(float(row[3]))

    ax3.plot(t, mag, 's', marker='x', color='k', markersize=3)
    ax3.set_ylabel("$|z_2|$", fontsize=11)
    # ax3.set_xlabel("$t_{3,4}/t_{1,2}$", fontsize=11)
    ax3.set_xlim([0, 4])
    ax3.axvspan(0, 1, alpha=0.5, color='slateblue')
    ax3.axvspan(1, 4, alpha=0.5, color='gold')

    ax4 = plt.subplot(gs[3], sharex=ax3)
    ax3.tick_params('x', direction='in', bottom=True)
    ax4.tick_params('x', direction='in', top=True)
    plt.setp(ax3.get_xticklabels(), visible=False)
    ax4.plot(t, ang, 's', marker='x', color='k', markersize=3)
    ax4.set_ylabel("$\\angle z_2 / (\pi/4)$", fontsize=11)
    ax4.set_xlabel("$t_{3,4}/t_{1,2}$", fontsize=11)
    ax4.set_xlim([0, 4])
    ax4.axvspan(0, 1, alpha=0.5, color='slateblue')
    ax4.axvspan(1, 4, alpha=0.5, color='gold')

    ####################################################################################################################

    # fig.text(0.04, 0.9, "(a)", fontsize=12)
    # fig.text(0.5, 0.9, "(b)", fontsize=12)
    # fig.text(0.04, 0.45, "(c)", fontsize=12)
    # fig.text(0.5, 0.45, "(d)", fontsize=12)

    plt.savefig(f"/home/bart/Documents/papers/BBH/figures/BBH2/{name}.png", bbox_inches='tight', dpi=300)
    plt.show()
