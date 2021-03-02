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
from scipy import stats
from fractions import Fraction


plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}')
# matplotlib.verbose.level = 'debug-annoying'

if __name__ == '__main__':

    fig = plt.figure()
    outer_grid = gridspec.GridSpec(2, 2, hspace=0.5, wspace=0.4)
    top_left_cell = outer_grid[0, 0]
    top_right_cell = outer_grid[0, 1]
    bottom_left_cell = outer_grid[1, 0]
    bottom_right_cell = outer_grid[1, 1]
    top_left_grid = gridspec.GridSpecFromSubplotSpec(2, 1, top_left_cell, hspace=0.1)
    top_right_grid = gridspec.GridSpecFromSubplotSpec(2, 1, top_right_cell, hspace=0.1)
    bottom_left_grid = gridspec.GridSpecFromSubplotSpec(2, 1, bottom_left_cell, hspace=0.1)
    bottom_right_grid = gridspec.GridSpecFromSubplotSpec(2, 1, bottom_right_cell, hspace=0.1)

    # define a list of easily-visible markers
    markers = [(3, 0, 0), (4, 0, 0), (5, 0, 0), (6, 0, 0), (4, 1, 0), (5, 1, 0), (6, 1, 0),
               (3, 2, 0), (4, 2, 0), (5, 2, 0), (6, 2, 0), 'X', 'x', 'd', 'D', 'P',
               '$a$', '$b$', '$c$', '$d$', '$e$', '$f$', '$g$', '$h$', '$i$', '$j$', '$k$', '$l$', '$m$', '$n$', '$o$',
               '$p$', '$q$', '$r$', '$s$', '$t$', '$u$', '$v$', '$w$', '$x$', '$y$', '$z$']

    ####################################################################################################################

    for filling in range(8):
        if filling == 0:
            ax1 = plt.subplot(top_left_grid[0])
            nphi = (1, 3)
            Ly_val = 7
        elif filling == 1:
            ax1 = plt.subplot(top_right_grid[0])
            nphi = (1, 4)
            Ly_val = 7
        elif filling == 2:
            ax1 = plt.subplot(top_left_grid[1])
            nphi = (1, 3)
            Ly_val = 7
        elif filling == 3:
            ax1 = plt.subplot(top_right_grid[1])
            nphi = (1, 4)
            Ly_val = 7
        elif filling == 4:
            ax1 = plt.subplot(bottom_left_grid[0])
            nphi = (1, 5)
            Ly_val = 7
        elif filling == 5:
            ax1 = plt.subplot(bottom_right_grid[0])
            nphi = (1, 6)
            Ly_val = 7
        elif filling == 6:
            ax1 = plt.subplot(bottom_left_grid[1])
            nphi = (1, 5)
            Ly_val = 7
        elif filling == 7:
            ax1 = plt.subplot(bottom_right_grid[1])
            nphi = (1, 6)
            Ly_val = 7

        for chi_num, chi in enumerate([50, 100, 150, 200, 250]):

            nu = (3, 7)
            n = Fraction(nu[0]*nphi[0], nu[1]*nphi[1]).limit_denominator(100)

            phiflow_dir = '/home/bart/PycharmProjects/infinite_cylinder/data/charge_pump/FerHofSqu1'
            phiflow_file = f'charge_pump_FerHofSqu1_chi_{chi}_t1_1_V_10_Coulomb_1_' \
                         f'n_{n.numerator}_{n.denominator}_nphi_{nphi[0]}_{nphi[1]}_LxMUC_1_Ly_{Ly_val}_phi_0_{nu[1]}_{nu[1]*10+1}.dat'
            phiflow_path = os.path.join(phiflow_dir, phiflow_file)

            # extract data from file
            with open(phiflow_path, 'r') as csvfile:
                plots = csv.reader(csvfile, delimiter='\t')
                phi = []
                charge = []
                for row in plots:
                    phi.append(float(row[0]))
                    charge.append(float(row[1]))

            charge = [i - charge[0] for i in charge]

            ax1.plot(phi, charge, '.', c=f'C{chi_num}', marker=markers[chi_num], fillstyle='none', markersize=2.5, label=f'${chi}$')

            if filling == 0:
                leg = ax1.legend(loc='upper center', handletextpad=0.3, handlelength=1, labelspacing=0.1, borderpad=0.3,
                                 framealpha=1,
                                 edgecolor='k', markerscale=2, fontsize=10, ncol=6, columnspacing=0.5,
                                 bbox_to_anchor=(1.2, 2.2), title='$\chi$', title_fontsize=11)
                leg.get_frame().set_linewidth(0.5)

            ax1.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
            ax1.set_xlim([0, nu[1]])
            ax1.set_xticks(np.arange(0, nu[1]+0.1, 1))
            ax1.set_yticks([-nu[0], 0])
            # ax1.set_ylim([-nu[0], 0])
            ax1.set_ylabel("$\\langle Q_\mathrm{L} \\rangle$", fontsize=11)
            if chi == 50:
                ax1.text(0.05*nu[1], -0.9*nu[0], f"$L_y={Ly_val}$", fontsize=11)
            if filling not in [2, 3, 6, 7]:
                plt.setp(ax1.get_xticklabels(), visible=False)
                ax1.set_title(f"$n_\phi={nphi[0]}/{nphi[1]}$")
                plt.tick_params(
                    axis='x',  # changes apply to the x-axis
                    which='both',  # both major and minor ticks are affected
                    bottom=False,  # ticks along the bottom edge are off
                    top=False,  # ticks along the top edge are off
                    labelbottom=False)
            else:
                ax1.set_xlabel("$\\Phi_x / 2\pi$", fontsize=11)

    ####################################################################################################################
    ####################################################################################################################

    fig.text(0.04, 0.875, "(a)", fontsize=12)
    fig.text(0.49, 0.875, "(b)", fontsize=12)
    fig.text(0.04, 0.41, "(c)", fontsize=12)
    fig.text(0.49, 0.41, "(d)", fontsize=12)

    plt.savefig("/home/bart/Documents/papers/FCI/bond_scaling.png", bbox_inches='tight', dpi=300)
    plt.show()
