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

    fig = plt.figure(figsize=(12, 12))
    gs = gridspec.GridSpec(6, 3, hspace=0.5)

    # define a list of easily-visible markers
    markers = [(3, 0, 0), (4, 0, 0), (5, 0, 0), (6, 0, 0), (4, 1, 0), (5, 1, 0), (6, 1, 0),
               (3, 2, 0), (4, 2, 0), (5, 2, 0), (6, 2, 0), 'X', 'x', 'd', 'D', 'P',
               '$a$', '$b$', '$c$', '$d$', '$e$', '$f$', '$g$', '$h$', '$i$', '$j$', '$k$', '$l$', '$m$', '$n$', '$o$',
               '$p$', '$q$', '$r$', '$s$', '$t$', '$u$', '$v$', '$w$', '$x$', '$y$', '$z$']

    ####################################################################################################################

    for filling in range(6):
        if filling == 0:
            ax1 = plt.subplot(gs[0])
            nu = (4, 7)
            Ly_val = 14
        elif filling == 1:
            ax1 = plt.subplot(gs[3])
            nu = (3, 5)
            Ly_val = 10
        elif filling == 2:
            ax1 = plt.subplot(gs[6])
            nu = (2, 3)
            Ly_val = 6
        elif filling == 3:
            ax1 = plt.subplot(gs[9])
            nu = (1, 3)
            Ly_val = 6
        elif filling == 4:
            ax1 = plt.subplot(gs[12])
            nu = (2, 5)
            Ly_val = 10
        elif filling == 5:
            ax1 = plt.subplot(gs[15])
            nu = (3, 7)
            Ly_val = 14

        for q in range(3, 9):

            nphi = (1, q)
            n = Fraction(nu[0]*nphi[0], nu[1]*nphi[1]).limit_denominator(100)

            Vflow_dir = '/home/bart/PycharmProjects/infinite_cylinder/data/corr_len_V_flow/FerHofSqu1'
            Vflow_file = f'corr_len_V_flow_FerHofSqu1_chi_250_t1_1_V_0_10_41_Coulomb_1_' \
                         f'n_{n.numerator}_{n.denominator}_nphi_{nphi[0]}_{nphi[1]}_LxMUC_1_Ly_{Ly_val}.dat'
            Vflow_path = os.path.join(Vflow_dir, Vflow_file)

            # extract data from file
            with open(Vflow_path, 'r') as csvfile:
                plots = csv.reader(csvfile, delimiter='\t')
                V = []
                xi = []
                for row in plots:
                    V.append(float(row[0]))
                    xi.append(float(row[1]))

            ax1.plot(V, xi, '.', c=f'C{q}', markersize=5, label=f'${nphi[0]}/{nphi[1]}$')

            if filling == 0:
                leg = ax1.legend(loc='upper center', handletextpad=0.3, handlelength=1, labelspacing=0.1, borderpad=0.3,
                                 framealpha=1,
                                 edgecolor='k', markerscale=1, fontsize=10, ncol=6, columnspacing=0.5,
                                 bbox_to_anchor=(1.7, 1.7), title='$n_\\phi$')
                leg.get_frame().set_linewidth(0.5)

            ax1.set_xlim([0, 10])
            ax1.set_xlabel("$V$", fontsize=11)
            ax1.set_ylabel("$\\xi$", fontsize=11)

    ####################################################################################################################

    for filling in range(6):
        if filling == 0:
            ax2 = plt.subplot(gs[1])
            nu = (4, 7)
            Ly_val = 14
        elif filling == 1:
            ax2 = plt.subplot(gs[4])
            nu = (3, 5)
            Ly_val = 10
        elif filling == 2:
            ax2 = plt.subplot(gs[7])
            nu = (2, 3)
            Ly_val = 6
        elif filling == 3:
            ax2 = plt.subplot(gs[10])
            nu = (1, 3)
            Ly_val = 6
        elif filling == 4:
            ax2 = plt.subplot(gs[13])
            nu = (2, 5)
            Ly_val = 10
        elif filling == 5:
            ax2 = plt.subplot(gs[16])
            nu = (3, 7)
            Ly_val = 14

        for q in range(3, 9):

            nphi = (1, q)
            n = Fraction(nu[0] * nphi[0], nu[1] * nphi[1]).limit_denominator(100)

            Vflow_dir = '/home/bart/PycharmProjects/infinite_cylinder/data/corr_len_V_flow/FerHofSqu1'
            Vflow_file = f'corr_len_V_flow_FerHofSqu1_chi_250_t1_1_V_0_1_41_Coulomb_1_' \
                         f'n_{n.numerator}_{n.denominator}_nphi_{nphi[0]}_{nphi[1]}_LxMUC_1_Ly_{Ly_val}.dat'
            Vflow_path = os.path.join(Vflow_dir, Vflow_file)

            # extract data from file
            with open(Vflow_path, 'r') as csvfile:
                plots = csv.reader(csvfile, delimiter='\t')
                V = []
                xi = []
                for row in plots:
                    V.append(float(row[0]))
                    xi.append(float(row[1]))

            ax2.plot(V, xi, '.', c=f'C{q}', markersize=5, label=f'${nphi[0]}/{nphi[1]}$')

            ax2.set_xlim([0, 1])
            ax2.set_xlabel("$V$", fontsize=11)
            ax2.set_ylabel("$\\xi$", fontsize=11)

    ####################################################################################################################

    for filling in range(6):
        if filling == 0:
            ax3 = plt.subplot(gs[2])
            nu = (4, 7)
            Ly_val = 14
        elif filling == 1:
            ax3 = plt.subplot(gs[5])
            nu = (3, 5)
            Ly_val = 10
        elif filling == 2:
            ax3 = plt.subplot(gs[8])
            nu = (2, 3)
            Ly_val = 6
        elif filling == 3:
            ax3 = plt.subplot(gs[11])
            nu = (1, 3)
            Ly_val = 6
        elif filling == 4:
            ax3 = plt.subplot(gs[14])
            nu = (2, 5)
            Ly_val = 10
        elif filling == 5:
            ax3 = plt.subplot(gs[17])
            nu = (3, 7)
            Ly_val = 14

        for q in range(3, 9):

            nphi = (1, q)
            n = Fraction(nu[0] * nphi[0], nu[1] * nphi[1]).limit_denominator(100)

            Vflow_dir = '/home/bart/PycharmProjects/infinite_cylinder/data/corr_len_V_flow/FerHofSqu1'
            Vflow_file = f'corr_len_V_flow_FerHofSqu1_chi_250_t1_1_V_0_10_41_Coulomb_1_' \
                         f'n_{n.numerator}_{n.denominator}_nphi_{nphi[0]}_{nphi[1]}_LxMUC_1_Ly_{Ly_val}.dat'
            Vflow_path = os.path.join(Vflow_dir, Vflow_file)

            # extract data from file
            with open(Vflow_path, 'r') as csvfile:
                plots = csv.reader(csvfile, delimiter='\t')
                V = []
                xi = []
                for row in plots:
                    V.append(float(row[0]))
                    xi.append(float(row[1]))

            ax3.plot(V, xi, '.', c=f'C{q}', markersize=5, label=f'${nphi[0]}/{nphi[1]}$')

            ax3.set_xlim([0, 0.5])
            ax3.set_xlabel("$V$", fontsize=11)
            ax3.set_ylabel("$\\xi$", fontsize=11)

    ####################################################################################################################
    ####################################################################################################################

    fig.text(0, 0.85, "(a) $\\nu=4/7$", fontsize=12)
    fig.text(0, 0.716, "(b) $\\nu=3/5$", fontsize=12)
    fig.text(0, 0.582, "(c) $\\nu=2/3$", fontsize=12)
    fig.text(0, 0.448, "(d) $\\nu=1/3$", fontsize=12)
    fig.text(0, 0.314, "(e) $\\nu=2/5$", fontsize=12)
    fig.text(0, 0.18, "(f) $\\nu=3/7$", fontsize=12)

    plt.savefig("/home/bart/Documents/papers/FCI/Vflow_analysis_complete.png", bbox_inches='tight', dpi=300)
    plt.show()
