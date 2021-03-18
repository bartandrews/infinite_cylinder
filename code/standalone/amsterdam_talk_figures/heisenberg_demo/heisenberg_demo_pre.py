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
    outer_grid = gridspec.GridSpec(2, 1, hspace=0.4, wspace=0.4, height_ratios=[2, 1])
    top_cell = outer_grid[0, 0]
    bottom_cell = outer_grid[1, 0]
    upper_grid = gridspec.GridSpecFromSubplotSpec(2, 1, top_cell, hspace=0)

    # define a list of easily-visible markers
    markers = [(3, 0, 0), (4, 0, 0), (5, 0, 0), (6, 0, 0), (4, 1, 0), (5, 1, 0), (6, 1, 0),
               (3, 2, 0), (4, 2, 0), (5, 2, 0), (6, 2, 0), 'X', 'x', 'd', 'D', 'P',
               '$a$', '$b$', '$c$', '$d$', '$e$', '$f$', '$g$', '$h$', '$i$', '$j$', '$k$', '$l$', '$m$', '$n$', '$o$',
               '$p$', '$q$', '$r$', '$s$', '$t$', '$u$', '$v$', '$w$', '$x$', '$y$', '$z$']

    ####################################################################################################################

    ax1 = plt.subplot(upper_grid[0])

    # extract data from file
    with open('heisenberg_entropy.dat', 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        D = []
        S = []
        for row in plots:
            D.append(float(row[0]))
            S.append(float(row[1]))

    ax1.plot(D, S, '.', c='k')
    ax1.set_xlim([-2, 2])
    ax1.set_ylabel("$S_\mathrm{vN}$", fontsize=12)
    plt.setp(ax1.get_xticklabels(), visible=False)
    plt.tick_params(
        axis='x',  # changes apply to the x-axis
        which='both',  # both major and minor ticks are affected
        bottom=False,  # ticks along the bottom edge are off
        top=False,  # ticks along the top edge are off
        labelbottom=False)
    # ax1.axvspan(-2, -1/3, alpha=0.5, color='slateblue')
    # ax1.axvspan(-1/3, 1, alpha=0.5, color='red')
    # ax1.axvspan(1, 2, alpha=0.5, color='green')
    # ax1.axvline(-1.5, linestyle='--', c='k', linewidth=0.5)
    # ax1.axvline(0, linestyle='--', c='k', linewidth=0.5)
    # ax1.axvline(1.5, linestyle='--', c='k', linewidth=0.5)

    ####################################################################################################################

    ax2 = plt.subplot(upper_grid[1])

    # extract data from file
    with open('heisenberg_O.dat', 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        D = []
        O = []
        for row in plots:
            D.append(float(row[0]))
            O.append(float(row[1]))

    ax2.plot(D, O, '.', c='k')
    ax2.set_xlim([-2, 2])
    ax2.xaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax2.set_ylabel("$\mathcal{O}_\mathcal{I}$", fontsize=12)
    ax2.set_xlabel("$D/J$", fontsize=12)
    # ax2.axvspan(-2, -1/3, alpha=0.5, color='slateblue')
    # ax2.axvspan(-1/3, 1, alpha=0.5, color='red')
    # ax2.axvspan(1, 2, alpha=0.5, color='green')
    # ax2.text(-1.4, 0.3, "trivial", fontsize=12)
    # ax2.text(0.1, 0, "Haldane", fontsize=12)
    # ax2.text(1.3, 0, "large $D$", fontsize=12)
    # ax2.axvline(-1.5, linestyle='--', c='k', linewidth=0.5)
    # ax2.axvline(0, linestyle='--', c='k', linewidth=0.5)
    # ax2.axvline(1.5, linestyle='--', c='k', linewidth=0.5)

    ####################################################################################################################

    # ax3 = plt.subplot(bottom_cell)
    #
    # # extract data from file
    # with open('heisenberg_On.dat', 'r') as csvfile:
    #     plots = csv.reader(csvfile, delimiter='\t')
    #     n = []
    #     overlap1 = []
    #     overlap2 = []
    #     overlap3 = []
    #     for row in plots:
    #         n.append(float(row[0]))
    #         overlap1.append(float(row[1]))
    #         overlap2.append(float(row[2]))
    #         overlap3.append(float(row[3]))
    #
    # ax3.plot(n, overlap3, '.', label="$1.5$", color='green', marker=markers[2])
    # ax3.plot(n, overlap1, '.', label="$-1.5$", color='slateblue', marker=markers[0])
    # ax3.plot(n, overlap2, '.', label="$0$", color='red', marker=markers[1])
    # ax3.set_xlim([0, 30])
    # ax3.set_ylabel("$\\langle\Psi|\mathcal{R}_\mathrm{part}|\Psi\\rangle$", fontsize=12)
    # ax3.set_xlabel("sites in partial reflection", fontsize=12)
    #
    # # leg = ax3.legend(loc='center left', handletextpad=0.3, handlelength=1, labelspacing=1, borderpad=0.3,
    # #                  framealpha=1,
    # #                  edgecolor='k', markerscale=2, fontsize=10, ncol=1, columnspacing=0.5,
    # #                  bbox_to_anchor=(0, 0.65), title='$D/J$', title_fontsize=11)
    # # leg.get_frame().set_linewidth(0.5)
    # ax3.text(0.2, 0.8, "$D/J=1.5$", fontsize=10)
    # ax3.text(0.2, 0.2, "$D/J=-1.5$", fontsize=10)
    # ax3.text(0.2, -0.6, "$D/J=0$", fontsize=10)


    ####################################################################################################################

    plt.savefig("/home/bart/Documents/presentations/2021_03_18/figures/heisenberg_demo_pre.png", bbox_inches='tight', dpi=300)
    plt.show()
