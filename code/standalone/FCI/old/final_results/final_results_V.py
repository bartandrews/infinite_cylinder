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


plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}')
# matplotlib.verbose.level = 'debug-annoying'

if __name__ == '__main__':

    fig = plt.figure(figsize=(6, 3))
    gs = gridspec.GridSpec(2, 2, hspace=0, wspace=0.4)

    # define a list of easily-visible markers
    markers = [(3, 0, 0), (4, 0, 0), (5, 0, 0), (6, 0, 0), (4, 1, 0), (5, 1, 0), (6, 1, 0),
               (3, 2, 0), (4, 2, 0), (5, 2, 0), (6, 2, 0), 'X', 'x', 'd', 'D', 'P',
               '$a$', '$b$', '$c$', '$d$', '$e$', '$f$', '$g$', '$h$', '$i$', '$j$', '$k$', '$l$', '$m$', '$n$', '$o$',
               '$p$', '$q$', '$r$', '$s$', '$t$', '$u$', '$v$', '$w$', '$x$', '$y$', '$z$']

    ####################################################################################################################

    ax1 = plt.subplot(gs[0], anchor=(0, 0.85))
    ax1.tick_params('x', direction='in', bottom=True)

    TFB_dir = '/standalone/FCI/final_results'
    TFB_file = 'TFB.dat'
    TFB_path = os.path.join(TFB_dir, TFB_file)

    # extract data from file
    with open(TFB_path, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        C = []
        Vcrit = [[] for i in range(6)]
        c = [[] for i in range(6)]
        for row in plots:
            C.append(int(row[0]))
            for i in range(6):
                Vcrit[i].append(float(row[2*i+1]))
                c[i].append(float(row[2*i+2]))

    lbl = [1, -1, 2, -2, 3, -3]
    for i in range(6):
        ax1.plot(C, Vcrit[i], marker=markers[i], color=f'C{i}', markersize=5, label=f'$r={lbl[i]}$')

    leg = ax1.legend(loc='upper center', handletextpad=0.3, handlelength=1, labelspacing=0.1, borderpad=0.3, framealpha=1,
                     edgecolor='k', markerscale=1, fontsize=10, ncol=3, columnspacing=0.5, bbox_to_anchor=(0.5, 1.44))
    leg.get_frame().set_linewidth(0.5)

    ax1.tick_params(axis="x", labelsize=10)
    ax1.set_xlim([1, 5])
    ax1.set_xticks(np.arange(1, 5.1, 1))
    ax1.tick_params(axis="y", labelsize=10)
    ax1.set_ylabel("$C$", fontsize=11)
    ax1.set_ylabel("$V_\\text{crit}$", fontsize=11)

    ####################################################################################################################

    ax2 = plt.subplot(gs[1])
    ax2.tick_params('x', direction='in', bottom=True)

    Hofstadter_dir = '/standalone/FCI/final_results'
    Hofstadter_file = 'Hofstadter.dat'
    Hofstadter_path = os.path.join(Hofstadter_dir, Hofstadter_file)

    # extract data from file
    with open(Hofstadter_path, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        C = []
        Vcrit = [[] for i in range(6)]
        c = [[] for i in range(6)]
        for row in plots:
            C.append(int(row[0]))
            for i in range(6):
                Vcrit[i].append(float(row[2 * i + 1]))
                c[i].append(float(row[2 * i + 2]))

    lbl = [1, -1, 2, -2, 3, -3]
    for i in range(6):
        ax2.plot(C, Vcrit[i], marker=markers[i], color=f'C{i}', markersize=5, label=f'$r={lbl[i]}$')

    leg = ax2.legend(loc='upper center', handletextpad=0.3, handlelength=1, labelspacing=0.1, borderpad=0.3,
                     framealpha=1,
                     edgecolor='k', markerscale=1, fontsize=10, ncol=3, columnspacing=0.5, bbox_to_anchor=(0.5, 1.44))
    leg.get_frame().set_linewidth(0.5)

    ax2.tick_params(axis="x", labelsize=10)
    ax2.set_xlim([1, 5])
    ax2.set_xticks(np.arange(1, 5.1, 1))
    ax2.tick_params(axis="y", labelsize=10)
    ax2.set_xlabel("$C$", fontsize=11)
    ax2.set_ylabel("$V_\\text{crit}$", fontsize=11)

    ####################################################################################################################

    ax3 = plt.subplot(gs[2])

    for i in range(6):
        ax3.plot(C, c[i], marker=markers[i], color=f'C{i}', markersize=5)

    ax3.tick_params(axis="x", labelsize=10)
    ax3.set_xlim([1, 5])
    ax3.set_xticks(np.arange(1, 5.1, 1))
    ax3.tick_params(axis="y", labelsize=10)
    ax3.set_xlabel("$C$", fontsize=11)
    ax3.set_ylabel("$c$", fontsize=11)

    ####################################################################################################################

    ax4 = plt.subplot(gs[3])

    for i in range(6):
        ax4.plot(C, c[i], marker=markers[i], color=f'C{i}', markersize=5)

    ax4.tick_params(axis="x", labelsize=10)
    ax4.set_xlim([1, 5])
    ax4.set_xticks(np.arange(1, 5.1, 1))
    ax4.tick_params(axis="y", labelsize=10)
    ax4.set_xlabel("$C$", fontsize=11)
    ax4.set_ylabel("$c$", fontsize=11)

    ####################################################################################################################
    ####################################################################################################################

    # fig.text(0.04, 0.9, "(a)", fontsize=12)
    # fig.text(0.5, 0.9, "(b)", fontsize=12)
    # fig.text(0.04, 0.45, "(c)", fontsize=12)
    # fig.text(0.5, 0.45, "(d)", fontsize=12)

    plt.setp(ax1.get_xticklabels(), visible=False)
    plt.setp(ax2.get_xticklabels(), visible=False)

    fig.text(0.22, 1.07, 'TFB model', fontsize=11)
    fig.text(0.65, 1.07, 'Hofstadter model', fontsize=11)

    for i in range(2, 5, 1):
        ax1.axvline(i, color='k', linewidth=0.5, ls='--')
        ax2.axvline(i, color='k', linewidth=0.5, ls='--')
        ax3.axvline(i, color='k', linewidth=0.5, ls='--')
        ax4.axvline(i, color='k', linewidth=0.5, ls='--')

    plt.savefig("/home/bart/Documents/papers/BT/final_results_V.png", bbox_inches='tight', dpi=300)
    plt.show()
