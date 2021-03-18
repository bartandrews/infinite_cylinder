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
plt.rc('text.latex', preamble=r'\usepackage{amsmath}\usepackage{amssymb}')
# matplotlib.verbose.level = 'debug-annoying'

if __name__ == '__main__':

    fig = plt.figure(figsize=(10, 4))
    gs = gridspec.GridSpec(2, 2, hspace=0, wspace=0.5)

    ####################################################################################################################

    ax1 = plt.subplot(gs[0])

    t = []
    mag = []
    ang = []

    with open('90.dat', 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        for row in plots:
            t.append(float(row[0]))
            mag.append(float(row[1]))
            ang.append(float(row[2]))

    ax1.plot(t, mag, 's', marker='x', color='k', markersize=3)
    # ax1.set_xlabel("$t_{3,4}/t_{1,2}$", fontsize=11)
    ax1.set_ylabel("magnitude", fontsize=11)
    ax1.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax1.set_xlim([0, 4])
    ax1.axvspan(0, 1, alpha=0.5, color='slateblue')
    ax1.axvspan(1, 4, alpha=0.5, color='gold')

    ax2 = plt.subplot(gs[2], sharex=ax1)
    ax1.tick_params('x', direction='in', bottom=True)
    ax2.tick_params('x', direction='in', top=True)
    plt.setp(ax1.get_xticklabels(), visible=False)
    ax2.plot(t, ang, 's', marker='x', color='k', markersize=3)
    ax2.set_ylabel("phase", fontsize=11)
    ax2.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax2.set_xlabel("$\gamma/\lambda$", fontsize=11)
    ax2.set_xlim([0, 4])
    ax2.axvspan(0, 1, alpha=0.5, color='slateblue')
    ax2.axvspan(1, 4, alpha=0.5, color='gold')
    ax1.set_title("$C^4_4=+1$ (SPT with $\mathbb{Z}_4$ classification?)")
    ax1.text(0.25, 0.1, "trivial", fontsize=12)
    ax1.text(2, 0.6, "topological", fontsize=12)
    ax2.text(0.1, 1, "$\langle \\tilde{C}_4 \\rangle \leq 1$", fontsize=12)
    ax2.text(1.9, 0.6, "$\langle \\tilde{C}_4 \\rangle = \\frac{1}{4}e^{\mathrm{i}\pi/2}$",
             fontsize=12)
    labels = [item.get_text() for item in ax2.get_yticklabels()]
    print(labels)
    labels[1] = '$0$'
    labels[2] = '$\pi/8$'
    labels[3] = '$\pi/4$'
    labels[4] = '$3\pi/8$'
    labels[5] = '$\pi/2$'
    print(labels)
    ax2.set_yticklabels(labels)

    ####################################################################################################################

    ax3 = plt.subplot(gs[1])

    t = []
    mag = []
    ang = []

    with open('180.dat', 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        for row in plots:
            t.append(float(row[0]))
            mag.append(float(row[1]))
            ang.append(float(row[2]))

    ax3.plot(t, mag, 's', marker='x', color='k', markersize=3)
    ax3.set_ylabel("magnitude", fontsize=11)
    ax3.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    # ax3.set_xlabel("$t_{3,4}/t_{1,2}$", fontsize=11)
    ax3.set_xlim([0, 4])
    # ax3.axvspan(0, 1, alpha=0.5, color='slateblue')
    # ax3.axvspan(1, 4, alpha=0.5, color='gold')

    ax4 = plt.subplot(gs[3], sharex=ax3)
    ax3.tick_params('x', direction='in', bottom=True)
    ax4.tick_params('x', direction='in', top=True)
    plt.setp(ax3.get_xticklabels(), visible=False)
    ax4.plot(t, ang, 's', marker='x', color='k', markersize=3)
    ax4.set_ylabel("phase", fontsize=11, labelpad=-5)
    ax4.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax4.set_xlabel("$\gamma/\lambda$", fontsize=11)
    ax4.set_xlim([0, 4])
    # ax4.axvspan(0, 1, alpha=0.5, color='slateblue')
    # ax4.axvspan(1, 4, alpha=0.5, color='gold')
    ax3.set_title("$C^2_2=+1$ (no SPT detected)")

    fig.text(0.04, 0.43, "$\langle \Psi | \\tilde{C}_4 |\Psi\\rangle$", fontsize=12, rotation='vertical')
    fig.text(0.495, 0.43, "$\langle \Psi | \\tilde{C}_2 |\Psi\\rangle$", fontsize=12, rotation='vertical')

    ####################################################################################################################

    plt.savefig(f"/home/bart/Documents/presentations/2021_03_18/figures/prelim_results.png", bbox_inches='tight', dpi=300)
    plt.show()
