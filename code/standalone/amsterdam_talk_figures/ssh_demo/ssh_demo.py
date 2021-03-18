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
plt.rc('text.latex', preamble=r'\usepackage{amsmath}\usepackage{amssymb}')
# matplotlib.verbose.level = 'debug-annoying'

if __name__ == '__main__':

    fig = plt.figure(figsize=(10, 3))
    gs = gridspec.GridSpec(2, 2, hspace=0, height_ratios=[1, 1])

    # define a list of easily-visible markers
    markers = [(3, 0, 0), (4, 0, 0), (5, 0, 0), (6, 0, 0), (4, 1, 0), (5, 1, 0), (6, 1, 0),
               (3, 2, 0), (4, 2, 0), (5, 2, 0), (6, 2, 0), 'X', 'x', 'd', 'D', 'P',
               '$a$', '$b$', '$c$', '$d$', '$e$', '$f$', '$g$', '$h$', '$i$', '$j$', '$k$', '$l$', '$m$', '$n$', '$o$',
               '$p$', '$q$', '$r$', '$s$', '$t$', '$u$', '$v$', '$w$', '$x$', '$y$', '$z$']

    ####################################################################################################################

    ax1 = plt.subplot(gs[0])

    # extract data from file
    with open('ssh_magnitude_plus.dat', 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        ratio = []
        magnitude = []
        for row in plots:
            ratio.append(float(row[0]))
            magnitude.append(float(row[1]))

    ax1.plot(ratio, magnitude, '.-', c='k')
    ax1.set_xlim([0, 2])
    # ax1.set_ylim([0.3, 1.2])
    ax1.set_ylabel("magnitude", fontsize=12)
    ax1.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    plt.setp(ax1.get_xticklabels(), visible=False)
    plt.tick_params(
        axis='x',  # changes apply to the x-axis
        which='both',  # both major and minor ticks are affected
        bottom=False,  # ticks along the bottom edge are off
        top=False,  # ticks along the top edge are off
        labelbottom=False)
    # ax1.axvspan(0, 1, alpha=0.5, color='slateblue')
    # ax1.axvspan(1, 2, alpha=0.5, color='gold')
    # ax1.text(0.4, 0.5, "trivial", fontsize=12)
    # ax1.text(1.4, 0.8, "topological", fontsize=12)
    ax1.set_title("$\mathcal{R}^2=+1$ (no SPT detected)")

    ####################################################################################################################

    ax2 = plt.subplot(gs[2])

    # extract data from file
    with open('ssh_phase_plus.dat', 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        ratio = []
        phase = []
        for row in plots:
            ratio.append(float(row[0]))
            phase.append(float(row[1]))

    ax2.plot(ratio, phase, '.-', c='k')
    ax2.set_xlim([0, 2])
    # ax2.set_ylim([-0.2, 2.2])
    ax2.xaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax2.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax2.set_ylabel("phase", fontsize=12, labelpad=-5)
    ax2.set_xlabel("$t_2/t_1$", fontsize=12)
    # ax2.axvspan(0, 1, alpha=0.5, color='slateblue')
    # ax2.axvspan(1, 2, alpha=0.5, color='gold')
    # ax2.text(0.3, 1, "$\langle \Psi | \mathcal{R}_\mathrm{part} | \Psi \\rangle \leq 1$", fontsize=12)
    # ax2.text(1.3, 0.7, "$\langle \Psi | \mathcal{R}_\mathrm{part} | \Psi \\rangle = \\frac{1}{2}e^{\mathrm{i}\pi/2}$",
    #          fontsize=12)

    # labels = [item.get_text() for item in ax2.get_yticklabels()]
    # print(labels)
    # labels[1] = '$0$'
    # labels[2] = '$\pi/4$'
    # labels[3] = '$\pi/2$'
    # print(labels)
    # ax2.set_yticklabels(labels)

    fig.text(0.03, 0.36, "$\langle \Psi | \mathcal{R}_\\text{part}|\Psi\\rangle$", fontsize=12, rotation='vertical')

    ####################################################################################################################

    ax3 = plt.subplot(gs[1])

    # extract data from file
    with open('ssh_magnitude_minus.dat', 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        ratio = []
        magnitude = []
        for row in plots:
            ratio.append(float(row[0]))
            magnitude.append(float(row[1]))

    ax3.plot(ratio, magnitude, '.-', c='k')
    ax3.set_xlim([0, 2])
    ax3.set_ylim([0.3, 1.2])
    ax3.set_ylabel("magnitude", fontsize=12)
    ax3.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    plt.setp(ax3.get_xticklabels(), visible=False)
    plt.tick_params(
        axis='x',  # changes apply to the x-axis
        which='both',  # both major and minor ticks are affected
        bottom=False,  # ticks along the bottom edge are off
        top=False,  # ticks along the top edge are off
        labelbottom=False)
    ax3.axvspan(0, 1, alpha=0.5, color='slateblue')
    ax3.axvspan(1, 2, alpha=0.5, color='gold')
    ax3.text(0.4, 0.5, "trivial", fontsize=12)
    ax3.text(1.3, 0.8, "topological", fontsize=12)
    ax3.set_title("$\mathcal{R}^2=-1$ (SPT with $\\mathbb{Z}_4$ classification)")

    ####################################################################################################################

    ax4 = plt.subplot(gs[3])

    # extract data from file
    with open('ssh_phase_minus.dat', 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        ratio = []
        phase = []
        for row in plots:
            ratio.append(float(row[0]))
            phase.append(float(row[1]))

    ax4.plot(ratio, phase, '.-', c='k')
    ax4.set_xlim([0, 2])
    ax4.set_ylim([-0.2, 2.2])
    ax4.xaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax4.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax4.set_ylabel("phase", fontsize=12)
    ax4.set_xlabel("$t_2/t_1$", fontsize=12)
    ax4.axvspan(0, 1, alpha=0.5, color='slateblue')
    ax4.axvspan(1, 2, alpha=0.5, color='gold')
    ax4.text(0.2, 1, "$\langle \Psi | \mathcal{R}_\mathrm{part} | \Psi \\rangle \leq 1$", fontsize=12)
    ax4.text(1.15, 0.7, "$\langle \Psi | \mathcal{R}_\mathrm{part} | \Psi \\rangle = \\frac{1}{2}e^{\mathrm{i}\pi/2}$", fontsize=12)

    labels = [item.get_text() for item in ax4.get_yticklabels()]
    print(labels)
    labels[1] = '$0$'
    labels[2] = '$\pi/4$'
    labels[3] = '$\pi/2$'
    print(labels)
    ax4.set_yticklabels(labels)

    fig.text(0.03, 0.36, "$\langle \Psi | \mathcal{R}_\\text{part}|\Psi\\rangle$", fontsize=12, rotation='vertical')
    # fig.text(0.46, 0.36, "$\langle \Psi | \mathcal{R}_\\text{part}|\Psi\\rangle$", fontsize=12, rotation='vertical')

    ####################################################################################################################

    plt.savefig("/home/bart/Documents/presentations/2021_03_18/figures/ssh_demo.png", bbox_inches='tight', dpi=300)
    plt.show()
