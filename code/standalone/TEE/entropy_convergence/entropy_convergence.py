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

plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}')
# matplotlib.verbose.level = 'debug-annoying'

if __name__ == '__main__':

    # extract data from file
    with open('FerHofSqu1.dat', 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        chi = []
        entanglement_entropy = []
        delta_entropy = []
        for row in plots:
            chi.append(float(row[0]))
            entanglement_entropy.append(float(row[1]))
            if len(row) > 2:
                delta_entropy.append(float(row[2]))

    print(chi)
    print(entanglement_entropy)
    print(delta_entropy)

    ####################################################################################################################

    fig = plt.figure(figsize=(11, 2))
    gs = gridspec.GridSpec(1, 3, width_ratios=[1, 1, 1], wspace=0.5)

    ####################################################################################################################

    ax1 = plt.subplot(gs[0])

    invchi = [1/chi_val for chi_val in chi]

    ax1.plot(invchi, entanglement_entropy, 's', marker='x', color='k', markersize=3)

    print(invchi)

    # lower limit
    plt.hlines(entanglement_entropy[-1:], xmin=0, xmax=invchi[0], color='k', linewidth=0.5, linestyles='dashed')
    #ax1.axhline(entanglement_entropy[-1:], xmin=0, xmax=invchi[0], color='k', linewidth=0.5, ls='--')
    # upper limit
    x = np.linspace(0, invchi[0], 100000)
    c, m = polyfit(invchi[-2:], entanglement_entropy[-2:], 1)
    ax1.plot(x, m * x + c, '--', c='k', linewidth=0.5)

    ax1.set_xlim(0)
    ax1.tick_params(axis="x", labelsize=8)
    ax1.tick_params(axis="y", labelsize=8)
    ax1.set_xlabel("$\\chi^{-1}$", fontsize=10)
    ax1.set_ylabel("$S_\mathrm{{vN}}$", fontsize=10)

    ####################################################################################################################

    ax2 = plt.subplot(gs[1])

    ax2.plot(invchi[-2:], entanglement_entropy[-2:], 's', marker='x', color='k', markersize=3)

    # lower limit
    plt.hlines(entanglement_entropy[-1:], xmin=0, xmax=invchi[-2:][1], color='k', linewidth=0.5, linestyles='dashed')
    # upper limit
    x = np.linspace(0, invchi[-2:][0], 100000)
    ax2.plot(x, m * x + c, '--', c='k', linewidth=0.5)
    fig.text(0.47, 0.7, "$S={gradient:.4f}\\chi^{{-1}}+{intercept:.4f}$".format(gradient=m, intercept=c), fontsize=8)
    # data point
    e = ax2.errorbar(0, np.average([entanglement_entropy[-1:], c]), marker='s', zorder=1, color='r', yerr=np.abs(entanglement_entropy[-1:]-c)/2, markersize=10, capsize=10, clip_on=False)
    for b in e[1]:  # unclip the error caps
        b.set_clip_on(False)
    for b in e[2]:  # unclip the error lines
        b.set_clip_on(False)
    # error percentage
    error = ((np.abs(entanglement_entropy[-1:]-c)/2) / np.average([entanglement_entropy[-1:], c]) * 100)[0]
    fig.text(0.425, 0.25, f"$\sim{error:.4f}\%$ error", fontsize=8)

    ax2.set_xlim(0)
    ax2.set_ylim(1.42655, 1.431)
    ax2.tick_params(axis="x", labelsize=8)
    ax2.tick_params(axis="y", labelsize=8)
    ax2.set_xlabel("$\\chi^{-1}$", fontsize=10)
    ax2.set_ylabel("$S_\mathrm{{vN}}$", fontsize=10)

    ####################################################################################################################

    ax3 = plt.subplot(gs[2])

    ln_delta_entropy = [np.log(delta_entropy_val) for delta_entropy_val in delta_entropy]

    ax3.plot(invchi[1:], ln_delta_entropy, 's', marker='x', color='k', markersize=3)

    ax3.set_xlim(0)
    ax3.tick_params(axis="x", labelsize=8)
    ax3.tick_params(axis="y", labelsize=8)
    ax3.set_xlabel("$\\chi^{-1}$", fontsize=10)
    ax3.set_ylabel("$\ln (\Delta S_\mathrm{{vN}})$", fontsize=10)

    ####################################################################################################################

    fig.text(0.07, 0.875, "(a)", fontsize=10)
    fig.text(0.35, 0.875, "(b)", fontsize=10)
    fig.text(0.66, 0.875, "(c)", fontsize=10)

    plt.savefig("/home/bart/Documents/papers/TEE/figures/entropy_convergence.png", bbox_inches='tight', dpi=300)
    plt.show()
