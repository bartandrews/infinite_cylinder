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

if __name__ == '__main__':

    fig = plt.figure(figsize=(11, 2))
    gs = gridspec.GridSpec(1, 4, width_ratios=[1, 1, 1, 1])
    ax0 = plt.subplot(gs[0])

    alpha = []
    chi_range = np.linspace(50, 300, 6)
    entanglement_energies = {}
    for chi in chi_range:
        entanglement_energies[chi] = []

    with open('FerHofHex1Hex5Orbital.dat', 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        for row in plots:
            alpha.append(float(row[0]))
            for i, chi in enumerate(chi_range):
                entanglement_energies[chi].append(float(row[i + 1]))

    for i, chi in enumerate(chi_range):
        ax0.plot(alpha[:int(chi)], entanglement_energies[chi][:int(chi)], 's', marker='o', color='C{}'.format(i), markersize=1, label='$\\chi={}$'.format(chi))

    ax0.set_xlim([0, 300])
    ax0.set_xticks(np.arange(0, 301, 100))
    ax0.set_ylim([0, 30])
    ax0.set_yticks(np.arange(0, 31, 10))
    ax0.tick_params(axis="x", labelsize=8)
    ax0.tick_params(axis="y", labelsize=8)
    ax0.set_xlabel("$\\alpha$", fontsize=10)
    ax0.set_ylabel("$\\epsilon_\\alpha$", fontsize=10)
    for chi in chi_range:
        ax0.axvline(chi, color='k', linewidth=0.5, ls='--')

    ax0.set_aspect('auto')

    ########################################################################################################################

    gs.update(wspace=0.5)

    ax1 = plt.subplot(gs[1])

    chi_range = np.linspace(50, 300, 6, dtype=int)
    entropy_elements = {}
    entanglement_entropy = np.zeros(len(chi_range))
    inv_chi = np.zeros(len(chi_range))
    for i, chi in enumerate(chi_range):
        inv_chi[i] = 1/chi

    for chi in chi_range:
        entropy_elements[chi] = []

    with open('FerHofHex1Hex5Orbital.dat', 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        for row in plots:
            for i, chi in enumerate(chi_range):
                entropy_elements[chi].append(float(row[i + 1])*np.exp(-float(row[i + 1])))

    for i, chi in enumerate(chi_range):
        entanglement_entropy[i] = sum(entropy_elements[chi])

    ax1.plot(inv_chi, entanglement_entropy, 's', marker='x', color='k', markersize=3)

    c, m = polyfit(inv_chi[-3:], entanglement_entropy[-3:], 1)
    ax1.plot(inv_chi, m * inv_chi + c, '-', c='C8', zorder=0)
    fig.text(0.3425, 0.15, "$S={gradient:.2f}\\chi^{{-1}}+{intercept:.2f}$".format(gradient=m, intercept=c), fontsize=8)

    #ax1.set_xlim([0, 500])
    #ax1.set_xticks(np.arange(0, 501, 100))
    #ax1.set_ylim([0, 30])
    #ax1.set_yticks(np.arange(0, 31, 10))
    ax1.tick_params(axis="x", labelsize=8)
    ax1.tick_params(axis="y", labelsize=8)
    ax1.set_xlabel("$\\chi^{-1}$", fontsize=10)
    ax1.set_ylabel("$S$", fontsize=10)

    ####################################################################################################################

    ax2 = plt.subplot(gs[2])

    alpha = []
    chi_range = np.linspace(50, 500, 10)
    entanglement_energies = {}
    for chi in chi_range:
        entanglement_energies[int(chi)] = []

    with open('FerHofHex1.dat', 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        for row in plots:
            alpha.append(float(row[0]))
            for i, chi in enumerate(chi_range):
                entanglement_energies[chi].append(float(row[i + 1]))

    for i, chi in enumerate(chi_range):
        ax2.plot(alpha[:int(chi)], entanglement_energies[chi][:int(chi)], 's', marker='o', color='C{}'.format(i),
                 markersize=1, label='$\\chi={}$'.format(chi))

    ax2.set_xlim([0, 500])
    ax2.set_xticks(np.arange(0, 501, 100))
    ax2.set_ylim([0, 30])
    ax2.set_yticks(np.arange(0, 31, 10))
    ax2.tick_params(axis="x", labelsize=8)
    ax2.tick_params(axis="y", labelsize=8)
    ax2.set_xlabel("$\\alpha$", fontsize=10)
    ax2.set_ylabel("$\\epsilon_\\alpha$", fontsize=10)
    for chi in chi_range:
        ax2.axvline(chi, color='k', linewidth=0.5, ls='--')

    ####################################################################################################################

    ax3 = plt.subplot(gs[3])

    chi_range = np.linspace(50, 500, 10, dtype=int)
    entropy_elements = {}
    entanglement_entropy = np.zeros(len(chi_range))
    inv_chi = np.zeros(len(chi_range))
    for i, chi in enumerate(chi_range):
        inv_chi[i] = 1 / chi

    for chi in chi_range:
        entropy_elements[chi] = []

    with open('FerHofHex1.dat', 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        for row in plots:
            for i, chi in enumerate(chi_range):
                entropy_elements[chi].append(float(row[i + 1]) * np.exp(-float(row[i + 1])))

    for i, chi in enumerate(chi_range):
        entanglement_entropy[i] = sum(entropy_elements[chi])

    ax3.plot(inv_chi, entanglement_entropy, 's', marker='x', color='k', markersize=3)

    c, m = polyfit(inv_chi[-5:], entanglement_entropy[-5:], 1)
    ax3.plot(inv_chi, m * inv_chi + c, '-', c='C8', zorder=0)
    fig.text(0.765, 0.15, "$S={gradient:.2f}\\chi^{{-1}}+{intercept:.2f}$".format(gradient=m, intercept=c), fontsize=8)

    # ax3.set_xlim([0, 500])
    # ax3.set_xticks(np.arange(0, 501, 100))
    # ax3.set_ylim([0, 30])
    # ax3.set_yticks(np.arange(0, 31, 10))
    ax3.tick_params(axis="x", labelsize=8)
    ax3.tick_params(axis="y", labelsize=8)
    ax3.set_xlabel("$\\chi^{-1}$", fontsize=10)
    ax3.set_ylabel("$S$", fontsize=10)

    ####################################################################################################################

    fig.text(0.07, 0.875, "(a)", fontsize=10)
    fig.text(0.275, 0.875, "(b)", fontsize=10)
    fig.text(0.49, 0.875, "(c)", fontsize=10)
    fig.text(0.7, 0.875, "(d)", fontsize=10)

    plt.savefig("/home/bart/Documents/papers/TBG_rev1/figures/schmidt_values.png", bbox_inches='tight', dpi=300)
    plt.show()
