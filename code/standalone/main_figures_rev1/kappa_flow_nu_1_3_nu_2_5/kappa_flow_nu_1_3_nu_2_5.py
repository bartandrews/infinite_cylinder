import numpy as np
import matplotlib
import csv
import matplotlib.pyplot as plt
import sys
from itertools import product
import matplotlib.gridspec as gridspec
from matplotlib.patches import Polygon
import matplotlib.patches as patches

if __name__ == '__main__':

    fig = plt.figure()

    gs = gridspec.GridSpec(3, 1, height_ratios=[1, 1, 3])

    ax0 = plt.subplot(gs[0])

    phi_IQH = []
    IQH_charge = []
    phi_FQH = []
    phi_FQH_out = []

    # an array of parameters, each of our curves depend on a specific
    # value of parameters
    chi_range = np.linspace(100, 150, 1)

    # norm is a class which, when called, can normalize data into the
    # [0.0, 1.0] interval.
    norm = matplotlib.colors.Normalize(
        vmin=np.min(chi_range),
        vmax=np.max(chi_range))

    # choose a colormap
    # c_m = matplotlib.cm.cool
    c_m = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white", "C9"])

    # create a ScalarMappable and initialize a data structure
    s_m = matplotlib.cm.ScalarMappable(cmap=c_m, norm=norm)
    s_m.set_array([])

    FQH_charge = {}
    for chi in chi_range:
        FQH_charge[chi] = []

    with open('DeltaQL.dat', 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        for row in plots:
            phi_FQH.append(float(row[0]))
            for i, chi in enumerate(chi_range):
                FQH_charge[chi].append(float(row[i + 1]))

    ax0.plot(phi_FQH, FQH_charge[chi], 's', marker='o', color='k',
             markersize=3, mec='k')

    # FQH_charge_out = {}
    # for chi in chi_range:
    #     FQH_charge_out[chi] = []
    #
    # with open('DeltaQL_invt2dash_flow.dat.hex1hex5orbital.outliers', 'r') as csvfile:
    #     plots = csv.reader(csvfile, delimiter='\t')
    #     for row in plots:
    #         phi_FQH_out.append(float(row[0]))
    #         for i, chi in enumerate(chi_range):
    #             FQH_charge_out[chi].append(float(row[i + 1]))
    #
    # ax0.plot(phi_FQH_out, FQH_charge_out[chi], 's', marker='o', color='r',
    #          markersize=3, mec='r')

    # for i, chi in enumerate(chi_range):
    #     if i == len(chi_range) - 1:
    #         ax0.plot(phi_FQH, FQH_charge[chi], marker='s', linestyle='--', color=s_m.to_rgba(chi),
    #                  markersize=5, mec='k')
    #     else:
    #         ax0.plot(phi_FQH, FQH_charge[chi], marker='s', linestyle='--', color=s_m.to_rgba(chi), markersize=5)

    # # cbax = plt.subplot()  # Place it where it should be.
    # # # --------------------------------------------------------
    # # cb = Colorbar(mappable=s_m, ax=cbax, orientation='horizontal', ticklocation='top')
    # # cb.set_label(r'Colorbar !', labelpad=10)

    # cbax = fig.add_axes([0.92, 0.625, 0.025, 0.255])
    # cb = plt.colorbar(s_m, cax=cbax, orientation='vertical', ticklocation='right', ticks=[50, 100, 150, 200, 250])
    # cb.ax.set_yticklabels(['$0.5$', '$1$', '$1.5$', '$2$', '$2.5$'])  # horizontal colorbar
    # cb.set_label("$\chi / 10^2$", fontsize=11, labelpad=7)

    metal0 = Polygon(((0, -10), (0.4, -10), (0.4, 10), (0, 10)), fc=(0, 0, 0, 0.2))
    ax0.add_artist(metal0)

    ax0.tick_params(axis="y", labelsize=10)

    ax0.set_xlim([0, 1])
    ax0.set_xticks(np.arange(0, 1, 11))
    ax0.set_ylim([-0.5, 1.5])
    ax0.set_yticks(np.arange(-0.5, 1.6, 0.5))
    ax0.set_ylabel("$\\langle \\Delta Q_\\mathrm{L}\\rangle$", fontsize=11)
    ax0.tick_params('x', direction='in', bottom=True)
    ax0.axhline(0, color='k', linewidth=0.5, ls='--')
    ax0.axhline(1, color='k', linewidth=0.5, ls='--')
    plt.setp(ax0.get_xticklabels(), visible=False)

    ########################################################################################################################

    ax1 = plt.subplot(gs[1], sharex=ax0)

    phi_IQH = []
    IQH_charge = []
    phi_FQH = []
    phi_FQH_out = []

    # an array of parameters, each of our curves depend on a specific
    # value of parameters
    chi_range = np.linspace(100, 150, 1)

    # norm is a class which, when called, can normalize data into the
    # [0.0, 1.0] interval.
    norm = matplotlib.colors.Normalize(
        vmin=np.min(chi_range),
        vmax=np.max(chi_range))
    norm2 = matplotlib.colors.Normalize(
        vmin=0,
        vmax=1)

    # choose a colormap
    # c_m = matplotlib.cm.cool
    c_m = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white", "C9"])
    c_m2 = matplotlib.colors.LinearSegmentedColormap.from_list("", ["black", "red"])

    # create a ScalarMappable and initialize a data structure
    s_m = matplotlib.cm.ScalarMappable(cmap=c_m, norm=norm)
    s_m2 = matplotlib.cm.ScalarMappable(cmap=c_m2, norm=norm2)
    s_m.set_array([])
    s_m2.set_array([])

    FQH_charge = {}
    for chi in chi_range:
        FQH_charge[chi] = []

    with open('corr_len.dat', 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        for row in plots:
            phi_FQH.append(float(row[0]))
            for i, chi in enumerate(chi_range):
                FQH_charge[chi].append(float(row[i + 1]))

    ax1.plot(phi_FQH, FQH_charge[chi], 's', marker='x', color='k', markersize=5, mec='k')

    # FQH_charge_out = {}
    # for chi in chi_range:
    #     FQH_charge_out[chi] = []
    #
    # with open('corr_func.dat.hex1hex5orbital.outliers', 'r') as csvfile:
    #     plots = csv.reader(csvfile, delimiter='\t')
    #     for row in plots:
    #         phi_FQH_out.append(float(row[0]))
    #         for i, chi in enumerate(chi_range):
    #             FQH_charge_out[chi].append(float(row[i + 1]))
    #
    # ax1.plot(phi_FQH_out, FQH_charge_out[chi], 's', marker='x', color='red',
    #          markersize=3)

    # for i, chi in enumerate(chi_range):
    #     if i == len(chi_range) - 1:
    #         ax0.plot(phi_FQH, FQH_charge[chi], marker='s', linestyle='--', color=s_m.to_rgba(chi),
    #                  markersize=5, mec='k')
    #     else:
    #         ax0.plot(phi_FQH, FQH_charge[chi], marker='s', linestyle='--', color=s_m.to_rgba(chi), markersize=5)

    # # cbax = plt.subplot()  # Place it where it should be.
    # # # --------------------------------------------------------
    # # cb = Colorbar(mappable=s_m, ax=cbax, orientation='horizontal', ticklocation='top')
    # # cb.set_label(r'Colorbar !', labelpad=10)

    # cbax = fig.add_axes([0.92, 0.625, 0.025, 0.255])
    # cb = plt.colorbar(s_m, cax=cbax, orientation='vertical', ticklocation='right', ticks=[50, 100, 150, 200, 250])
    # cb.ax.set_yticklabels(['$0.5$', '$1$', '$1.5$', '$2$', '$2.5$'])  # horizontal colorbar
    # cb.set_label("$\chi / 10^2$", fontsize=11, labelpad=7)

    metal0 = Polygon(((0, -10), (0.4, -10), (0.4, 10), (0, 10)), fc=(0, 0, 0, 0.2))
    ax1.add_artist(metal0)
    fig.text(0.225, 0.9, 'FQH state', fontsize=11)
    fig.text(0.6, 0.9, 'top. trivial', fontsize=11)

    ax1.axhline(2.87640603731211, color='k', linewidth=0.5, ls='--')

    ax1.tick_params(axis="y", labelsize=10)

    ax1.set_xlim([0, 1])
    ax1.set_xticks(np.arange(0, 1, 11))
    ax1.set_ylim([2.6, 3])
    ax1.set_yticks(np.arange(2.6, 3, 0.1))
    ax1.set_ylabel("$\\ln \\, \\xi$", fontsize=11)
    ax1.tick_params('x', direction='in', bottom=True)
    plt.setp(ax1.get_xticklabels(), visible=False)

########################################################################################################################

    gs.update(hspace=0)

    ax2 = plt.subplot(gs[2], sharex=ax1)

    x = []
    y = []
    z = []

    with open('ent_spec.dat', 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        for row in plots:
            x.append(float(row[1]))
            y.append(float(row[2]))
            z.append(int(row[0]))

    for value in np.linspace(4, -4, 9, dtype=int):
        xvalue = []
        yvalue = []
        for i in range(len(x)):
            if z[i] == value:
                xvalue.append(x[i])
                yvalue.append(y[i])
        ax2.scatter(xvalue, yvalue, marker='_', c='C{}'.format((value + 4) % 10), label='{}'.format(value))

    ax2.legend(loc='center left', handletextpad=0, borderpad=0.2, framealpha=1, edgecolor='k', markerscale=2,
               fontsize=12, ncol=1, bbox_to_anchor=(1.01, 0.5))

    metal0 = Polygon(((0, 0), (0.4, 0), (0.4, 20), (0, 20)), fc=(0, 0, 0, 0.2))
    ax2.add_artist(metal0)

    ax2.set_xlim([0, 1])
    ax2.set_xticks(np.arange(0, 1.1, 0.1))
    ax2.set_ylim([0, 16])
    ax2.set_yticks(np.arange(0, 16, 2))
    ax2.set_ylabel("$\epsilon_\\alpha$", fontsize=11)
    ax2.set_xlabel("$\kappa$", fontsize=11)
    plt.setp(ax2.get_xticklabels(), visible=True)

    ax2.axhline(0.19703914636968, color='k', linewidth=0.5, ls='--')
    ax2.axhline(2.7, color='k', linewidth=0.5, ls='--')
    ax2.axhline(3.2, color='k', linewidth=0.5, ls='--')
    ax2.axhline(3.31168612016415, color='k', linewidth=0.5, ls='--')
    ax2.axhline(4.29508476900778, color='k', linewidth=0.5, ls='--')
    ax2.axhline(4.73812358078833, color='k', linewidth=0.5, ls='--')

    ax2.tick_params(axis="x", labelsize=10)
    ax2.tick_params(axis="y", labelsize=10)

    fig.text(0.01, 0.87, "(a)", fontsize=12)
    fig.text(0.01, 0.71, "(b)", fontsize=12)
    fig.text(0.01, 0.56, "(c)", fontsize=12)

    plt.savefig("/home/bart/Documents/papers/TBG_rev1/figures/kappa_flow.png", bbox_inches='tight', dpi=300)
    plt.show()
