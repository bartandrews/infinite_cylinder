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

    gs = gridspec.GridSpec(2, 1, height_ratios=[1, 2])

    ax0 = plt.subplot(gs[0])

    phi_CI = []
    CI_charge = []
    phi_FCI = []

    # an array of parameters, each of our curves depend on a specific
    # value of parameters
    chi_range = np.linspace(100, 500, 5)

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

    FCI_charge = {}
    for chi in chi_range:
        FCI_charge[chi] = []

    with open('corr_func_example.dat', 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        for row in plots:
            phi_FCI.append(float(row[0]))
            for i, chi in enumerate(chi_range):
                FCI_charge[chi].append(float(row[i + 1]))

    for i, chi in enumerate(chi_range):
        if i == len(chi_range) - 1:
            ax0.plot(phi_FCI, FCI_charge[chi], marker='s', linestyle='--', color=s_m.to_rgba(chi),
                     markersize=5, mec='k')
        else:
            ax0.plot(phi_FCI, FCI_charge[chi], marker='s', linestyle='--', color=s_m.to_rgba(chi), markersize=5)

    # cbax = plt.subplot()  # Place it where it should be.
    # # --------------------------------------------------------
    # cb = Colorbar(mappable=s_m, ax=cbax, orientation='horizontal', ticklocation='top')
    # cb.set_label(r'Colorbar !', labelpad=10)

    cbax = fig.add_axes([0.92, 0.625, 0.025, 0.255])
    cb = plt.colorbar(s_m, cax=cbax, orientation='vertical', ticklocation='right', ticks=[100, 200, 300, 400, 500])
    cb.ax.set_yticklabels(['$1$', '$2$', '$3$', '$4$', '$5$'])  # horizontal colorbar
    cb.set_label("$\chi / 10^2$", fontsize=11, labelpad=7)

    metal0 = Polygon(((0, -1.5), (1.5, -1.5), (1.5, 1), (0, 1)), fc=(0, 0, 0, 0.1))
    ax0.add_artist(metal0)
    fig.text(0.275, 0.9, 'metal', fontsize=11)
    fig.text(0.69, 0.9, 'FCI', fontsize=11)

    ax0.tick_params(axis="y", labelsize=10)

    ax0.set_xlim([0, 3])
    ax0.set_xticks(np.arange(0, 3, 1))
    ax0.set_ylim([-1.5, 1])
    ax0.set_yticks(np.arange(-1.5, 1.5, 0.5))
    ax0.set_ylabel("$\\xi$", fontsize=11)
    ax0.tick_params('x', direction='in', bottom=True)
    plt.setp(ax0.get_xticklabels(), visible=False)

########################################################################################################################

    gs.update(hspace=0)

    ax1 = plt.subplot(gs[1], sharex=ax0)

    x = []
    y = []
    z = []

    with open('ent_spec_V_flow_example.dat', 'r') as csvfile:
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
        ax1.scatter(xvalue, yvalue, marker='_', c='C{}'.format((value + 4) % 10), label='{}'.format(value))

    ax1.legend(loc='center left', handletextpad=0, borderpad=0.2, framealpha=1, edgecolor='k', markerscale=2,
               fontsize=12, ncol=1, bbox_to_anchor=(1.01, 0.5))

    metal1 = Polygon(((0, 0), (1.5, 0), (1.5, 10), (0, 10)), fc=(0, 0, 0, 0.1))
    ax1.add_artist(metal1)

    ax1.set_xlim([0, 3])
    ax1.set_xticks(np.arange(0, 4, 1))
    ax1.set_ylim([0, 10])
    ax1.set_yticks(np.arange(0, 10, 1))
    ax1.set_ylabel("$\epsilon_\\alpha$", fontsize=11)
    ax1.set_xlabel("$V_1$", fontsize=11)

    ax1.tick_params(axis="x", labelsize=10)
    ax1.tick_params(axis="y", labelsize=10)

    plt.text(-0.45, 15, "(a)", fontsize=12)
    plt.text(-0.45, 9.5, "(b)", fontsize=12)

    plt.savefig("/home/bart/Documents/papers/TBG/figures/V_flow.png", bbox_inches='tight', dpi=300)
    plt.show()
