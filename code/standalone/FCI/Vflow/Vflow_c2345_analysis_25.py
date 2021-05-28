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

    fig = plt.figure(figsize=(6, 4))
    outer_grid = gridspec.GridSpec(1, 2, wspace=0.6, hspace=0.3)
    left_cell = outer_grid[0, 0]
    right_cell = outer_grid[0, 1]
    # lower_left_cell = outer_grid[1, 0]
    # lower_right_cell = outer_grid[1, 1]
    left_inner_grid = gridspec.GridSpecFromSubplotSpec(3, 1, left_cell, hspace=0, height_ratios=[1, 1, 2])
    right_inner_grid = gridspec.GridSpecFromSubplotSpec(3, 1, right_cell, hspace=0, height_ratios=[1, 1, 2])

    # define a list of easily-visible markers
    markers = [(3, 0, 0), (4, 0, 0), (5, 0, 0), (6, 0, 0), (4, 1, 0), (5, 1, 0), (6, 1, 0),
               (3, 2, 0), (4, 2, 0), (5, 2, 0), (6, 2, 0), 'X', 'x', 'd', 'D', 'P',
               '$a$', '$b$', '$c$', '$d$', '$e$', '$f$', '$g$', '$h$', '$i$', '$j$', '$k$', '$l$', '$m$', '$n$', '$o$',
               '$p$', '$q$', '$r$', '$s$', '$t$', '$u$', '$v$', '$w$', '$x$', '$y$', '$z$']

    ####################################################################################################################

    # ax1 = plt.subplot(upper_left_cell)
    #
    # for i, chi_val in enumerate(np.arange(100, 400, 100)):
    #     corr_len_dir = '/home/bart/PycharmProjects/infinite_cylinder/data/corr_len_V_flow/FerHofSqu1'
    #     corr_len_file = f'corr_len_V_flow_FerHofSqu1_chi_{chi_val}_t1_1_V_0_10_41_Coulomb_1_n_1_55_nphi_6_11_LxMUC_1_Ly_10.dat'
    #     corr_len_path = os.path.join(corr_len_dir, corr_len_file)
    #
    #     # extract data from file
    #     with open(corr_len_path, 'r') as csvfile:
    #         plots = csv.reader(csvfile, delimiter='\t')
    #         V = []
    #         xi = []
    #         for row in plots:
    #             V.append(float(row[0]))
    #             xi.append(float(row[1]))
    #
    #     ax1.plot(V, xi, c=f"C{i}", marker=markers[i], markersize=2, label=f"$\chi={chi_val}$")
    #
    # ax1.legend(loc='upper right', handletextpad=0, handlelength=1, borderpad=0, framealpha=1, edgecolor='w',
    #            markerscale=1, fontsize=10, ncol=2)
    #
    # ax1.set_xlim([0, 10])
    # ax1.set_xlabel("$V$", fontsize=11)
    # ax1.xaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    # ax1.set_ylabel("$\\xi$", fontsize=11)
    # ax1.set_title("(a) $\\nu=1/5$, $n_\phi=6/11$", fontsize=12)
    # ax1.axvspan(0, 2, alpha=0.5, color='grey')
    # line1 = plt.Line2D((.125, .125), (.645, .748), color="k", linewidth=1, linestyle='--', alpha=0.5)
    # fig.add_artist(line1)
    # line2 = plt.Line2D((.435, .139), (.645, .748), color="k", linewidth=1, linestyle='--', alpha=0.5)
    # fig.add_artist(line2)

    ####################################################################################################################

    ax2 = plt.subplot(left_inner_grid[0])

    for i, chi_val in enumerate(np.arange(100, 600, 100)):
        corr_len_dir = '/home/bart/PycharmProjects/infinite_cylinder/data/corr_len_V_flow/FerHofSqu1'
        corr_len_file = f'corr_len_V_flow_FerHofSqu1_chi_{chi_val}_t1_1_V_0_3_25_Coulomb_1_n_1_55_nphi_6_11_LxMUC_1_Ly_10.dat'
        corr_len_path = os.path.join(corr_len_dir, corr_len_file)

        # extract data from file
        with open(corr_len_path, 'r') as csvfile:
            plots = csv.reader(csvfile, delimiter='\t')
            V = []
            xi = []
            for row in plots:
                V.append(float(row[0]))
                xi.append(float(row[1]))

        ax2.plot(V, xi, c=f"C{i}", marker=markers[i], markersize=2, label=f"${chi_val/100:g}$")

    ax2.legend(loc='center left', handletextpad=0.5, handlelength=1, borderpad=0.2, framealpha=1, edgecolor=None,
               markerscale=2, fontsize=10, ncol=1, bbox_to_anchor=(1, 0), title="$\chi/10^2$")

    ax2.axvline(0.00678, c='b', linestyle=':', zorder=-2)
    ax2.axvline(0.300, c='r', linestyle=':', zorder=-2)
    ax2.grid(color='k', linestyle='-', linewidth=0.3, axis='x')
    ax2.set_xlim([0, 3])
    ax2.set_xticks(np.arange(0, 3.1, 1))
    ax2.set_yticks(np.arange(0, 450, 200))
    ax2.set_ylabel("$\\xi$", fontsize=11)
    plt.setp(ax2.get_xticklabels(), visible=False)
    ax2.tick_params(
        axis='x',  # changes apply to the x-axis
        which='both',  # both major and minor ticks are affected
        bottom=False,  # ticks along the bottom edge are off
        top=False,  # ticks along the top edge are off
        labelbottom=False)
    ax2.set_title("(a) $\\nu=1/5$, $n_\phi=6/11$", fontsize=12)
    
    ####################################################################################################################

    ax3 = plt.subplot(left_inner_grid[1])

    for i, chi_val in enumerate(np.arange(100, 600, 100)):
        ent_dir = '/home/bart/PycharmProjects/infinite_cylinder/data/ent_V_flow/FerHofSqu1'
        ent_file = f'ent_V_flow_FerHofSqu1_chi_{chi_val}_t1_1_V_0_3_25_Coulomb_1_n_1_55_nphi_6_11_LxMUC_1_Ly_10.dat'
        ent_path = os.path.join(ent_dir, ent_file)

        # extract data from file
        with open(ent_path, 'r') as csvfile:
            plots = csv.reader(csvfile, delimiter='\t')
            V = []
            xi = []
            for row in plots:
                V.append(float(row[0]))
                xi.append(float(row[1]))

        ax3.plot(V, xi, c=f"C{i}", marker=markers[i], markersize=2, label=f"$\chi={chi_val}$")

    ax3.axvline(0.00678, c='b', linestyle=':', zorder=-2)
    ax3.axvline(0.300, c='r', linestyle=':', zorder=-2)
    ax3.grid(color='k', linestyle='-', linewidth=0.3, axis='x')
    ax3.set_xlim([0, 3])
    ax3.set_xticks(np.arange(0, 3.1, 1))
    ax3.set_yticks(np.arange(0, 3, 1))
    ax3.set_ylabel("$S_\mathrm{vN}$", fontsize=11)
    plt.setp(ax3.get_xticklabels(), visible=False)
    ax3.tick_params(
        axis='x',  # changes apply to the x-axis
        which='both',  # both major and minor ticks are affected
        bottom=False,  # ticks along the bottom edge are off
        top=False,  # ticks along the top edge are off
        labelbottom=False)

    ####################################################################################################################

    ax4 = plt.subplot(left_inner_grid[2])

    ent_spec_dir = '/home/bart/PycharmProjects/infinite_cylinder/data/ent_spec_V_flow/FerHofSqu1'
    ent_spec_file = f'ent_spec_V_flow_FerHofSqu1_chi_500_t1_1_V_0_3_25_Coulomb_1_n_1_55_nphi_6_11_LxMUC_1_Ly_10.dat'
    ent_spec_path = os.path.join(ent_spec_dir, ent_spec_file)

    x = []
    y = []
    z = []

    with open(ent_spec_path, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        for row in plots:
            x.append(float(row[1]))
            y.append(float(row[2]))
            z.append(int(row[0]))

    for value in np.linspace(3, -3, 7, dtype=int):
        xvalue = []
        yvalue = []
        for i in range(len(x)):
            if z[i] == value:
                xvalue.append(x[i])
                yvalue.append(y[i])
        if value != 0:
            ax4.scatter(xvalue, yvalue, marker='_', c='C{}'.format((value + 4) % 10), label='${}$'.format(value), s=20)
        else:
            ax4.scatter(xvalue, yvalue, marker='x', c='C{}'.format((value + 4) % 10), label='${}$'.format(value), s=20)

    ax4.legend(loc='center left', handletextpad=0, handlelength=1, borderpad=0.2, framealpha=1, edgecolor=None,
               markerscale=1,
               fontsize=10, ncol=1, bbox_to_anchor=(1, 0.5))

    ax4.tick_params(axis="x", labelsize=10)
    ax4.tick_params(axis="y", labelsize=10)

    ax4.axvline(0.00678, c='b', linestyle=':', zorder=-2)
    ax4.axvline(0.300, c='r', linestyle=':', zorder=-2)
    ax4.grid(color='k', linestyle='-', linewidth=0.3, axis='x')
    ax4.set_xlim([0, 3])
    ax4.set_xticks(np.arange(0, 3.1, 1))
    ax4.set_ylim([0, 10])
    ax4.set_yticks(np.arange(0, 10, 2))
    ax4.xaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax4.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax4.set_xlabel("$V$", fontsize=11)
    ax4.set_ylabel("$\\epsilon_\\alpha$", fontsize=11)
    ax4.text(2.02, 9.2, '$\chi=500$', fontsize=10, verticalalignment='top',
             bbox=dict(boxstyle='round', facecolor='white', alpha=1))

    ####################################################################################################################
    ####################################################################################################################

    # ax5 = plt.subplot(upper_right_cell)
    #
    # for i, chi_val in enumerate(np.arange(100, 400, 100)):
    #     corr_len_dir = '/home/bart/PycharmProjects/infinite_cylinder/data/corr_len_V_flow/FerHofSqu1'
    #     corr_len_file = f'corr_len_V_flow_FerHofSqu1_chi_{chi_val}_t1_1_V_0_10_41_Coulomb_1_n_2_81_nphi_5_9_LxMUC_1_Ly_9.dat'
    #     corr_len_path = os.path.join(corr_len_dir, corr_len_file)
    #
    #     # extract data from file
    #     with open(corr_len_path, 'r') as csvfile:
    #         plots = csv.reader(csvfile, delimiter='\t')
    #         V = []
    #         xi = []
    #         for row in plots:
    #             V.append(float(row[0]))
    #             xi.append(float(row[1]))
    #
    #     ax5.plot(V, xi, c=f"C{i}", marker=markers[i], markersize=2, label=f"$\chi={chi_val}$")
    #
    # ax5.legend(loc='upper right', handletextpad=0, handlelength=1, borderpad=0, framealpha=1, edgecolor='w',
    #            markerscale=1, fontsize=10, ncol=2)
    #
    # ax5.set_xlim([0, 10])
    # ax5.set_xlabel("$V$", fontsize=11)
    # ax5.xaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    # ax5.set_ylabel("$\\xi$", fontsize=11)
    # ax5.set_title("(b) $\\nu=2/9$, $n_\phi=5/9$", fontsize=12)
    # ax5.axvspan(0, 0.4, alpha=0.5, color='grey')
    # line1 = plt.Line2D((.59, .59), (.645, .748), color="k", linewidth=1, linestyle='--', alpha=0.5)
    # fig.add_artist(line1)
    # line2 = plt.Line2D((.9, .602), (.645, .748), color="k", linewidth=1, linestyle='--', alpha=0.5)
    # fig.add_artist(line2)

    ####################################################################################################################

    ax6 = plt.subplot(right_inner_grid[0])

    for i, chi_val in enumerate(np.arange(100, 600, 100)):
        corr_len_dir = '/home/bart/PycharmProjects/infinite_cylinder/data/corr_len_V_flow/FerHofSqu1'
        corr_len_file = f'corr_len_V_flow_FerHofSqu1_chi_{chi_val}_t1_1_V_0_3_25_Coulomb_1_n_2_81_nphi_5_9_LxMUC_1_Ly_9.dat'
        corr_len_path = os.path.join(corr_len_dir, corr_len_file)

        # extract data from file
        with open(corr_len_path, 'r') as csvfile:
            plots = csv.reader(csvfile, delimiter='\t')
            V = []
            xi = []
            for row in plots:
                V.append(float(row[0]))
                xi.append(float(row[1]))

        ax6.plot(V, xi, c=f"C{i}", marker=markers[i], markersize=2, label=f"${chi_val/100:g}$")

    ax6.legend(loc='center left', handletextpad=0.5, handlelength=1, borderpad=0.2, framealpha=1, edgecolor=None,
               markerscale=2, fontsize=10, ncol=1, bbox_to_anchor=(1, 0), title="$\chi/10^2$")

    ax6.axvline(0.0230, c='b', linestyle=':', zorder=-2)
    ax6.axvline(0.309, c='r', linestyle=':', zorder=-2)
    ax6.grid(color='k', linestyle='-', linewidth=0.3, axis='x')
    ax6.set_xlim([0, 3])
    ax6.set_xticks(np.arange(0, 3.1, 1))
    ax6.set_yticks(np.arange(0, 350, 100))
    ax6.set_ylim([0, 350])
    ax6.set_ylabel("$\\xi$", fontsize=11)
    plt.setp(ax6.get_xticklabels(), visible=False)
    ax6.tick_params(
        axis='x',  # changes apply to the x-axis
        which='both',  # both major and minor ticks are affected
        bottom=False,  # ticks along the bottom edge are off
        top=False,  # ticks along the top edge are off
        labelbottom=False)
    ax6.set_title("(b) $\\nu=2/9$, $n_\phi=5/9$", fontsize=12)

    ####################################################################################################################

    ax7 = plt.subplot(right_inner_grid[1])

    for i, chi_val in enumerate(np.arange(100, 600, 100)):
        ent_dir = '/home/bart/PycharmProjects/infinite_cylinder/data/ent_V_flow/FerHofSqu1'
        ent_file = f'ent_V_flow_FerHofSqu1_chi_{chi_val}_t1_1_V_0_3_25_Coulomb_1_n_2_81_nphi_5_9_LxMUC_1_Ly_9.dat'
        ent_path = os.path.join(ent_dir, ent_file)

        # extract data from file
        with open(ent_path, 'r') as csvfile:
            plots = csv.reader(csvfile, delimiter='\t')
            V = []
            xi = []
            for row in plots:
                V.append(float(row[0]))
                xi.append(float(row[1]))

        ax7.plot(V, xi, c=f"C{i}", marker=markers[i], markersize=2, label=f"$\chi={chi_val}$")

    ax7.axvline(0.0230, c='b', linestyle=':', zorder=-2)
    ax7.axvline(0.309, c='r', linestyle=':', zorder=-2)
    ax7.grid(color='k', linestyle='-', linewidth=0.3, axis='x')
    ax7.set_xlim([0, 3])
    ax7.set_xticks(np.arange(0, 3.1, 1))
    ax7.set_yticks(np.arange(0, 3, 1))
    ax7.set_ylabel("$S_\mathrm{vN}$", fontsize=11)
    plt.setp(ax7.get_xticklabels(), visible=False)
    ax7.tick_params(
        axis='x',  # changes apply to the x-axis
        which='both',  # both major and minor ticks are affected
        bottom=False,  # ticks along the bottom edge are off
        top=False,  # ticks along the top edge are off
        labelbottom=False)

    ####################################################################################################################

    ax8 = plt.subplot(right_inner_grid[2])

    ent_spec_dir = '/home/bart/PycharmProjects/infinite_cylinder/data/ent_spec_V_flow/FerHofSqu1'
    ent_spec_file = f'ent_spec_V_flow_FerHofSqu1_chi_500_t1_1_V_0_3_25_Coulomb_1_n_2_81_nphi_5_9_LxMUC_1_Ly_9.dat'
    ent_spec_path = os.path.join(ent_spec_dir, ent_spec_file)

    x = []
    y = []
    z = []

    with open(ent_spec_path, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        for row in plots:
            x.append(float(row[1]))
            y.append(float(row[2]))
            z.append(int(row[0]))

    for value in np.linspace(3, -3, 7, dtype=int):
        xvalue = []
        yvalue = []
        for i in range(len(x)):
            if z[i] == value:
                xvalue.append(x[i])
                yvalue.append(y[i])
        if value != 0:
            ax8.scatter(xvalue, yvalue, marker='_', c='C{}'.format((value + 4) % 10), label='${}$'.format(value), s=20)
        else:
            ax8.scatter(xvalue, yvalue, marker='x', c='C{}'.format((value + 4) % 10), label='${}$'.format(value), s=20)

    ax8.legend(loc='center left', handletextpad=0, handlelength=1, borderpad=0.2, framealpha=1, edgecolor=None,
               markerscale=1,
               fontsize=10, ncol=1, bbox_to_anchor=(1, 0.5))

    ax8.tick_params(axis="x", labelsize=10)
    ax8.tick_params(axis="y", labelsize=10)

    ax8.axvline(0.0230, c='b', linestyle=':', zorder=-2)
    ax8.axvline(0.309, c='r', linestyle=':', zorder=-2)
    ax8.grid(color='k', linestyle='-', linewidth=0.3, axis='x')
    ax8.set_xlim([0, 3])
    ax8.set_xticks(np.arange(0, 3.1, 1))
    ax8.set_ylim([0, 10])
    ax8.set_yticks(np.arange(0, 10, 2))
    ax8.xaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax8.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax8.set_xlabel("$V$", fontsize=11)
    ax8.set_ylabel("$\\epsilon_\\alpha$", fontsize=11)
    ax8.text(2.02, 9.2, '$\chi=500$', fontsize=10, verticalalignment='top',
             bbox=dict(boxstyle='round', facecolor='white', alpha=1))
    
    ####################################################################################################################

    # fig.text(0.03, 0.87, "(a)", fontsize=12)
    # fig.text(0.5, 0.87, "(b)", fontsize=12)

    plt.savefig("/home/bart/Documents/papers/FCI/Vflow_c2345_analysis_25.png", bbox_inches='tight', dpi=300)
    plt.show()
