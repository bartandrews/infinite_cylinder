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

    fig = plt.figure(figsize=(6, 5))
    outer_grid = gridspec.GridSpec(2, 2, wspace=0.5, hspace=0.3, height_ratios=[1, 4])
    upper_left_cell = outer_grid[0, 0]
    upper_right_cell = outer_grid[0, 1]
    lower_left_cell = outer_grid[1, 0]
    lower_right_cell = outer_grid[1, 1]
    left_inner_grid = gridspec.GridSpecFromSubplotSpec(3, 1, lower_left_cell, hspace=0, height_ratios=[1, 1, 2])
    right_inner_grid = gridspec.GridSpecFromSubplotSpec(3, 1, lower_right_cell, hspace=0, height_ratios=[1, 1, 2])

    # define a list of easily-visible markers
    markers = [(3, 0, 0), (4, 0, 0), (5, 0, 0), (6, 0, 0), (4, 1, 0), (5, 1, 0), (6, 1, 0),
               (3, 2, 0), (4, 2, 0), (5, 2, 0), (6, 2, 0), 'X', 'x', 'd', 'D', 'P',
               '$a$', '$b$', '$c$', '$d$', '$e$', '$f$', '$g$', '$h$', '$i$', '$j$', '$k$', '$l$', '$m$', '$n$', '$o$',
               '$p$', '$q$', '$r$', '$s$', '$t$', '$u$', '$v$', '$w$', '$x$', '$y$', '$z$']

    ####################################################################################################################

    ax1 = plt.subplot(upper_left_cell)

    corr_len_dir = '/home/bart/PycharmProjects/infinite_cylinder/data/corr_len_V_flow/FerHofSqu1'
    corr_len_file = f'corr_len_V_flow_FerHofSqu1_chi_25_t1_1_V_0_10_41_Coulomb_1_n_1_33_nphi_6_11_LxMUC_1_Ly_6.dat'
    corr_len_path = os.path.join(corr_len_dir, corr_len_file)

    # extract data from file
    with open(corr_len_path, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        V = []
        xi = []
        for row in plots:
            V.append(float(row[0]))
            xi.append(float(row[1]))

    ax1.plot(V, xi, 'x', c='k', markersize=2)

    ax1.set_xlim([0, 10])
    ax1.set_xlabel("$V$", fontsize=11)
    ax1.xaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax1.set_ylabel("$\\xi$", fontsize=11)
    ax1.set_title("(a) $C=2$, $\\nu=1/3$", fontsize=12)
    ax1.axvspan(0, 0.4, alpha=0.5, color='grey')
    line1 = plt.Line2D((.125, .125), (.645, .748), color="k", linewidth=1, linestyle='--', alpha=0.5)
    fig.add_artist(line1)
    line2 = plt.Line2D((.435, .139), (.645, .748), color="k", linewidth=1, linestyle='--', alpha=0.5)
    fig.add_artist(line2)

    ####################################################################################################################

    ax2 = plt.subplot(left_inner_grid[0])

    corr_len_dir = '/home/bart/PycharmProjects/infinite_cylinder/data/corr_len_V_flow/FerHofSqu1'
    corr_len_file = f'corr_len_V_flow_FerHofSqu1_chi_25_t1_1_V_0_10_41_Coulomb_1_n_1_33_nphi_6_11_LxMUC_1_Ly_6.dat'
    corr_len_path = os.path.join(corr_len_dir, corr_len_file)

    # extract data from file
    with open(corr_len_path, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        V = []
        xi = []
        for row in plots:
            V.append(float(row[0]))
            xi.append(float(row[1]))

    ax2.plot(V, xi, 'x', c='k', markersize=2)

    ax2.set_xlim([0, 10])
    ax2.set_ylabel("$\\xi$", fontsize=11)
    plt.setp(ax2.get_xticklabels(), visible=False)
    ax2.tick_params(
        axis='x',  # changes apply to the x-axis
        which='both',  # both major and minor ticks are affected
        bottom=False,  # ticks along the bottom edge are off
        top=False,  # ticks along the top edge are off
        labelbottom=False)
    
    ####################################################################################################################

    ax3 = plt.subplot(left_inner_grid[1])

    ent_dir = '/home/bart/PycharmProjects/infinite_cylinder/data/ent_V_flow/FerHofSqu1'
    ent_file = f'ent_V_flow_FerHofSqu1_chi_25_t1_1_V_0_10_41_Coulomb_1_n_1_33_nphi_6_11_LxMUC_1_Ly_6.dat'
    ent_path = os.path.join(ent_dir, ent_file)

    # extract data from file
    with open(ent_path, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        V = []
        xi = []
        for row in plots:
            V.append(float(row[0]))
            xi.append(float(row[1]))

    ax3.plot(V, xi, 'x', c='k', markersize=2)

    ax3.set_xlim([0, 10])
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
    ent_spec_file = f'ent_spec_V_flow_FerHofSqu1_chi_25_t1_1_V_0_10_41_Coulomb_1_n_1_33_nphi_6_11_LxMUC_1_Ly_6.dat'
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
            ax4.scatter(xvalue, yvalue, marker='_', c='C{}'.format((value + 4) % 10), label='{}'.format(value), s=20)
        else:
            ax4.scatter(xvalue, yvalue, marker='x', c='C{}'.format((value + 4) % 10), label='{}'.format(value), s=20)

    ax4.legend(loc='center left', handletextpad=0, handlelength=1, borderpad=0.2, framealpha=1, edgecolor='k',
               markerscale=1,
               fontsize=10, ncol=1, bbox_to_anchor=(1, 0.5))

    ax4.tick_params(axis="x", labelsize=10)
    ax4.tick_params(axis="y", labelsize=10)

    ax4.set_xlim([0, 10])
    ax4.set_ylim([0, 10])
    ax4.xaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax4.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax4.set_xlabel("$V$", fontsize=11)
    ax4.set_ylabel("$\\epsilon_\\alpha$", fontsize=11)

    ####################################################################################################################
    ####################################################################################################################

    ax5 = plt.subplot(upper_right_cell)

    corr_len_dir = '/home/bart/PycharmProjects/infinite_cylinder/data/corr_len_V_flow/FerHofSqu1'
    corr_len_file = f'corr_len_V_flow_FerHofSqu1_chi_50_t1_1_V_0_10_41_Coulomb_1_n_1_55_nphi_6_11_LxMUC_1_Ly_10.dat'
    corr_len_path = os.path.join(corr_len_dir, corr_len_file)

    # extract data from file
    with open(corr_len_path, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        V = []
        xi = []
        for row in plots:
            V.append(float(row[0]))
            xi.append(float(row[1]))

    ax5.plot(V, xi, 'x', c='k', markersize=2)

    ax5.set_xlim([0, 10])
    ax5.set_xlabel("$V$", fontsize=11)
    ax5.xaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax5.set_ylabel("$\\xi$", fontsize=11)
    ax5.set_title("(b) $C=2$, $\\nu=1/5$", fontsize=12)
    ax5.axvspan(0, 0.4, alpha=0.5, color='grey')
    line1 = plt.Line2D((.59, .59), (.645, .748), color="k", linewidth=1, linestyle='--', alpha=0.5)
    fig.add_artist(line1)
    line2 = plt.Line2D((.9, .602), (.645, .748), color="k", linewidth=1, linestyle='--', alpha=0.5)
    fig.add_artist(line2)

    ####################################################################################################################

    ax6 = plt.subplot(right_inner_grid[0])

    corr_len_dir = '/home/bart/PycharmProjects/infinite_cylinder/data/corr_len_V_flow/FerHofSqu1'
    corr_len_file = f'corr_len_V_flow_FerHofSqu1_chi_50_t1_1_V_0_10_41_Coulomb_1_n_1_55_nphi_6_11_LxMUC_1_Ly_10.dat'
    corr_len_path = os.path.join(corr_len_dir, corr_len_file)

    # extract data from file
    with open(corr_len_path, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        V = []
        xi = []
        for row in plots:
            V.append(float(row[0]))
            xi.append(float(row[1]))

    ax6.plot(V, xi, 'x', c='k', markersize=2)

    ax6.set_xlim([0, 10])
    # ax6.set_ylim([5, 10])
    ax6.set_ylabel("$\\xi$", fontsize=11)
    plt.setp(ax6.get_xticklabels(), visible=False)
    ax6.tick_params(
        axis='x',  # changes apply to the x-axis
        which='both',  # both major and minor ticks are affected
        bottom=False,  # ticks along the bottom edge are off
        top=False,  # ticks along the top edge are off
        labelbottom=False)

    ####################################################################################################################

    ax7 = plt.subplot(right_inner_grid[1])

    ent_dir = '/home/bart/PycharmProjects/infinite_cylinder/data/ent_V_flow/FerHofSqu1'
    ent_file = f'ent_V_flow_FerHofSqu1_chi_50_t1_1_V_0_10_41_Coulomb_1_n_1_55_nphi_6_11_LxMUC_1_Ly_10.dat'
    ent_path = os.path.join(ent_dir, ent_file)

    # extract data from file
    with open(ent_path, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        V = []
        xi = []
        for row in plots:
            V.append(float(row[0]))
            xi.append(float(row[1]))

    ax7.plot(V, xi, 'x', c='k', markersize=2)

    ax7.set_xlim([0, 10])
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
    ent_spec_file = f'ent_spec_V_flow_FerHofSqu1_chi_50_t1_1_V_0_10_41_Coulomb_1_n_1_55_nphi_6_11_LxMUC_1_Ly_10.dat'
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
            ax8.scatter(xvalue, yvalue, marker='_', c='C{}'.format((value + 4) % 10), label='{}'.format(value), s=20)
        else:
            ax8.scatter(xvalue, yvalue, marker='x', c='C{}'.format((value + 4) % 10), label='{}'.format(value), s=20)

    ax8.legend(loc='center left', handletextpad=0, handlelength=1, borderpad=0.2, framealpha=1, edgecolor='k',
               markerscale=1,
               fontsize=10, ncol=1, bbox_to_anchor=(1, 0.5))

    ax8.tick_params(axis="x", labelsize=10)
    ax8.tick_params(axis="y", labelsize=10)

    ax8.set_xlim([0, 10])
    ax8.set_ylim([0, 10])
    ax8.xaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax8.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax8.set_xlabel("$V$", fontsize=11)
    ax8.set_ylabel("$\\epsilon_\\alpha$", fontsize=11)
    
    ####################################################################################################################

    # fig.text(0.03, 0.87, "(a)", fontsize=12)
    # fig.text(0.5, 0.87, "(b)", fontsize=12)

    plt.savefig("/home/bart/Documents/papers/FCI/Vflow_c2345_analysis.png", bbox_inches='tight', dpi=300)
    plt.show()
