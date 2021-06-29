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

    fig = plt.figure(figsize=(6, 9))
    outer_grid = gridspec.GridSpec(3, 2, wspace=0.5, hspace=0.3, height_ratios=[1, 5, 2])
    upper_left_cell = outer_grid[0, 0]
    upper_right_cell = outer_grid[0, 1]
    middle_left_cell = outer_grid[1, 0]
    middle_right_cell = outer_grid[1, 1]
    lower_left_cell = outer_grid[2, 0]
    lower_right_cell = outer_grid[2, 1]
    middle_left_inner_grid = gridspec.GridSpecFromSubplotSpec(4, 1, middle_left_cell, hspace=0, height_ratios=[1, 1, 1, 2])
    middle_right_inner_grid = gridspec.GridSpecFromSubplotSpec(4, 1, middle_right_cell, hspace=0, height_ratios=[1, 1, 1, 2])

    # define a list of easily-visible markers
    markers = [(3, 0, 0), (4, 0, 0), (5, 0, 0), (6, 0, 0), (4, 1, 0), (5, 1, 0), (6, 1, 0),
               (3, 2, 0), (4, 2, 0), (5, 2, 0), (6, 2, 0), 'X', 'x', 'd', 'D', 'P',
               '$a$', '$b$', '$c$', '$d$', '$e$', '$f$', '$g$', '$h$', '$i$', '$j$', '$k$', '$l$', '$m$', '$n$', '$o$',
               '$p$', '$q$', '$r$', '$s$', '$t$', '$u$', '$v$', '$w$', '$x$', '$y$', '$z$']

    ####################################################################################################################

    ax1 = plt.subplot(upper_left_cell)

    for i, chi_val in enumerate(np.arange(200, 1100, 200)):
        corr_len_dir = '/home/bart/PycharmProjects/infinite_cylinder/data/corr_len_V_flow/FerHofSqu1'
        corr_len_file = f'corr_len_V_flow_FerHofSqu1_chi_{chi_val}_t1_1_V_0_10_41_Coulomb_1_n_1_12_nphi_1_4_LxMUC_1_Ly_6.dat'
        corr_len_path = os.path.join(corr_len_dir, corr_len_file)

        # extract data from file
        with open(corr_len_path, 'r') as csvfile:
            plots = csv.reader(csvfile, delimiter='\t')
            V = []
            xi = []
            for row in plots:
                V.append(float(row[0]))
                xi.append(float(row[1]))

        ax1.plot(V, xi, c=f"C{2*i+1}", marker=markers[i], markersize=2, label=f"${chi_val/100:g}$")

    ax1.legend(loc='upper right', handletextpad=0, handlelength=1, borderpad=0, framealpha=1, edgecolor='w',
               markerscale=2, fontsize=10, ncol=5, columnspacing=0.5, labelspacing=0.25, title="$\chi/10^2$")

    ax1.axvline(0.215, c='b', linestyle=':', zorder=-2)
    ax1.axvline(1.53, c='r', linestyle=':', zorder=-2)
    ax1.set_xlim([0, 10])
    ax1.set_xticks(np.arange(0, 10.1, 2))
    ax1.set_xlabel("$V$", fontsize=11)
    ax1.xaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax1.set_ylabel("$\\xi$", fontsize=11)
    ax1.set_title("(a) $n_\phi=1/4$", fontsize=12)
    ax1.axvspan(0, 0.6, alpha=0.5, color='grey')
    line1 = plt.Line2D((.125, .125), (.735, .8), color="k", linewidth=1, linestyle='--', alpha=0.5)
    fig.add_artist(line1)
    line2 = plt.Line2D((.435, .144), (.735, .8), color="k", linewidth=1, linestyle='--', alpha=0.5)
    fig.add_artist(line2)

    ####################################################################################################################

    ax2 = plt.subplot(middle_left_inner_grid[0])

    for i, chi_val in enumerate(np.arange(200, 1100, 200)):
        corr_len_dir = '/home/bart/PycharmProjects/infinite_cylinder/data/corr_len_V_flow/FerHofSqu1'
        corr_len_file = f'corr_len_V_flow_FerHofSqu1_chi_{chi_val}_t1_1_V_0_1_41_Coulomb_1_n_1_12_nphi_1_4_LxMUC_1_Ly_6.dat'
        corr_len_path = os.path.join(corr_len_dir, corr_len_file)

        # extract data from file
        with open(corr_len_path, 'r') as csvfile:
            plots = csv.reader(csvfile, delimiter='\t')
            V = []
            xi = []
            for row in plots:
                V.append(float(row[0]))
                xi.append(float(row[1]))

        ax2.plot(V, xi, c=f"C{2*i+1}", marker=markers[i], markersize=2, label=f"$\chi={chi_val}$")

    ax2.axvline(0.215, c='b', linestyle=':', zorder=-2)
    ax2.grid(color='k', linestyle='-', linewidth=0.3, axis='x')
    ax2.set_xlim([0, 0.6])
    ax2.set_xticks(np.arange(0, 0.61, 0.2))
    ax2.set_yticks(np.arange(0, 200.1, 50))
    ax2.set_ylabel("$\\xi$", fontsize=11)
    plt.setp(ax2.get_xticklabels(), visible=False)
    ax2.tick_params(
        axis='x',  # changes apply to the x-axis
        which='both',  # both major and minor ticks are affected
        bottom=False,  # ticks along the bottom edge are off
        top=False,  # ticks along the top edge are off
        labelbottom=False)
    
    ####################################################################################################################

    ax3 = plt.subplot(middle_left_inner_grid[1])

    for i, chi_val in enumerate(np.arange(200, 1100, 200)):
        ent_dir = '/home/bart/PycharmProjects/infinite_cylinder/data/ent_V_flow/FerHofSqu1'
        ent_file = f'ent_V_flow_FerHofSqu1_chi_{chi_val}_t1_1_V_0_1_41_Coulomb_1_n_1_12_nphi_1_4_LxMUC_1_Ly_6.dat'
        ent_path = os.path.join(ent_dir, ent_file)

        # extract data from file
        with open(ent_path, 'r') as csvfile:
            plots = csv.reader(csvfile, delimiter='\t')
            V = []
            SvN = []
            for row in plots:
                V.append(float(row[0]))
                SvN.append(float(row[1]))

        ax3.plot(V, SvN, c=f"C{2*i+1}", marker=markers[i], markersize=2, label=f"$\chi={chi_val}$")

    ax3.axvline(0.215, c='b', linestyle=':', zorder=-2)
    ax3.grid(color='k', linestyle='-', linewidth=0.3, axis='x')
    ax3.set_xlim([0, 0.6])
    ax3.set_xticks(np.arange(0, 0.61, 0.2))
    ax3.set_yticks(np.arange(0, 4, 1))
    ax3.set_ylabel("$S_\mathrm{vN}$", fontsize=11)
    plt.setp(ax3.get_xticklabels(), visible=False)
    ax3.tick_params(
        axis='x',  # changes apply to the x-axis
        which='both',  # both major and minor ticks are affected
        bottom=False,  # ticks along the bottom edge are off
        top=False,  # ticks along the top edge are off
        labelbottom=False)

    ####################################################################################################################

    ax4 = plt.subplot(middle_left_inner_grid[2])

    for i, chi_val in enumerate(np.arange(200, 1100, 200)):

        if chi_val == 200:
            V_min = 0
            V_samples = 25
        else:
            V_min = 0.025
            V_samples = 24

        energy_dir = '/home/bart/PycharmProjects/infinite_cylinder/data/energy_V_flow/FerHofSqu1'
        energy_file = f'energy_V_flow_FerHofSqu1_chi_{chi_val}_t1_1_V_{V_min}_0.6_{V_samples}_Coulomb_1_n_1_12_nphi_1_4_LxMUC_1_Ly_6.dat'
        energy_path = os.path.join(energy_dir, energy_file)

        # extract data from file
        with open(energy_path, 'r') as csvfile:
            plots = csv.reader(csvfile, delimiter='\t')
            V = []
            kin = []
            pot = []
            tot = []
            for row in plots:
                V.append(float(row[0]))
                kin.append(float(row[1]))
                pot.append(float(row[2]))
                tot.append(float(row[3]))

        if V_min == 0:
            ax4.plot(V[1:], np.log(pot[1:]), c=f"C{2 * i + 1}", marker=markers[i], markersize=2,
                     label=f"${chi_val / 100:g}$")
        else:
            ax4.plot(V, np.log(pot), c=f"C{2 * i + 1}", marker=markers[i], markersize=2,
                     label=f"${chi_val / 100:g}$")

    ax4.axvline(0.215, c='b', linestyle=':', zorder=-2)
    ax4.grid(color='k', linestyle='-', linewidth=0.3, axis='x')
    ax4.set_xlim([0, 0.6])
    ax4.set_xticks(np.arange(0, 0.61, 0.2))
    # ax4.set_yticks(np.arange(0, 200.1, 50))
    ax4.set_ylabel("$\ln \langle \hat{V} \\rangle$", fontsize=11)
    plt.setp(ax4.get_xticklabels(), visible=False)
    ax4.tick_params(
        axis='x',  # changes apply to the x-axis
        which='both',  # both major and minor ticks are affected
        bottom=False,  # ticks along the bottom edge are off
        top=False,  # ticks along the top edge are off
        labelbottom=False)
    
    ####################################################################################################################

    ax5 = plt.subplot(middle_left_inner_grid[3])

    ent_spec_dir = '/home/bart/PycharmProjects/infinite_cylinder/data/ent_spec_V_flow/FerHofSqu1'
    ent_spec_file = f'ent_spec_V_flow_FerHofSqu1_chi_1000_t1_1_V_0_1_41_Coulomb_1_n_1_12_nphi_1_4_LxMUC_1_Ly_6.dat'
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
            ax5.scatter(xvalue, yvalue, marker='_', c='C{}'.format((value + 4) % 10), label='${}$'.format(value), s=20)
        else:
            ax5.scatter(xvalue, yvalue, marker='x', c='C{}'.format((value + 4) % 10), label='${}$'.format(value), s=20)

    ax5.legend(loc='center left', handletextpad=0, handlelength=1, borderpad=0.2, framealpha=1, edgecolor=None,
               markerscale=1,
               fontsize=10, ncol=1, bbox_to_anchor=(1, 0.5))

    ax5.tick_params(axis="x", labelsize=10)
    ax5.tick_params(axis="y", labelsize=10)

    ax5.axvline(0.215, c='b', linestyle=':', zorder=-2)
    ax5.grid(color='k', linestyle='-', linewidth=0.3, axis='x')
    ax5.set_xlim([0, 0.6])
    ax5.set_xticks(np.arange(0, 0.61, 0.2))
    ax5.set_yticks(np.arange(0, 10, 2))
    ax5.set_ylim([0, 10])
    ax5.xaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax5.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax5.set_xlabel("$V$", fontsize=11)
    ax5.set_ylabel("$\\epsilon_\\alpha$", fontsize=11)
    ax5.text(0.415, 9.2, '$\chi=10^3$', fontsize=10, verticalalignment='top',
             bbox=dict(boxstyle='round', facecolor='white', alpha=1))

    ####################################################################################################################

    ax6 = plt.subplot(lower_left_cell)

    ent_spec_mom_dir = '/home/bart/PycharmProjects/infinite_cylinder/data/ent_spec_mom/FerHofSqu1'
    ent_spec_mom_file = 'ent_spec_mom_FerHofSqu1_chi_1000_chiK_1000_t1_1_V_0.6_Coulomb_1_n_1_12_nphi_1_4_LxMUC_1_Ly_6.dat'
    ent_spec_mom_path = os.path.join(ent_spec_mom_dir, ent_spec_mom_file)

    x = []
    y = []
    z = []

    with open(ent_spec_mom_path, 'r') as csvfile:
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
            ax6.scatter(xvalue, yvalue, marker='_', c='C{}'.format((value + 4) % 10), label='${}$'.format(value))
        else:
            ax6.scatter(xvalue, yvalue, marker='x', c='C{}'.format((value + 4) % 10), label='${}$'.format(value))

    ax6.set_yticks(np.arange(0, 15.1, 5))
    ax6.set_ylim([0, 15])

    # ax6.legend(loc='center left', handletextpad=0, handlelength=1, borderpad=0.2, framealpha=1, edgecolor=None,
    #            markerscale=1,
    #            fontsize=10, ncol=1, bbox_to_anchor=(1, 0.5))
    ax6.set_xlim([-np.pi / 3, np.pi / 3])
    ax6.set_xlabel("$k_a / \pi$", fontsize=11)
    ax6.set_ylabel("$\epsilon_{\\alpha}$", fontsize=11)

    # ax6.plot([-1, 0.18], [0, 14.2], color='k', linestyle=':', linewidth=1)
    # ax6.plot([-0.73, 0.18], [0, 6], color='k', linestyle=':', linewidth=1)

    ax6.tick_params(axis="x", labelsize=10)
    ax6.tick_params(axis="y", labelsize=10)

    ####################################################################################################################
    ####################################################################################################################

    ax7 = plt.subplot(upper_right_cell)

    for i, chi_val in enumerate(np.arange(200, 1100, 200)):
        corr_len_dir = '/home/bart/PycharmProjects/infinite_cylinder/data/corr_len_V_flow/FerHofSqu1'
        corr_len_file = f'corr_len_V_flow_FerHofSqu1_chi_{chi_val}_t1_1_V_0_10_41_Coulomb_1_n_1_18_nphi_1_6_LxMUC_1_Ly_6.dat'
        corr_len_path = os.path.join(corr_len_dir, corr_len_file)

        # extract data from file
        with open(corr_len_path, 'r') as csvfile:
            plots = csv.reader(csvfile, delimiter='\t')
            V = []
            xi = []
            for row in plots:
                V.append(float(row[0]))
                xi.append(float(row[1]))

        ax7.plot(V, xi, c=f"C{2*i+1}", marker=markers[i], markersize=2, label=f"${chi_val/100:g}$")

    ax7.legend(loc='upper right', handletextpad=0, handlelength=1, borderpad=0, framealpha=1, edgecolor='w',
               markerscale=2, fontsize=10, ncol=5, columnspacing=0.5, labelspacing=0.25, title="$\chi/10^2$")

    ax7.axvline(0.0192, c='b', linestyle=':', zorder=-2)
    ax7.axvline(1.48, c='r', linestyle=':', zorder=-2)
    ax7.set_xlim([0, 10])
    ax7.set_xticks(np.arange(0, 10.1, 2))
    ax7.set_xlabel("$V$", fontsize=11)
    ax7.xaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax7.set_ylabel("$\\xi$", fontsize=11)
    ax7.set_title("(b) $n_\phi=1/6$", fontsize=12)
    ax7.axvspan(0, 0.6, alpha=0.5, color='grey')
    line1 = plt.Line2D((.59, .59), (.735, .8), color="k", linewidth=1, linestyle='--', alpha=0.5)
    fig.add_artist(line1)
    line2 = plt.Line2D((.9, .61), (.735, .8), color="k", linewidth=1, linestyle='--', alpha=0.5)
    fig.add_artist(line2)

    ####################################################################################################################

    ax8 = plt.subplot(middle_right_inner_grid[0])

    for i, chi_val in enumerate(np.arange(200, 1100, 200)):
        corr_len_dir = '/home/bart/PycharmProjects/infinite_cylinder/data/corr_len_V_flow/FerHofSqu1'
        corr_len_file = f'corr_len_V_flow_FerHofSqu1_chi_{chi_val}_t1_1_V_0_1_41_Coulomb_1_n_1_18_nphi_1_6_LxMUC_1_Ly_6.dat'
        corr_len_path = os.path.join(corr_len_dir, corr_len_file)

        # extract data from file
        with open(corr_len_path, 'r') as csvfile:
            plots = csv.reader(csvfile, delimiter='\t')
            V = []
            xi = []
            for row in plots:
                V.append(float(row[0]))
                xi.append(float(row[1]))

        ax8.plot(V, xi, c=f"C{2*i+1}", marker=markers[i], markersize=2, label=f"$\chi={chi_val}$")

    ax8.axvline(0.0192, c='b', linestyle=':', zorder=-2)
    ax8.grid(color='k', linestyle='-', linewidth=0.3, axis='x')
    ax8.set_xlim([0, 0.6])
    ax8.set_xticks(np.arange(0, 0.61, 0.2))
    ax8.set_yticks(np.arange(0, 200.1, 50))
    # ax8.set_ylim([5, 10])
    ax8.set_ylabel("$\\xi$", fontsize=11)
    plt.setp(ax8.get_xticklabels(), visible=False)
    ax8.tick_params(
        axis='x',  # changes apply to the x-axis
        which='both',  # both major and minor ticks are affected
        bottom=False,  # ticks along the bottom edge are off
        top=False,  # ticks along the top edge are off
        labelbottom=False)

    ####################################################################################################################

    ax9 = plt.subplot(middle_right_inner_grid[1])

    for i, chi_val in enumerate(np.arange(200, 1100, 200)):
        ent_dir = '/home/bart/PycharmProjects/infinite_cylinder/data/ent_V_flow/FerHofSqu1'
        ent_file = f'ent_V_flow_FerHofSqu1_chi_{chi_val}_t1_1_V_0_1_41_Coulomb_1_n_1_18_nphi_1_6_LxMUC_1_Ly_6.dat'
        ent_path = os.path.join(ent_dir, ent_file)

        # extract data from file
        with open(ent_path, 'r') as csvfile:
            plots = csv.reader(csvfile, delimiter='\t')
            V = []
            SvN = []
            for row in plots:
                V.append(float(row[0]))
                SvN.append(float(row[1]))

        ax9.plot(V, SvN, c=f"C{2*i+1}", marker=markers[i], markersize=2, label=f"$\chi={chi_val}$")

    ax9.axvline(0.0192, c='b', linestyle=':', zorder=-2)
    ax9.grid(color='k', linestyle='-', linewidth=0.3, axis='x')
    ax9.set_xlim([0, 0.6])
    ax9.set_xticks(np.arange(0, 0.61, 0.2))
    ax9.set_yticks(np.arange(0, 4, 1))
    ax9.set_ylabel("$S_\mathrm{vN}$", fontsize=11)
    plt.setp(ax9.get_xticklabels(), visible=False)
    ax9.tick_params(
        axis='x',  # changes apply to the x-axis
        which='both',  # both major and minor ticks are affected
        bottom=False,  # ticks along the bottom edge are off
        top=False,  # ticks along the top edge are off
        labelbottom=False)

    ####################################################################################################################

    ax10 = plt.subplot(middle_right_inner_grid[2])

    for i, chi_val in enumerate(np.arange(200, 1100, 200)):

        if chi_val == 200:
            V_min = 0
            V_samples = 25
        else:
            V_min = 0.025
            V_samples = 24

        energy_dir = '/home/bart/PycharmProjects/infinite_cylinder/data/energy_V_flow/FerHofSqu1'
        energy_file = f'energy_V_flow_FerHofSqu1_chi_{chi_val}_t1_1_V_{V_min}_0.6_{V_samples}_Coulomb_1_n_1_18_nphi_1_6_LxMUC_1_Ly_6.dat'
        energy_path = os.path.join(energy_dir, energy_file)

        # extract data from file
        with open(energy_path, 'r') as csvfile:
            plots = csv.reader(csvfile, delimiter='\t')
            V = []
            kin = []
            pot = []
            tot = []
            for row in plots:
                V.append(float(row[0]))
                kin.append(float(row[1]))
                pot.append(float(row[2]))
                tot.append(float(row[3]))

        if V_min == 0:
            ax10.plot(V[1:], np.log(pot[1:]), c=f"C{2 * i + 1}", marker=markers[i], markersize=2,
                     label=f"${chi_val / 100:g}$")
        else:
            ax10.plot(V, np.log(pot), c=f"C{2 * i + 1}", marker=markers[i], markersize=2,
                      label=f"${chi_val / 100:g}$")

    ax10.axvline(0.0192, c='b', linestyle=':', zorder=-2)
    ax10.grid(color='k', linestyle='-', linewidth=0.3, axis='x')
    ax10.set_xlim([0, 0.6])
    ax10.set_xticks(np.arange(0, 0.61, 0.2))
    # ax10.set_yticks(np.arange(0, 200.1, 50))
    ax10.set_ylabel("$\ln \langle \hat{V} \\rangle$", fontsize=11)
    plt.setp(ax10.get_xticklabels(), visible=False)
    ax10.tick_params(
        axis='x',  # changes apply to the x-axis
        which='both',  # both major and minor ticks are affected
        bottom=False,  # ticks along the bottom edge are off
        top=False,  # ticks along the top edge are off
        labelbottom=False)
    
    ####################################################################################################################

    ax11 = plt.subplot(middle_right_inner_grid[3])

    ent_spec_dir = '/home/bart/PycharmProjects/infinite_cylinder/data/ent_spec_V_flow/FerHofSqu1'
    ent_spec_file = f'ent_spec_V_flow_FerHofSqu1_chi_1000_t1_1_V_0_1_41_Coulomb_1_n_1_18_nphi_1_6_LxMUC_1_Ly_6.dat'
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
            ax11.scatter(xvalue, yvalue, marker='_', c='C{}'.format((value + 4) % 10), label='${}$'.format(value), s=20)
        else:
            ax11.scatter(xvalue, yvalue, marker='x', c='C{}'.format((value + 4) % 10), label='${}$'.format(value), s=20)

    ax11.legend(loc='center left', handletextpad=0, handlelength=1, borderpad=0.2, framealpha=1, edgecolor=None,
               markerscale=1,
               fontsize=10, ncol=1, bbox_to_anchor=(1, 0.5))

    ax11.tick_params(axis="x", labelsize=10)
    ax11.tick_params(axis="y", labelsize=10)

    ax11.axvline(0.0192, c='b', linestyle=':', zorder=-2)
    ax11.grid(color='k', linestyle='-', linewidth=0.3, axis='x')
    ax11.set_xlim([0, 0.6])
    ax11.set_xticks(np.arange(0, 0.61, 0.2))
    ax11.set_yticks(np.arange(0, 10, 2))
    ax11.set_ylim([0, 10])
    ax11.xaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax11.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax11.set_xlabel("$V$", fontsize=11)
    ax11.set_ylabel("$\\epsilon_\\alpha$", fontsize=11)
    ax11.text(0.415, 9.2, '$\chi=10^3$', fontsize=10, verticalalignment='top',
             bbox=dict(boxstyle='round', facecolor='white', alpha=1))
    
    ####################################################################################################################

    ax12 = plt.subplot(lower_right_cell)

    ent_spec_mom_dir = '/home/bart/PycharmProjects/infinite_cylinder/data/ent_spec_mom/FerHofSqu1'
    ent_spec_mom_file = 'ent_spec_mom_FerHofSqu1_chi_1000_chiK_1000_t1_1_V_0.6_Coulomb_1_n_1_18_nphi_1_6_LxMUC_1_Ly_6.dat.rotated'
    ent_spec_mom_path = os.path.join(ent_spec_mom_dir, ent_spec_mom_file)

    x = []
    y = []
    z = []

    with open(ent_spec_mom_path, 'r') as csvfile:
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
            ax12.scatter(xvalue, yvalue, marker='_', c='C{}'.format((value + 4) % 10), label='${}$'.format(value))
        else:
            ax12.scatter(xvalue, yvalue, marker='x', c='C{}'.format((value + 4) % 10), label='${}$'.format(value))

    ax12.set_yticks(np.arange(0, 20.1, 5))
    ax12.set_ylim([0, 20])

    # ax12.legend(loc='center left', handletextpad=0, handlelength=1, borderpad=0.2, framealpha=1, edgecolor=None,
    #            markerscale=1,
    #            fontsize=10, ncol=1, bbox_to_anchor=(1, 0.5))
    ax12.set_xlim([-np.pi / 3, np.pi / 3])
    ax12.set_xlabel("$\\tilde{k}_a / \pi$", fontsize=11)
    ax12.set_ylabel("$\epsilon_{\\alpha}$", fontsize=11)

    # ax12.plot([-1.05, 0.3], [5, 18.2], color='k', linestyle=':', linewidth=1)
    # ax12.plot([-0.3, 0.3], [0, 5], color='k', linestyle=':', linewidth=1)

    ax12.tick_params(axis="x", labelsize=10)
    ax12.tick_params(axis="y", labelsize=10)

    ####################################################################################################################


    left_con = ConnectionPatch(xyA=(0, 15), xyB=(0.6, 0), coordsA="data", coordsB="data",
                                axesA=ax6, axesB=ax5, connectionstyle="angle3,angleA=40,angleB=200", arrowstyle='->',
                                facecolor='k', edgecolor='k')
    right_con = ConnectionPatch(xyA=(0, 20), xyB=(0.6, 0), coordsA="data", coordsB="data",
                               axesA=ax12, axesB=ax11, connectionstyle="angle3,angleA=40,angleB=200", arrowstyle='->',
                               facecolor='k', edgecolor='k')

    ax6.add_artist(left_con)
    ax12.add_artist(right_con)

    # fig.text(0.03, 0.87, "(a)", fontsize=12)
    # fig.text(0.5, 0.87, "(b)", fontsize=12)

    plt.savefig("/home/bart/Documents/papers/FCI/Vflow_analysis_full.png", bbox_inches='tight', dpi=300)
    plt.show()
