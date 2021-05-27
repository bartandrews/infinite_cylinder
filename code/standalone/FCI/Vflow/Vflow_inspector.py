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
from pathlib import Path


plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}')
# matplotlib.verbose.level = 'debug-annoying'

if __name__ == '__main__':

    # Configuration to inspect #########################################################################################
    n = (1, 12)
    nphi = (1, 4)
    Ly = 6
    ####################################################################################################################

    fig = plt.figure(figsize=(20, 10))
    outer_grid = gridspec.GridSpec(1, 2)
    left_cell = outer_grid[0, 0]
    right_cell = outer_grid[0, 1]
    left_inner_grid = gridspec.GridSpecFromSubplotSpec(2, 2, left_cell)
    right_inner_grid = gridspec.GridSpecFromSubplotSpec(2, 2, right_cell)

    # define a list of easily-visible markers
    markers = [(3, 0, 0), (4, 0, 0), (5, 0, 0), (6, 0, 0), (4, 1, 0), (5, 1, 0), (6, 1, 0),
               (3, 2, 0), (4, 2, 0), (5, 2, 0), (6, 2, 0), 'X', 'x', 'd', 'D', 'P',
               '$a$', '$b$', '$c$', '$d$', '$e$', '$f$', '$g$', '$h$', '$i$', '$j$', '$k$', '$l$', '$m$', '$n$', '$o$',
               '$p$', '$q$', '$r$', '$s$', '$t$', '$u$', '$v$', '$w$', '$x$', '$y$', '$z$']

    # V_max = 10 #######################################################################################################

    ax1 = plt.subplot(left_inner_grid[0, 0])

    for i, chi_val in enumerate(np.arange(100, 1050, 100)):
        corr_len_dir = '/home/bart/PycharmProjects/infinite_cylinder/data/corr_len_V_flow/FerHofSqu1'
        corr_len_file = f'corr_len_V_flow_FerHofSqu1_chi_{chi_val}_t1_1_V_0_10_41_Coulomb_1_n_{n[0]}_{n[1]}_nphi_{nphi[0]}_{nphi[1]}_LxMUC_1_Ly_{Ly}.dat'
        corr_len_path = os.path.join(corr_len_dir, corr_len_file)

        if Path(corr_len_path).is_file():
            # extract data from file
            with open(corr_len_path, 'r') as csvfile:
                plots = csv.reader(csvfile, delimiter='\t')
                V = []
                xi = []
                for row in plots:
                    V.append(float(row[0]))
                    xi.append(float(row[1]))

            ax1.plot(V, xi, c=f"C{i}", marker=markers[i], markersize=5, label=f"$\chi={chi_val}$")

    ax1.legend(loc='upper right', handletextpad=0, handlelength=1, borderpad=0.2, framealpha=1, edgecolor='k',
               markerscale=1, fontsize=10, ncol=1)

    ax1.set_xlim([0, 10])
    ax1.set_xlabel("$V$", fontsize=11)
    ax1.xaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax1.set_ylabel("$\\xi$", fontsize=11)

    ####################################################################################################################

    ax2 = plt.subplot(left_inner_grid[0, 1])

    for i, chi_val in enumerate(np.arange(100, 1050, 100)):
        corr_len_dir = '/home/bart/PycharmProjects/infinite_cylinder/data/ent_V_flow/FerHofSqu1'
        corr_len_file = f'ent_V_flow_FerHofSqu1_chi_{chi_val}_t1_1_V_0_10_41_Coulomb_1_n_{n[0]}_{n[1]}_nphi_{nphi[0]}_{nphi[1]}_LxMUC_1_Ly_{Ly}.dat'
        corr_len_path = os.path.join(corr_len_dir, corr_len_file)

        if Path(corr_len_path).is_file():
            # extract data from file
            with open(corr_len_path, 'r') as csvfile:
                plots = csv.reader(csvfile, delimiter='\t')
                V = []
                SvN = []
                for row in plots:
                    V.append(float(row[0]))
                    SvN.append(float(row[1]))

            ax2.plot(V, SvN, c=f"C{i}", marker=markers[i], markersize=5, label=f"$\chi={chi_val}$")

            chi_val_max = chi_val-100  # for next two subplots

    ax2.legend(loc='upper right', handletextpad=0, handlelength=1, borderpad=0.2, framealpha=1, edgecolor='k',
               markerscale=1, fontsize=10, ncol=1)

    ax2.set_xlim([0, 10])
    ax2.set_xlabel("$V$", fontsize=11)
    ax2.xaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax2.set_ylabel("$S_\\mathrm{vN}$", fontsize=11)
    
    ####################################################################################################################

    ax3 = plt.subplot(left_inner_grid[1, 0])

    corr_len_dir = '/home/bart/PycharmProjects/infinite_cylinder/data/ent_corr_len/FerHofSqu1'
    corr_len_file = f'ent_corr_len_FerHofSqu1_chi_{chi_val_max}_t1_1_V_0_10_41_Coulomb_1_n_{n[0]}_{n[1]}_nphi_{nphi[0]}_{nphi[1]}_LxMUC_1_Ly_{Ly}.dat'
    corr_len_path = os.path.join(corr_len_dir, corr_len_file)

    # extract data from file
    with open(corr_len_path, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        logxi = []
        SvN = []
        V = []
        for row in plots:
            logxi.append(float(row[0]))
            SvN.append(float(row[1]))
            V.append(float(row[2]))

    # ax3.plot(logxi, SvN, 'x-', c='k', markersize=5)
    im = ax3.scatter(logxi, SvN, c=V, cmap='gist_rainbow')
    plt.colorbar(im, label='$V$')

    ax3.set_xlabel("$\ln\\xi$", fontsize=11)
    ax3.xaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax3.set_ylabel("$S_\\mathrm{vN}$", fontsize=11)
    ax3.set_title(f"$\chi_\mathrm{{max}}={chi_val_max}$")

    ####################################################################################################################

    ax4 = plt.subplot(left_inner_grid[1, 1])

    ent_spec_dir = '/home/bart/PycharmProjects/infinite_cylinder/data/ent_spec_V_flow/FerHofSqu1'
    ent_spec_file = f'ent_spec_V_flow_FerHofSqu1_chi_{chi_val_max}_t1_1_V_0_10_41_Coulomb_1_n_{n[0]}_{n[1]}_nphi_{nphi[0]}_{nphi[1]}_LxMUC_1_Ly_{Ly}.dat'
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
    ax4.set_title(f"$\chi_\mathrm{{max}}={chi_val_max}$")

    # V_max = 1 ########################################################################################################

    ax5 = plt.subplot(right_inner_grid[0, 0])

    for i, chi_val in enumerate(np.arange(100, 1050, 100)):
        corr_len_dir = '/home/bart/PycharmProjects/infinite_cylinder/data/corr_len_V_flow/FerHofSqu1'
        corr_len_file = f'corr_len_V_flow_FerHofSqu1_chi_{chi_val}_t1_1_V_0_1_41_Coulomb_1_n_{n[0]}_{n[1]}_nphi_{nphi[0]}_{nphi[1]}_LxMUC_1_Ly_{Ly}.dat'
        corr_len_path = os.path.join(corr_len_dir, corr_len_file)

        if Path(corr_len_path).is_file():
            # extract data from file
            with open(corr_len_path, 'r') as csvfile:
                plots = csv.reader(csvfile, delimiter='\t')
                V = []
                xi = []
                for row in plots:
                    V.append(float(row[0]))
                    xi.append(float(row[1]))

            ax5.plot(V, xi, c=f"C{i}", marker=markers[i], markersize=5, label=f"$\chi={chi_val}$")

    ax5.legend(loc='upper right', handletextpad=0, handlelength=1, borderpad=0.2, framealpha=1, edgecolor='k',
               markerscale=1, fontsize=10, ncol=1)

    ax5.set_xlim([0, 1])
    ax5.set_xlabel("$V$", fontsize=11)
    ax5.xaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax5.set_ylabel("$\\xi$", fontsize=11)

    ####################################################################################################################

    ax6 = plt.subplot(right_inner_grid[0, 1])

    for i, chi_val in enumerate(np.arange(100, 1050, 100)):
        corr_len_dir = '/home/bart/PycharmProjects/infinite_cylinder/data/ent_V_flow/FerHofSqu1'
        corr_len_file = f'ent_V_flow_FerHofSqu1_chi_{chi_val}_t1_1_V_0_1_41_Coulomb_1_n_{n[0]}_{n[1]}_nphi_{nphi[0]}_{nphi[1]}_LxMUC_1_Ly_{Ly}.dat'
        corr_len_path = os.path.join(corr_len_dir, corr_len_file)

        if Path(corr_len_path).is_file():
            # extract data from file
            with open(corr_len_path, 'r') as csvfile:
                plots = csv.reader(csvfile, delimiter='\t')
                V = []
                SvN = []
                for row in plots:
                    V.append(float(row[0]))
                    SvN.append(float(row[1]))

            ax6.plot(V, SvN, c=f"C{i}", marker=markers[i], markersize=5, label=f"$\chi={chi_val}$")

            chi_val_max = chi_val  # for next two subplots

    ax6.legend(loc='upper right', handletextpad=0, handlelength=1, borderpad=0.2, framealpha=1, edgecolor='k',
               markerscale=1, fontsize=10, ncol=1)

    ax6.set_xlim([0, 1])
    ax6.set_xlabel("$V$", fontsize=11)
    ax6.xaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax6.set_ylabel("$S_\\mathrm{vN}$", fontsize=11)

    ####################################################################################################################

    ax7 = plt.subplot(right_inner_grid[1, 0])

    corr_len_dir = '/home/bart/PycharmProjects/infinite_cylinder/data/ent_corr_len/FerHofSqu1'
    corr_len_file = f'ent_corr_len_FerHofSqu1_chi_{chi_val_max}_t1_1_V_0_1_41_Coulomb_1_n_{n[0]}_{n[1]}_nphi_{nphi[0]}_{nphi[1]}_LxMUC_1_Ly_{Ly}.dat'
    corr_len_path = os.path.join(corr_len_dir, corr_len_file)

    # extract data from file
    with open(corr_len_path, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        logxi = []
        SvN = []
        V = []
        for row in plots:
            logxi.append(float(row[0]))
            SvN.append(float(row[1]))
            V.append(float(row[2]))

    # ax3.plot(logxi, SvN, 'x-', c='k', markersize=5)
    im = ax7.scatter(logxi, SvN, c=V, cmap='gist_rainbow')
    plt.colorbar(im, label='$V$')

    ax7.set_xlabel("$\ln\\xi$", fontsize=11)
    ax7.xaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax7.set_ylabel("$S_\\mathrm{vN}$", fontsize=11)
    ax7.set_title(f"$\chi_\mathrm{{max}}={chi_val_max}$")

    ####################################################################################################################

    ax8 = plt.subplot(right_inner_grid[1, 1])

    ent_spec_dir = '/home/bart/PycharmProjects/infinite_cylinder/data/ent_spec_V_flow/FerHofSqu1'
    ent_spec_file = f'ent_spec_V_flow_FerHofSqu1_chi_{chi_val_max}_t1_1_V_0_1_41_Coulomb_1_n_{n[0]}_{n[1]}_nphi_{nphi[0]}_{nphi[1]}_LxMUC_1_Ly_{Ly}.dat'
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

    ax8.set_xlim([0, 1])
    ax8.set_ylim([0, 10])
    ax8.xaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax8.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax8.set_xlabel("$V$", fontsize=11)
    ax8.set_ylabel("$\\epsilon_\\alpha$", fontsize=11)
    ax8.set_title(f"$\chi_\mathrm{{max}}={chi_val_max}$")

    ####################################################################################################################

    fig.suptitle(f"$n={n[0]}/{n[1]},n_\phi={nphi[0]}/{nphi[1]},L_y={Ly}$")
    plt.figtext(0.285, 0.91, "$V_\mathrm{max}=10$", fontsize=11)
    plt.figtext(0.71, 0.91, "$V_\mathrm{max}=1$", fontsize=11)
    plt.show()
