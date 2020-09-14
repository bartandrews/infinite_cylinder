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


plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}')
# matplotlib.verbose.level = 'debug-annoying'

if __name__ == '__main__':

    fig = plt.figure()

    outer_grid = gridspec.GridSpec(2, 2, hspace=0.5, wspace=0.4, height_ratios=[2, 1])
    top_left_cell = outer_grid[0, 0]
    top_right_cell = outer_grid[0, 1]
    bottom_left_cell = outer_grid[1, 0]
    bottom_right_cell = outer_grid[1, 1]

    top_left_grid = gridspec.GridSpecFromSubplotSpec(3, 1, top_left_cell, hspace=0, height_ratios=[2, 1, 1])
    top_right_grid = gridspec.GridSpecFromSubplotSpec(3, 1, top_right_cell, hspace=0, height_ratios=[2, 1, 1])

    ####################################################################################################################

    ax1 = plt.subplot(top_left_grid[0], anchor=(0, 0.85))
    ax1.tick_params('x', direction='in', bottom=False)

    ent_spec_V_flow_dir = '/home/bart/PycharmProjects/infinite_cylinder/data/ent_spec_V_flow/FerHalHexC1'
    ent_spec_V_flow_file = 'ent_spec_V_flow_FerHalHexC1_chi_100_t1_1_V_0_10_41_Coulomb_1_n_1_2_LxMUC_1_Ly_3.dat'
    ent_spec_V_flow_path = os.path.join(ent_spec_V_flow_dir, ent_spec_V_flow_file)

    x = []
    y = []
    z = []

    with open(ent_spec_V_flow_path, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        for row in plots:
            x.append(float(row[1]))
            y.append(float(row[2]))
            z.append(int(row[0]))

    for value in np.linspace(-3, 3, 7, dtype=int):
        xvalue = []
        yvalue = []
        for i in range(len(x)):
            if z[i] == value:
                xvalue.append(x[i])
                yvalue.append(y[i])
        if value != 0:
            ax1.scatter(xvalue, yvalue, marker='_', c=f'C{value+4}', label=f'${value}$')
        else:
            ax1.scatter(xvalue, yvalue, marker='x', c=f'C{value+4}', label=f'${value}$')

    # Shrink current axis by 20%
    # box = ax1.get_position()
    # ax1.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    leg = ax1.legend(loc='upper center', handletextpad=0, handlelength=1, borderpad=0.2, framealpha=1, edgecolor='k',
               markerscale=1,
               fontsize=10, ncol=10, columnspacing=0.5, bbox_to_anchor=(0.5, 1.35))

    leg.get_frame().set_linewidth(0.5)

    ax1.tick_params(axis="x", labelsize=10)
    ax1.set_xlim([0, 10])
    ax1.set_xticks(np.arange(0, 10.1, 2))
    ax1.tick_params(axis="y", labelsize=10)
    ax1.set_ylim([0, 8])
    ax1.set_xlabel("$V$", fontsize=11)
    ax1.set_ylabel("$\epsilon_{\\alpha}$", fontsize=11)

    ####################################################################################################################

    ax2 = plt.subplot(top_left_grid[1], sharex=ax1)
    ax2.tick_params('x', direction='in', bottom=False)

    ent_V_flow_dir = '/home/bart/PycharmProjects/infinite_cylinder/data/ent_V_flow/FerHalHexC1'
    ent_V_flow_file = 'ent_V_flow_FerHalHexC1_chi_100_t1_1_V_0_10_41_Coulomb_1_n_1_2_LxMUC_1_Ly_3.dat'
    ent_V_flow_path = os.path.join(ent_V_flow_dir, ent_V_flow_file)

    # extract data from file
    with open(ent_V_flow_path, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        V = []
        SvN = []
        for row in plots:
            V.append(float(row[0]))
            SvN.append(float(row[1]))

    ax2.plot(V, SvN, 's', marker='x', color='k', markersize=3)

    ax2.tick_params(axis="x", labelsize=10)
    ax2.set_xlim([0, 10])
    ax2.set_xticks(np.arange(0, 10.1, 2))
    ax2.tick_params(axis="y", labelsize=10)
    ax2.set_ylim([0, 1.5])
    ax2.set_xlabel("$V$", fontsize=11)
    ax2.set_ylabel("$S_\\text{vN}$", fontsize=11)

    ####################################################################################################################

    ax3 = plt.subplot(top_left_grid[2], sharex=ax2)

    corr_len_V_flow_dir = '/home/bart/PycharmProjects/infinite_cylinder/data/corr_len_V_flow/FerHalHexC1'
    corr_len_V_flow_file = 'corr_len_V_flow_FerHalHexC1_chi_100_t1_1_V_0_10_41_Coulomb_1_n_1_2_LxMUC_1_Ly_3.dat'
    corr_len_V_flow_path = os.path.join(corr_len_V_flow_dir, corr_len_V_flow_file)

    # extract data from file
    with open(corr_len_V_flow_path, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        V = []
        xi = []
        for row in plots:
            V.append(float(row[0]))
            xi.append(float(row[1]))

    ax3.plot(V, xi, 's', marker='x', color='k', markersize=3)

    ax3.axvline(V[xi.index(max(xi))], color='k', linewidth=0.5, ls='--')
    ax2.axvline(V[xi.index(max(xi))], color='k', linewidth=0.5, ls='--')
    ax1.axvline(V[xi.index(max(xi))], color='k', linewidth=0.5, ls='--')

    ax3.tick_params(axis="x", labelsize=10)
    ax3.set_xlim([0, 10])
    ax3.set_xticks(np.arange(0, 10.1, 2))
    ax3.tick_params(axis="y", labelsize=10)
    ax3.set_ylim([0, 15])
    ax3.set_xlabel("$V$", fontsize=11)
    ax3.set_ylabel("$\\xi$", fontsize=11)

    ####################################################################################################################

    ax4 = plt.subplot(bottom_left_cell)

    ent_corr_len_dir = '/home/bart/PycharmProjects/infinite_cylinder/data/ent_corr_len/FerHalHexC1'
    ent_corr_len_file = 'ent_corr_len_FerHalHexC1_chi_100_t1_1_V_0_10_41_Coulomb_1_n_1_2_LxMUC_1_Ly_3.dat'
    ent_corr_len_path = os.path.join(ent_corr_len_dir, ent_corr_len_file)

    lower_V_limit = 2
    upper_V_limit = 3

    # extract data from file
    with open(ent_corr_len_path, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        ln_xi = []
        SvN = []
        for row in plots:
            if lower_V_limit <= float(row[2]) <= upper_V_limit:
                ln_xi.append(float(row[0]))
                SvN.append(float(row[1]))

    # clockwise, starting from bottom left
    metal0 = Polygon(((lower_V_limit, 0), (lower_V_limit, 2), (upper_V_limit, 2), (upper_V_limit, 0)), fc=(0, 0, 0, 0.2))
    metal1 = Polygon(((lower_V_limit, 0), (lower_V_limit, 15), (upper_V_limit, 15), (upper_V_limit, 0)), fc=(0, 0, 0, 0.2))
    ax2.add_artist(metal0)
    ax3.add_artist(metal1)

    ax3.annotate('', xy=(0.2, 0), xycoords='axes fraction', xytext=(0.004, -1.54), arrowprops=dict(arrowstyle="-", color='k', linewidth=0.5))
    ax3.annotate('', xy=(0.29, 0), xycoords='axes fraction', xytext=(1.01, -1.54), arrowprops=dict(arrowstyle="-", color='k', linewidth=0.5))

    ax4.plot(ln_xi, SvN, 's', marker='x', color='k', markersize=3)

    parameters, cov = np.polyfit(ln_xi, SvN, 1, cov=True)
    _, _, r_value, _, _ = stats.linregress(ln_xi, SvN)
    m, m_err, c, c_err = parameters[0], np.sqrt(cov[0][0]), parameters[1], np.sqrt(cov[1][1])
    r2_value = r_value * r_value
    print("SvN = m*ln(xi) + c")
    print(f"(m, m_err, c, c_err) = ({m:.3f}, {m_err:.3f}, {c:.3f}, {c_err:.3f})")
    xvalues = np.linspace(min(ln_xi), max(ln_xi))
    ax4.plot(xvalues, m * xvalues + c, '-', c='k', zorder=0)
    ax4.text(0.1, 0.8, f"$S_\mathrm{{vN}}\sim({m:.3f}\pm{m_err:.3f})\ln\\xi$\n$R^2={r2_value:.5f}$", fontsize=10)

    ax4.tick_params(axis="x", labelsize=10)
    ax4.set_xlim([0, 3])
    # ax4.set_xticks(np.arange(0, 10.1, 2))
    ax4.tick_params(axis="y", labelsize=10)
    ax4.set_ylim([0, 1.2])
    ax4.set_xlabel("$\ln \\xi$", fontsize=11)
    ax4.set_ylabel("$S_\\text{vN}$", fontsize=11)

    ####################################################################################################################
    ####################################################################################################################

    ax5 = plt.subplot(top_right_grid[0], anchor=(0, 0.85))
    ax5.tick_params('x', direction='in', bottom=False)

    ent_spec_V_flow_dir = '/home/bart/PycharmProjects/infinite_cylinder/data/ent_spec_V_flow/FerHalHexC1'
    ent_spec_V_flow_file = 'ent_spec_V_flow_FerHalHexC1_chi_100_t1_1_V_0_10_41_Coulomb_1_n_1_2_LxMUC_1_Ly_3.dat'
    ent_spec_V_flow_path = os.path.join(ent_spec_V_flow_dir, ent_spec_V_flow_file)

    x = []
    y = []
    z = []

    with open(ent_spec_V_flow_path, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        for row in plots:
            x.append(float(row[1]))
            y.append(float(row[2]))
            z.append(int(row[0]))

    for value in np.linspace(-3, 3, 7, dtype=int):
        xvalue = []
        yvalue = []
        for i in range(len(x)):
            if z[i] == value:
                xvalue.append(x[i])
                yvalue.append(y[i])
        if value != 0:
            ax5.scatter(xvalue, yvalue, marker='_', c=f'C{value+4}', label=f'${value}$')
        else:
            ax5.scatter(xvalue, yvalue, marker='x', c=f'C{value+4}', label=f'${value}$')

    # Shrink current axis by 20%
    # box = ax1.get_position()
    # ax1.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    leg = ax5.legend(loc='upper center', handletextpad=0, handlelength=1, borderpad=0.2, framealpha=1, edgecolor='k',
                     markerscale=1,
                     fontsize=10, ncol=10, columnspacing=0.5, bbox_to_anchor=(0.5, 1.35))

    leg.get_frame().set_linewidth(0.5)

    ax5.tick_params(axis="x", labelsize=10)
    ax5.set_xlim([0, 10])
    ax5.set_xticks(np.arange(0, 10.1, 2))
    ax5.tick_params(axis="y", labelsize=10)
    ax5.set_ylim([0, 8])
    ax5.set_xlabel("$V$", fontsize=11)
    ax5.set_ylabel("$\epsilon_{\\alpha}$", fontsize=11)

    ####################################################################################################################

    ax6 = plt.subplot(top_right_grid[1], sharex=ax5)
    ax6.tick_params('x', direction='in', bottom=False)

    ent_V_flow_dir = '/home/bart/PycharmProjects/infinite_cylinder/data/ent_V_flow/FerHalHexC1'
    ent_V_flow_file = 'ent_V_flow_FerHalHexC1_chi_100_t1_1_V_0_10_41_Coulomb_1_n_1_2_LxMUC_1_Ly_3.dat'
    ent_V_flow_path = os.path.join(ent_V_flow_dir, ent_V_flow_file)

    # extract data from file
    with open(ent_V_flow_path, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        V = []
        SvN = []
        for row in plots:
            V.append(float(row[0]))
            SvN.append(float(row[1]))

    ax6.plot(V, SvN, 's', marker='x', color='k', markersize=3)

    ax6.tick_params(axis="x", labelsize=10)
    ax6.set_xlim([0, 10])
    ax6.set_xticks(np.arange(0, 10.1, 2))
    ax6.tick_params(axis="y", labelsize=10)
    ax6.set_ylim([0, 1.5])
    ax6.set_xlabel("$V$", fontsize=11)
    ax6.set_ylabel("$S_\\text{vN}$", fontsize=11)

    ####################################################################################################################

    ax7 = plt.subplot(top_right_grid[2], sharex=ax6)

    corr_len_V_flow_dir = '/home/bart/PycharmProjects/infinite_cylinder/data/corr_len_V_flow/FerHalHexC1'
    corr_len_V_flow_file = 'corr_len_V_flow_FerHalHexC1_chi_100_t1_1_V_0_10_41_Coulomb_1_n_1_2_LxMUC_1_Ly_3.dat'
    corr_len_V_flow_path = os.path.join(corr_len_V_flow_dir, corr_len_V_flow_file)

    # extract data from file
    with open(corr_len_V_flow_path, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        V = []
        xi = []
        for row in plots:
            V.append(float(row[0]))
            xi.append(float(row[1]))

    ax7.plot(V, xi, 's', marker='x', color='k', markersize=3)

    ax7.axvline(V[xi.index(max(xi))], color='k', linewidth=0.5, ls='--')
    ax6.axvline(V[xi.index(max(xi))], color='k', linewidth=0.5, ls='--')
    ax5.axvline(V[xi.index(max(xi))], color='k', linewidth=0.5, ls='--')

    ax7.tick_params(axis="x", labelsize=10)
    ax7.set_xlim([0, 10])
    ax7.set_xticks(np.arange(0, 10.1, 2))
    ax7.tick_params(axis="y", labelsize=10)
    ax7.set_ylim([0, 15])
    ax7.set_xlabel("$V$", fontsize=11)
    ax7.set_ylabel("$\\xi$", fontsize=11)

    ####################################################################################################################

    ax8 = plt.subplot(bottom_right_cell)

    ent_corr_len_dir = '/home/bart/PycharmProjects/infinite_cylinder/data/ent_corr_len/FerHalHexC1'
    ent_corr_len_file = 'ent_corr_len_FerHalHexC1_chi_100_t1_1_V_0_10_41_Coulomb_1_n_1_2_LxMUC_1_Ly_3.dat'
    ent_corr_len_path = os.path.join(ent_corr_len_dir, ent_corr_len_file)

    lower_V_limit = 2
    upper_V_limit = 3

    # extract data from file
    with open(ent_corr_len_path, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        ln_xi = []
        SvN = []
        for row in plots:
            if lower_V_limit <= float(row[2]) <= upper_V_limit:
                ln_xi.append(float(row[0]))
                SvN.append(float(row[1]))

    # clockwise, starting from bottom left
    metal0 = Polygon(((lower_V_limit, 0), (lower_V_limit, 2), (upper_V_limit, 2), (upper_V_limit, 0)),
                     fc=(0, 0, 0, 0.2))
    metal1 = Polygon(((lower_V_limit, 0), (lower_V_limit, 15), (upper_V_limit, 15), (upper_V_limit, 0)),
                     fc=(0, 0, 0, 0.2))
    ax6.add_artist(metal0)
    ax7.add_artist(metal1)

    ax7.annotate('', xy=(0.2, 0), xycoords='axes fraction', xytext=(0.004, -1.54),
                 arrowprops=dict(arrowstyle="-", color='k', linewidth=0.5))
    ax7.annotate('', xy=(0.29, 0), xycoords='axes fraction', xytext=(1.01, -1.54),
                 arrowprops=dict(arrowstyle="-", color='k', linewidth=0.5))

    ax8.plot(ln_xi, SvN, 's', marker='x', color='k', markersize=3)

    parameters, cov = np.polyfit(ln_xi, SvN, 1, cov=True)
    _, _, r_value, _, _ = stats.linregress(ln_xi, SvN)
    m, m_err, c, c_err = parameters[0], np.sqrt(cov[0][0]), parameters[1], np.sqrt(cov[1][1])
    r2_value = r_value * r_value
    print("SvN = m*ln(xi) + c")
    print(f"(m, m_err, c, c_err) = ({m:.3f}, {m_err:.3f}, {c:.3f}, {c_err:.3f})")
    xvalues = np.linspace(min(ln_xi), max(ln_xi))
    ax8.plot(xvalues, m * xvalues + c, '-', c='k', zorder=0)
    ax8.text(0.1, 0.8, f"$S_\mathrm{{vN}}\sim({m:.3f}\pm{m_err:.3f})\ln\\xi$\n$R^2={r2_value:.5f}$", fontsize=10)

    ax8.tick_params(axis="x", labelsize=10)
    ax8.set_xlim([0, 3])
    # ax4.set_xticks(np.arange(0, 10.1, 2))
    ax8.tick_params(axis="y", labelsize=10)
    ax8.set_ylim([0, 1.2])
    ax8.set_xlabel("$\ln \\xi$", fontsize=11)
    ax8.set_ylabel("$S_\\text{vN}$", fontsize=11)

    ####################################################################################################################
    ####################################################################################################################

    fig.text(0.04, 0.92, "(a)", fontsize=12)
    fig.text(0.04, 0.67, "(b)", fontsize=12)
    fig.text(0.04, 0.56, "(c)", fontsize=12)
    fig.text(0.04, 0.32, "(d)", fontsize=12)

    fig.text(0.49, 0.92, "(e)", fontsize=12)
    fig.text(0.49, 0.67, "(f)", fontsize=12)
    fig.text(0.49, 0.56, "(g)", fontsize=12)
    fig.text(0.49, 0.32, "(h)", fontsize=12)

    plt.setp(ax1.get_xticklabels(), visible=False)
    plt.setp(ax2.get_xticklabels(), visible=False)

    plt.setp(ax5.get_xticklabels(), visible=False)
    plt.setp(ax6.get_xticklabels(), visible=False)

    fig.text(0.22, 0.97, 'TFB model', fontsize=12)
    fig.text(0.65, 0.97, 'Hofstadter model', fontsize=12)

    plt.savefig("/home/bart/Documents/papers/BT/case_study_V.png", bbox_inches='tight', dpi=300)
    plt.show()
