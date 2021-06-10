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
from matplotlib.patches import Patch
from matplotlib.lines import Line2D


plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}')
# matplotlib.verbose.level = 'debug-annoying'


def line_of_best_fit(x_list, y_list):

    parameters, cov = np.polyfit(x_list, y_list, 1, cov=True)
    _, _, r_value, _, _ = stats.linregress(x_list, y_list)
    m, m_err, c, c_err = parameters[0], np.sqrt(cov[0][0]), parameters[1], np.sqrt(cov[1][1])
    r2_value = r_value*r_value

    return m, m_err, c, c_err, r2_value


if __name__ == '__main__':

    corr_len_scale = True

    fig = plt.figure(figsize=(6, 2.75))
    gs = gridspec.GridSpec(2, 2, hspace=0.6, wspace=0.6)

    # define a list of easily-visible markers
    markers = [(3, 0, 0), (4, 0, 0), (5, 0, 0), (6, 0, 0), (4, 1, 0), (5, 1, 0), (6, 1, 0),
               (3, 2, 0), (4, 2, 0), (5, 2, 0), (6, 2, 0), 'X', 'x', 'd', 'D', 'P',
               '$a$', '$b$', '$c$', '$d$', '$e$', '$f$', '$g$', '$h$', '$i$', '$j$', '$k$', '$l$', '$m$', '$n$', '$o$',
               '$p$', '$q$', '$r$', '$s$', '$t$', '$u$', '$v$', '$w$', '$x$', '$y$', '$z$']

    ####################################################################################################################

    corrfunc_dir = '/home/bart/PycharmProjects/infinite_cylinder/data/corr_func_ext/FerHofSqu1'

    ax = plt.subplot(gs[0])  # 071829 #################################################################################
    nu = (1, 7)

    corrfunc_file = f'corr_func_ext_FerHofSqu1_chi_50_t1_1_V_10_Coulomb_1_n_1_98_nphi_5_14_LxMUC_1_Ly_14.dat'
    corrfunc_path = os.path.join(corrfunc_dir, corrfunc_file)
    with open(corrfunc_path, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        sites = np.arange(100*14)  # Landau gauge
        corr_func = []
        for i, row in enumerate(plots):
            corr_func.append(row[0])
    corr_func = [float(i) for i in corr_func]
    min_val = min(corr_func)
    corr_func = [float(i) - min_val for i in corr_func]

    # offset
    corr_func_offset = []
    (m1, m1_err, c, c_err, r2_value1) = line_of_best_fit(sites[-200:], corr_func[-200:])
    for i in corr_func:
        corr_func_offset.append(i - c)

    ax.plot(sites[1:], corr_func_offset[1:], '.-', c=f'C1', marker=markers[1], fillstyle='none', markersize=5)
    if corr_len_scale:
        ax.axvline(47.15388054186111, c=f'C1', ls='--', zorder=-3)
    # ax.plot(sites_cont, corr_func_cont, '--', c=f'C1')

    corrfunc_file = f'corr_func_ext_FerHofSqu1_chi_50_t1_1_V_10_Coulomb_1_n_1_140_nphi_7_20_LxMUC_1_Ly_14.dat'
    corrfunc_path = os.path.join(corrfunc_dir, corrfunc_file)
    with open(corrfunc_path, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        sites = np.arange(100*20)  # Landau gauge
        corr_func = []
        for i, row in enumerate(plots):
            corr_func.append(row[0])
    corr_func = [float(i) for i in corr_func]
    min_val = min(corr_func)
    corr_func = [float(i) - min_val for i in corr_func]

    # offset
    corr_func_offset = []
    (m1, m1_err, c, c_err, r2_value1) = line_of_best_fit(sites[-200:], corr_func[-200:])
    for i in corr_func:
        corr_func_offset.append(i - c)

    ax.plot(sites[1:], corr_func_offset[1:], '.-', c=f'C1', marker=markers[2], fillstyle='none', markersize=5)
    if corr_len_scale:
        ax.axvline(54.381396991666556, c=f'C1', ls='--', zorder=-3)
    # ax.plot(sites_cont, corr_func_cont, '--', c=f'C1')

    corrfunc_file = f'corr_func_ext_FerHofSqu1_chi_100_t1_1_V_10_Coulomb_1_n_1_98_nphi_5_14_LxMUC_1_Ly_14.dat'
    corrfunc_path = os.path.join(corrfunc_dir, corrfunc_file)
    with open(corrfunc_path, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        sites = np.arange(100*14)  # Landau gauge
        corr_func = []
        for i, row in enumerate(plots):
            corr_func.append(row[0])
    corr_func = [float(i) for i in corr_func]
    min_val = min(corr_func)
    corr_func = [float(i) - min_val for i in corr_func]

    # offset
    corr_func_offset = []
    (m1, m1_err, c, c_err, r2_value1) = line_of_best_fit(sites[-200:], corr_func[-200:])
    for i in corr_func:
        corr_func_offset.append(i - c)

    ax.plot(sites[1:], corr_func_offset[1:], '.-', c=f'C2', marker=markers[1], fillstyle='none', markersize=5)
    if corr_len_scale:
        ax.axvline(41.898911948440734, c=f'C2', ls='--', zorder=-3)
    # ax.plot(sites_cont, corr_func_cont, '--', c=f'C2')

    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    if corr_len_scale:
        ax.set_xlim([0, 60])
        ax.set_xticks(np.arange(0, 61, 20))
    else:
        ax.set_xlim([0, 20])
        ax.set_xticks(np.arange(0, 21, 10))
    # ax.set_xlim([0, 100])
    # ax.set_xticks(np.arange(0, 14 + 0.1, 2))
    # ax.set_ylim(0)
    ax.set_xlabel("$x$", fontsize=11)
    ax.set_ylabel("$g(x)$", fontsize=11)
    # ax.text(0.35 * 14, 0.0053, f"$\\nu={nu[0]}/{nu[1]}$", fontsize=11)
    ax.annotate(f"$\\nu={nu[0]}/{nu[1]}$", xy=(0.69, 0.87), xycoords='axes fraction', fontsize=10,
                 verticalalignment='top',
                 bbox=dict(boxstyle='round', facecolor='white', alpha=1))

    nphi_legend_elements = [
        Line2D([0], [0], linestyle='none', marker=markers[0], color='k', label='$3/8$', fillstyle='none', markersize=5),
        Line2D([0], [0], linestyle='none', marker=markers[7], color='k', label='$4/11$', fillstyle='none', markersize=5),
        Line2D([0], [0], linestyle='none', marker=markers[1], color='k', label='$5/14$', fillstyle='none', markersize=5),
        Line2D([0], [0], linestyle='none', marker=markers[8], color='k', label='$6/17$', fillstyle='none', markersize=5),
        Line2D([0], [0], linestyle='none', marker=markers[2], color='k', label='$7/20$', fillstyle='none', markersize=5),
        Line2D([0], [0], linestyle='none', marker=markers[9], color='k', label='$8/23$', fillstyle='none', markersize=5)]
    leg1 = ax.legend(handles=nphi_legend_elements, loc='center', handletextpad=0.3, handlelength=1, labelspacing=0.1,
                      borderpad=0.3, framealpha=1, edgecolor='k', markerscale=2, fontsize=10, ncol=2, columnspacing=0.5,
                      bbox_to_anchor=(0.15, -1.2), title='$n_\phi$', title_fontsize=11)
    leg1.get_frame().set_linewidth(0.5)

    ax1 = plt.subplot(gs[1])  # 071829 ##################################################################################
    nu = (1, 5)

    corrfunc_file = f'corr_func_ext_FerHofSqu1_chi_50_t1_1_V_10_Coulomb_1_n_1_85_nphi_6_17_LxMUC_1_Ly_10.dat'
    corrfunc_path = os.path.join(corrfunc_dir, corrfunc_file)
    with open(corrfunc_path, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        sites = np.arange(100*17)  # Landau gauge
        corr_func = []
        for i, row in enumerate(plots):
            corr_func.append(row[0])
    corr_func = [float(i) for i in corr_func]
    min_val = min(corr_func)
    corr_func = [float(i) - min_val for i in corr_func]

    # offset
    corr_func_offset = []
    (m1, m1_err, c, c_err, r2_value1) = line_of_best_fit(sites[-200:], corr_func[-200:])
    for i in corr_func:
        corr_func_offset.append(i - c)

    ax1.plot(sites[1:], corr_func_offset[1:], '.-', c=f'C1', marker=markers[8], fillstyle='none', markersize=5)
    if corr_len_scale:
        ax1.axvline(54.457403049059415, c=f'C1', ls='--', zorder=-3)
    # ax1.plot(sites_cont, corr_func_cont, '--', c=f'C1')

    corrfunc_file = f'corr_func_ext_FerHofSqu1_chi_50_t1_1_V_10_Coulomb_1_n_1_100_nphi_7_20_LxMUC_1_Ly_10.dat'
    corrfunc_path = os.path.join(corrfunc_dir, corrfunc_file)
    with open(corrfunc_path, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        sites = np.arange(100*20)  # Landau gauge
        corr_func = []
        for i, row in enumerate(plots):
            corr_func.append(row[0])
    corr_func = [float(i) for i in corr_func]
    min_val = min(corr_func)
    corr_func = [float(i) - min_val for i in corr_func]

    # offset
    corr_func_offset = []
    (m1, m1_err, c, c_err, r2_value1) = line_of_best_fit(sites[-200:], corr_func[-200:])
    for i in corr_func:
        corr_func_offset.append(i - c)

    ax1.plot(sites[1:], corr_func_offset[1:], '.-', c=f'C1', marker=markers[2], fillstyle='none', markersize=5)
    if corr_len_scale:
        ax1.axvline(52.34473066303418, c=f'C1', ls='--', zorder=-3)
    # ax1.plot(sites_cont, corr_func_cont, '--', c=f'C1')

    ax1.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    if corr_len_scale:
        ax1.set_xlim([0, 60])
        ax1.set_xticks(np.arange(0, 61, 20))
    else:
        ax1.set_xlim([0, 20])
        ax1.set_xticks(np.arange(0, 21, 10))
    # ax1.set_xlim([0, 100])
    # ax1.set_xticks(np.arange(0, 10+0.1, 1))
    # ax1.set_ylim(0)
    ax1.set_xlabel("$x$", fontsize=11)
    ax1.set_ylabel("$g(x)$", fontsize=11)
    # ax1.text(0.35*10, 0.02, f"$\\nu={nu[0]}/{nu[1]}$", fontsize=11)
    ax1.annotate(f"$\\nu={nu[0]}/{nu[1]}$", xy=(0.69, 0.87), xycoords='axes fraction', fontsize=10,
                 verticalalignment='top',
                 bbox=dict(boxstyle='round', facecolor='white', alpha=1))

    ax2 = plt.subplot(gs[3])  # 071829 #################################################################################
    nu = (2, 11)

    corrfunc_file = f'corr_func_ext_FerHofSqu1_chi_50_t1_1_V_10_Coulomb_1_n_2_121_nphi_4_11_LxMUC_1_Ly_11.dat'
    corrfunc_path = os.path.join(corrfunc_dir, corrfunc_file)
    with open(corrfunc_path, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        sites = np.arange(100*11)  # Landau gauge
        corr_func = []
        for i, row in enumerate(plots):
            corr_func.append(row[0])
    corr_func = [float(i) for i in corr_func]
    min_val = min(corr_func)
    corr_func = [float(i) - min_val for i in corr_func]

    # offset
    corr_func_offset = []
    (m1, m1_err, c, c_err, r2_value1) = line_of_best_fit(sites[-200:], corr_func[-200:])
    for i in corr_func:
        corr_func_offset.append(i - c)

    ax2.plot(sites[1:], corr_func_offset[1:], '.-', c=f'C1', marker=markers[7], fillstyle='none', markersize=5)
    if corr_len_scale:
        ax2.axvline(33.153975928046435, c=f'C1', ls='--', zorder=-3)
    # ax2.plot(sites_cont, corr_func_cont, '--', c=f'C1')

    corrfunc_file = f'corr_func_ext_FerHofSqu1_chi_100_t1_1_V_10_Coulomb_1_n_2_121_nphi_4_11_LxMUC_1_Ly_11.dat'
    corrfunc_path = os.path.join(corrfunc_dir, corrfunc_file)
    with open(corrfunc_path, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        sites = np.arange(100*11)  # Landau gauge
        corr_func = []
        for i, row in enumerate(plots):
            corr_func.append(row[0])
    corr_func = [float(i) for i in corr_func]
    min_val = min(corr_func)
    corr_func = [float(i) - min_val for i in corr_func]

    # offset
    corr_func_offset = []
    (m1, m1_err, c, c_err, r2_value1) = line_of_best_fit(sites[-200:], corr_func[-200:])
    for i in corr_func:
        corr_func_offset.append(i - c)

    ax2.plot(sites[1:], corr_func_offset[1:], '.-', c=f'C2', marker=markers[7], fillstyle='none', markersize=5)
    if corr_len_scale:
        ax2.axvline(47.79041868796284, c=f'C2', ls='--', zorder=-3)
    # ax2.plot(sites_cont, corr_func_cont, '--', c=f'C2')

    ax2.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    if corr_len_scale:
        ax2.set_xlim([0, 50])
        ax2.set_xticks(np.arange(0, 51, 25))
    else:
        ax2.set_xlim([0, 20])
        ax2.set_xticks(np.arange(0, 21, 10))
    # ax2.set_xlim([0, 100])
    # ax2.set_xticks(np.arange(0, 11 + 0.1, 1))
    # ax2.set_ylim(0)
    ax2.set_xlabel("$x$", fontsize=11)
    ax2.set_ylabel("$g(x)$", fontsize=11)
    # ax2.text(0.35 * 11, 0.0205, f"$\\nu={nu[0]}/{nu[1]}$", fontsize=11)
    ax2.annotate(f"$\\nu={nu[0]}/{nu[1]}$", xy=(0.65, 0.87), xycoords='axes fraction', fontsize=10,
                 verticalalignment='top',
                 bbox=dict(boxstyle='round', facecolor='white', alpha=1))

    chi_legend_elements = [Patch(facecolor='C1', label='$50$'), Patch(facecolor='C2', label='$100$')]
    leg2 = ax2.legend(handles=chi_legend_elements, loc='center', handletextpad=0.3, handlelength=1, labelspacing=0.1,
                      borderpad=0.3, framealpha=1, edgecolor='k', markerscale=2, fontsize=10, ncol=2, columnspacing=0.5,
                      bbox_to_anchor=(-0.8, 0.4), title='$\chi$', title_fontsize=11)
    leg2.get_frame().set_linewidth(0.5)

    ####################################################################################################################
    ####################################################################################################################

    # fig.text(0, 0.865, "(f)", fontsize=12)
    # fig.text(0.48, 0.865, "(g)", fontsize=12)
    # fig.text(0.48, 0.39, "(h)", fontsize=12)
    fig.text(0-0.01, 0.885, "(f)", fontsize=12)
    fig.text(0.48-0.01, 0.885, "(g)", fontsize=12)
    fig.text(0.48-0.01, 0.41, "(h)", fontsize=12)

    plt.savefig("/home/bart/Documents/papers/FCI/corrfuncext_c3_analysis.png", bbox_inches='tight', dpi=300)
    plt.show()
