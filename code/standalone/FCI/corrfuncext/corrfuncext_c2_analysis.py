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

    fig = plt.figure(figsize=(6, 4.5))
    gs = gridspec.GridSpec(3, 2, hspace=0.6, wspace=0.6)

    # define a list of easily-visible markers
    markers = [(3, 0, 0), (4, 0, 0), (5, 0, 0), (6, 0, 0), (4, 1, 0), (5, 1, 0), (6, 1, 0),
               (3, 2, 0), (4, 2, 0), (5, 2, 0), (6, 2, 0), 'X', 'x', 'd', 'D', 'P',
               '$a$', '$b$', '$c$', '$d$', '$e$', '$f$', '$g$', '$h$', '$i$', '$j$', '$k$', '$l$', '$m$', '$n$', '$o$',
               '$p$', '$q$', '$r$', '$s$', '$t$', '$u$', '$v$', '$w$', '$x$', '$y$', '$z$']

    ####################################################################################################################

    corrfunc_dir = '/home/bart/PycharmProjects/infinite_cylinder/data/corr_func_ext/FerHofSqu1'

    ax = plt.subplot(gs[0])  # 071829 #################################################################################
    nu = (1, 5)

    corrfunc_file = f'corr_func_ext_FerHofSqu1_chi_50_t1_1_V_10_Coulomb_1_n_1_55_nphi_6_11_LxMUC_1_Ly_10.dat'
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

    ax.plot(sites[1:], corr_func_offset[1:], '.-', c=f'C1', marker=markers[8], fillstyle='none', markersize=5)
    if corr_len_scale:
        ax.axvline(24.50763761530804, c=f'C1', ls='--', zorder=-3)

    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    if corr_len_scale:
        ax.set_xlim([0, 30])
        ax.set_xticks(np.arange(0, 31, 10))
    else:
        ax.set_xlim([0, 20])
        ax.set_xticks(np.arange(0, 21, 10))
    # ax.set_xlim([0, 110])
    # ax.set_xticks(np.arange(0, 10 + 0.1, 1))
    # ax.set_ylim(bottom=0)
    ax.set_xlabel("$x$", fontsize=11)
    ax.set_ylabel("$g(x)$", fontsize=11)
    # ax.text(0.05*10, 0, f"$\\nu={nu[0]}/{nu[1]}$", fontsize=11)  # 0.011425
    ax.annotate(f"$\\nu={nu[0]}/{nu[1]}$", xy=(0.69, 0.87), xycoords='axes fraction', fontsize=10, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='white', alpha=1))
    # ax.set_title(f"$\\nu={nu[0]}/{nu[1]}$", fontsize=11)

    ax1 = plt.subplot(gs[1])  # 071829 ##################################################################################
    nu = (1, 3)

    corrfunc_file = f'corr_func_ext_FerHofSqu1_chi_25_t1_1_V_10_Coulomb_1_n_1_27_nphi_5_9_LxMUC_1_Ly_6.dat'
    corrfunc_path = os.path.join(corrfunc_dir, corrfunc_file)
    with open(corrfunc_path, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        sites = np.arange(100*9)  # Landau gauge
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

    ax1.plot(sites[1:], corr_func_offset[1:], '.-', c=f'C0', marker=markers[1], fillstyle='none', markersize=5)
    if corr_len_scale:
        ax1.axvline(12.381443957055637, c=f'C0', ls='--', zorder=-3)
    # ax1.plot(sites_cont, corr_func_cont, '--', c=f'C0')

    corrfunc_file = f'corr_func_ext_FerHofSqu1_chi_25_t1_1_V_10_Coulomb_1_n_1_33_nphi_6_11_LxMUC_1_Ly_6.dat'
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

    ax1.plot(sites[1:], corr_func_offset[1:], '.-', c=f'C0', marker=markers[8], fillstyle='none', markersize=5)
    if corr_len_scale:
        ax1.axvline(16.978643738691847, c=f'C0', ls='--', zorder=-3)
    # ax1.plot(sites_cont, corr_func_cont, '--', c=f'C0')

    ax1.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax1.set_xlim([0, 20])
    ax1.set_xticks(np.arange(0, 21, 10))
    # ax1.set_xticks(np.arange(0, 6+0.1, 1))
    # ax1.set_ylim(0)
    ax1.set_xlabel("$x$", fontsize=11)
    ax1.set_ylabel("$g(x)$", fontsize=11)
    # ax1.text(0.3*6, 0, f"$\\nu={nu[0]}/{nu[1]}$", fontsize=11)  # 0.09
    ax1.annotate(f"$\\nu={nu[0]}/{nu[1]}$", xy=(0.69, 0.87), xycoords='axes fraction', fontsize=10,
                verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='white', alpha=1))

    ax2 = plt.subplot(gs[2])  # 071829 #################################################################################
    nu = (2, 9)

    corrfunc_file = f'corr_func_ext_FerHofSqu1_chi_50_t1_1_V_10_Coulomb_1_n_2_81_nphi_5_9_LxMUC_1_Ly_9.dat'
    corrfunc_path = os.path.join(corrfunc_dir, corrfunc_file)
    with open(corrfunc_path, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        sites = np.arange(100*9)  # Landau gauge
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

    ax2.plot(sites[1:], corr_func_offset[1:], '.-', c=f'C1', marker=markers[1], fillstyle='none', markersize=5)
    if corr_len_scale:
        ax2.axvline(63.80363247769049, c=f'C1', ls='--', zorder=-3)
    # ax2.plot(sites_cont, corr_func_cont, '--', c=f'C1')

    ax2.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    if corr_len_scale:
        ax2.set_xlim([0, 70])
        ax2.set_xticks(np.arange(0, 71, 35))
    else:
        ax2.set_xlim([0, 20])
        ax2.set_xticks(np.arange(0, 21, 10))
    # ax2.set_xlim([0, 90])
    # ax2.set_xticks(np.arange(0, 9 + 0.1, 1))
    # ax2.set_ylim(0)
    ax2.set_xlabel("$x$", fontsize=11)
    ax2.set_ylabel("$g(x)$", fontsize=11)
    # ax2.text(0.35*9, 0, f"$\\nu={nu[0]}/{nu[1]}$", fontsize=11)  # 0.0308
    ax2.annotate(f"$\\nu={nu[0]}/{nu[1]}$", xy=(0.69, 0.87), xycoords='axes fraction', fontsize=10,
                 verticalalignment='top',
                 bbox=dict(boxstyle='round', facecolor='white', alpha=1))

    nphi_legend_elements = [
        Line2D([0], [0], linestyle='none', marker=markers[0], color='k', label='$3/5$', fillstyle='none', markersize=5),
        Line2D([0], [0], linestyle='none', marker=markers[7], color='k', label='$4/7$', fillstyle='none', markersize=5),
        Line2D([0], [0], linestyle='none', marker=markers[1], color='k', label='$5/9$', fillstyle='none', markersize=5),
        Line2D([0], [0], linestyle='none', marker=markers[8], color='k', label='$6/11$', fillstyle='none', markersize=5),
        Line2D([0], [0], linestyle='none', marker=markers[2], color='k', label='$7/13$', fillstyle='none', markersize=5),
        Line2D([0], [0], linestyle='none', marker=markers[9], color='k', label='$8/15$', fillstyle='none', markersize=5)]
    leg1 = ax2.legend(handles=nphi_legend_elements, loc='center', handletextpad=0.3, handlelength=1, labelspacing=0.1,
                      borderpad=0.3, framealpha=1, edgecolor='k', markerscale=2, fontsize=10, ncol=2, columnspacing=0.5,
                      bbox_to_anchor=(0.15, -1.2), title='$n_\phi$', title_fontsize=11)
    leg1.get_frame().set_linewidth(0.5)

    ax3 = plt.subplot(gs[3])  # 071829 #################################################################################
    nu = (2, 7)

    corrfunc_file = f'corr_func_ext_FerHofSqu1_chi_100_t1_1_V_10_Coulomb_1_n_2_105_nphi_8_15_LxMUC_1_Ly_14.dat'
    corrfunc_path = os.path.join(corrfunc_dir, corrfunc_file)
    with open(corrfunc_path, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        sites = np.arange(100*15)  # Landau gauge
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

    ax3.plot(sites[1:], corr_func_offset[1:], '.-', c=f'C2', marker=markers[9], fillstyle='none', markersize=5)
    if corr_len_scale:
        ax3.axvline(54.80752722646903, c=f'C2', ls='--', zorder=-3)
    # ax3.plot(sites_cont, corr_func_cont, '--', c=f'C2')

    corrfunc_file = f'corr_func_ext_FerHofSqu1_chi_150_t1_1_V_10_Coulomb_1_n_2_105_nphi_8_15_LxMUC_1_Ly_14.dat'
    corrfunc_path = os.path.join(corrfunc_dir, corrfunc_file)
    with open(corrfunc_path, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        sites = np.arange(100*15)  # Landau gauge
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

    ax3.plot(sites[1:], corr_func_offset[1:], '.-', c=f'C3', marker=markers[9], fillstyle='none', markersize=5)
    if corr_len_scale:
        ax3.axvline(69.50673627086873, c=f'C3', ls='--', zorder=-3)
    # ax3.plot(sites_cont, corr_func_cont, '--', c=f'C3')

    ax3.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    if corr_len_scale:
        ax3.set_xlim([0, 70])
        ax3.set_xticks(np.arange(0, 71, 35))
    else:
        ax3.set_xlim([0, 20])
        ax3.set_xticks(np.arange(0, 21, 10))
    # ax3.set_xlim([0, 150])
    # ax3.set_xticks(np.arange(0, 14 + 0.1, 2))
    # ax3.set_ylim(0)
    ax3.set_xlabel("$x$", fontsize=11)
    ax3.set_ylabel("$g(x)$", fontsize=11)
    # ax3.text(0.35*14, 0, f"$\\nu={nu[0]}/{nu[1]}$", fontsize=11)  # 0.01265
    ax3.annotate(f"$\\nu={nu[0]}/{nu[1]}$", xy=(0.69, 0.87), xycoords='axes fraction', fontsize=10,
                 verticalalignment='top',
                 bbox=dict(boxstyle='round', facecolor='white', alpha=1))

    ax4 = plt.subplot(gs[5])  # 071829 #################################################################################
    nu = (3, 11)

    corrfunc_file = f'corr_func_ext_FerHofSqu1_chi_50_t1_1_V_10_Coulomb_1_n_1_55_nphi_8_15_LxMUC_1_Ly_11.dat'
    corrfunc_path = os.path.join(corrfunc_dir, corrfunc_file)
    with open(corrfunc_path, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        sites = np.arange(100*15)  # Landau gauge
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

    ax4.plot(sites[1:], corr_func_offset[1:], '.-', c=f'C1', marker=markers[9], fillstyle='none', markersize=5)
    if corr_len_scale:
        ax4.axvline(51.53909790997497, c=f'C1', ls='--', zorder=-3)
    # ax4.plot(sites_cont, corr_func_cont, '--', c=f'C1')

    corrfunc_file = f'corr_func_ext_FerHofSqu1_chi_100_t1_1_V_10_Coulomb_1_n_3_143_nphi_7_13_LxMUC_1_Ly_11.dat'
    corrfunc_path = os.path.join(corrfunc_dir, corrfunc_file)
    with open(corrfunc_path, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        sites = np.arange(100*13)  # Landau gauge
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

    ax4.plot(sites[1:], corr_func_offset[1:], '.-', c=f'C2', marker=markers[2], fillstyle='none', markersize=5)
    if corr_len_scale:
        ax4.axvline(70.21306138185336, c=f'C2', ls='--', zorder=-3)
    # ax4.plot(sites_cont, corr_func_cont, '--', c=f'C2')

    corrfunc_file = f'corr_func_ext_FerHofSqu1_chi_150_t1_1_V_10_Coulomb_1_n_3_121_nphi_6_11_LxMUC_1_Ly_11.dat'
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

    ax4.plot(sites[1:], corr_func_offset[1:], '.-', c=f'C3', marker=markers[8], fillstyle='none', markersize=5)
    if corr_len_scale:
        ax4.axvline(71.91924134595075, c=f'C3', ls='--', zorder=-3)
    # ax4.plot(sites_cont, corr_func_cont, '--', c=f'C3')

    ax4.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    if corr_len_scale:
        ax4.set_xlim([0, 80])
        ax4.set_xticks(np.arange(0, 81, 20))
    else:
        ax4.set_xlim([0, 20])
        ax4.set_xticks(np.arange(0, 21, 10))
    # ax4.set_xlim([0, 110])
    # ax4.set_xticks(np.arange(0, 11 + 0.1, 1))
    # ax4.set_ylim(0)
    ax4.set_xlabel("$x$", fontsize=11)
    ax4.set_ylabel("$g(x)$", fontsize=11)
    # ax4.text(0.35*11, 0, f"$\\nu={nu[0]}/{nu[1]}$", fontsize=11)  # 0.0185
    ax4.annotate(f"$\\nu={nu[0]}/{nu[1]}$", xy=(0.65, 0.87), xycoords='axes fraction', fontsize=10,
                 verticalalignment='top',
                 bbox=dict(boxstyle='round', facecolor='white', alpha=1))

    chi_legend_elements = [Patch(facecolor='C0', label='$25$'), Patch(facecolor='C1', label='$50$'),
                           Patch(facecolor='C2', label='$100$'), Patch(facecolor='C3', label='$150$')]
    leg2 = ax4.legend(handles=chi_legend_elements, loc='center', handletextpad=0.3, handlelength=1, labelspacing=0.1,
                      borderpad=0.3, framealpha=1, edgecolor='k', markerscale=2, fontsize=10, ncol=2, columnspacing=0.5,
                      bbox_to_anchor=(-0.8, 0.4), title='$\chi$', title_fontsize=11)
    leg2.get_frame().set_linewidth(0.5)

    ####################################################################################################################
    ####################################################################################################################

    # fig.text(0, 0.87, "(a)", fontsize=12)
    # fig.text(0.48, 0.87, "(b)", fontsize=12)
    # fig.text(0, 0.575, "(c)", fontsize=12)
    # fig.text(0.48, 0.575, "(d)", fontsize=12)
    # fig.text(0.48, 0.28, "(e)", fontsize=12)
    fig.text(0-0.02, 0.89, "(a)", fontsize=12)
    fig.text(0.48-0.02, 0.89, "(b)", fontsize=12)
    fig.text(0-0.02, 0.595, "(c)", fontsize=12)
    fig.text(0.48-0.02, 0.595, "(d)", fontsize=12)
    fig.text(0.48-0.02, 0.30, "(e)", fontsize=12)

    plt.savefig("/home/bart/Documents/papers/FCI/corrfuncext_c2_analysis.png", bbox_inches='tight', dpi=300)
    plt.show()
