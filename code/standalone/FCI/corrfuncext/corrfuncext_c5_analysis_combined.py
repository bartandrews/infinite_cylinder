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

    corr_len_scale = False

    fig = plt.figure(figsize=(6, 1.1))
    outer_grid = gridspec.GridSpec(1, 2, hspace=0.6, wspace=0.8)
    left = outer_grid[0, 0]
    right = outer_grid[0, 1]
    left_inner_grid = gridspec.GridSpecFromSubplotSpec(1, 2, left, wspace=0)
    right_inner_grid = gridspec.GridSpecFromSubplotSpec(1, 2, right, wspace=0)

    # define a list of easily-visible markers
    markers = [(3, 0, 0), (4, 0, 0), (5, 0, 0), (6, 0, 0), (4, 1, 0), (5, 1, 0), (6, 1, 0),
               (3, 2, 0), (4, 2, 0), (5, 2, 0), (6, 2, 0), 'X', 'x', 'd', 'D', 'P',
               '$a$', '$b$', '$c$', '$d$', '$e$', '$f$', '$g$', '$h$', '$i$', '$j$', '$k$', '$l$', '$m$', '$n$', '$o$',
               '$p$', '$q$', '$r$', '$s$', '$t$', '$u$', '$v$', '$w$', '$x$', '$y$', '$z$']

    ####################################################################################################################

    corrfunc_dir = '/home/bart/PycharmProjects/infinite_cylinder/data/corr_func/FerHofSqu1'
    corrfuncext_dir = '/home/bart/PycharmProjects/infinite_cylinder/data/corr_func_ext/FerHofSqu1'

    ax2 = plt.subplot(left_inner_grid[0, 0])  # 071829 #################################################################################
    nu = (1, 9)

    corrfuncext_file = f'corr_func_ext_FerHofSqu1_chi_25_t1_1_V_10_Coulomb_1_n_1_216_nphi_5_24_LxMUC_1_Ly_18.dat'
    corrfuncext_path = os.path.join(corrfuncext_dir, corrfuncext_file)
    with open(corrfuncext_path, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        sites = np.arange(100*24)  # Landau gauge
        corr_func_ext = []
        for i, row in enumerate(plots):
            corr_func_ext.append(row[0])
    corr_func_ext = [float(i) for i in corr_func_ext]
    min_val = min(corr_func_ext)
    corr_func_ext = [float(i) - min_val for i in corr_func_ext]
    # offset
    corr_func_offset = []
    (m1, m1_err, c, c_err, r2_value1) = line_of_best_fit(sites[-200:], corr_func_ext[-200:])
    for i in corr_func_ext:
        corr_func_offset.append(i - c)
    ax2.plot(sites[1:], corr_func_offset[1:], '.-', c=f'C0', marker=markers[1], fillstyle='none', markersize=5)
    if corr_len_scale:
        ax2.axvline(65.36765041432689, c=f'C0', ls='--', zorder=-3)

    ax2.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    if corr_len_scale:
        ax2.set_xlim([0, 70])
        ax2.set_xticks(np.arange(0, 70, 35))
    else:
        ax2.set_xlim([0, 20])
        ax2.set_xticks(np.arange(0, 20, 10))
    # ax2.set_ylim(0)
    ax2.set_xlabel("$x$", fontsize=11)
    # ax2.set_ylabel("$\langle :\mathrel{\\rho_{0,0} \\rho_{x,0}}: \\rangle$", fontsize=11)
    # ax2.text(0.675 * 18, 0.0059451, f"$\\nu={nu[0]}/{nu[1]}$", fontsize=11)

    ax2_2 = plt.subplot(left_inner_grid[0, 1])  # 071829 #################################################################################
    nu = (1, 9)

    corrfunc_file = f'corr_func_FerHofSqu1_chi_25_t1_1_V_10_Coulomb_1_n_1_216_nphi_5_24_LxMUC_1_Ly_18.dat'
    corrfunc_path = os.path.join(corrfunc_dir, corrfunc_file)
    with open(corrfunc_path, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        sites = np.arange(18)
        corr_func = []
        for i, row in enumerate(plots):
            if i == 0:
                corr_func = row
                break
    corr_func = [float(i) for i in corr_func]
    min_val = min(corr_func)
    corr_func = [float(i) - min_val for i in corr_func]
    sites_cont = [sites[-1], sites[-1] + 1]
    corr_func_cont = [corr_func[-1], corr_func[0]]
    ax2_2.plot(sites[1:], corr_func[1:], '.-', c=f'C0', marker=markers[1], fillstyle='none', markersize=5)
    # ax2.plot(sites_cont, corr_func_cont, '--', c=f'C0')

    ax2_2.yaxis.tick_right()
    ax2_2.yaxis.set_label_position("right")
    ax2_2.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax2_2.set_xlim([0, 18])
    ax2_2.set_xticks(np.arange(0, 18 + 0.1, 6))
    # ax2.set_ylim(0)
    ax2_2.set_xlabel("$y$", fontsize=11)
    # ax2_2.set_ylabel("$\langle :\mathrel{\\rho_{0,0} \\rho_{0,y}}: \\rangle$", fontsize=11)
    ax2_2.text(18/5, 0.0059459, f"$\\nu={nu[0]}/{nu[1]}$", fontsize=11)

    ax3 = plt.subplot(right_inner_grid[0, 0])  # 071829 #################################################################################
    nu = (2, 19)

    corrfuncext_file = f'corr_func_ext_FerHofSqu1_chi_25_t1_1_V_10_Coulomb_1_n_2_361_nphi_4_19_LxMUC_1_Ly_19.dat'
    corrfuncext_path = os.path.join(corrfuncext_dir, corrfuncext_file)
    with open(corrfuncext_path, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        sites = np.arange(100*19)  # Landau gauge
        corr_func_ext = []
        for i, row in enumerate(plots):
            corr_func_ext.append(row[0])
    corr_func_ext = [float(i) for i in corr_func_ext]
    min_val = min(corr_func_ext)
    corr_func_ext = [float(i) - min_val for i in corr_func_ext]
    # offset
    corr_func_offset = []
    (m1, m1_err, c, c_err, r2_value1) = line_of_best_fit(sites[-200:], corr_func_ext[-200:])
    for i in corr_func_ext:
        corr_func_offset.append(i - c)
    ax3.plot(sites[1:], corr_func_offset[1:], '.-', c=f'C0', marker=markers[7], fillstyle='none', markersize=5)
    if corr_len_scale:
        ax3.axvline(56.758693816403486, c=f'C0', ls='--', zorder=-3)

    ax3.xaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax3.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    if corr_len_scale:
        ax3.set_xlim([0, 60])
        ax3.set_xticks(np.arange(0, 60, 20))
    else:
        ax3.set_xlim([0, 20])
        ax3.set_xticks(np.arange(0, 20, 10))
    # ax3.set_ylim(0)
    ax3.set_xlabel("$x$", fontsize=11)
    # ax3.set_ylabel("$\langle :\mathrel{\\rho_{0,0} \\rho_{x,0}}: \\rangle$", fontsize=11)
    # ax3.text(0.55 * 19, 0.021185, f"$\\nu={nu[0]}/{nu[1]}$", fontsize=11)

    ax3_2 = plt.subplot(right_inner_grid[0, 1])  # 071829 #################################################################################
    nu = (2, 19)

    corrfunc_file = f'corr_func_FerHofSqu1_chi_25_t1_1_V_10_Coulomb_1_n_2_361_nphi_4_19_LxMUC_1_Ly_19.dat'
    corrfunc_path = os.path.join(corrfunc_dir, corrfunc_file)
    with open(corrfunc_path, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        sites = np.arange(19)
        corr_func = []
        for i, row in enumerate(plots):
            if i == 0:
                corr_func = row
                break
    corr_func = [float(i) for i in corr_func]
    min_val = min(corr_func)
    corr_func = [float(i) - min_val for i in corr_func]
    sites_cont = [sites[-1], sites[-1] + 1]
    corr_func_cont = [corr_func[-1], corr_func[0]]
    ax3_2.plot(sites[1:], corr_func[1:], '.-', c=f'C0', marker=markers[7], fillstyle='none', markersize=5)
    # ax3.plot(sites_cont, corr_func_cont, '--', c=f'C0')

    ax3_2.yaxis.tick_right()
    ax3_2.yaxis.set_label_position("right")
    ax3_2.xaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax3_2.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax3_2.set_xlim([0, 19])
    ax3_2.set_xticks(np.arange(0, 19 + 0.1, 4.75))
    # ax3.set_ylim(0)
    ax3_2.set_xlabel("$y$", fontsize=11)
    # ax3_2.set_ylabel("$\langle :\mathrel{\\rho_{0,0} \\rho_{0,y}}: \\rangle$", fontsize=11)
    ax3_2.text(19/5, 0.021185, f"$\\nu={nu[0]}/{nu[1]}$", fontsize=11)

    ####################################################################################################################
    ####################################################################################################################

    fig.text(0.02, 0.925, "(c)", fontsize=12)
    fig.text(0.5, 0.925, "(d)", fontsize=12)

    plt.savefig("/home/bart/Documents/papers/FCI/corrfuncext_c5_analysis_combined.png", bbox_inches='tight', dpi=300)
    plt.show()
