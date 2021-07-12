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

    fig = plt.figure(figsize=(6, 1.1))
    gs = gridspec.GridSpec(1, 2, hspace=0.6, wspace=0.6)

    # define a list of easily-visible markers
    markers = [(3, 0, 0), (4, 0, 0), (5, 0, 0), (6, 0, 0), (4, 1, 0), (5, 1, 0), (6, 1, 0),
               (3, 2, 0), (4, 2, 0), (5, 2, 0), (6, 2, 0), 'X', 'x', 'd', 'D', 'P',
               '$a$', '$b$', '$c$', '$d$', '$e$', '$f$', '$g$', '$h$', '$i$', '$j$', '$k$', '$l$', '$m$', '$n$', '$o$',
               '$p$', '$q$', '$r$', '$s$', '$t$', '$u$', '$v$', '$w$', '$x$', '$y$', '$z$']

    ####################################################################################################################

    corrfunc_dir = '/home/bart/PycharmProjects/infinite_cylinder/data/corr_func_ext/FerHofSqu1'

    ax1 = plt.subplot(gs[0])  # 071829 #################################################################################
    nu = (1, 3)

    for count, Ly_val in enumerate([6, 9, 12]):
        corrfunc_file = f'corr_func_ext_FerHofSqu1_chi_1000_t1_1_V_10_Coulomb_1_n_1_12_nphi_1_4_LxMUC_1_Ly_{Ly_val}.dat'
        corrfunc_path = os.path.join(corrfunc_dir, corrfunc_file)
        with open(corrfunc_path, 'r') as csvfile:
            plots = csv.reader(csvfile, delimiter='\t')
            sites = np.arange(100*4)  # Landau gauge
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

        ax1.plot(sites[1:], corr_func_offset[1:], '.-', c=f'C{count}', marker=markers[count], fillstyle='none', markersize=5, label=f'{int(Ly_val)}')

    if corr_len_scale:
        ax1.axvline(9.307335857576474, c=f'C0', ls='--', zorder=-3)
    # ax1.plot(sites_cont, corr_func_cont, '--', c=f'C0')

    ax1.xaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax1.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    if corr_len_scale:
        ax1.set_xlim([0, 10])
        ax1.set_xticks(np.arange(0, 11, 1))
    else:
        ax1.set_xlim([0, 20])
        ax1.set_xticks(np.arange(0, 21, 10))
    # ax2.set_xlim([0, 100])
    # ax2.set_xticks(np.arange(0, 18 + 0.1, 2))
    # ax2.set_ylim(0)
    ax1.set_xlabel("$x$", fontsize=11)
    ax1.set_ylabel("$g(x)$", fontsize=11)
    # ax2.text(0.675 * 18, 0.0059451, f"$\\nu={nu[0]}/{nu[1]}$", fontsize=11)
    ax1.annotate(f"$\\nu={nu[0]}/{nu[1]}$", xy=(0.69, 0.87), xycoords='axes fraction', fontsize=10,
                verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='white', alpha=1))

    leg = ax1.legend(loc='upper center', handletextpad=0.5, handlelength=0, labelspacing=0.2, borderpad=0.35,
                     framealpha=1,
                     edgecolor='k', markerscale=0.8, fontsize=10, ncol=10, columnspacing=1, bbox_to_anchor=(1.19, 1.8),
                     title='$L_y$')

    ax2 = plt.subplot(gs[1])  # 071829 #################################################################################
    nu = (1, 3)

    for count, Ly_val in enumerate([6, 9, 12]):
        corrfunc_file = f'corr_func_ext_FerHofSqu1_chi_1000_t1_1_V_10_Coulomb_1_n_1_12_nphi_1_4_LxMUC_1_Ly_{Ly_val}.dat'
        corrfunc_path = os.path.join(corrfunc_dir, corrfunc_file)
        with open(corrfunc_path, 'r') as csvfile:
            plots = csv.reader(csvfile, delimiter='\t')
            sites = np.arange(100*4)  # Landau gauge
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

        ax2.plot(sites[1:], corr_func_offset[1:], '.-', c=f'C{count}', marker=markers[count], fillstyle='none', markersize=5)

    if corr_len_scale:
        ax2.axvline(23.947622591804954, c='C0', ls='--', zorder=-3)
    # ax2.plot(sites_cont, corr_func_cont, '--', c=f'C0')

    ax2.xaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax2.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    if corr_len_scale:
        ax2.set_xlim([0, 300])
        ax2.set_xticks(np.arange(0, 301, 50))
    else:
        ax2.set_xlim([0, 20])
        ax2.set_xticks(np.arange(0, 21, 10))
    # ax2.set_xlim([0, 100])
    # ax2.set_xticks(np.arange(0, 19 + 0.1, 4.75))
    # ax2.set_ylim(0)
    ax2.set_xlabel("$x$", fontsize=11)
    ax2.set_ylabel("$g(x)$", fontsize=11)
    # ax2.text(0.55 * 19, 0.021185, f"$\\nu={nu[0]}/{nu[1]}$", fontsize=11)
    ax2.annotate(f"$\\nu={nu[0]}/{nu[1]}$", xy=(0.69, 0.87), xycoords='axes fraction', fontsize=10,
                 verticalalignment='top',
                 bbox=dict(boxstyle='round', facecolor='white', alpha=1))

    ####################################################################################################################
    ####################################################################################################################

    # fig.text(-0.01, 0.925, "(a)", fontsize=12)
    # fig.text(0.46-0.005, 0.925, "(b)", fontsize=12)

    plt.savefig("/home/bart/Documents/papers/FCI/corrfuncext_Ly.png", bbox_inches='tight', dpi=300)
    plt.show()
