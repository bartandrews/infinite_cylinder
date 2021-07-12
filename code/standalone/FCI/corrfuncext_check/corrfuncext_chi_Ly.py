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
plt.rc('text.latex', preamble=r'\usepackage{amsmath}\usepackage{nicefrac}')
# matplotlib.verbose.level = 'debug-annoying'


def line_of_best_fit(x_list, y_list):

    parameters, cov = np.polyfit(x_list, y_list, 1, cov=True)
    _, _, r_value, _, _ = stats.linregress(x_list, y_list)
    m, m_err, c, c_err = parameters[0], np.sqrt(cov[0][0]), parameters[1], np.sqrt(cov[1][1])
    r2_value = r_value*r_value

    return m, m_err, c, c_err, r2_value


if __name__ == '__main__':

    corr_len_scale = True

    fig = plt.figure(figsize=(6, 4))
    gs = gridspec.GridSpec(2, 2, hspace=1.4, wspace=0.6)

    # define a list of easily-visible markers
    markers = [(3, 0, 0), (4, 0, 0), (5, 0, 0), (6, 0, 0), (4, 1, 0), (5, 1, 0), (6, 1, 0),
               (3, 2, 0), (4, 2, 0), (5, 2, 0), (6, 2, 0), 'X', 'x', 'd', 'D', 'P',
               '$a$', '$b$', '$c$', '$d$', '$e$', '$f$', '$g$', '$h$', '$i$', '$j$', '$k$', '$l$', '$m$', '$n$', '$o$',
               '$p$', '$q$', '$r$', '$s$', '$t$', '$u$', '$v$', '$w$', '$x$', '$y$', '$z$']

    ####################################################################################################################

    corrfunc_dir = '/home/bart/PycharmProjects/infinite_cylinder/data/corr_func_ext/FerHofSqu1'

    ax1 = plt.subplot(gs[0])  # 071829 #################################################################################
    nu = (1, 3)

    for count, chi_val in enumerate([100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]):
        corrfunc_file = f'corr_func_ext_FerHofSqu1_chi_{chi_val}_t1_1_V_10_Coulomb_1_n_1_12_nphi_1_4_LxMUC_1_Ly_6.dat'
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

        ax1.plot(sites[1:], corr_func_offset[1:], '.-', c=f'C{count}', marker=markers[count], fillstyle='none', markersize=5, label=f'${int(chi_val/100)}$')

    if corr_len_scale:
        ax1.axvline(8.407415400658548, c=f'C0', ls='--', zorder=-3)  # chi=1000
        ax1.axvline(8.643074063066836, c=f'C1', ls='--', zorder=-3)  # chi=1000
        ax1.axvline(8.847635048461914, c=f'C2', ls='--', zorder=-3)  # chi=1000
        ax1.axvline(8.960241594318733, c=f'C3', ls='--', zorder=-3)  # chi=1000
        ax1.axvline(9.055314847464626, c=f'C4', ls='--', zorder=-3)  # chi=1000
        ax1.axvline(9.11603903541259, c=f'C5', ls='--', zorder=-3)  # chi=1000
        ax1.axvline(9.172761261043753, c=f'C6', ls='--', zorder=-3)  # chi=1000
        ax1.axvline(9.225647027419644, c=f'C7', ls='--', zorder=-3)  # chi=1000
        ax1.axvline(9.27112133610095, c=f'C8', ls='--', zorder=-3)  # chi=1000
        ax1.axvline(9.307335857576474, c=f'C9', ls='--', zorder=-3)  # chi=1000
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
    ax1.set_title("$\{C,\\nu,n_\phi,L_y\}=\{1,\\nicefrac{1}{3},\\nicefrac{1}{4},6\}$", fontsize=11)
    # ax2.text(0.675 * 18, 0.0059451, f"$\\nu={nu[0]}/{nu[1]}$", fontsize=11)
    # ax1.annotate(f"$\\nu={nu[0]}/{nu[1]}$", xy=(0.69, 0.87), xycoords='axes fraction', fontsize=10,
    #             verticalalignment='top',
    #             bbox=dict(boxstyle='round', facecolor='white', alpha=1))

    leg = ax1.legend(loc='upper center', handletextpad=0.5, handlelength=0, labelspacing=0.2, borderpad=0.35,
                     framealpha=1,
                     edgecolor='k', markerscale=0.8, fontsize=10, ncol=10, columnspacing=1, bbox_to_anchor=(1.3, 1.95),
                     title='$\chi/10^2$')
    leg.get_frame().set_linewidth(0.5)

    ax2 = plt.subplot(gs[1])  # 071829 #################################################################################
    nu = (1, 5)

    for count, chi_val in enumerate([100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]):
        corrfunc_file = f'corr_func_ext_FerHofSqu1_chi_{chi_val}_t1_1_V_10_Coulomb_1_n_1_55_nphi_6_11_LxMUC_1_Ly_10.dat'
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

        ax2.plot(sites[1:], corr_func_offset[1:], '.-', c=f'C{count}', marker=markers[count], fillstyle='none', markersize=5)

    if corr_len_scale:
        ax2.axvline(24.379929534532238, c='C0', ls='--', zorder=-3)  # chi=1000
        ax2.axvline(22.98406770307601, c='C1', ls='--', zorder=-3)  # chi=1000
        ax2.axvline(23.436468230949316, c='C2', ls='--', zorder=-3)  # chi=1000
        ax2.axvline(23.75000503480662, c='C3', ls='--', zorder=-3)  # chi=1000
        ax2.axvline(23.720028947117555, c='C4', ls='--', zorder=-3)  # chi=1000
        ax2.axvline(23.746961653980712, c='C5', ls='--', zorder=-3)  # chi=1000
        ax2.axvline(23.828870159959166, c='C6', ls='--', zorder=-3)  # chi=1000
        ax2.axvline(23.725361032665475, c='C7', ls='--', zorder=-3)  # chi=1000
        ax2.axvline(23.705015235245472, c='C8', ls='--', zorder=-3)  # chi=1000
        ax2.axvline(23.947622591804954, c='C9', ls='--', zorder=-3)  # chi=1000
    # ax2.plot(sites_cont, corr_func_cont, '--', c=f'C0')

    ax2.xaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax2.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    if corr_len_scale:
        ax2.set_xlim([0, 30])
        ax2.set_xticks(np.arange(0, 31, 5))
    else:
        ax2.set_xlim([0, 20])
        ax2.set_xticks(np.arange(0, 21, 10))
    # ax2.set_xlim([0, 100])
    # ax2.set_xticks(np.arange(0, 19 + 0.1, 4.75))
    # ax2.set_ylim(0)
    ax2.set_xlabel("$x$", fontsize=11)
    ax2.set_ylabel("$g(x)$", fontsize=11)
    ax2.set_title("$\{C,\\nu,n_\phi,L_y\}=\{2,\\nicefrac{1}{5},\\nicefrac{6}{11},10\}$", fontsize=11)
    # ax2.text(0.55 * 19, 0.021185, f"$\\nu={nu[0]}/{nu[1]}$", fontsize=11)
    # ax2.annotate(f"$\\nu={nu[0]}/{nu[1]}$", xy=(0.69, 0.87), xycoords='axes fraction', fontsize=10,
    #              verticalalignment='top',
    #              bbox=dict(boxstyle='round', facecolor='white', alpha=1))

    ax3 = plt.subplot(gs[2])  # 071829 #################################################################################
    nu = (1, 3)

    for count, Ly_val in enumerate([6, 9, 12]):
        corrfunc_file = f'corr_func_ext_FerHofSqu1_chi_1000_t1_1_V_10_Coulomb_1_n_1_12_nphi_1_4_LxMUC_1_Ly_{Ly_val}.dat'
        corrfunc_path = os.path.join(corrfunc_dir, corrfunc_file)
        with open(corrfunc_path, 'r') as csvfile:
            plots = csv.reader(csvfile, delimiter='\t')
            sites = np.arange(100 * 4)  # Landau gauge
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

        ax3.plot(sites[1:], corr_func_offset[1:], '.-', c=f'C{count}', marker=markers[count], fillstyle='none',
                 markersize=5, label=f'${int(Ly_val)}$')

    if corr_len_scale:
        ax3.axvline(9.307335857576474, c=f'C0', ls='--', zorder=-3)  # chi=1000
        ax3.axvline(11.536121582309628, c=f'C1', ls='--', zorder=-3)  # chi=1000
        ax3.axvline(13.072935466607136, c=f'C2', ls='--', zorder=-3)  # chi=1000
    # ax3.plot(sites_cont, corr_func_cont, '--', c=f'C0')

    ax3.xaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax3.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    if corr_len_scale:
        ax3.set_xlim([0, 14])
        ax3.set_xticks(np.arange(0, 15, 2))
    else:
        ax3.set_xlim([0, 20])
        ax3.set_xticks(np.arange(0, 21, 10))
    # ax3.set_xlim([0, 100])
    # ax3.set_xticks(np.arange(0, 18 + 0.1, 2))
    # ax3.set_ylim(0)
    ax3.set_xlabel("$x$", fontsize=11)
    ax3.set_ylabel("$g(x)$", fontsize=11)
    ax3.set_title("$\{C,\\nu,n_\phi,\chi\}=\{1,\\nicefrac{1}{3},\\nicefrac{1}{4},10^3\}$", fontsize=11)
    # ax3.text(0.675 * 18, 0.0059451, f"$\\nu={nu[0]}/{nu[1]}$", fontsize=11)
    # ax3.annotate(f"$\\nu={nu[0]}/{nu[1]}$", xy=(0.69, 0.87), xycoords='axes fraction', fontsize=10,
    #              verticalalignment='top',
    #              bbox=dict(boxstyle='round', facecolor='white', alpha=1))

    leg = ax3.legend(loc='upper center', handletextpad=0.5, handlelength=0, labelspacing=0.2, borderpad=0.35,
                     framealpha=1,
                     edgecolor='k', markerscale=0.8, fontsize=10, ncol=10, columnspacing=1, bbox_to_anchor=(0.5, 1.9),
                     title='$L_y$')
    leg.get_frame().set_linewidth(0.5)

    ax4 = plt.subplot(gs[3])  # 071829 #################################################################################
    nu = (1, 5)

    for count, Ly_val in enumerate([10, 15, 20]):
        corrfunc_file = f'corr_func_ext_FerHofSqu1_chi_100_t1_1_V_10_Coulomb_1_n_1_55_nphi_6_11_LxMUC_1_Ly_{Ly_val}.dat'
        corrfunc_path = os.path.join(corrfunc_dir, corrfunc_file)
        with open(corrfunc_path, 'r') as csvfile:
            plots = csv.reader(csvfile, delimiter='\t')
            sites = np.arange(100 * 11)  # Landau gauge
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

        ax4.plot(sites[1:], corr_func_offset[1:], '.-', c=f'C{count}', marker=markers[count], fillstyle='none',
                 markersize=5, label=f'${int(Ly_val)}$')

    if corr_len_scale:
        ax4.axvline(24.379929534532238, c=f'C0', ls='--', zorder=-3)  # chi=100
        ax4.axvline(38.52972334472212, c=f'C1', ls='--', zorder=-3)  # chi=100
        ax4.axvline(33.59220089136867, c=f'C2', ls='--', zorder=-3)  # chi=100
    # ax4.plot(sites_cont, corr_func_cont, '--', c=f'C0')

    ax4.xaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax4.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    if corr_len_scale:
        ax4.set_xlim([0, 40])
        ax4.set_xticks(np.arange(0, 41, 5))
    else:
        ax4.set_xlim([0, 20])
        ax4.set_xticks(np.arange(0, 21, 10))
    # ax4.set_xlim([0, 100])
    # ax4.set_xticks(np.arange(0, 19 + 0.1, 4.75))
    # ax4.set_ylim(0)
    ax4.set_xlabel("$x$", fontsize=11)
    ax4.set_ylabel("$g(x)$", fontsize=11)
    ax4.set_title("$\{C,\\nu,n_\phi,\chi\}=\{2,\\nicefrac{1}{5},\\nicefrac{6}{11},10^3\}$", fontsize=11)
    # ax4.text(0.55 * 19, 0.021185, f"$\\nu={nu[0]}/{nu[1]}$", fontsize=11)
    # ax4.annotate(f"$\\nu={nu[0]}/{nu[1]}$", xy=(0.69, 0.87), xycoords='axes fraction', fontsize=10,
    #              verticalalignment='top',
    #              bbox=dict(boxstyle='round', facecolor='white', alpha=1))

    leg = ax4.legend(loc='upper center', handletextpad=0.5, handlelength=0, labelspacing=0.2, borderpad=0.35,
                     framealpha=1,
                     edgecolor='k', markerscale=0.8, fontsize=10, ncol=10, columnspacing=1, bbox_to_anchor=(0.5, 1.9),
                     title='$L_y$')
    leg.get_frame().set_linewidth(0.5)

    ####################################################################################################################
    ####################################################################################################################

    fig.text(0.02, 1.02, "(a)", fontsize=12)
    fig.text(0.02, 0.47, "(b)", fontsize=12)

    plt.savefig("/home/bart/Documents/papers/FCI/corrfuncext_chi_Ly.png", bbox_inches='tight', dpi=300)
    plt.show()
