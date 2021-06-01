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
import math


plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}')
# matplotlib.verbose.level = 'debug-annoying'


def line_of_best_fit(x_list, y_list):

    parameters, cov = np.polyfit(x_list, y_list, 1, cov=True)
    _, _, r_value, _, _ = stats.linregress(x_list, y_list)
    m, m_err, c, c_err = parameters[0], np.sqrt(cov[0][0]), parameters[1], np.sqrt(cov[1][1])
    r2_value = r_value*r_value

    return m, m_err, c, c_err, r2_value


def roundup(x):
    return int(math.ceil(x / 10.0)) * 10


if __name__ == '__main__':

    corr_len_scale = False

    fig = plt.figure(figsize=(6, 6.5))  # 6, 6
    outer_grid = gridspec.GridSpec(4, 2, hspace=0.6, wspace=0.5)  # 0.6, 0.6
    very_top_left = outer_grid[0, 0]
    very_top_right = outer_grid[0, 1]
    top_left = outer_grid[1, 0]
    top_right = outer_grid[1, 1]
    bottom_left = outer_grid[2, 0]
    bottom_right = outer_grid[2, 1]
    very_bottom_left = outer_grid[3, 0]
    very_bottom_right = outer_grid[3, 1]
    very_top_left_inner_grid = gridspec.GridSpecFromSubplotSpec(1, 2, very_top_left, wspace=0)
    very_top_right_inner_grid = gridspec.GridSpecFromSubplotSpec(1, 2, very_top_right, wspace=0)
    top_left_inner_grid = gridspec.GridSpecFromSubplotSpec(1, 2, top_left, wspace=0)
    top_right_inner_grid = gridspec.GridSpecFromSubplotSpec(1, 2, top_right, wspace=0)
    bottom_left_inner_grid = gridspec.GridSpecFromSubplotSpec(1, 2, bottom_left, wspace=0)
    bottom_right_inner_grid = gridspec.GridSpecFromSubplotSpec(1, 2, bottom_right, wspace=0)
    very_bottom_left_inner_grid = gridspec.GridSpecFromSubplotSpec(1, 2, very_bottom_left, wspace=0)
    very_bottom_right_inner_grid = gridspec.GridSpecFromSubplotSpec(1, 2, very_bottom_right, wspace=0)

    # define a list of easily-visible markers
    markers = [(3, 0, 0), (4, 0, 0), (5, 0, 0), (6, 0, 0), (4, 1, 0), (5, 1, 0), (6, 1, 0),
               (3, 2, 0), (4, 2, 0), (5, 2, 0), (6, 2, 0), 'X', 'x', 'd', 'D', 'P',
               '$a$', '$b$', '$c$', '$d$', '$e$', '$f$', '$g$', '$h$', '$i$', '$j$', '$k$', '$l$', '$m$', '$n$', '$o$',
               '$p$', '$q$', '$r$', '$s$', '$t$', '$u$', '$v$', '$w$', '$x$', '$y$', '$z$']

    ####################################################################################################################

    for filling in range(8):
        if filling == 0:
            ax1 = plt.subplot(very_top_left_inner_grid[0, 0])
            ax2 = plt.subplot(very_top_left_inner_grid[0, 1])
            nu = (1, 3)
            Ly_val = 6
            ax2.text(Ly_val/4.5, 0.05, f"$\\nu={nu[0]}/{nu[1]}$", fontsize=11)
            q_list = [3, 4, 5, 6, 7, 8]
        elif filling == 1:
            ax1 = plt.subplot(very_top_right_inner_grid[0, 0])
            ax2 = plt.subplot(very_top_right_inner_grid[0, 1])
            nu = (2, 3)
            Ly_val = 6
            ax2.text(Ly_val/4.5, 0.05, f"$\\nu={nu[0]}/{nu[1]}$", fontsize=11)
            q_list = [5, 6, 7, 8]
        elif filling == 2:
            ax1 = plt.subplot(top_left_inner_grid[0, 0])
            ax2 = plt.subplot(top_left_inner_grid[0, 1])
            nu = (2, 5)
            Ly_val = 10
            ax2.text(Ly_val / 4.5, 0.025, f"$\\nu={nu[0]}/{nu[1]}$", fontsize=11)
            q_list = [3, 4, 5, 6, 7, 8]
        elif filling == 3:
            ax1 = plt.subplot(top_right_inner_grid[0, 0])
            ax2 = plt.subplot(top_right_inner_grid[0, 1])
            nu = (3, 5)
            Ly_val = 10
            ax2.text(Ly_val / 4.5, 0.025, f"$\\nu={nu[0]}/{nu[1]}$", fontsize=11)
            q_list = [5, 6, 7, 8]
        elif filling == 4:
            ax1 = plt.subplot(bottom_left_inner_grid[0, 0])
            ax2 = plt.subplot(bottom_left_inner_grid[0, 1])
            nu = (3, 7)
            Ly_val = 14
            ax2.text(Ly_val / 4.5, 0.025, f"$\\nu={nu[0]}/{nu[1]}$", fontsize=11)
            q_list = [5, 6, 7, 8]
        elif filling == 5:
            ax1 = plt.subplot(bottom_right_inner_grid[0, 0])
            ax2 = plt.subplot(bottom_right_inner_grid[0, 1])
            nu = (4, 7)
            Ly_val = 14
            ax2.text(Ly_val / 4.5, 0.07, f"$\\nu={nu[0]}/{nu[1]}$", fontsize=11)
            q_list = [6, 7, 8]
        elif filling == 6:
            ax1 = plt.subplot(very_bottom_left_inner_grid[0, 0])
            ax2 = plt.subplot(very_bottom_left_inner_grid[0, 1])
            nu = (4, 9)
            Ly_val = 18
            ax2.text(Ly_val/4.5, 0.016, f"$\\nu={nu[0]}/{nu[1]}$", fontsize=11)
            q_list = [8]
        elif filling == 7:
            ax1 = plt.subplot(very_bottom_right_inner_grid[0, 0])
            ax2 = plt.subplot(very_bottom_right_inner_grid[0, 1])
            nu = (5, 9)
            Ly_val = 18
            ax2.text(Ly_val / 4.5, 0.04, f"$\\nu={nu[0]}/{nu[1]}$", fontsize=11)
            q_list = [8]

        corr_len = []
        for q in q_list:

            nphi = (1, q)
            n = Fraction(nu[0]*nphi[0], nu[1]*nphi[1]).limit_denominator(100)

            corrfuncext_dir = '/home/bart/PycharmProjects/infinite_cylinder/data/corr_func_ext/FerHofSqu1'
            corrfuncext_file = f'corr_func_ext_FerHofSqu1_chi_250_t1_1_V_10_Coulomb_1_' \
                               f'n_{n.numerator}_{n.denominator}_nphi_{nphi[0]}_{nphi[1]}_LxMUC_1_Ly_{Ly_val}.dat'
            corrfuncext_path = os.path.join(corrfuncext_dir, corrfuncext_file)

            # extract corrfuncext data from file
            with open(corrfuncext_path, 'r') as csvfile:
                plots = csv.reader(csvfile, delimiter='\t')
                sites = np.arange(100*nphi[1])  # Landau gauge
                corr_func_ext = []
                for i, row in enumerate(plots):
                    corr_func_ext.append(row[0])

            corr_func_ext = [float(i) for i in corr_func_ext]
            min_val = min(corr_func_ext)
            corr_func_ext = [float(i)-min_val for i in corr_func_ext]

            # offset
            corr_func_offset = []
            (m1, m1_err, c, c_err, r2_value1) = line_of_best_fit(sites[-200:], corr_func_ext[-200:])
            for i in corr_func_ext:
                corr_func_offset.append(i-c)

            ax1.plot(sites[1:], corr_func_offset[1:], '.-', c=f'C{q-3}', marker=markers[q-3], fillstyle='none', markersize=2.5, label=f'${nphi[0]}/{nphi[1]}$')

            if filling == 0:
                leg = ax1.legend(loc='upper center', handletextpad=0.3, handlelength=1, labelspacing=0.1, borderpad=0.3,
                                 framealpha=1,
                                 edgecolor='k', markerscale=2, fontsize=10, ncol=6, columnspacing=0.5,
                                 bbox_to_anchor=(2.5, 2), title='$n_\\phi$', title_fontsize=11)
                leg.get_frame().set_linewidth(0.5)

            ax1.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
            # if Ly_val < 14:
            #     ax1.set_xticks(np.arange(0, Ly_val+0.1, 1))
            # else:
            #     ax1.set_xticks(np.arange(0, Ly_val + 0.1, 2))
            # ax1.set_yticks([-nu[0], 0])
            # if q == 8:
            #     ax1.set_ylim(bottom=0)
            ax1.set_xlabel("$x$", fontsize=11)
            # ax1.set_ylabel("$\langle :\mathrel{\\rho_{0,0} \\rho_{x,0}}: \\rangle$", fontsize=11)
            # if q == 3:
            #     ax1.text(0.05*nu[1], -0.9*nu[0], f"$\\nu={nu[0]}/{nu[1]}$", fontsize=11)

            if corr_len_scale:
                corrlen_dir = '/home/bart/PycharmProjects/infinite_cylinder/logs/observables/FerHofSqu1'
                corrlen_file = f'log_observables_FerHofSqu1_chi_250_t1_1_V_10_Coulomb_1_' \
                               f'n_{n.numerator}_{n.denominator}_nphi_{nphi[0]}_{nphi[1]}_LxMUC_1_Ly_{Ly_val}.dat'
                corrlen_path = os.path.join(corrlen_dir, corrlen_file)

                # extract corrlen from file
                for line in list(open(corrlen_path)):  # find the max SnV
                    if "xi" in line:
                        xi_line = line.rstrip()
                        xi = float(xi_line.split()[-1])
                        ax1.axvline(xi, c=f'C{q-3}', ls='--', zorder=-3)
                        corr_len.append(xi)
                        break

                if q==8:
                    ax1.set_xlim([0, roundup(np.max(corr_len))])
                    ax1.set_xticks(np.arange(0, roundup(np.max(corr_len)), 10))
            else:
                ax1.set_xlim([0, 20])
                ax1.set_xticks(np.arange(0, 20, 10))

            corrfunc_dir = '/home/bart/PycharmProjects/infinite_cylinder/data/corr_func/FerHofSqu1'
            corrfunc_file = f'corr_func_FerHofSqu1_chi_250_t1_1_V_10_Coulomb_1_' \
                            f'n_{n.numerator}_{n.denominator}_nphi_{nphi[0]}_{nphi[1]}_LxMUC_1_Ly_{Ly_val}.dat'
            corrfunc_path = os.path.join(corrfunc_dir, corrfunc_file)

            # extract corrfunc data from file
            with open(corrfunc_path, 'r') as csvfile:
                plots = csv.reader(csvfile, delimiter='\t')
                sites = np.arange(Ly_val)
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

            ax2.plot(sites, corr_func, '.-', c=f'C{q - 3}', marker=markers[q - 3], fillstyle='none', markersize=2.5,
                     label=f'${nphi[0]}/{nphi[1]}$')
            ax2.plot(sites_cont, corr_func_cont, '--', c=f'C{q - 3}')

            ax2.yaxis.tick_right()
            ax2.yaxis.set_label_position("right")
            ax2.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
            ax2.set_xlim([0, Ly_val])
            if Ly_val < 14:
                ax2.set_xticks(np.arange(0, Ly_val + 0.1, 2))
            elif Ly_val == 14:
                ax2.set_xticks(np.arange(0, Ly_val + 0.1, 7))
            elif Ly_val == 18:
                ax2.set_xticks(np.arange(0, Ly_val + 0.1, 6))
            # ax1.set_yticks([-nu[0], 0])
            if q == 8:
                ax2.set_ylim(bottom=0)
            ax2.set_xlabel("$y$", fontsize=11)
            # ax2.set_ylabel("$\langle :\mathrel{\\rho_{0,0} \\rho_{0,y}}: \\rangle$", fontsize=11)

    ####################################################################################################################
    ####################################################################################################################

    fig.text(0.02, 0.88, "(a)", fontsize=12)
    fig.text(0.495, 0.88, "(b)", fontsize=12)
    fig.text(0.02, 0.665, "(c)", fontsize=12)
    fig.text(0.495, 0.665, "(d)", fontsize=12)
    fig.text(0.02, 0.455, "(e)", fontsize=12)
    fig.text(0.495, 0.455, "(f)", fontsize=12)
    fig.text(0.02, 0.241, "(g)", fontsize=12)
    fig.text(0.495, 0.241, "(h)", fontsize=12)

    plt.savefig("/home/bart/Documents/papers/FCI/corrfuncext_analysis_combined.png", bbox_inches='tight', dpi=300)
    plt.show()
