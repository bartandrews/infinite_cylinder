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

    fig = plt.figure(figsize=(6, 6))
    gs = gridspec.GridSpec(4, 2, hspace=0.6, wspace=0.6)

    # define a list of easily-visible markers
    markers = [(3, 0, 0), (4, 0, 0), (5, 0, 0), (6, 0, 0), (4, 1, 0), (5, 1, 0), (6, 1, 0),
               (3, 2, 0), (4, 2, 0), (5, 2, 0), (6, 2, 0), 'X', 'x', 'd', 'D', 'P',
               '$a$', '$b$', '$c$', '$d$', '$e$', '$f$', '$g$', '$h$', '$i$', '$j$', '$k$', '$l$', '$m$', '$n$', '$o$',
               '$p$', '$q$', '$r$', '$s$', '$t$', '$u$', '$v$', '$w$', '$x$', '$y$', '$z$']

    ####################################################################################################################

    for filling in range(8):
        if filling == 0:
            ax1 = plt.subplot(gs[0])
            nu = (1, 3)
            Ly_val = 6
            ax1.text(Ly_val/3, 0.05, f"$\\nu={nu[0]}/{nu[1]}$", fontsize=11)
            q_list = [3, 4, 5, 6, 7, 8]
        elif filling == 1:
            ax1 = plt.subplot(gs[1])
            nu = (2, 3)
            Ly_val = 6
            ax1.text(Ly_val/3, 0.05, f"$\\nu={nu[0]}/{nu[1]}$", fontsize=11)
            q_list = [5, 6, 7, 8]
        elif filling == 2:
            ax1 = plt.subplot(gs[2])
            nu = (2, 5)
            Ly_val = 10
            ax1.text(Ly_val / 3, 0.025, f"$\\nu={nu[0]}/{nu[1]}$", fontsize=11)
            q_list = [3, 4, 5, 6, 7, 8]
        elif filling == 3:
            ax1 = plt.subplot(gs[3])
            nu = (3, 5)
            Ly_val = 10
            ax1.text(Ly_val / 3, 0.025, f"$\\nu={nu[0]}/{nu[1]}$", fontsize=11)
            q_list = [5, 6, 7, 8]
        elif filling == 4:
            ax1 = plt.subplot(gs[4])
            nu = (3, 7)
            Ly_val = 14
            ax1.text(Ly_val / 3, 0.025, f"$\\nu={nu[0]}/{nu[1]}$", fontsize=11)
            q_list = [5, 6, 7, 8]
        elif filling == 5:
            ax1 = plt.subplot(gs[5])
            nu = (4, 7)
            Ly_val = 14
            ax1.text(Ly_val / 3, 0.07, f"$\\nu={nu[0]}/{nu[1]}$", fontsize=11)
            q_list = [6, 7, 8]
        elif filling == 6:
            ax1 = plt.subplot(gs[6])
            nu = (4, 9)
            Ly_val = 18
            ax1.text(Ly_val/3, 0.016, f"$\\nu={nu[0]}/{nu[1]}$", fontsize=11)
            q_list = [8]
        elif filling == 7:
            ax1 = plt.subplot(gs[7])
            nu = (5, 9)
            Ly_val = 18
            ax1.text(Ly_val / 3, 0.04, f"$\\nu={nu[0]}/{nu[1]}$", fontsize=11)
            q_list = [8]

        for q in q_list:

            nphi = (1, q)
            n = Fraction(nu[0]*nphi[0], nu[1]*nphi[1]).limit_denominator(100)

            corrfunc_dir = '/home/bart/PycharmProjects/infinite_cylinder/data/corr_func/FerHofSqu1'
            corrfunc_file = f'corr_func_FerHofSqu1_chi_250_t1_1_V_10_Coulomb_1_' \
                         f'n_{n.numerator}_{n.denominator}_nphi_{nphi[0]}_{nphi[1]}_LxMUC_1_Ly_{Ly_val}.dat'
            corrfunc_path = os.path.join(corrfunc_dir, corrfunc_file)

            # extract data from file
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
            corr_func = [float(i)-min_val for i in corr_func]

            sites_cont = [sites[-1], sites[-1]+1]
            corr_func_cont = [corr_func[-1], corr_func[0]]

            print(sites)
            print(corr_func)

            ax1.plot(sites, corr_func, '.-', c=f'C{q-3}', marker=markers[q-3], fillstyle='none', markersize=2.5, label=f'${nphi[0]}/{nphi[1]}$')
            ax1.plot(sites_cont, corr_func_cont, '--', c=f'C{q-3}')

            if filling == 0:
                leg = ax1.legend(loc='upper center', handletextpad=0.3, handlelength=1, labelspacing=0.1, borderpad=0.3,
                                 framealpha=1,
                                 edgecolor='k', markerscale=2, fontsize=10, ncol=6, columnspacing=0.5,
                                 bbox_to_anchor=(1.2, 1.9), title='$n_\\phi$', title_fontsize=11)
                leg.get_frame().set_linewidth(0.5)

            ax1.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
            ax1.set_xlim([0, Ly_val])
            if Ly_val < 14:
                ax1.set_xticks(np.arange(0, Ly_val+0.1, 1))
            else:
                ax1.set_xticks(np.arange(0, Ly_val + 0.1, 2))
            # ax1.set_yticks([-nu[0], 0])
            if q == 8:
                ax1.set_ylim(bottom=0)
            ax1.set_xlabel("$y$", fontsize=11)
            ax1.set_ylabel("$\langle :\mathrel{\\rho_{0,0} \\rho_{0,y}}: \\rangle$", fontsize=11)
            # if q == 3:
            #     ax1.text(0.05*nu[1], -0.9*nu[0], f"$\\nu={nu[0]}/{nu[1]}$", fontsize=11)

    ####################################################################################################################
    ####################################################################################################################

    fig.text(0, 0.87, "(a)", fontsize=12)
    fig.text(0.48, 0.87, "(b)", fontsize=12)
    fig.text(0, 0.655, "(c)", fontsize=12)
    fig.text(0.48, 0.655, "(d)", fontsize=12)
    fig.text(0, 0.445, "(e)", fontsize=12)
    fig.text(0.48, 0.445, "(f)", fontsize=12)
    fig.text(0, 0.231, "(g)", fontsize=12)
    fig.text(0.48, 0.231, "(h)", fontsize=12)

    plt.savefig("/home/bart/Documents/papers/FCI/corrfunc_analysis.png", bbox_inches='tight', dpi=300)
    plt.show()
