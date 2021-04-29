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

if __name__ == '__main__':

    fig = plt.figure(figsize=(6, 1.1))
    gs = gridspec.GridSpec(1, 2, hspace=0.6, wspace=0.6)

    # define a list of easily-visible markers
    markers = [(3, 0, 0), (4, 0, 0), (5, 0, 0), (6, 0, 0), (4, 1, 0), (5, 1, 0), (6, 1, 0),
               (3, 2, 0), (4, 2, 0), (5, 2, 0), (6, 2, 0), 'X', 'x', 'd', 'D', 'P',
               '$a$', '$b$', '$c$', '$d$', '$e$', '$f$', '$g$', '$h$', '$i$', '$j$', '$k$', '$l$', '$m$', '$n$', '$o$',
               '$p$', '$q$', '$r$', '$s$', '$t$', '$u$', '$v$', '$w$', '$x$', '$y$', '$z$']

    ####################################################################################################################

    phiflow_dir = '/home/bart/PycharmProjects/infinite_cylinder/data/charge_pump/FerHofSqu1'
    corrfunc_dir = '/home/bart/PycharmProjects/infinite_cylinder/data/corr_func/FerHofSqu1'

    ax2 = plt.subplot(gs[0])  # 071829 #################################################################################
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
    ax2.plot(sites[1:], corr_func[1:], '.-', c=f'C0', marker=markers[1], fillstyle='none', markersize=5)
    # ax2.plot(sites_cont, corr_func_cont, '--', c=f'C0')

    ax2.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax2.set_xlim([0, 18])
    ax2.set_xticks(np.arange(0, 18 + 0.1, 2))
    # ax2.set_ylim(0)
    ax2.set_xlabel("$y$", fontsize=11)
    ax2.set_ylabel("$\langle :\mathrel{\\rho_{0,0} \\rho_{0,y}}: \\rangle$", fontsize=11)
    ax2.text(0.675 * 18, 0.0059451, f"$\\nu={nu[0]}/{nu[1]}$", fontsize=11)

    ax3 = plt.subplot(gs[1])  # 071829 #################################################################################
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
    ax3.plot(sites[1:], corr_func[1:], '.-', c=f'C0', marker=markers[7], fillstyle='none', markersize=5)
    # ax3.plot(sites_cont, corr_func_cont, '--', c=f'C0')

    ax3.xaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax3.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax3.set_xlim([0, 19])
    ax3.set_xticks(np.arange(0, 19 + 0.1, 4.75))
    # ax3.set_ylim(0)
    ax3.set_xlabel("$y$", fontsize=11)
    ax3.set_ylabel("$\langle :\mathrel{\\rho_{0,0} \\rho_{0,y}}: \\rangle$", fontsize=11)
    ax3.text(0.55 * 19, 0.021185, f"$\\nu={nu[0]}/{nu[1]}$", fontsize=11)

    # chi_legend_elements = [Patch(facecolor='C1', label='$50$'), Patch(facecolor='C2', label='$100$')]
    # leg2 = ax3.legend(handles=chi_legend_elements, loc='center', handletextpad=0.3, handlelength=1, labelspacing=0.1,
    #                   borderpad=0.3, framealpha=1, edgecolor='k', markerscale=2, fontsize=10, ncol=2, columnspacing=0.5,
    #                   bbox_to_anchor=(-0.65, 0.4), title='$\chi$', title_fontsize=11)
    # leg2.get_frame().set_linewidth(0.5)

    ####################################################################################################################
    ####################################################################################################################

    fig.text(-0.02, 0.925, "(c)", fontsize=12)
    fig.text(0.46, 0.925, "(d)", fontsize=12)

    plt.savefig("/home/bart/Documents/papers/FCI/corrfunc_c5_analysis.png", bbox_inches='tight', dpi=300)
    plt.show()
