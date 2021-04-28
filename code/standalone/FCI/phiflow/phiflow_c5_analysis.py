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
    gs = gridspec.GridSpec(1, 2, hspace=0.6, wspace=0.4)

    # define a list of easily-visible markers
    markers = [(3, 0, 0), (4, 0, 0), (5, 0, 0), (6, 0, 0), (4, 1, 0), (5, 1, 0), (6, 1, 0),
               (3, 2, 0), (4, 2, 0), (5, 2, 0), (6, 2, 0), 'X', 'x', 'd', 'D', 'P',
               '$a$', '$b$', '$c$', '$d$', '$e$', '$f$', '$g$', '$h$', '$i$', '$j$', '$k$', '$l$', '$m$', '$n$', '$o$',
               '$p$', '$q$', '$r$', '$s$', '$t$', '$u$', '$v$', '$w$', '$x$', '$y$', '$z$']

    ####################################################################################################################

    phiflow_dir = '/home/bart/PycharmProjects/infinite_cylinder/data/charge_pump/FerHofSqu1'

    ax2 = plt.subplot(gs[0])  # 071829 #################################################################################
    nu = (1, 9)

    phiflow_file = f'charge_pump_FerHofSqu1_chi_25_t1_1_V_10_Coulomb_1_n_1_216_nphi_5_24_LxMUC_1_Ly_18_phi_0_9_91.dat.modified'
    phiflow_path = os.path.join(phiflow_dir, phiflow_file)
    with open(phiflow_path, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        phi = []
        charge = []
        for row in plots:
            phi.append(float(row[0]))
            charge.append(float(row[1]))
    charge = [i - charge[0] for i in charge]
    ax2.plot(phi, charge, '.', c=f'C0', marker=markers[1], fillstyle='none', markersize=5, markeredgewidth=0.2)

    ax2.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax2.set_xlim([0, nu[1]])
    ax2.set_xticks(np.arange(0, nu[1] + 0.1, 1))
    ax2.set_yticks([-5 * nu[0], 0])
    ax2.set_xlabel("$\\Phi_x / 2\pi$", fontsize=11)
    ax2.set_ylabel("$\\langle Q_\mathrm{L} \\rangle$", fontsize=11)
    ax2.text(0.05 * nu[1], -0.9 * 5 * nu[0], f"$\\nu={nu[0]}/{nu[1]}$", fontsize=11)

    ax3 = plt.subplot(gs[1])  # 071829 #################################################################################
    nu = (2, 19)

    phiflow_file = f'charge_pump_FerHofSqu1_chi_25_t1_1_V_10_Coulomb_1_n_2_361_nphi_4_19_LxMUC_1_Ly_19_phi_0_19_191.dat'
    phiflow_path = os.path.join(phiflow_dir, phiflow_file)
    with open(phiflow_path, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        phi = []
        charge = []
        for row in plots:
            phi.append(float(row[0]))
            charge.append(float(row[1]))
    charge = [i - charge[0] for i in charge]
    ax3.plot(phi, charge, '.', c=f'C0', marker=markers[7], fillstyle='none', markersize=5, markeredgewidth=0.2)

    ax3.xaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax3.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax3.set_xlim([0, nu[1]])
    ax3.set_xticks(np.arange(0, nu[1] + 0.1, 4.75))
    ax3.set_yticks([-5 * nu[0], 0])
    ax3.set_xlabel("$\\Phi_x / 2\pi$", fontsize=11)
    ax3.set_ylabel("$\\langle Q_\mathrm{L} \\rangle$", fontsize=11)
    ax3.text(0.05 * nu[1], -0.9 * 5 * nu[0], f"$\\nu={nu[0]}/{nu[1]}$", fontsize=11)

    # chi_legend_elements = [Patch(facecolor='C1', label='$50$'), Patch(facecolor='C2', label='$100$')]
    # leg2 = ax3.legend(handles=chi_legend_elements, loc='center', handletextpad=0.3, handlelength=1, labelspacing=0.1,
    #                   borderpad=0.3, framealpha=1, edgecolor='k', markerscale=2, fontsize=10, ncol=2, columnspacing=0.5,
    #                   bbox_to_anchor=(-0.65, 0.4), title='$\chi$', title_fontsize=11)
    # leg2.get_frame().set_linewidth(0.5)

    ####################################################################################################################
    ####################################################################################################################

    fig.text(0.035, 0.865, "(c)", fontsize=12)
    fig.text(0.48, 0.865, "(d)", fontsize=12)

    plt.savefig("/home/bart/Documents/papers/FCI/phiflow_c5_analysis.png", bbox_inches='tight', dpi=300)
    plt.show()
