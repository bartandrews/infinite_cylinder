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

    fig = plt.figure(figsize=(6, 2))  # figsize=(6, 1.1)
    gs = gridspec.GridSpec(1, 2, wspace=0.4)

    # define a list of easily-visible markers
    markers = [(3, 0, 0), (4, 0, 0), (5, 0, 0), (6, 0, 0), (4, 1, 0), (5, 1, 0), (6, 1, 0),
               (3, 2, 0), (4, 2, 0), (5, 2, 0), (6, 2, 0), 'X', 'x', 'd', 'D', 'P',
               '$a$', '$b$', '$c$', '$d$', '$e$', '$f$', '$g$', '$h$', '$i$', '$j$', '$k$', '$l$', '$m$', '$n$', '$o$',
               '$p$', '$q$', '$r$', '$s$', '$t$', '$u$', '$v$', '$w$', '$x$', '$y$', '$z$']

    ####################################################################################################################

    ax1 = plt.subplot(gs[0])  # 071829 #################################################################################

    with open('code/standalone/FCI/Vcritical/Vcritical_C1.dat', 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        p = []
        VcW_rm4, VcW_rm4_err = [], []
        VcW_rm3, VcW_rm3_err = [], []
        VcW_rm2, VcW_rm2_err = [], []
        VcW_r1, VcW_r1_err = [], []
        VcW_r2, VcW_r2_err = [], []
        VcW_r3, VcW_r3_err = [], []
        for i, row in enumerate(plots):
            if int(row[0]) == -4:
                p.append(float(row[1]))
                VcW_rm4.append(float(row[2])/float(row[4]))
                VcW_rm4_err.append(float(row[3]) / float(row[4]))
            elif int(row[0]) == -3:
                VcW_rm3.append(float(row[2])/float(row[4]))
                VcW_rm3_err.append(float(row[3]) / float(row[4]))
            elif int(row[0]) == -2:
                VcW_rm2.append(float(row[2])/float(row[4]))
                VcW_rm2_err.append(float(row[3]) / float(row[4]))
            elif int(row[0]) == 1:
                VcW_r1.append(float(row[2])/float(row[4]))
                VcW_r1_err.append(float(row[3]) / float(row[4]))
            elif int(row[0]) == 2:
                VcW_r2.append(float(row[2])/float(row[4]))
                VcW_r2_err.append(float(row[3]) / float(row[4]))
            elif int(row[0]) == 3:
                VcW_r3.append(float(row[2])/float(row[4]))
                VcW_r3_err.append(float(row[3]) / float(row[4]))

    # ax1.plot(p, VcW_rm4, '-', c='C0', marker=markers[0], fillstyle='none', markersize=5)
    ax1.errorbar(p, VcW_rm4, yerr=VcW_rm4_err, c='C0', marker=markers[0], fillstyle='none', markersize=5, label='$-4$', capsize=3)
    ax1.errorbar(p, VcW_rm3, yerr=VcW_rm3_err, c='C1', marker=markers[1], fillstyle='none', markersize=5, label='$-3$', capsize=3)
    ax1.errorbar(p, VcW_rm2, yerr=VcW_rm2_err, c='C2', marker=markers[2], fillstyle='none', markersize=5, label='$-2$', capsize=3)
    ax1.errorbar(np.NaN, np.NaN, yerr=np.NaN, c='C3', marker=markers[3], fillstyle='none', markersize=5, label='$-1$', capsize=3)
    ax1.errorbar(p, VcW_r1, yerr=VcW_r1_err, c='C4', marker=markers[4], fillstyle='none', markersize=5, label='$1$', capsize=3)
    ax1.errorbar(p, VcW_r2, yerr=VcW_r2_err, c='C5', marker=markers[5], fillstyle='none', markersize=5, label='$2$', capsize=3)
    ax1.errorbar(p, VcW_r3, yerr=VcW_r3_err, c='C6', marker=markers[6], fillstyle='none', markersize=5, label='$3$', capsize=3)

    left, bottom, width, height = [0.195, 0.55, 0.15, 0.3]
    ax1sub = fig.add_axes([left, bottom, width, height])
    with open('code/standalone/FCI/Vcritical/Vcritical_C2.dat', 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        for i, row in enumerate(plots):
            p = int(row[1])
            VcW = float(row[2])/float(row[4])
            VcW_err = float(row[3]) / float(row[4])
            ax1sub.errorbar(p, VcW, yerr=VcW_err, c=f'C{int(row[0])+4}', marker=markers[int(row[0])+4], fillstyle='none', markersize=5, capsize=3)
    ax1sub.set_xticks(np.arange(5, 8 + 0.1, 1))

    ax1.set_xlabel("$p$", fontsize=11)
    ax1.set_xticks(np.arange(4, 9 + 0.1, 1))
    ax1.set_ylabel("$V_\mathrm{crit} / W$", fontsize=11)

    leg = ax1.legend(loc='upper center', handletextpad=0.5, handlelength=0, labelspacing=0.2, borderpad=0.35,
                     framealpha=1,
                     edgecolor='k', markerscale=0.8, fontsize=10, ncol=8, columnspacing=1, bbox_to_anchor=(1.19, 1.49),
                     title='$r$')
    leg.get_frame().set_linewidth(0.5)

    ax2 = plt.subplot(gs[1])  # 071829 #################################################################################

    with open('code/standalone/FCI/Vcritical/Vcritical_C1.dat', 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        p = []
        VcD_rm4, VcD_rm4_err = [], []
        VcD_rm3, VcD_rm3_err = [], []
        VcD_rm2, VcD_rm2_err = [], []
        VcD_r1, VcD_r1_err = [], []
        VcD_r2, VcD_r2_err = [], []
        VcD_r3, VcD_r3_err = [], []
        for i, row in enumerate(plots):
            if int(row[0]) == -4:
                p.append(float(row[1]))
                VcD_rm4.append(float(row[2]) / float(row[5]))
                VcD_rm4_err.append(float(row[3]) / float(row[5]))
            elif int(row[0]) == -3:
                VcD_rm3.append(float(row[2]) / float(row[5]))
                VcD_rm3_err.append(float(row[3]) / float(row[5]))
            elif int(row[0]) == -2:
                VcD_rm2.append(float(row[2]) / float(row[5]))
                VcD_rm2_err.append(float(row[3]) / float(row[5]))
            elif int(row[0]) == 1:
                VcD_r1.append(float(row[2]) / float(row[5]))
                VcD_r1_err.append(float(row[3]) / float(row[5]))
            elif int(row[0]) == 2:
                VcD_r2.append(float(row[2]) / float(row[5]))
                VcD_r2_err.append(float(row[3]) / float(row[5]))
            elif int(row[0]) == 3:
                VcD_r3.append(float(row[2]) / float(row[5]))
                VcD_r3_err.append(float(row[3]) / float(row[5]))

    ax2.errorbar(p, VcD_rm4, yerr=VcD_rm4_err, c='C0', marker=markers[0], fillstyle='none', markersize=5, label='$-4$', capsize=3)
    ax2.errorbar(p, VcD_rm3, yerr=VcD_rm3_err, c='C1', marker=markers[1], fillstyle='none', markersize=5, label='$-3$', capsize=3)
    ax2.errorbar(p, VcD_rm2, yerr=VcD_rm2_err, c='C2', marker=markers[2], fillstyle='none', markersize=5, label='$-2$', capsize=3)
    # ax2.plot(4, 20, '.-', c='C3', marker=markers[3], fillstyle='none', markersize=5, label='$-1$')
    ax2.errorbar(p, VcD_r1, yerr=VcD_r1_err, c='C4', marker=markers[4], fillstyle='none', markersize=5, label='$1$', capsize=3)
    ax2.errorbar(p, VcD_r2, yerr=VcD_r2_err, c='C5', marker=markers[5], fillstyle='none', markersize=5, label='$2$', capsize=3)
    ax2.errorbar(p, VcD_r3, yerr=VcD_r3_err, c='C6', marker=markers[6], fillstyle='none', markersize=5, label='$3$', capsize=3)

    left, bottom, width, height = [0.72, 0.55, 0.15, 0.3]
    ax2sub = fig.add_axes([left, bottom, width, height])
    with open('code/standalone/FCI/Vcritical/Vcritical_C2.dat', 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        for i, row in enumerate(plots):
            p = int(row[1])
            VcD = float(row[2])/float(row[5])
            VcD_err = float(row[3])/float(row[5])
            ax2sub.errorbar(p, VcD, yerr=VcD_err, c=f'C{int(row[0])+4}', marker=markers[int(row[0])+4], fillstyle='none', markersize=5, capsize=3)
    ax2sub.set_xticks(np.arange(5, 8 + 0.1, 1))

    ax2.set_xlabel("$p$", fontsize=11)
    ax2.set_xticks(np.arange(4, 9 + 0.1, 1))
    ax2.set_ylabel("$V_\mathrm{crit} / \Delta$", fontsize=11)

    ####################################################################################################################
    ####################################################################################################################

    fig.text(0.03, 0.87, "(a)", fontsize=12)
    fig.text(0.5, 0.87, "(b)", fontsize=12)

    plt.savefig("/home/bart/Documents/papers/FCI/Vcritical.png", bbox_inches='tight', dpi=300)
    plt.show()
