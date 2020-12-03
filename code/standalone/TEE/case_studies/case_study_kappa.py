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

plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}')
# matplotlib.verbose.level = 'debug-annoying'

if __name__ == '__main__':

    fig = plt.figure()
    gs = gridspec.GridSpec(2, 2, hspace=0.4, wspace=0.4)

    # nu=2/5 charge_pump ###############################################################################################

    ax1 = plt.subplot(gs[0], anchor=(0, 0.85))

    charge_pump_dir = '/home/bart/PycharmProjects/infinite_cylinder/data/charge_pump/FerHofSqu1'
    charge_pump_file = 'charge_pump_FerHofSqu1_chi_400_t1_1_V_10_Coulomb_1_n_2_35_nphi_1_7_LxMUC_1_Ly_10_phi_0_5_51.dat'
    charge_pump_path = os.path.join(charge_pump_dir, charge_pump_file)

    # extract data from file
    with open(charge_pump_path, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        phi = []
        charge = []
        for row in plots:
            phi.append(float(row[0]))
            charge.append(float(row[1]))

    charge = [i - charge[0] for i in charge]

    ax1.plot(phi, charge, 's', marker='x', color='k', markersize=3, label='$\kappa=1$')

    charge_pump_file = 'charge_pump_FerHofSqu1_chi_500_t1_1_V_10.555_Coulomb_2_n_2_35_nphi_1_7_LxMUC_1_Ly_10_phi_0_5_51.dat.reflected'
    charge_pump_path = os.path.join(charge_pump_dir, charge_pump_file)

    # extract data from file
    with open(charge_pump_path, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        phi = []
        charge = []
        for row in plots:
            phi.append(float(row[0]))
            charge.append(float(row[1]))

    charge = [i - charge[0] for i in charge]

    ax1.plot(phi, charge, 's', marker='x', color='g', markersize=3, label='$\kappa=2$')

    charge_pump_file = 'charge_pump_FerHofSqu1_chi_500_t1_1_V_10.555_Coulomb_3_n_2_35_nphi_1_7_LxMUC_1_Ly_10_phi_0_5_51.dat.reflected'
    charge_pump_path = os.path.join(charge_pump_dir, charge_pump_file)

    # extract data from file
    with open(charge_pump_path, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        phi = []
        charge = []
        for row in plots:
            phi.append(float(row[0]))
            charge.append(float(row[1]))

    charge = [i - charge[0] for i in charge]

    ax1.plot(phi, charge, 's', marker='x', color='r', markersize=3, label='$\kappa=3$')

    ax1.axhline(charge[0], color='k', linewidth=0.5, ls='--')
    ax1.axhline(charge[-1], color='k', linewidth=0.5, ls='--')
    ax1.axvline(1, color='k', linewidth=0.5, ls='--')
    ax1.axvline(2, color='k', linewidth=0.5, ls='--')
    ax1.axvline(3, color='k', linewidth=0.5, ls='--')
    ax1.axvline(4, color='k', linewidth=0.5, ls='--')
    ax1.axhline(1, color='k', linewidth=0.5, ls='--')
    ax1.text(0.45, 1.4, "$\displaystyle{\sigma_\\text{H}=\\frac{2}{5}\\frac{e^2}{h}}$",
             bbox=dict(facecolor='white', alpha=1, lw=0.5, ls='-'))

    ax1.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax1.tick_params(axis="x", labelsize=10)
    ax1.set_xlim([0, 5])
    ax1.set_xticks(np.arange(0, 5.1, 1))
    ax1.tick_params(axis="y", labelsize=10)
    ax1.set_ylim([0, 2])
    ax1.set_xlabel("$\\Phi_x / 2\pi$", fontsize=11)
    ax1.set_ylabel("$\langle Q_\mathrm{L} \\rangle$", fontsize=11)

    ax1.legend(loc='upper center', bbox_to_anchor=(1.175, 1.4), handletextpad=0, borderpad=0.4,
              framealpha=1, edgecolor='k', markerscale=2, fontsize=10, ncol=8, labelspacing=0.3, columnspacing=0)

    # nu=2/5 ent_spec_mom ##############################################################################################

    ax2 = plt.subplot(gs[1])

    ent_spec_mom_dir = '/home/bart/PycharmProjects/infinite_cylinder/data/ent_spec_mom/FerHofSqu1'
    ent_spec_mom_file = 'ent_spec_mom_FerHofSqu1_chi_800_chiK_800_t1_1_V_10_Coulomb_1_n_2_35_nphi_1_7_LxMUC_1_Ly_10.dat.rotated'
    ent_spec_mom_path = os.path.join(ent_spec_mom_dir, ent_spec_mom_file)

    x = []
    y = []
    z = []

    with open(ent_spec_mom_path, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        for row in plots:
            x.append(float(row[1]))
            y.append(float(row[2]))
            z.append(int(row[0]))

    for value in np.linspace(3, -3, 7, dtype=int):
        xvalue = []
        yvalue = []
        for i in range(len(x)):
            if z[i] == value:
                xvalue.append(x[i])
                yvalue.append(y[i])
        if value != 0:
            ax2.scatter(xvalue, yvalue, marker='_', c='k', label='{}'.format(value))
        else:
            ax2.scatter(xvalue, yvalue, marker='_', c='k', label='{}'.format(value))

    ent_spec_mom_file = 'ent_spec_mom_FerHofSqu1_chi_2000_chiK_2000_t1_1_V_10.555_Coulomb_2_n_2_35_nphi_1_7_LxMUC_1_Ly_10.dat.rotated'
    ent_spec_mom_path = os.path.join(ent_spec_mom_dir, ent_spec_mom_file)

    x = []
    y = []
    z = []

    with open(ent_spec_mom_path, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        for row in plots:
            x.append(float(row[1])+0.02879)
            y.append(float(row[2]))
            z.append(int(row[0]))

    for value in np.linspace(3, -3, 7, dtype=int):
        xvalue = []
        yvalue = []
        for i in range(len(x)):
            if z[i] == value:
                xvalue.append(x[i])
                yvalue.append(y[i])
        if value != 0:
            ax2.scatter(xvalue, yvalue, marker='_', c='g', label='{}'.format(value))
        else:
            ax2.scatter(xvalue, yvalue, marker='_', c='g', label='{}'.format(value))

    ent_spec_mom_file = 'ent_spec_mom_FerHofSqu1_chi_2000_chiK_2000_t1_1_V_10.555_Coulomb_3_n_2_35_nphi_1_7_LxMUC_1_Ly_10.dat.rotated'
    ent_spec_mom_path = os.path.join(ent_spec_mom_dir, ent_spec_mom_file)

    x = []
    y = []
    z = []

    with open(ent_spec_mom_path, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        for row in plots:
            x.append(float(row[1])-0.04702)
            y.append(float(row[2]))
            z.append(int(row[0]))

    for value in np.linspace(3, -3, 7, dtype=int):
        xvalue = []
        yvalue = []
        for i in range(len(x)):
            if z[i] == value:
                xvalue.append(x[i])
                yvalue.append(y[i])
        if value != 0:
            ax2.scatter(xvalue, yvalue, marker='_', c='r', label='{}'.format(value))
        else:
            ax2.scatter(xvalue, yvalue, marker='_', c='r', label='{}'.format(value))

    ax2.set_yticks(np.arange(0, 15.1, 5))
    ax2.set_ylim([0, 15])

    ax2.set_xlim([-np.pi / 3, np.pi / 3])
    ax2.set_xlabel("$\\tilde{k}_a / \pi$", fontsize=11)
    ax2.set_ylabel("$\epsilon_{\\alpha}$", fontsize=11)

    ax2.tick_params(axis="x", labelsize=10)
    ax2.tick_params(axis="y", labelsize=10)

    # ax2.plot([-0.25, 1 * np.pi / 3], [0, 11], color='k', linestyle=':', linewidth=1)
    # ax2.plot([0.65, 1 * np.pi / 3], [0, 3.5], color='k', linestyle=':', linewidth=1)

    # nu=3/7 charge_pump ###############################################################################################

    ax3 = plt.subplot(gs[2], anchor=(0, 0.85))

    charge_pump_dir = '/home/bart/PycharmProjects/infinite_cylinder/data/charge_pump/FerHofSqu1'
    charge_pump_file = 'charge_pump_FerHofSqu1_chi_500_t1_1_V_10_Coulomb_1_n_3_70_nphi_1_10_LxMUC_1_Ly_14_phi_0_7_71.dat'
    charge_pump_path = os.path.join(charge_pump_dir, charge_pump_file)

    # extract data from file
    with open(charge_pump_path, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        phi = []
        charge = []
        for row in plots:
            phi.append(float(row[0]))
            charge.append(float(row[1]))

    charge = [i - charge[0] for i in charge]

    ax3.plot(phi, charge, 's', marker='x', color='k', markersize=3)

    charge_pump_file = 'charge_pump_FerHofSqu1_chi_500_t1_1_V_12.6157_Coulomb_3_n_3_70_nphi_1_10_LxMUC_1_Ly_14_phi_0_7_71.dat.reflected'
    charge_pump_path = os.path.join(charge_pump_dir, charge_pump_file)

    # extract data from file
    with open(charge_pump_path, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        phi = []
        charge = []
        for row in plots:
            phi.append(float(row[0]))
            charge.append(float(row[1]))

    charge = [i - charge[0] for i in charge]

    ax3.plot(phi, charge, 's', marker='x', color='g', markersize=3)

    # charge_pump_file = 'charge_pump_FerHofSqu1_chi_500_t1_1_V_12.6157_Coulomb_3_n_3_70_nphi_1_10_LxMUC_1_Ly_14_phi_0_7_71.dat.reflected'
    # charge_pump_path = os.path.join(charge_pump_dir, charge_pump_file)
    #
    # # extract data from file
    # with open(charge_pump_path, 'r') as csvfile:
    #     plots = csv.reader(csvfile, delimiter='\t')
    #     phi = []
    #     charge = []
    #     for row in plots:
    #         phi.append(float(row[0]))
    #         charge.append(float(row[1]))
    #
    # charge = [i - charge[0] for i in charge]
    #
    # ax3.plot(phi, charge, 's', marker='x', color='r', markersize=3)

    # ax1.axhline(charge[0], color='k', linewidth=0.5, ls='--')
    # ax1.axhline(charge[-1], color='k', linewidth=0.5, ls='--')
    ax3.axvline(1, color='k', linewidth=0.5, ls='--')
    ax3.axvline(2, color='k', linewidth=0.5, ls='--')
    ax3.axvline(3, color='k', linewidth=0.5, ls='--')
    ax3.axvline(4, color='k', linewidth=0.5, ls='--')
    ax3.axvline(5, color='k', linewidth=0.5, ls='--')
    ax3.axvline(6, color='k', linewidth=0.5, ls='--')
    ax3.axhline(1, color='k', linewidth=0.5, ls='--')
    ax3.axhline(2, color='k', linewidth=0.5, ls='--')
    ax3.text(0.6, 2.1, "$\displaystyle{\sigma_\\text{H}=\\frac{3}{7}\\frac{e^2}{h}}$",
             bbox=dict(facecolor='white', alpha=1, lw=0.5, ls='-'))

    ax3.tick_params(axis="x", labelsize=10)
    ax3.set_xlim([0, 7])
    ax3.set_xticks(np.arange(0, 7.1, 1))
    ax3.tick_params(axis="y", labelsize=10)
    ax3.set_ylim([0, 3])
    ax3.set_xlabel("$\\Phi_x / 2\pi$", fontsize=11)
    ax3.set_ylabel("$\langle Q_\mathrm{L} \\rangle$", fontsize=11)

    # nu=3/7 ent_spec_mom ##############################################################################################

    ax4 = plt.subplot(gs[3])

    ent_spec_mom_dir = '/home/bart/PycharmProjects/infinite_cylinder/data/ent_spec_mom/FerHofSqu1'
    ent_spec_mom_file = 'ent_spec_mom_FerHofSqu1_chi_2000_chiK_2000_t1_1_V_10_Coulomb_1_n_3_70_nphi_1_10_LxMUC_1_Ly_14.dat'
    ent_spec_mom_path = os.path.join(ent_spec_mom_dir, ent_spec_mom_file)

    x = []
    y = []
    z = []

    with open(ent_spec_mom_path, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        for row in plots:
            x.append(float(row[1]))
            y.append(float(row[2]))
            z.append(int(row[0]))

    for value in np.linspace(3, -3, 7, dtype=int):
        xvalue = []
        yvalue = []
        for i in range(len(x)):
            if z[i] == value:
                xvalue.append(x[i])
                yvalue.append(y[i])
        if value != 0:
            ax4.scatter(xvalue, yvalue, marker='_', c='k', label='{}'.format(value))
        else:
            ax4.scatter(xvalue, yvalue, marker='_', c='k', label='{}'.format(value))

    ent_spec_mom_file = 'ent_spec_mom_FerHofSqu1_chi_2000_chiK_2000_t1_1_V_12.6157_Coulomb_2_n_3_70_nphi_1_10_LxMUC_1_Ly_14.dat.rotated'
    ent_spec_mom_path = os.path.join(ent_spec_mom_dir, ent_spec_mom_file)

    x = []
    y = []
    z = []

    with open(ent_spec_mom_path, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        for row in plots:
            x.append(float(row[1])+0.08085)
            y.append(float(row[2]))
            z.append(int(row[0]))

    for value in np.linspace(3, -3, 7, dtype=int):
        xvalue = []
        yvalue = []
        for i in range(len(x)):
            if z[i] == value:
                xvalue.append(x[i])
                yvalue.append(y[i])
        if value != 0:
            ax4.scatter(xvalue, yvalue, marker='_', c='g', label='{}'.format(value))
        else:
            ax4.scatter(xvalue, yvalue, marker='_', c='g', label='{}'.format(value))

    # ent_spec_mom_file = 'ent_spec_mom_FerHofSqu1_chi_2000_chiK_2000_t1_1_V_12.6157_Coulomb_3_n_3_70_nphi_1_10_LxMUC_1_Ly_14.dat.rotated'
    # ent_spec_mom_path = os.path.join(ent_spec_mom_dir, ent_spec_mom_file)
    #
    # x = []
    # y = []
    # z = []
    #
    # with open(ent_spec_mom_path, 'r') as csvfile:
    #     plots = csv.reader(csvfile, delimiter='\t')
    #     for row in plots:
    #         x.append(float(row[1])-0.02216)
    #         y.append(float(row[2]))
    #         z.append(int(row[0]))
    #
    # for value in np.linspace(3, -3, 7, dtype=int):
    #     xvalue = []
    #     yvalue = []
    #     for i in range(len(x)):
    #         if z[i] == value:
    #             xvalue.append(x[i])
    #             yvalue.append(y[i])
    #     if value != 0:
    #         ax4.scatter(xvalue, yvalue, marker='_', c='r', label='{}'.format(value))
    #     else:
    #         ax4.scatter(xvalue, yvalue, marker='_', c='r', label='{}'.format(value))

    ax4.set_yticks(np.arange(0, 15.1, 5))
    ax4.set_ylim([0, 15])

    ax4.set_xlim([-np.pi / 3, np.pi / 3])
    ax4.set_xlabel("$\\tilde{k}_a / \pi$", fontsize=11)
    ax4.set_ylabel("$\epsilon_{\\alpha}$", fontsize=11)

    ax4.tick_params(axis="x", labelsize=10)
    ax4.tick_params(axis="y", labelsize=10)

    # ax4.plot([-0.25, 1 * np.pi / 3], [0, 11], color='k', linestyle=':', linewidth=1)
    # ax4.plot([0.65, 1 * np.pi / 3], [0, 3.5], color='k', linestyle=':', linewidth=1)

    # complete plot ####################################################################################################

    fig.text(0.04, 0.88, "(a)", fontsize=12)
    fig.text(0.5, 0.88, "(b)", fontsize=12)
    fig.text(0.04, 0.43, "(c)", fontsize=12)
    fig.text(0.5, 0.43, "(d)", fontsize=12)

    plt.savefig("/home/bart/Documents/papers/TEE/case_study_kappa.png", bbox_inches='tight', dpi=300)
    plt.show()
