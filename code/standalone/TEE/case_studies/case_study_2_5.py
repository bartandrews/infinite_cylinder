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

if __name__ == '__main__':

    fig = plt.figure()
    gs = gridspec.GridSpec(2, 2, hspace=0.5, wspace=0.4)

    ####################################################################################################################

    ax1 = plt.subplot(gs[0], anchor=(0, 0.85))

    charge_pump_dir = '/home/bart/PycharmProjects/infinite_cylinder/data/charge_pump/BosHofSqu1'
    charge_pump_file = 'charge_pump_BosHofSqu1_chi_50_t1_1_n_1_8_nphi_1_4_LxMUC_1_Ly_4_phi_0_2_21.dat'
    charge_pump_path = os.path.join(charge_pump_dir, charge_pump_file)

    # extract data from file
    with open(charge_pump_path, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        phi = []
        charge = []
        for row in plots:
            phi.append(float(row[0]))
            charge.append(float(row[1]))

    ax1.plot(phi, charge, 's', marker='x', color='k', markersize=3)

    ax1.axhline(charge[0], color='k', linewidth=0.5, ls='--')
    ax1.axhline(charge[-1], color='k', linewidth=0.5, ls='--')
    ax1.axvline(1, color='k', linewidth=0.5, ls='--')
    ax1.axvline(2, color='k', linewidth=0.5, ls='--')

    ax1.set_xlim(0)
    ax1.tick_params(axis="x", labelsize=10)
    ax1.tick_params(axis="y", labelsize=10)
    ax1.set_xlabel("$\\Phi_x / 2\pi$", fontsize=11)
    ax1.set_ylabel("$\langle Q_\mathrm{L} \\rangle$", fontsize=11)

    ############################################

    left, bottom, width, height = [0.11, 0.215 + 0.12, 0.25, 0.25]  # left = 0.12
    ax5 = fig.add_axes([left, bottom, width, height], projection='3d')

    r = 5

    end = 6
    point_size = 2

    for y in range(0, end, 2):
        theta = np.linspace(0, 2 * np.pi, 11)
        z = r * np.cos(theta)
        x = r * np.sin(theta)
        ax5.scatter(x, y, z, label='parametric curve', c='b', s=point_size)

        theta = np.linspace((1 / 20) * 2 * np.pi, (1 + 1 / 20) * 2 * np.pi, 11)
        y = y + 1 / 3
        z = r * np.cos(theta)
        x = r * np.sin(theta)
        ax5.scatter(x, y, z, label='parametric curve', c='r', s=point_size)

        theta = np.linspace((1 / 20) * 2 * np.pi, (1 + 1 / 20) * 2 * np.pi, 11)
        y = y + 2 / 3
        z = r * np.cos(theta)
        x = r * np.sin(theta)
        ax5.scatter(x, y, z, label='parametric curve', c='b', s=point_size)

        theta = np.linspace(0, 2 * np.pi, 11)
        y = y + 1 / 3
        z = r * np.cos(theta)
        x = r * np.sin(theta)
        ax5.scatter(x, y, z, label='parametric curve', c='r', s=point_size)

    # cut
    xx, zz = np.meshgrid(np.linspace(-7, 7, 10), np.linspace(-7, 7, 10))
    yy = end / 2
    ax5.plot_surface(xx, yy, zz, alpha=0.25)

    # arrow
    ax5.quiver(0, end, 0, 0, -(end + end / 4), 0, arrow_length_ratio=0.06, color='k', lw=0.5)
    ax5.text2D(-0.078, 0.006, '$\Phi_x$', fontsize=11)

    # cylinder
    angle = np.linspace(0, 2 * np.pi, 21)
    length = np.linspace(0, end - 2 / 3, 13)
    theta, y = np.meshgrid(angle, length)
    z = r * np.cos(theta)
    x = r * np.sin(theta)
    ax5.plot_surface(x, y, z, alpha=0.25)

    # ring
    theta = np.linspace(0, 2 * np.pi, 21)
    y = np.repeat(end, 21)
    z = r * np.cos(theta)
    x = r * np.sin(theta)
    ax5.plot(x, y, z, label='parametric curve', c='k', lw=0.5)

    theta = np.linspace(0, 2 * np.pi, 1)
    y = np.repeat(end, 1)
    z = r * np.cos(theta)
    x = r * np.sin(theta)
    ax5.scatter(x, y, z, label='parametric curve', c='k', marker='v', s=5)
    ax5.text2D(0.05, 0.03, '$L_y$', fontsize=11)

    ax5.view_init(elev=10., azim=-15)
    ax5.axis('off')
    ax5.patch.set_alpha(0)

    ####################################################################################################################

    ax2 = plt.subplot(gs[1])

    ent_spec_mom_dir = '/home/bart/PycharmProjects/infinite_cylinder/data/ent_spec_mom/FerHofSqu1'
    ent_spec_mom_file = 'ent_spec_mom_FerHofSqu1_chi_500_chiK_1000_t1_1_V_10_Coulomb_1_n_1_9_nphi_1_3_LxMUC_1_Ly_9.dat'
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
            ax2.scatter(xvalue, yvalue, marker='_', c='C{}'.format((value + 4) % 10), label='{}'.format(value))
        else:
            ax2.scatter(xvalue, yvalue, marker='x', c='C{}'.format((value + 4) % 10), label='{}'.format(value))

    ax2.set_yticks(np.arange(0, 15.1, 5))
    ax2.set_ylim([0, 15])

    # Shrink current axis by 20%
    box = ax2.get_position()
    ax2.set_position([box.x0, box.y0, box.width * 0.8, box.height])

    ax2.legend(loc='center left', handletextpad=0, handlelength=1, borderpad=0.2, framealpha=1, edgecolor='k', markerscale=1,
               fontsize=10, ncol=1, bbox_to_anchor=(1, 0.5))
    ax2.set_xlim([-np.pi / 3, np.pi / 3])
    ax2.set_xlabel("$k_\\alpha / \pi$", fontsize=11)
    ax2.set_ylabel("$\epsilon_{\\alpha}$", fontsize=11)

    ax2.tick_params(axis="x", labelsize=10)
    ax2.tick_params(axis="y", labelsize=10)

    # ax2.plot([-0.25, 1 * np.pi / 3], [0, 11], color='k', linestyle=':', linewidth=1)
    # ax2.plot([0.65, 1 * np.pi / 3], [0, 3.5], color='k', linestyle=':', linewidth=1)

    ####################################################################################################################

    ax3 = plt.subplot(gs[2])

    density_dir = '/home/bart/PycharmProjects/infinite_cylinder/data/density/FerHofSqu1'
    density_file = 'density_FerHofSqu1_chi_500_t1_1_V_10_Coulomb_1_n_1_9_nphi_1_3_LxMUC_1_Ly_9.dat'
    density_path = os.path.join(density_dir, density_file)

    # extract data from file
    with open(density_path, 'r') as file:
        mat = [[float(num) for num in line.split('\t')] for line in file]

    matrix = np.array(mat).transpose()

    divider = make_axes_locatable(ax3)
    cax = divider.append_axes('right', size='10%', pad=0.10)
    ax3.set_xlabel("$L_x$", fontsize=11)
    ax3.set_ylabel("$L_y$", fontsize=11)
    im = ax3.imshow(matrix, origin='lower')

    # def fmt(x, pos):
    #     a, b = '{:.2e}'.format(x).split('e')
    #     b = int(b)
    #     return r'${} \times 10^{{{}}}$'.format(a, b)
    #
    # fig.colorbar(im, cax=cax, orientation='vertical', label='$\langle \\rho_i \\rangle$', format=ticker.FuncFormatter(fmt))

    fig.colorbar(im, cax=cax, orientation='vertical', label='$\langle \\rho_i \\rangle$')

    ax3.set_position([0.1, 0.1, 0.35, 0.35])

    ####################################################################################################################

    ax4 = plt.subplot(gs[3])

    corr_func_dir = '/home/bart/PycharmProjects/infinite_cylinder/data/corr_func/FerHofSqu1'
    corr_func_file = 'corr_func_FerHofSqu1_chi_500_t1_1_V_10_Coulomb_1_n_1_9_nphi_1_3_LxMUC_1_Ly_9.dat'
    corr_func_path = os.path.join(corr_func_dir, corr_func_file)

    # extract data from file
    with open(corr_func_path, 'r') as file:
        mat = [[float(num) for num in line.split('\t')] for line in file]

    matrix = np.array(mat).transpose()

    divider = make_axes_locatable(ax4)
    cax = divider.append_axes('right', size='10%', pad=0.10)
    ax4.set_xlabel("$L_x$", fontsize=11)
    ax4.set_ylabel("$L_y$", fontsize=11)
    im = ax4.imshow(matrix, origin='lower')
    fig.colorbar(im, cax=cax, orientation='vertical', label='$\langle :\\rho_0 \\rho_i : \\rangle$')

    ax4.set_position([0.54, 0.1, 0.35, 0.35])

    ####################################################################################################################

    fig.text(0.05, 0.9, "(a)", fontsize=12)
    fig.text(0.49, 0.9, "(b)", fontsize=12)
    fig.text(0.05, 0.45, "(c)", fontsize=12)
    fig.text(0.49, 0.45, "(d)", fontsize=12)

    plt.savefig("/home/bart/Documents/papers/TEE/figures/case_study_2_5.png", bbox_inches='tight', dpi=300)
    plt.show()
