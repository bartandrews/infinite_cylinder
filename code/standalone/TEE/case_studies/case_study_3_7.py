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
    gs = gridspec.GridSpec(2, 2, hspace=0.5, wspace=0.4)

    ####################################################################################################################

    ax1 = plt.subplot(gs[0], anchor=(0, 0.85))

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

    ax1.plot(phi, charge, 's', marker='x', color='k', markersize=3)

    #ax1.axhline(charge[0], color='k', linewidth=0.5, ls='--')
    #ax1.axhline(charge[-1], color='k', linewidth=0.5, ls='--')
    ax1.axvline(1, color='k', linewidth=0.5, ls='--')
    ax1.axvline(2, color='k', linewidth=0.5, ls='--')
    ax1.axvline(3, color='k', linewidth=0.5, ls='--')
    ax1.axvline(4, color='k', linewidth=0.5, ls='--')
    ax1.axvline(5, color='k', linewidth=0.5, ls='--')
    ax1.axvline(6, color='k', linewidth=0.5, ls='--')
    ax1.axhline(1, color='k', linewidth=0.5, ls='--')
    ax1.axhline(2, color='k', linewidth=0.5, ls='--')
    ax1.text(4.6, 0.33, "$\displaystyle{\sigma_\\text{H}=\\frac{3}{7}\\frac{e^2}{h}}$", bbox=dict(facecolor='white', alpha=1, lw=0.5, ls='-'))

    ax1.tick_params(axis="x", labelsize=10)
    ax1.set_xlim([0, 7])
    ax1.set_xticks(np.arange(0, 7.1, 1))
    ax1.tick_params(axis="y", labelsize=10)
    ax1.set_ylim([0, 3])
    ax1.set_xlabel("$\\Phi_x / 2\pi$", fontsize=11)
    ax1.set_ylabel("$\langle Q_\mathrm{L} \\rangle$", fontsize=11)

    ############################################

    # left, bottom, width, height = [0.215, 0.535, 0.25, 0.25]
    left, bottom, width, height = [0.11, 0.665, 0.25, 0.25]  # left = 0.12
    ax5 = fig.add_axes([left, bottom, width, height], projection='3d')

    r = 5

    end = 15
    point_size = 1.5

    for y in range(0, end, 2):
        theta = np.linspace(0, 2 * np.pi, 16)
        z = r * np.cos(theta)
        x = r * np.sin(theta)
        ax5.scatter(x, y, z, label='parametric curve', c='b', s=point_size)

    # cut
    xx, zz = np.meshgrid(np.linspace(-7, 7, 10), np.linspace(-7, 7, 10))
    yy = (end-1) / 2
    ax5.plot_surface(xx, yy, zz, alpha=0.25)

    # arrow
    ax5.quiver(0, end, 0, 0, -(end + end / 4), 0, arrow_length_ratio=0.06, color='k', lw=0.5)
    ax5.text2D(-0.078, 0.006, '$\Phi_x$', fontsize=11)

    # cylinder
    angle = np.linspace(0, 2 * np.pi, 21)
    length = np.linspace(0, end-1, 13)
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
    ax2.set_xlabel("$k_a / \pi$", fontsize=11)
    ax2.set_ylabel("$\epsilon_{\\alpha}$", fontsize=11)

    ax2.tick_params(axis="x", labelsize=10)
    ax2.tick_params(axis="y", labelsize=10)

    # ax2.plot([-0.25, 1 * np.pi / 3], [0, 11], color='k', linestyle=':', linewidth=1)
    # ax2.plot([0.65, 1 * np.pi / 3], [0, 3.5], color='k', linestyle=':', linewidth=1)

    ####################################################################################################################

    ax3 = plt.subplot(gs[2])

    density_dir = '/home/bart/PycharmProjects/infinite_cylinder/data/density/FerHofSqu1'
    density_file = 'density_FerHofSqu1_chi_2000_t1_1_V_10_Coulomb_1_n_3_70_nphi_1_10_LxMUC_1_Ly_14.dat'
    density_path = os.path.join(density_dir, density_file)

    # for the colorbar scale
    av_rho = 3/70

    # for tick labelling (in units of lattice constant)
    Lx, Ly = 10, 14

    # extract data from file
    with open(density_path, 'r') as file:
        mat = [[float(num)-av_rho for num in line.split('\t')] for line in file]

    matrix = np.array(mat).transpose()

    divider = make_axes_locatable(ax3)
    cax = divider.append_axes('right', size='10%', pad=0.10)
    ax3.set_xlabel("$L_x / a$", fontsize=11)
    ax3.set_xticks(np.arange(Lx, step=1))
    ax3.set_xticklabels([i + 1 for i in range(Lx)])
    ax3.set_ylabel("$L_y / a$", fontsize=11)
    ax3.set_yticks(np.arange(Ly, step=1))
    ax3.set_yticklabels([i + 1 for i in range(Ly)])
    im = ax3.imshow(matrix, origin='lower')

    def fmt(x, pos):
        a, b = '{:.0e}'.format(x).split('e')
        b = int(b)
        if int(a) == 0:
            return r'$0$'
        else:
            return r'${}$'.format(a)

    fig.colorbar(im, cax=cax, orientation='vertical', label='$(\langle \\rho_i \\rangle  - n)/10^{-3}$', format=ticker.FuncFormatter(fmt))

    # fig.colorbar(im, cax=cax, orientation='vertical', label='$\langle \\rho_i \\rangle - \\bar{\\rho}$')

    ax3.set_position([0.1, 0.1, 0.35, 0.35])

    ####################################################################################################################

    ax4 = plt.subplot(gs[3])

    corr_func_dir = '/home/bart/PycharmProjects/infinite_cylinder/data/corr_func/FerHofSqu1'
    corr_func_file = 'corr_func_FerHofSqu1_chi_2000_t1_1_V_10_Coulomb_1_n_3_70_nphi_1_10_LxMUC_1_Ly_14.dat'
    corr_func_path = os.path.join(corr_func_dir, corr_func_file)

    # for tick labelling (in units of lattice constant)
    Lx, Ly = 10, 14

    # extract data from file
    with open(corr_func_path, 'r') as file:
        mat = [[float(num) for num in line.split('\t')] for line in file]

    matrix = np.array(mat).transpose()

    divider = make_axes_locatable(ax4)
    cax = divider.append_axes('right', size='10%', pad=0.10)
    ax4.set_xlabel("$L_x / a$", fontsize=11)
    ax4.set_xticks(np.arange(Lx, step=1))
    ax4.set_xticklabels([i+1 for i in range(Lx)])
    ax4.set_ylabel("$L_y / a$", fontsize=11)
    ax4.set_yticks(np.arange(Ly, step=1))
    ax4.set_yticklabels([i + 1 for i in range(Ly)])
    im = ax4.imshow(matrix, origin='lower', vmin=0)
    fig.colorbar(im, cax=cax, orientation='vertical', label='$\langle :\mathrel{\\rho_0 \\rho_i}: \\rangle$')

    ax4.set_position([0.54, 0.1, 0.35, 0.35])

    ####################################################################################################################

    fig.text(0.04, 0.9, "(e)", fontsize=12)
    fig.text(0.5, 0.9, "(f)", fontsize=12)
    fig.text(0.04, 0.45, "(g)", fontsize=12)
    fig.text(0.5, 0.45, "(h)", fontsize=12)

    plt.savefig("/home/bart/Documents/papers/TEE/figures/case_study_3_7.png", bbox_inches='tight', dpi=300)
    plt.show()
