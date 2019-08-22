import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import csv
from numpy.polynomial.polynomial import polyfit
from mpl_toolkits.mplot3d import Axes3D

if __name__ == '__main__':

    fig = plt.figure(figsize=(7.5, 3))  # height = 3.5
    gs = gridspec.GridSpec(2, 2, width_ratios=[1, 1], height_ratios=[1, 1])

    ####################################################################################################################

    ax1 = plt.subplot(gs[:, 0])

    phi_CI = []
    CI_charge = []
    phi_FCI = []
    FCI_charge = []

    with open('charge_pump_CI.dat', 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        for row in plots:
            phi_CI.append(float(row[0]))
            CI_charge.append(float(row[1]))

    with open('charge_pump_FCI.dat', 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        for row in plots:
            phi_FCI.append(float(row[0]))
            FCI_charge.append(float(row[1]))

    ax1.plot(phi_CI, CI_charge, marker='o', linestyle='--', c='C8', label='CI')
    ax1.plot(phi_FCI, FCI_charge, marker='s', linestyle='--', c='C9', label='FCI')
    ax1.legend(loc='lower right', handletextpad=0.2, borderpad=0.4, framealpha=1, edgecolor='k', markerscale=1,
               fontsize=10)

    ax1.set_xticks(np.arange(0, 2.1, 0.2))
    ax1.set_xlim([0, 2])
    ax1.set_xlabel("$\Phi_y / 2\pi$", fontsize=10)
    ax1.set_yticks(np.arange(0, 2.1, 0.2))
    ax1.set_ylim([0, 2])
    ax1.set_ylabel("$\langle Q^{L}(\Phi_y) \\rangle$", fontsize=10)

    ax1.axhline(1, color='k', linewidth=0.5, ls='--')
    ax1.axvline(1, color='k', linewidth=0.5, ls='--')

    ax1.tick_params(axis="x", labelsize=8)
    ax1.tick_params(axis="y", labelsize=8)
    ax1.tick_params(axis="z", labelsize=8)

    ############################################

    left, bottom, width, height = [0.115, 0.47, 0.2, 0.4]  # left = 0.12
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
    yy = end/2
    ax5.plot_surface(xx, yy, zz, alpha=0.25)

    # arrow
    ax5.quiver(0, end, 0, 0, -(end + end/4), 0, arrow_length_ratio=0.06, color='k', lw=0.5)
    ax5.text2D(-0.078, 0.006, '$\Phi_y$', fontsize=10)

    # cylinder
    angle = np.linspace(0, 2 * np.pi, 21)
    length = np.linspace(0, end-2/3, 13)
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
    ax5.text2D(0.05, 0.03, '$L$', fontsize=10)

    ax5.view_init(elev=10., azim=-15)
    ax5.axis('off')
    ax5.patch.set_alpha(0)

    ####################################################################################################################

    ax2 = plt.subplot(gs[0, 1])

    x = []
    y = []
    z = []

    with open('ent_spec_flow_CI.dat', 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        for row in plots:
            x.append(float(row[1]))
            y.append(float(row[2]))
            z.append(int(row[0]))

    for value in np.linspace(-5, 5, 11, dtype=int):
        xvalue = []
        yvalue = []
        for i in range(len(x)):
            if z[i] == value:
                xvalue.append(x[i])
                yvalue.append(y[i])
        ax2.scatter(xvalue, yvalue, marker='_', c='C{}'.format((value+5)%10), label='{}'.format(value))

    ax2.set_ylim([0, 10])
    # ax2.legend(loc='upper right', handletextpad=0, borderpad=0.2, framealpha=1, edgecolor='k', markerscale=2, fontsize=10)
    ax2.set_xlim([0, 2])
    ax2.set_xlabel("$\Phi_y / 2\pi$", fontsize=10)
    ax2.set_ylabel("$\epsilon_{\\alpha}$", fontsize=10)

    ax2.tick_params(axis="x", labelsize=8)
    ax2.tick_params(axis="y", labelsize=8)
    ax2.tick_params(axis="z", labelsize=8)

    ax2.axvline(1, color='k', linewidth=0.5, ls='--')
    # ax2.axhline(0.153205653777937, color='k', linewidth=0.5, ls=':')
    # ax2.axhline(2.073942791479203, color='k', linewidth=0.5, ls=':')
    # ax2.axhline(4.410048930142911, color='k', linewidth=0.5, ls=':')
    # ax2.axhline(5.700361870617707, color='k', linewidth=0.5, ls=':')
    # ax2.axhline(7.870452206203998, color='k', linewidth=0.5, ls=':')
    # ax2.axhline(8.400756721645187, color='k', linewidth=0.5, ls=':')
    # ax2.axhline(8.802885632428440, color='k', linewidth=0.5, ls=':')

    fig.text(0.57, 0.85, 'CI', fontsize=10, verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=1))

    gs.update(wspace=0.25, hspace=0.2)

    ####################################################################################################################

    ax3 = plt.subplot(gs[1, 1])

    x = []
    y = []
    z = []

    with open('ent_spec_flow_FCI.dat', 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        for row in plots:
            x.append(float(row[1]))
            y.append(float(row[2]))
            z.append(int(row[0]))

    for value in np.linspace(-5, 5, 11, dtype=int):
        xvalue = []
        yvalue = []
        for i in range(len(x)):
            if z[i] == value:
                xvalue.append(x[i])
                yvalue.append(y[i])
        ax3.scatter(xvalue, yvalue, marker='_', c='C{}'.format((value + 5)%10), label='{}'.format(value))

    ax3.set_yticks(np.arange(0, 12, 3))
    ax3.set_ylim([0, 12])
    ax3.legend(loc='center left', handletextpad=0, borderpad=0.2, framealpha=1, edgecolor='k', markerscale=2,
                fontsize=10, ncol=1, bbox_to_anchor=(1, 1))
    ax3.set_xlim([0, 2])
    ax3.set_xlabel("$\Phi_y / 2\pi$", fontsize=10)
    ax3.set_ylabel("$\epsilon_{\\alpha}$", fontsize=10)

    ax3.tick_params(axis="x", labelsize=8)
    ax3.tick_params(axis="y", labelsize=8)
    ax3.tick_params(axis="z", labelsize=8)

    ax3.axvline(1, color='k', linewidth=0.5, ls='--')

    fig.text(0.57, 0.46, 'FCI', fontsize=10, verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=1))

    gs.update(wspace=0.25, hspace=0)
    plt.setp(ax2.get_xticklabels(), visible=False)

    ####################################################################################################################

    fig.text(0.05, 0.875, "(a)", color="k", fontsize=12)
    fig.text(0.49, 0.875, "(b)", color="k", fontsize=12)
    fig.text(0.49, 0.49, "(c)", color="k", fontsize=12)

    plt.savefig("/home/bart/Documents/papers/TBG/figures/CI_1x2_detail_poster.png", bbox_inches='tight', dpi=300)
    plt.show()