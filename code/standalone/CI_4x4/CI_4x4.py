import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import csv
from numpy.polynomial.polynomial import polyfit
from mpl_toolkits.mplot3d import Axes3D

if __name__ == '__main__':

    fig = plt.figure(figsize=(8.5, 8))
    gs = gridspec.GridSpec(2, 2, width_ratios=[1, 1], height_ratios=[1, 1])

    ####################################################################################################################

    ax1 = plt.subplot(gs[0])

    x = []
    y = []

    with open('charge_pump_example.dat', 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        for row in plots:
            x.append(float(row[0]))
            y.append(float(row[1]))

    ax1.scatter(x, y, s=3, c='k')

    ax1.set_xlim([0, 1])
    ax1.set_xlabel("$\Phi / 2\pi$", fontsize=14)
    ax1.set_ylim([-1, 0])
    ax1.set_ylabel("$\langle Q^{L} \\rangle$", fontsize=14)

    ax1.text(0.6, -0.25, '$|C|=1$', fontsize=14)

    ax1.tick_params(axis="x", labelsize=14)
    ax1.tick_params(axis="y", labelsize=14)
    ax1.tick_params(axis="z", labelsize=14)

    ############################################

    left, bottom, width, height = [0.12, 0.5, 0.2, 0.2]
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
    ax5.text2D(-0.078, 0.006, '$\Phi_y$', fontsize=12)

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
    ax5.text2D(0.05, 0.03, '$L$', fontsize=12)

    ax5.view_init(elev=10., azim=-15)
    ax5.axis('off')
    ax5.patch.set_alpha(0)

    ####################################################################################################################

    ax2 = plt.subplot(gs[1])

    x = []
    y = []
    z = []

    with open('ent_spec_flow_example.dat', 'r') as csvfile:
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
        ax2.scatter(xvalue, yvalue, marker='_', c='C{}'.format(value+3), label='{}'.format(value))

    ax2.set_ylim([0, 8])
    ax2.legend(loc='upper right', handletextpad=0, borderpad=0.2, framealpha=1, edgecolor='k', markerscale=2, fontsize=12)
    ax2.set_xlim([0, 1])
    ax2.set_xlabel("$\Phi / 2\pi$", fontsize=14)
    ax2.set_ylabel("$\epsilon_{\\alpha}$", fontsize=14)

    ax2.tick_params(axis="x", labelsize=14)
    ax2.tick_params(axis="y", labelsize=14)
    ax2.tick_params(axis="z", labelsize=14)

    gs.update(wspace=0.25, hspace=0.2)
    ####################################################################################################################

    ax3 = plt.subplot(gs[2])

    x = []
    y = []
    z = []

    with open('ent_spec_mom_example.dat', 'r') as csvfile:
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
        if value == 0:
            ax3.scatter(xvalue, yvalue, marker='_', c='C{}'.format(value+3), label='{}'.format(value), zorder=2)
        else:
            ax3.scatter(xvalue, yvalue, marker='.', c='C{}'.format(value+3), label='{}'.format(value))

    # ax3.scatter(x, y, marker='_', c=z, label=z)

    ax3.legend(loc='upper left', handletextpad=0, borderpad=0.2, framealpha=1, edgecolor='k', markerscale=2, fontsize=12)
    ax3.set_xlim([-1, 1])
    ax3.set_xlabel("$k_{\\alpha} / \pi$", fontsize=14)
    ax3.set_ylim([0, 10])
    ax3.set_ylabel("$\epsilon_{\\alpha}$", fontsize=14)

    ax3.tick_params(axis="x", labelsize=14)
    ax3.tick_params(axis="y", labelsize=14)
    ax3.tick_params(axis="z", labelsize=14)

    ####################################################################################################################

    ax4 = plt.subplot(gs[3])

    x = []
    y1 = []
    y2 = []

    with open('ent_scal_example.dat', 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        for row in plots:
            x.append(int(row[0]))
            y1.append(float(row[1]))
            y2.append(float(row[2]))

    ax4.plot(x, y1, '.', marker='o', markersize=6, label='$S_{vN}$', c='C0')
    c, m = polyfit(x, y1, 1)
    xvalues = np.arange(9)
    ax4.plot(xvalues, m*xvalues + c, '-', c='C0')

    ax4.plot(x, y2, '.', marker='s', markersize=6, label='$S_{\infty}$', c='C1')
    c, m = polyfit(x, y2, 1)
    xvalues = np.arange(9)
    ax4.plot(xvalues, m * xvalues + c, '-', c='C1')

    ax4.axhline(-0.5, color='k', linewidth=1, ls='--')

    ax4.legend(loc='upper left', handletextpad=0, borderpad=0.2, framealpha=1, edgecolor='k', fontsize=12)
    ax4.set_xlim([0, 8])
    ax4.set_xlabel("$L$", fontsize=14)
    ax4.set_ylabel("$S$", fontsize=14)

    ax4.tick_params(axis="x", labelsize=14)
    ax4.tick_params(axis="y", labelsize=14)
    ax4.tick_params(axis="z", labelsize=14)

    ####################################################################################################################

    fig.text(0.04, 0.875, "(a)", color="k", fontsize=14)
    fig.text(0.04, 0.455, "(c)", color="k", fontsize=14)
    fig.text(0.49, 0.875, "(b)", color="k", fontsize=14)
    fig.text(0.49, 0.455, "(d)", color="k", fontsize=14)

    plt.savefig("/home/bart/Documents/papers/TBG/figures/CI_4x4.png", bbox_inches='tight', dpi=300)
    plt.show()