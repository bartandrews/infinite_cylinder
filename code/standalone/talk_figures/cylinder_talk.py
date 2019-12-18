import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import csv
from matplotlib.patches import ConnectionPatch
from numpy.polynomial.polynomial import polyfit
from mpl_toolkits.mplot3d import Axes3D

if __name__ == '__main__':

    fig = plt.figure()

    ax5 = plt.subplot(111, projection='3d')

    r = 5

    end = 6
    point_size = 10

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

    # # arrow
    # ax5.quiver(0, end, 0, 0, -(end + end/4), 0, arrow_length_ratio=0.06, color='k', lw=0.5)
    # ax5.text2D(-0.078, 0.006, '$\Phi_y$', fontsize=12)

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
    ax5.text2D(0.05, 0.03, '$L_y$', fontsize=14)

    ax5.view_init(elev=10., azim=-15)
    ax5.axis('off')
    ax5.patch.set_alpha(0)

    fig.text(0.35, 0.75, "$|\\alpha_\mathrm{L}\\rangle$", fontsize=14)
    fig.text(0.6, 0.75, "$|\\alpha_\mathrm{R}\\rangle$", fontsize=14)

    fig.subplots_adjust(top=1, bottom=0, right=1, left=0, hspace=0, wspace=0)

    plt.savefig("/home/bart/Documents/papers/TBG_talk/figures/cylinder_talk.png", bbox_inches='tight', dpi=300)
    plt.show()