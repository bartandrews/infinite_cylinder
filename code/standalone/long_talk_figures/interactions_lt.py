import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import matplotlib.gridspec as gridspec
from matplotlib.path import Path
import matplotlib.patches as patches

# lattice vectors (normalized)

a1 = (1/2) * np.array([1, np.sqrt(3)])
a2 = (1/2) * np.array([-1, np.sqrt(3)])
avec = np.vstack((a1, a2))

a=1


def draw_orbitals_at(center_x, center_y):
    verts_py1 = [
        (center_x, center_y),  # P0
        (center_x-0.2, center_y+0.2),  # P1
        (center_x+0.2, center_y+0.2),  # P2
        (center_x, center_y),  # P3
    ]

    verts_py2 = [
        (center_x, center_y),  # P0
        (center_x-0.2, center_y-0.2),  # P1
        (center_x+0.2, center_y-0.2),  # P2
        (center_x, center_y),  # P3
    ]

    verts_px1 = [
        (center_x, center_y),  # P0
        (center_x+0.2, center_y+0.2),  # P1
        (center_x+0.2, center_y-0.2),  # P2
        (center_x, center_y),  # P3
    ]

    verts_px2 = [
        (center_x, center_y),  # P0
        (center_x-0.2, center_y+0.2),  # P1
        (center_x-0.2, center_y-0.2),  # P2
        (center_x, center_y),  # P3
    ]

    codes = [
        Path.MOVETO,
        Path.CURVE4,
        Path.CURVE4,
        Path.CURVE4,
    ]

    path_py1 = Path(verts_py1, codes)
    path_py2 = Path(verts_py2, codes)
    path_px1 = Path(verts_px1, codes)
    path_px2 = Path(verts_px2, codes)

    patch_py1 = patches.PathPatch(path_py1, facecolor='none', lw=2)
    patch_py2 = patches.PathPatch(path_py2, facecolor='none', lw=2)
    patch_px1 = patches.PathPatch(path_px1, facecolor='none', lw=2)
    patch_px2 = patches.PathPatch(path_px2, facecolor='none', lw=2)

    ax.add_patch(patch_py1)
    ax.add_patch(patch_py2)
    ax.add_patch(patch_px1)
    ax.add_patch(patch_px2)


if __name__ == '__main__':

    fig = plt.figure()

    ax = plt.subplot(111)

    delta = np.zeros((3, 2))
    delta[0, :] = (-1 / 3) * avec[0, :] + (-1 / 3) * avec[1, :]
    delta[1, :] = (2 / 3) * avec[0, :] + (-1 / 3) * avec[1, :]
    delta[2, :] = (-1 / 3) * avec[0, :] + (2 / 3) * avec[1, :]

    for m in range(-6, 6):
        for n in range(-6, 6):
            xcoord = 0 + 2*n*avec[0, 0] + m*avec[0, 0]
            ycoord = 0 + m*avec[0, 1]
            for i in range(3):
                ax.arrow(xcoord, ycoord, delta[i, 0], delta[i, 1], color='gray', width=0.0001, ls=':', alpha=0.5)

    # for m in range(-6, 6):
    #     for n in range(-6, 6):
    #         xcoord = 0 + 2*n*avec[0, 0] + m*avec[0, 0]
    #         ycoord = 0 + m*avec[0, 1]
    #         ax.arrow(xcoord, ycoord, avec[0, 0], avec[0, 1], color='black', width=0.0001, ls='--', alpha=0.5)
    #         ax.arrow(xcoord, ycoord, 2*avec[0, 0], 0, color='black', width=0.0001, ls='--', alpha=0.5)
    #         ax.arrow(xcoord, ycoord, avec[0, 0], -avec[0, 1], color='black', width=0.0001, ls='--', alpha=0.5)


    ax.scatter(delta[:, 0], delta[:, 1], color='blue', zorder=2)
    ax.set_xlabel('$m$', fontsize=11)
    ax.set_ylabel('$n$', fontsize=11)

    def xformat(value, tick_number):
        return "{}".format(int(value / (a/2)))

    def yformat(value, tick_number):
        return "{}".format(int(round(value/(np.sqrt(3)/6))))

    ax.set_xlim([-2*avec[0, 0], 2*avec[0, 0]])
    ax.set_xticks(np.arange(-2 * avec[0, 0], 2 * avec[0, 0]+0.1, step=avec[0, 0]))
    ax.set_ylim([-1*avec[0, 1], 1*avec[0, 1]])
    ax.set_yticks(np.arange(-1 * avec[0, 1], 1 * avec[0, 1]+0.1, step=(1/3)*avec[0, 1]))

    ax.xaxis.set_major_formatter(plt.FuncFormatter(xformat))
    ax.yaxis.set_major_formatter(plt.FuncFormatter(yformat))

    ax.tick_params(axis='both', which='major', labelsize=10)

    draw_orbitals_at(0, 0)
    draw_orbitals_at(0.5, 0.29)
    draw_orbitals_at(-0.5, 0.29)
    draw_orbitals_at(0, -0.58)
    plt.annotate(s='', xy=(-0.05, -0.2), xytext=(-0.2, -0.05), arrowprops=dict(arrowstyle='<->', color='green'))
    plt.annotate(s='', xy=(0.33, 0.23), xytext=(0.1, 0.1), arrowprops=dict(arrowstyle='<->', color='brown'))
    ax.text(-0.35,-0.275,"$U=10^2$", fontsize=12, color='green')
    ax.text(0.05, 0.25, "$V=10$", fontsize=12, color='brown')
    ax.text(0.47, 0.5, "$p_y$", fontsize=12)
    ax.text(0.69, 0.28, "$p_x$", fontsize=12)

    plt.text(-0.95, 0.9, "not to scale", fontsize=12)


    ax.set_aspect('equal', adjustable='box')

    plt.savefig("/home/bart/Documents/presentations/2020_02_21/figures/interactions_lt.png", bbox_inches='tight', dpi=300)
    plt.show()