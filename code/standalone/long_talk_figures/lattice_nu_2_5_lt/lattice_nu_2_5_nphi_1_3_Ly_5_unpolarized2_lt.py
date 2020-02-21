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

    radius = 0.3

    verts_py1 = [
        (center_x, center_y),  # P0
        (center_x-radius, center_y+radius),  # P1
        (center_x+radius, center_y+radius),  # P2
        (center_x, center_y),  # P3
    ]

    verts_py2 = [
        (center_x, center_y),  # P0
        (center_x-radius, center_y-radius),  # P1
        (center_x+radius, center_y-radius),  # P2
        (center_x, center_y),  # P3
    ]

    verts_px1 = [
        (center_x, center_y),  # P0
        (center_x+radius, center_y+radius),  # P1
        (center_x+radius, center_y-radius),  # P2
        (center_x, center_y),  # P3
    ]

    verts_px2 = [
        (center_x, center_y),  # P0
        (center_x-radius, center_y+radius),  # P1
        (center_x-radius, center_y-radius),  # P2
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

    patch_py1 = patches.PathPatch(path_py1, facecolor='none', lw=1, zorder=3)
    patch_py2 = patches.PathPatch(path_py2, facecolor='none', lw=1, zorder=3)
    patch_px1 = patches.PathPatch(path_px1, facecolor='none', lw=1, zorder=3)
    patch_px2 = patches.PathPatch(path_px2, facecolor='none', lw=1, zorder=3)

    ax.add_patch(patch_py1)
    ax.add_patch(patch_py2)
    ax.add_patch(patch_px1)
    ax.add_patch(patch_px2)


if __name__ == '__main__':

    fig = plt.figure()

    # gs = gridspec.GridSpec(1, 2, width_ratios=[1.4, 1])

    ax = plt.subplot(111)

    delta = np.zeros((3, 2))
    delta[0, :] = (-1 / 3) * avec[0, :] + (-1 / 3) * avec[1, :]
    delta[1, :] = (2 / 3) * avec[0, :] + (-1 / 3) * avec[1, :]
    delta[2, :] = (-1 / 3) * avec[0, :] + (2 / 3) * avec[1, :]

    secondNN = np.zeros((6, 2))
    # positive direction for A sites / negative direction for B sites
    secondNN[0, :] = avec[0, :]
    secondNN[1, :] = -avec[0, :] + avec[1, :]
    secondNN[2, :] = -avec[1, :]
    # negative direction for A sites / positive direction for B sites
    secondNN[3, :] = avec[1, :]
    secondNN[4, :] = -avec[0, :]
    secondNN[5, :] = avec[0, :] - avec[1, :]

    fifthNN = np.zeros((6, 2))
    fifthNN[0, :] = -avec[0, :] - avec[1, :]
    fifthNN[1, :] = -avec[0, :] + 2 * avec[1, :]
    fifthNN[2, :] = 2 * avec[0, :] - avec[1, :]
    fifthNN[3, :] = avec[0, :] + avec[1, :]
    fifthNN[4, :] = -2 * avec[0, :] + avec[1, :]
    fifthNN[5, :] = avec[0, :] - 2 * avec[1, :]

    for m in range(-6, 6):
        for n in range(-6, 6):
            xcoord = 0 + 2*n*avec[0, 0] + m*avec[0, 0]
            ycoord = 0 + m*avec[0, 1]
            for i in range(3):
                ax.arrow(xcoord, ycoord, delta[i, 0], delta[i, 1], color='gray', width=0.0001, ls=':', alpha=0.5)

    for m in range(-6, 6):
        for n in range(-6, 6):
            xcoord = 0 + 2*n*avec[0, 0] + m*avec[0, 0]
            ycoord = 0 + m*avec[0, 1]
            ax.arrow(xcoord, ycoord, avec[0, 0], avec[0, 1], color='black', width=0.0001, ls='--', alpha=0.5)
            ax.arrow(xcoord, ycoord, 2*avec[0, 0], 0, color='black', width=0.0001, ls='--', alpha=0.5)
            ax.arrow(xcoord, ycoord, avec[0, 0], -avec[0, 1], color='black', width=0.0001, ls='--', alpha=0.5)


    # ax.scatter(delta[:, 0], delta[:, 1], color='blue', zorder=2)
    # ax.scatter(secondNN[0:3, 0], secondNN[0:3, 1], color='green')
    # ax.scatter(secondNN[3:6, 0], secondNN[3:6, 1], color='green', marker='x')
    # ax.scatter(fifthNN[0:3, 0], fifthNN[0:3, 1], color='red', zorder=2)
    # ax.scatter(fifthNN[3:6, 0], fifthNN[3:6, 1], color='red', marker='x', zorder=2)
    ax.set_xlabel('$m$', fontsize=11)
    ax.set_ylabel('$n$', fontsize=11)

    # ax.arrow(0, 0, avec[0, 0], avec[0, 1], color='black', head_width=0.1, length_includes_head=True)
    # ax.arrow(0, 0, avec[1, 0], avec[1, 1], color='black', head_width=0.1, length_includes_head=True)
    # ax.text(avec[0, 0]/3, 0.4, "$\mathbf{a}_1$")
    # ax.text(avec[0, 0]/3, -0.5, "$\mathbf{a}_2$")

    ax.annotate("", xy=(-1*avec[0, 0], -3*avec[0, 1]), xycoords='data', xytext=(4*avec[0, 0], 2*avec[0, 1]), textcoords='data',
               arrowprops=dict(arrowstyle="<->", connectionstyle="arc3", color='k', lw=2))
    ax.text(1.8 * avec[0, 0], -0.6 * avec[0, 1], '$L_y=5$', fontsize=14)

    ax.annotate("", xy=(-2 * avec[0, 0], 2 * avec[0, 1]), xycoords='data', xytext=(4 * avec[0, 0], 2 * avec[0, 1]),
                textcoords='data',
                arrowprops=dict(arrowstyle="<->", connectionstyle="arc3", color='k', lw=2))
    ax.text(0.3 * avec[0, 0], 2.2 * avec[0, 1], '$L_x=1$', fontsize=14)

    # ax.annotate("", xy=(0, 0), xycoords='data', xytext=(avec[0, 0], 0), textcoords='data',
    #             arrowprops=dict(arrowstyle="<->", connectionstyle="arc3", color='k', lw=2))
    # ax.text(0.35 * avec[0, 0], -0.4 * avec[0, 1], 'b', fontsize=11)
    #
    # ax.annotate("", xy=(avec[0, 0], 0), xycoords='data', xytext=(avec[0, 0], (1 / 3) * avec[0, 1]), textcoords='data',
    #             arrowprops=dict(arrowstyle="<->", connectionstyle="arc3", color='k', lw=2))
    # ax.text(1.2 * avec[0, 0], 0.25 * (1 / 3) * avec[0, 1], 'c', fontsize=11)


    def xformat(value, tick_number):
        return "{}".format(int(value / (a/2)))

    def yformat(value, tick_number):
        return "{}".format(int(round(value/(np.sqrt(3)/6))))

    ax.set_xlim([-7*avec[0, 0], 4*avec[0, 0]])
    ax.set_xticks(np.arange(-7 * avec[0, 0], 4 * avec[0, 0]+0.1, step=2*avec[0, 0]))
    ax.set_ylim([-3*avec[0, 1], 2*avec[0, 1]])
    ax.set_yticks(np.arange(-3 * avec[0, 1], 2 * avec[0, 1]+0.1, step=avec[0, 1]))

    ax.xaxis.set_major_formatter(plt.FuncFormatter(xformat))
    ax.yaxis.set_major_formatter(plt.FuncFormatter(yformat))

    ax.tick_params(axis='both', which='major', labelsize=10)

    muc = Polygon(((-7*avec[0, 0], -3*avec[0, 1]), (-1*avec[0, 0], 3*avec[0, 1]), (5*avec[0, 0], 3*avec[0, 1]), (-7*avec[0, 0]+3*2*avec[0, 0], -3*avec[0, 1]+0)), fc=(1, 1, 0, 0.2))

    for y in range(5):
        for x in range(3):
            if x == 0 and y == 0:
                ax.scatter(-7*avec[0, 0]+2*x*avec[0, 0]+y*avec[0, 0]+0, -3*avec[0, 1]+0+y*avec[0, 1], color='k', zorder=2, facecolors='white')
            elif x == 1 and y == 2:
                ax.scatter(-7*avec[0, 0]+2*x*avec[0, 0]+y*avec[0, 0]+0, -3*avec[0, 1]+0+y*avec[0, 1], color='k', zorder=2, facecolors='white')
            else:
                ax.scatter(-7 * avec[0, 0] + 2 * x * avec[0, 0] + y * avec[0, 0] + 0, -3 * avec[0, 1] + 0 + y * avec[0, 1], color='k', zorder=2, facecolors='white')
                ax.scatter(-7*avec[0, 0]+delta[1, 0]+2*x*avec[0, 0]+y*avec[0, 0], -3*avec[0, 1]+delta[1, 1]+y*avec[0, 1], color='k', zorder=2, facecolors='white')

    draw_orbitals_at(-7*avec[0, 0]+delta[1, 0]+2*0*avec[0, 0]+0*avec[0, 0], -3*avec[0, 1]+delta[1, 1]+0*avec[0, 1])
    draw_orbitals_at(-7*avec[0, 0]+delta[1, 0]+2*1*avec[0, 0]+2*avec[0, 0], -3*avec[0, 1]+delta[1, 1]+2*avec[0, 1])

    # 1st up
    ax.scatter(-7*avec[0, 0]+delta[1, 0]+2*0*avec[0, 0]+0*avec[0, 0], -3*avec[0, 1]+delta[1, 1]+0*avec[0, 1] +0.22, color='k', zorder=4, facecolors='w')
    # 1st right
    ax.scatter(-7 * avec[0, 0] + delta[1, 0] + 2 * 0 * avec[0, 0] + 0 * avec[0, 0] +0.21, -3 * avec[0, 1] + delta[1, 1] + 0 * avec[0, 1]+0.01, color='k', zorder=4, facecolors='k')

    # 2nd up
    ax.scatter(-7*avec[0, 0]+delta[1, 0]+2*1*avec[0, 0]+2*avec[0, 0], -3*avec[0, 1]+delta[1, 1]+2*avec[0, 1]+0.21, color='k', zorder=4, facecolors='k')
    # 2nd right
    ax.scatter(-7*avec[0, 0]+delta[1, 0]+2*1*avec[0, 0]+2*avec[0, 0]+0.21, -3*avec[0, 1]+delta[1, 1]+2*avec[0, 1]+0.01, color='k', zorder=4, facecolors='w')

    #uc = Polygon(((0, 0), (avec[0, 0], avec[0, 1]), (3 * avec[0, 0], avec[0, 1]), (2 * avec[0, 0], 0)), fill=False, ec='g')

    fig.text(0.2, 0.9, "$n_\phi=1/3 \;\;\;\;\; \odot\mathbf{B}$", fontsize=14)

    ax.add_artist(muc)
    #ax.add_artist(uc)

    ax.set_aspect('equal', adjustable='box')

########################################################################################################################

    # ax1 = plt.subplot(gs[1])
    #
    # box = ax1.get_position()
    # box.y0 = box.y0 - 0.015
    # box.y1 = box.y1 - 0.015
    # ax1.set_position(box)
    #
    # ax1.grid(True)
    # ax1.spines['left'].set_position('zero')
    # ax1.spines['right'].set_color('none')
    # ax1.spines['bottom'].set_position('zero')
    # ax1.spines['top'].set_color('none')
    #
    # ax1.arrow(0, 0, 3 * np.pi + 1, 0, color='k', head_width=0.1, length_includes_head=True)
    # ax1.arrow(0, 0, 0, 4 * np.pi + 1, color='k', head_width=0.1, length_includes_head=True)
    #
    # ax1.arrow(0, 0, np.pi, -np.pi, color='r', head_width=0.1, length_includes_head=True, ls="--", zorder=3)
    # ax1.arrow(0, 0, 2 * np.pi, 2 * np.pi, color='r', head_width=0.1, length_includes_head=True, ls="--", zorder=3)
    # ax1.plot([np.pi, 3 * np.pi], [-np.pi, np.pi], color='b', linestyle=':', linewidth=1)
    # ax1.plot([2 * np.pi, 3 * np.pi], [2 * np.pi, np.pi], color='b', linestyle=':', linewidth=1)
    #
    # ax1.set_xlim([0 - 1, 3 * np.pi + 1])
    # ax1.set_ylim([-np.pi - 1, 4 * np.pi + 1])
    #
    # muc_r = Polygon(((0, 0), (0, 4 * np.pi), (0.5*np.pi, 4 * np.pi), (0.5*np.pi, 0)), fc=(1, 1, 0, 0.2))
    #
    # ax1.set_xticks(np.arange(0, 4 * np.pi, step=np.pi))
    # ax1.set_yticks(np.arange(-np.pi, 5 * np.pi, step=np.pi))
    #
    #
    # def custom(value, tick_number):
    #
    #     if int(round(value / (np.pi))) == 0:
    #         return 0
    #     elif int(round(value / (np.pi))) == 1:
    #         return "$\pi$"
    #     elif int(round(value / (np.pi))) == -1:
    #         return "$-\pi$"
    #     else:
    #         return "{}$\pi$".format(int(round(value / (np.pi))))
    #
    #
    # ax1.xaxis.set_major_formatter(plt.FuncFormatter(custom))
    # ax1.yaxis.set_major_formatter(plt.FuncFormatter(custom))
    #
    # ax1.tick_params(axis='both', which='major', labelsize=10)
    #
    # ax1.text(3 * np.pi + 1, 0.5, '$x$', fontsize=11)
    # ax1.text(0.5, 4 * np.pi + 1, '$y$', fontsize=11)
    #
    # ax1.text(np.pi / 2 - 1.15, -np.pi / 2 - 0.75, '$\mathbf{b}_1$', fontsize=11)
    # ax1.text(3 * np.pi / 2 - 1.35, 3 * np.pi / 2 + 0.2, '$\mathbf{b}_2$', fontsize=11)
    #
    # # ax.set_xlabel("x")
    # # ax.set_ylabel("y")
    #
    # ax1.add_artist(muc_r)
    #
    # uc_r = Polygon(((0, 0), (0, 4 * np.pi), (np.pi, 4 * np.pi), (np.pi, 0)), fill=False, ec='g', zorder=2)
    # ax1.add_artist(uc_r)
    #
    # ax1.set_aspect('equal', adjustable='box')

    # plt.text(-23, 14, "(a)", fontsize=12)
    # plt.text(-3, 14, "(b)", fontsize=12)

    plt.savefig("/home/bart/Documents/presentations/2020_02_21/figures/lattice_nu_2_5_nphi_1_3_Ly_5_unpolarized2_lt.png", bbox_inches='tight', dpi=300)
    plt.show()