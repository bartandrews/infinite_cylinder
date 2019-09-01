import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

# lattice vectors (normalized)

a1 = (1/2) * np.array([1, np.sqrt(3)])
a2 = (1/2) * np.array([-1, np.sqrt(3)])
avec = np.vstack((a1, a2))

a=1

if __name__ == '__main__':

    fig = plt.figure()

    ax = fig.add_subplot(1,1,1)

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


    ax.scatter(delta[:, 0], delta[:, 1], color='blue')
    # ax.scatter(secondNN[0:3, 0], secondNN[0:3, 1], color='green')
    # ax.scatter(secondNN[3:6, 0], secondNN[3:6, 1], color='green', marker='x')
    ax.scatter(fifthNN[0:3, 0], fifthNN[0:3, 1], color='red')
    ax.scatter(fifthNN[3:6, 0], fifthNN[3:6, 1], color='red', marker='x')
    ax.set_xlabel('$m$')
    ax.set_ylabel('$n$')

    # ax.arrow(0, 0, avec[0, 0], avec[0, 1], color='black', head_width=0.1, length_includes_head=True)
    # ax.arrow(0, 0, avec[1, 0], avec[1, 1], color='black', head_width=0.1, length_includes_head=True)
    # ax.text(avec[0, 0]/3, 0.4, "$\mathbf{a}_1$")
    # ax.text(avec[0, 0]/3, -0.5, "$\mathbf{a}_2$")

    ax.annotate("", xy=(0, 0), xycoords='data', xytext=(avec[0, 0], avec[0, 1]), textcoords='data',
                arrowprops=dict(arrowstyle="<->", connectionstyle="arc3", color='k', lw=2))
    plt.text(0.1 * avec[0, 0], 0.5 * avec[0, 1], 'a', fontsize=14)

    ax.annotate("", xy=(0, 0), xycoords='data', xytext=(avec[0, 0], 0), textcoords='data',
                arrowprops=dict(arrowstyle="<->", connectionstyle="arc3", color='k', lw=2))
    plt.text(0.35 * avec[0, 0], -0.25 * avec[0, 1], 'b', fontsize=14)

    ax.annotate("", xy=(avec[0, 0], 0), xycoords='data', xytext=(avec[0, 0], (1 / 3) * avec[0, 1]), textcoords='data',
                arrowprops=dict(arrowstyle="<->", connectionstyle="arc3", color='k', lw=2))
    plt.text(1.2 * avec[0, 0], 0.25 * (1 / 3) * avec[0, 1], 'c', fontsize=14)


    def xformat(value, tick_number):
        return "{}".format(int(value / (a/2)))

    def yformat(value, tick_number):
        return "{}".format(int(round(value/(np.sqrt(3)/6))))

    ax.set_xlim([-6*avec[0, 0], 6*avec[0, 0]])
    plt.xticks(np.arange(-6 * avec[0, 0], 6 * avec[0, 0]+0.1, step=2*avec[0, 0]))
    ax.set_ylim([-3*avec[0, 1], 3*avec[0, 1]])
    plt.yticks(np.arange(-3 * avec[0, 1], 3 * avec[0, 1], step=avec[0, 1]))

    ax.xaxis.set_major_formatter(plt.FuncFormatter(xformat))
    ax.yaxis.set_major_formatter(plt.FuncFormatter(yformat))

    triangle1 = Polygon(((0, 0), (avec[0, 0], avec[0, 1]), (5*avec[0, 0], avec[0, 1]), (4*avec[0, 0], 0)),
                        fc=(1, 1, 0, 0.2))

    ax.add_artist(triangle1)

    ax.set_aspect('equal', adjustable='box')

    # fig.canvas.manager.full_screen_toggle()
    plt.show()