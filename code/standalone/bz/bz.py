import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

if __name__ == '__main__':

    fig = plt.figure()

    ax = fig.add_subplot(1, 1, 1)

    ax.grid(True)
    ax.spines['left'].set_position('zero')
    ax.spines['right'].set_color('none')
    ax.spines['bottom'].set_position('zero')
    ax.spines['top'].set_color('none')

    ax.arrow(0, 0, 3*np.pi+1, 0, color='k', head_width=0.1, length_includes_head=True)
    ax.arrow(0, 0, 0, 4 * np.pi+1, color='k', head_width=0.1, length_includes_head=True)

    ax.arrow(0, 0, np.pi, -np.pi, color='r', head_width=0.1, length_includes_head=True, ls="--")
    ax.arrow(0, 0, 2*np.pi, 2*np.pi, color='r', head_width=0.1, length_includes_head=True, ls="--")
    ax.plot([np.pi, 3*np.pi], [-np.pi, np.pi], color='b', linestyle=':', linewidth=1)
    ax.plot([2*np.pi, 3*np.pi], [2*np.pi, np.pi], color='b', linestyle=':', linewidth=1)

    ax.set_xlim([0 - 1, 3*np.pi+1])
    ax.set_ylim([-np.pi-1, 4*np.pi+1])

    box = Polygon(((0, 0), (0, 4*np.pi), (np.pi, 4*np.pi), (np.pi, 0)),
                        fc=(1, 1, 0, 0.2))

    plt.xticks(np.arange(0, 4*np.pi, step=np.pi))
    plt.yticks(np.arange(-np.pi, 5*np.pi, step=np.pi))

    def custom(value, tick_number):

        if int(round(value / (np.pi))) == 0:
            return 0
        elif int(round(value / (np.pi))) == 1:
            return "$\pi$"
        elif int(round(value / (np.pi))) == -1:
            return "$-\pi$"
        else:
            return "{}$\pi$".format(int(round(value / (np.pi))))

    ax.xaxis.set_major_formatter(plt.FuncFormatter(custom))
    ax.yaxis.set_major_formatter(plt.FuncFormatter(custom))

    ax.text(3*np.pi+1, 0.5, '$x$')
    ax.text(0.5, 4*np.pi+1, '$y$')

    ax.text(np.pi/2-0.75, -np.pi/2-0.75, '$\mathbf{b}_1$')
    ax.text(3*np.pi/2-1, 3*np.pi/2+0.2, '$\mathbf{b}_2$')

    # ax.set_xlabel("x")
    # ax.set_ylabel("y")

    ax.add_artist(box)

    ax.set_aspect('equal', adjustable='box')

    plt.savefig('bz_tight.pdf', bbox_inches='tight', pad_inches=0)
    plt.show()
