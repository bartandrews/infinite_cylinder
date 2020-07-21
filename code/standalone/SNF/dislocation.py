import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

if __name__ == '__main__':

    fig = plt.figure()

    ax = fig.add_subplot(1,1,1)

    for m in range(-5, 5):
        for n in range(-5, 5):
            if m % 2 == 0:
                clr='red'
            else:
                clr='blue'
            if n == 0 and -2 <= m <= 2:
                ax.arrow(m, n, 1, 1, color='black', width=0.01, ls='--', alpha=0.5)
                ax.arrow(m, n, 1, 0, color='black', width=0.01, ls='--', alpha=0.5)
            else:
                if n!=4 and n!=-3:
                    ax.arrow(m, n, 1, 0, color='black', width=0.01, ls='--', alpha=0.5)
                if m==3 and n==0:
                    continue
                else:
                    if m!=-5:
                        ax.arrow(m, n, 0, 1, color=clr, width=0.1, ls='--', alpha=0.5, head_length=0, head_width=0)

    ax.set_xticks(np.arange(-5, 6, step=1))
    ax.set_yticks(np.arange(-3, 5, step=1))
    ax.set_aspect('equal', adjustable='box')

    ax.text(-2.85, 0.15, "X", fontsize=36)
    ax.text(3.15, 0.15, "X", fontsize=36)

    #path
    ax.arrow(-4, 2, 3, 0, color='black', width=0.1, ls='--', alpha=1, head_length=0, head_width=0)
    ax.arrow(-4, -1, 3, 0, color='black', width=0.1, ls='--', alpha=1, head_length=0, head_width=0)
    ax.arrow(-4, -1, 0, 3, color='black', width=0.1, ls='--', alpha=1, head_length=0, head_width=0)
    ax.arrow(-1, -1, 0, 1, color='black', width=0.1, ls='--', alpha=1, head_length=0, head_width=0)
    ax.arrow(-1, 0, 1, 1, color='black', width=0.1, ls='--', alpha=1, head_length=0.5, head_width=0.5, length_includes_head=True)

    ax.xaxis.set_ticks([])
    ax.yaxis.set_ticks([])
    ax.axis('off')

    # fig.canvas.manager.full_screen_toggle()
    plt.savefig("/home/bart/Documents/SNF2020/research_plan/figures/1d.png", bbox_inches='tight', dpi=1000)
    plt.show()