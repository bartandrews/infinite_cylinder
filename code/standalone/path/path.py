import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

if __name__ == '__main__':

    fig = plt.figure()

    ax = fig.add_subplot(1, 1, 1)

    hexagon = Polygon(((-1, 0), (-0.5, np.sqrt(3)/2), (0.5, np.sqrt(3)/2), (1, 0), (0.5, -np.sqrt(3)/2), (-0.5, -np.sqrt(3)/2)), fc='khaki', ec="k", ls="--")

    ax.add_artist(hexagon)

    ax.set_xlim([-1, 1])
    ax.set_ylim([-1, 1])

    ax.arrow(0.5, np.sqrt(3)/2, -0.5, -np.sqrt(3)/2, color='k', head_width=0.05, length_includes_head=True, lw=2)
    ax.arrow(0, 0, 3/4, np.sqrt(3) / 4, color='k', head_width=0.05, length_includes_head=True, lw=2)
    ax.arrow(3/4, np.sqrt(3) / 4, 0.25, -np.sqrt(3) / 4, color='k', head_width=0.05, length_includes_head=True, lw=2)

    ax.text(-0.1, 0, "Γ", fontsize=36)
    ax.text(1/2, np.sqrt(3)/2, "K", fontsize=36)
    ax.text(3/4, np.sqrt(3)/4, "M", fontsize=36)
    ax.text(1, 0, "K'", fontsize=36)

    ax.set_aspect('equal', adjustable='box')
    ax.axis('off')

    plt.savefig('/home/bart/Documents/papers/TBG/figures/path_tight.pdf', bbox_inches='tight', pad_inches=0, transparent=True)

    plt.show()