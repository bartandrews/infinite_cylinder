import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import Polygon

if __name__ == '__main__':

    plt.rcParams['legend.fontsize'] = 10

    fig = plt.figure()
    ax = fig.gca(projection='3d')

    r = 5

    for y in range(0, 12, 2):

        theta = np.linspace(0, 2 * np.pi, 21)
        z = r * np.cos(theta)
        x = r * np.sin(theta)
        ax.scatter(x, y, z, label='parametric curve', c='b')

        theta = np.linspace((1/40)*2*np.pi, (1+1/40)*2*np.pi, 21)
        y = y + 1/3
        z = r * np.cos(theta)
        x = r * np.sin(theta)
        ax.scatter(x, y, z, label='parametric curve', c='r')

        theta = np.linspace((1 / 40) * 2 * np.pi, (1 + 1 / 40) * 2 * np.pi, 21)
        y = y + 2/3
        z = r * np.cos(theta)
        x = r * np.sin(theta)
        ax.scatter(x, y, z, label='parametric curve', c='b')

        theta = np.linspace(0, 2 * np.pi, 21)
        y = y + 1/3
        z = r * np.cos(theta)
        x = r * np.sin(theta)
        ax.scatter(x, y, z, label='parametric curve', c='r')

    # cut
    xx, zz = np.meshgrid(np.linspace(-7, 7, 10), np.linspace(-7, 7, 10))
    yy = 6
    ax.plot_surface(xx, yy, zz, alpha=0.25)

    # arrow
    ax.quiver(0, 12, 0, 0, -15, 0, arrow_length_ratio=0.06, color='k')
    ax.text2D(-0.068, 0.005, '$\Phi_y$', fontsize=14)

    # cylinder
    angle = np.linspace(0, 2 * np.pi, 21)
    length = np.linspace(0, 34/3, 13)
    theta, y = np.meshgrid(angle, length)
    z = r * np.cos(theta)
    x = r * np.sin(theta)
    ax.plot_surface(x, y, z, alpha=0.25)

    # ring
    theta = np.linspace(0, 2 * np.pi, 21)
    y = np.repeat(12, 21)
    z = r * np.cos(theta)
    x = r * np.sin(theta)
    ax.plot(x, y, z, label='parametric curve', c='k')

    theta = np.linspace(0, 2 * np.pi, 1)
    y = np.repeat(12, 1)
    z = r * np.cos(theta)
    x = r * np.sin(theta)
    ax.scatter(x, y, z, label='parametric curve', c='k', marker='v', s=30)
    ax.text2D(0.05, 0.03, '$L$', fontsize=14)

    ax.view_init(elev=10., azim=-15)
    ax.axis('off')

    plt.savefig("/home/bart/Documents/papers/TBG/figures/tube.png", bbox_inches='tight', dpi=300)
    plt.show()
