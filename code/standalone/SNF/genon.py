import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as patches
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

def plot_implicit(axis, fn, bbox=(-3.3,3.3)):
    ''' create a plot of an implicit function
    fn  ...implicit function (plot where fn==0)
    bbox ..the x,y,and z limits of plotted interval'''
    xmin, xmax, ymin, ymax, zmin, zmax = bbox*3
    A = np.linspace(xmin, xmax, 1000) # resolution of the contour
    B = np.linspace(xmin, xmax, 100) # number of slices
    A1,A2 = np.meshgrid(A,A) # grid on which the contour is plotted

    for z in B: # plot contours in the XY plane
        X,Y = A1,A2
        Z = fn(X,Y,z)
        cset = axis.contour(X, Y, Z+z, [z], zdir='z')
        # [z] defines the only level to plot for this contour for this value of z

    for y in B: # plot contours in the XZ plane
        X,Z = A1,A2
        Y = fn(X,y,Z)
        cset = axis.contour(X, Y+y, Z, [y], zdir='y')

    for x in B: # plot contours in the YZ plane
        Y,Z = A1,A2
        X = fn(x,Y,Z)
        cset = axis.contour(X+x, Y, Z, [x], zdir='x')

    # must set plot limits because the contour will likely extend
    # way beyond the displayed level.  Otherwise matplotlib extends the plot limits
    # to encompass all values in the contour.
    axis.set_zlim3d(zmin,zmax)
    axis.set_xlim3d(xmin,xmax)
    axis.set_ylim3d(ymin,ymax)

if __name__ == '__main__':

    fig = plt.figure(figsize=(13.75, 2.5))
    gs = gridspec.GridSpec(1, 4, width_ratios=[1, 1, 1, 1], wspace=0.25)

    ####################################################################################################################

    ax1 = plt.subplot(gs[0])
    # xy, width, height
    layer1 = patches.Rectangle((0.5, 0), 2, 1, edgecolor=None, facecolor="red", zorder=2)
    layer2 = patches.Rectangle((0.25, 0.25), 2, 1, edgecolor=None, facecolor="blue", zorder=1)
    ax1.add_patch(layer1)
    ax1.add_patch(layer2)
    ax1.set_xlim([0, 2.75])
    ax1.set_ylim([-0.25, 1.5])
    ax1.set_aspect('equal')
    #ax1.xaxis.set_ticks([])
    #ax1.yaxis.set_ticks([])
    #ax1.axis('off')

    # enlarge current axis by 20%
    box = ax1.get_position()
    ax1.set_position([box.x0-0.08, box.y0, box.width * 1.7, box.height])

    ####################################################################################################################

    ax2 = plt.subplot(gs[1], projection='3d')

    # Generate torus mesh
    angle = np.linspace(0, 2 * np.pi, 32)
    theta, phi = np.meshgrid(angle, angle)
    r, R = .25, 1.
    X = (R + r * np.cos(phi)) * np.cos(theta)
    Y = (R + r * np.cos(phi)) * np.sin(theta)
    Z = r * np.sin(phi)

    ax2.set_xlim3d(-1, 1)
    ax2.set_ylim3d(-1, 1)
    ax2.set_zlim3d(-1, 1)
    ax2.plot_surface(X, Y, Z, color='w', rstride=1, cstride=1)

    ####################################################################################################################

    ax3 = plt.subplot(gs[2], projection='3d')

    def goursat_tangle(x, y, z):
        a, b, c = 0.0, -5.0, 11.8
        return x ** 4 + y ** 4 + z ** 4 + a * (x ** 2 + y ** 2 + z ** 2) ** 2 + b * (x ** 2 + y ** 2 + z ** 2) + c

    def double_torus(x, y, z):
        R, r, a = 1.2, 0.5, np.sqrt(2)
        return -a**2 + ((-r**2 + R**2)**2 - 2*(r**2 + R**2)*((-r - R + x)**2 + y**2) + 2*(-r**2 + R**2)*z**2 + ((-r - R + x)**2 + y**2 + z**2)**2)*((-r**2 + R**2)**2 - 2*(r**2 + R**2)*((r + R + x)**2 + y**2) + 2*(-r**2 + R**2)*z**2 + ((r + R + x)**2 + y**2 + z**2)**2)

    plot_implicit(ax3, double_torus)

    ####################################################################################################################

    ax4 = plt.subplot(gs[3])

    ax4.plot([1, 2, 3], [1, 2, 3])
    
    ####################################################################################################################

    fig.text(0.08, 0.9, "(a)", fontsize=12)
    fig.text(0.29, 0.9, "(b)", fontsize=12)
    fig.text(0.5, 0.9, "(c)", fontsize=12)
    fig.text(0.71, 0.9, "(d)", fontsize=12)

    plt.savefig("/home/bart/Documents/SNF2020/research_plan/figures/genon.png", bbox_inches='tight', dpi=300)
    plt.show()
