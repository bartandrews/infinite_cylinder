import numpy as np
import matplotlib.pyplot as plt
import sys
from itertools import product
import matplotlib.gridspec as gridspec
from matplotlib.patches import Polygon
#from numpy import arange
#from pylab import meshgrid,cm,imshow,contour,clabel,colorbar,axis,title,show
from mpl_toolkits.mplot3d import axes3d, Axes3D
from pylab import cm, imshow, contour, colorbar
import matplotlib.patches as patches

# lattice vectors (normalized)

a1 = (1/2) * np.array([np.sqrt(3), 1])
a2 = (1/2) * np.array([np.sqrt(3), -1])
avec = np.vstack((a1, a2))

# reciprocal lattice vectors

b1 = (2.*np.pi) * np.array([1/np.sqrt(3), 1])
b2 = (2.*np.pi) * np.array([1/np.sqrt(3), -1])
bvec = np.vstack((b1, b2))

# high symmetry points

K1 = np.array([2/3, 1/3])
K2 = np.array([1/3, 2/3])
GA = np.array([0., 0.])
MM = np.array([0.5, 0.5])


def hamiltonian(k, num_bands_1):

    t1, t2, t2dash = 1, -0.025, 1/18

    delta = np.zeros((3, 2))
    delta[0, :] = (1 / 3) * avec[0, :] + (1 / 3) * avec[1, :]
    delta[1, :] = (-2 / 3) * avec[0, :] + (1 / 3) * avec[1, :]
    delta[2, :] = (1 / 3) * avec[0, :] + (-2 / 3) * avec[1, :]

    fifthNN = np.zeros((6, 2))
    fifthNN[0, :] = -avec[0, :] - avec[1, :]
    fifthNN[1, :] = -avec[0, :] + 2*avec[1, :]
    fifthNN[2, :] = 2*avec[0, :] - avec[1, :]
    fifthNN[3, :] = avec[0, :] + avec[1, :]
    fifthNN[4, :] = -2*avec[0, :] + avec[1, :]
    fifthNN[5, :] = avec[0, :] - 2*avec[1, :]

    Hamiltonian = np.zeros((num_bands_1, num_bands_1), dtype=np.complex128)

    f = 0
    for i in range(3):
        f += t1 * np.exp(1j * k.dot(delta[i, :]))
    f2 = 0
    for i in [0, 1, 2]:  # range(3)
        f2 += (t2 - 1j * t2dash/2) * np.exp(1j * k.dot(fifthNN[i, :]))
    xi = 0
    for i in [0, 1, 2]:  # range(3)
        xi += (t2dash/2) * np.exp(1j * k.dot(fifthNN[i, :]))

    # spin up block

    Hamiltonian[0][0] = f2 + np.conj(f2)
    Hamiltonian[0][1] = f
    Hamiltonian[1][0] = np.conj(f)
    Hamiltonian[1][1] = f2 + np.conj(f2)

    Hamiltonian[0][2] = xi  #- np.conj(xi)
    Hamiltonian[1][3] = xi  #- np.conj(xi)

    Hamiltonian[2][0] = -np.conj(xi)  # np.conj(Hamiltonian[0][2])
    Hamiltonian[3][1] = -np.conj(xi)  # np.conj(Hamiltonian[1][3])

    Hamiltonian[2][2] = Hamiltonian[0][0]
    Hamiltonian[2][3] = Hamiltonian[0][1]
    Hamiltonian[3][2] = Hamiltonian[1][0]
    Hamiltonian[3][3] = Hamiltonian[1][1]

    # spin down block (the minimal model is spin degenerate)

    Hamiltonian[4][4] = Hamiltonian[0][0]
    Hamiltonian[4][5] = Hamiltonian[0][1]
    Hamiltonian[5][4] = Hamiltonian[1][0]
    Hamiltonian[5][5] = Hamiltonian[1][1]

    Hamiltonian[4][6] = Hamiltonian[0][2]
    Hamiltonian[5][7] = Hamiltonian[1][3]

    Hamiltonian[6][4] = Hamiltonian[2][0]
    Hamiltonian[7][5] = Hamiltonian[3][1]

    Hamiltonian[6][6] = Hamiltonian[2][2]
    Hamiltonian[6][7] = Hamiltonian[2][3]
    Hamiltonian[7][6] = Hamiltonian[3][2]
    Hamiltonian[7][7] = Hamiltonian[3][3]

    return Hamiltonian


if __name__ == '__main__':

    ####################
    # 1) Minimal model #
    ####################

    num_bands_1 = 8
    max_kx, max_ky = 101, 101

    eigval_bands_1 = np.zeros((num_bands_1, max_kx, max_ky))

    for idx_kx in range(max_kx):
        kx = -np.pi + idx_kx * (2*np.pi / max_kx)
        for idx_ky in range(max_ky):
            ky = -np.pi + idx_ky * (2 * np.pi / max_ky)
            eigvals = np.linalg.eigvals(hamiltonian(np.array([kx, ky]), num_bands_1))
            idx = np.argsort(eigvals)
            for band in range(num_bands_1):
                eigval_bands_1[band, idx_kx, idx_ky] = np.real(eigvals[idx[band]])

    ##############
    # Final Plot #
    ##############

    fig = plt.figure()
    ax = plt.subplot(1, 1, 1)

    kx = range(max_kx)
    ky = range(max_ky)

    KX, KY = np.meshgrid(kx, ky)

    E1 = np.zeros(num_bands_1, dtype=object)

    i = 7

    E1[i] = eigval_bands_1[i, KX, KY]
    hex = patches.RegularPolygon((50, 50), 6, radius=50, orientation=0, facecolor='none')
    ax.add_patch(hex)
    im = ax.imshow(E1[i], cmap=plt.get_cmap('RdBu_r'), clip_path=hex, clip_on=True)
    cset = ax.contour(E1[i], np.arange(2, 4, 0.1), linewidths=1, colors='black')
    for collection in cset.collections:
        collection.set_clip_path(hex)

    plt.axis('off')

    radius = 50
    center = 50

    plt.arrow(center+radius*np.cos(np.pi/6), center, 0, radius*np.sin(np.pi/6), color='k', head_width=0.05, length_includes_head=True, lw=2)
    plt.arrow(center, center, radius*np.cos(np.pi/6), -radius*np.sin(np.pi/6), color='k', head_width=0.05, length_includes_head=True, lw=2)
    plt.arrow(center, center, radius*np.cos(np.pi/6), 0, color='k', head_width=0.05, length_includes_head=True, lw=2)

    plt.text(center, center, "Γ", fontsize=36)
    plt.text(center+radius*np.cos(np.pi/6), center-radius*np.sin(np.pi/6), "K", fontsize=36)
    plt.text(center+radius*np.cos(np.pi/6), center, "M", fontsize=36)
    plt.text(center+radius*np.cos(np.pi/6), center+radius*np.sin(np.pi/6), "K'", fontsize=36)

    plt.show()