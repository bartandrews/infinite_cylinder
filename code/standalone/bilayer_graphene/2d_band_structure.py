# coding=utf-8
import numpy as np
import matplotlib.pyplot as plt
import sys
from itertools import product
import matplotlib.gridspec as gridspec
from matplotlib.patches import Polygon
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

    t, t0 = 3, 0.3

    delta = np.zeros((3, 2))
    delta[0, :] = (1 / 3) * avec[0, :] + (1 / 3) * avec[1, :]
    delta[1, :] = (-2 / 3) * avec[0, :] + (1 / 3) * avec[1, :]
    delta[2, :] = (1 / 3) * avec[0, :] + (-2 / 3) * avec[1, :]

    Hamiltonian = np.zeros((num_bands_1, num_bands_1), dtype=np.complex128)

    f = 0
    for i in range(3):
        f += np.exp(1j * k.dot(delta[i, :]))

    Hamiltonian[0][0] = 0
    Hamiltonian[0][1] = 0
    Hamiltonian[0][2] = -t*f
    Hamiltonian[0][3] = 0

    Hamiltonian[1][1] = 0
    Hamiltonian[1][2] = t0
    Hamiltonian[1][3] = -t*f

    Hamiltonian[2][2] = 0
    Hamiltonian[2][3] = 0

    Hamiltonian[3][3] = 0

    # H.c.

    Hamiltonian[1][0] = np.conj(Hamiltonian[0][1])

    Hamiltonian[2][0] = np.conj(Hamiltonian[0][2])
    Hamiltonian[2][1] = np.conj(Hamiltonian[1][2])

    Hamiltonian[3][0] = np.conj(Hamiltonian[0][3])
    Hamiltonian[3][1] = np.conj(Hamiltonian[1][3])
    Hamiltonian[3][2] = np.conj(Hamiltonian[2][3])

    return Hamiltonian


if __name__ == '__main__':

    ####################
    # 1) Minimal model #
    ####################

    num_bands_1 = 4

    ####################################################################################################################

    count, nk = 0, 100

    eigval_bands_1 = np.zeros((num_bands_1, 3*nk))

    for i in range(nk):
        k = K1 + (GA - K1) * float(i) / float(nk - 1)
        k = np.matmul(k, bvec)
        eigvals = np.linalg.eigvals(hamiltonian(k, num_bands_1))
        idx = np.argsort(eigvals)
        for band in range(num_bands_1):
            eigval_bands_1[band, count] = np.real(eigvals[idx[band]])
        count += 1

    for i in range(nk):
        k = GA + (MM - GA) * float(i) / float(nk - 1)
        k = np.matmul(k, bvec)
        eigvals = np.linalg.eigvals(hamiltonian(k, num_bands_1))
        idx = np.argsort(eigvals)
        for band in range(num_bands_1):
            eigval_bands_1[band, count] = np.real(eigvals[idx[band]])
        count += 1

    for i in range(nk):
        k = MM + (K2 - MM) * float(i) / float(nk - 1)
        k = np.matmul(k, bvec)
        eigvals = np.linalg.eigvals(hamiltonian(k, num_bands_1))
        idx = np.argsort(eigvals)
        for band in range(num_bands_1):
            eigval_bands_1[band, count] = np.real(eigvals[idx[band]])
        count += 1

    ####################################################################################################################

    samples_x, samples_y = 201, 201
    max_idx_x, max_idx_y = samples_x - 1, samples_y - 1

    energy_matrix_1 = np.zeros((num_bands_1, samples_x, samples_y))
    u_matrix_1 = np.zeros((num_bands_1, num_bands_1, samples_x, samples_y), dtype=np.complex128)

    for idx_x in range(samples_x):
        frac_kx = idx_x / max_idx_x
        for idx_y in range(samples_y):
            frac_ky = idx_y/max_idx_y

            # # wavevector used for u (depends on boundary conditions)
            # if idx_x == max_idx_x and idx_y == max_idx_y:
            #     k_u = np.matmul(np.array([0, 0]), bvec)
            # elif idx_x == max_idx_x:
            #     k_u = np.matmul(np.array([0, frac_ky]), bvec)
            # elif idx_y == max_idx_y:
            #     k_u = np.matmul(np.array([frac_kx, 0]), bvec)
            # else:
            #     k_u = np.matmul(np.array([frac_kx, frac_ky]), bvec)

            k_u = np.matmul(np.array([frac_kx, frac_ky]), bvec)

            eigvals, eigvecs = np.linalg.eig(hamiltonian(k_u, num_bands_1))
            idx = np.argsort(eigvals)
            for band in range(num_bands_1):
                energy_matrix_1[band, idx_x, idx_y] = np.real(eigvals[idx[band]])
                u_matrix_1[:, band, idx_x, idx_y] = eigvecs[:, idx[band]]

    ##############
    # Final Plot #
    ##############

    fig = plt.figure()  # figsize=(6, 4)

    gs = gridspec.GridSpec(1, 2, width_ratios=[2, 1])

    ax0 = plt.subplot(gs[0])

    # color index
    # colors = plt.cm.RdBu_r(np.linspace(0, 1, 2))
    cidx = ['b', 'g', 'r', 'c']

    for nb in range(num_bands_1):
        plt.scatter(np.linspace(0, 3*nk-1, 3*nk), eigval_bands_1[nb, :], c=cidx[nb], s=0.5)
    ax0.set_ylabel('$E$ / eV', fontsize=11)
    ax0.axvline(nk, color='k', linewidth=0.5)
    ax0.axvline(2*nk, color='k', linewidth=0.5)
    ax0.axhline(0, color='k', linewidth=0.5, ls='--')
    plt.xlim((0, 3*nk-1))
    plt.setp(ax0.get_xticklabels(), visible=True)
    ax0.tick_params(axis='both', which='major', labelsize=10)
    plt.xticks([0, nk, 2*nk, 3*nk-1], ["K", "Γ", "M", "K'"], fontsize=11)

    ####################################################################################################################

    ax1 = plt.subplot(gs[1], sharey=ax0)
    n, bins, patches = ax1.hist(np.ndarray.flatten(energy_matrix_1), 200,
                                density=False, orientation='horizontal', histtype='step', color='k',  lw=0.5)
    plt.setp(ax1.get_yticklabels(), visible=False)
    ax1.axhline(0, color='k', linewidth=0.5, ls='--')
    plt.xticks([1000, 2000], ["1000", "2000"])
    plt.setp(ax1.get_xticklabels(), visible=False)
    ax1.set_xlabel('DOS', fontsize=11)
    plt.tick_params(
        axis='x',  # changes apply to the x-axis
        which='both',  # both major and minor ticks are affected
        bottom=False,  # ticks along the bottom edge are off
        labelbottom=False)

    gs.update(wspace=0)

    plt.savefig("/home/bart/Documents/FSK2020/proposal/figures/bilayer_graphene/2d_band_structure.png", bbox_inches='tight', dpi=300)
    plt.show()
