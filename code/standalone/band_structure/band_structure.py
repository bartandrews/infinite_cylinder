import numpy as np
import matplotlib.pyplot as plt
import sys
from itertools import product
import matplotlib.gridspec as gridspec
from matplotlib.patches import Polygon

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

    t1, t2, t2dash = 0.331, 0.01, 0

    delta = np.zeros((3, 2))
    delta[0, :] = (1 / 3) * avec[0, :] + (1 / 3) * avec[1, :]
    delta[1, :] = (-2 / 3) * avec[0, :] + (1 / 3) * avec[1, :]
    delta[2, :] = (1 / 3) * avec[0, :] + (-2 / 3) * avec[1, :]

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
    for i in range(0, 3):
        f2 += t2 * np.exp(1j * k.dot(fifthNN[i, :]))
    xi = 0
    for i in range(3):
        xi += t2dash * np.exp(1j * k.dot(fifthNN[i, :]))

    # spin up block

    Hamiltonian[0][0] = f2 + np.conj(f2)
    Hamiltonian[0][1] = f
    Hamiltonian[1][0] = np.conj(f)
    Hamiltonian[1][1] = np.conj(f2) + f2

    Hamiltonian[0][2] = xi - np.conj(xi)
    Hamiltonian[1][3] = xi - np.conj(xi)

    Hamiltonian[2][0] = np.conj(Hamiltonian[0][2])
    Hamiltonian[3][1] = np.conj(Hamiltonian[1][3])

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


def hamiltonian2(k, M):

    t1, t2 = 0.331, -0.01
    a = 1
    c = np.sqrt(3) * a / 6
    eta = 3 * k[0] * M * a / 2
    alpha = float(p / q)

    matrix = np.zeros(shape=(M, M), dtype=np.complex128)

    def A(phi, m):
        return t2 * (2 * np.cos(4 * np.pi * phi * m + 12 * k[1] * c)
                     + 2 * np.cos(2 * np.pi * phi * (m - 3 / 2) + 6 * k[1] * c)
                     + 2 * np.cos(2 * np.pi * phi * (m + 3 / 2) + 6 * k[1] * c) + 3)

    def B(phi, m):
        return t1 * (2 * np.exp(1j * np.pi * phi / 3) * np.cos(np.pi * phi * (m + 1 / 2) + 3 * k[1] * c))

    def C(phi):
        return t1 * np.exp(1j * np.pi * phi / 3)

    def D(phi, m):
        return t2 * (2 * np.cos(np.pi * phi * (m - 5 / 2) + 3 * k[1] * c)
                     + 2 * np.cos(np.pi * phi * (m + 7 / 2) + 3 * k[1] * c)
                     + 2 * np.cos(3 * np.pi * phi * (m - 1 / 2) + 9 * k[1] * c)
                     + 2 * np.cos(3 * np.pi * phi * (m + 3 / 2) + 9 * k[1] * c))

    def G(phi, m):
        return t2 * (2 * np.cos(2 * np.pi * phi * (m + 1) + 6 * k[1] * c)
                     + 2 * np.cos(3 * np.pi * phi))

    for i in range(M):
        matrix[i, i] = A(alpha, i + 1)
    for i in range(M - 1):
        matrix[i, i + 1] = B(alpha, i + 1)
        matrix[i + 1, i] = np.conj(B(alpha, i + 1))
    for i in range(M - 2):
        matrix[i, i + 2] = np.conj(C(alpha))
        matrix[i + 2, i] = C(alpha)
    for i in range(M - 3):
        matrix[i, i + 3] = D(alpha, i + 2)
        matrix[i + 3, i] = D(alpha, i + 2)
    for i in range(M - 6):
        matrix[i, i + 6] = G(alpha, i + 3)
        matrix[i + 6, i] = G(alpha, i + 3)

    # top-right
    matrix[0][M - 1] = np.conj(B(alpha, M))

    matrix[0][M - 2] = C(alpha)
    matrix[1][M - 1] = C(alpha)

    matrix[0][M - 3] = D(alpha, M - 1) * np.exp(-1j * eta)
    matrix[1][M - 2] = D(alpha, M) * np.exp(-1j * eta)
    matrix[2][M - 1] = D(alpha, 1) * np.exp(-1j * eta)

    matrix[0][M - 6] = G(alpha, M - 3)
    matrix[1][M - 5] = G(alpha, M - 2)
    matrix[2][M - 4] = G(alpha, M - 1)
    matrix[3][M - 3] = G(alpha, M)
    matrix[4][M - 2] = G(alpha, 1)
    matrix[5][M - 1] = G(alpha, 2)

    # bottom-left
    matrix[M - 1][0] = B(alpha, M) * np.exp(1j * eta)

    matrix[M - 2][0] = np.conj(C(alpha))
    matrix[M - 1][1] = np.conj(C(alpha))

    matrix[M - 3][0] = D(alpha, M - 1) * np.exp(1j * eta)
    matrix[M - 2][1] = D(alpha, M) * np.exp(1j * eta)
    matrix[M - 1][2] = D(alpha, 1) * np.exp(1j * eta)

    matrix[M - 6][0] = G(alpha, M - 3) * np.exp(1j * eta)
    matrix[M - 5][1] = G(alpha, M - 2) * np.exp(1j * eta)
    matrix[M - 4][2] = G(alpha, M - 1) * np.exp(1j * eta)
    matrix[M - 3][3] = G(alpha, M) * np.exp(1j * eta)
    matrix[M - 2][4] = G(alpha, 1)
    matrix[M - 1][5] = G(alpha, 2)

    return matrix


def berry_curv(ev, ev_alpha, ev_beta, ev_alpha_beta):

    bc = - np.imag(np.log(np.conj(ev).dot(ev_alpha) * np.conj(ev_alpha).dot(ev_alpha_beta)
                          * np.conj(ev_alpha_beta).dot(ev_beta) * np.conj(ev_beta).dot(ev)))

    return bc


if __name__ == '__main__':

    ####################
    # 1) Minimal model #
    ####################

    num_bands_1 = 8

    ####################################################################################################################

    count, nk = 0, 30

    eigval_bands_1 = np.zeros((num_bands_1, 3*nk))

    for i in range(nk):
        k = K1 - K1 * float(i) / float(nk - 1)
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

    samples_x, samples_y = 101, 101
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

    berry_flux_matrix_1 = np.zeros((num_bands_1, samples_x-1, samples_y-1))

    for band in range(num_bands_1):
        for idx_x in range(max_idx_x):
            for idx_y in range(max_idx_y):
                berry_flux_matrix_1[band, idx_x, idx_y] = berry_curv(u_matrix_1[:, band, idx_x, idx_y],
                                                                     u_matrix_1[:, band, idx_x + 1, idx_y],
                                                                     u_matrix_1[:, band, idx_x, idx_y + 1],
                                                                     u_matrix_1[:, band, idx_x + 1, idx_y + 1])

    chern_numbers_1 = np.zeros(num_bands_1)
    for band in range(num_bands_1):
        chern_numbers_1[band] = np.sum(berry_flux_matrix_1[band, :, :]) / (2 * np.pi)

    ########################################
    # 2) Minimal model with magnetic field #
    ########################################

    p, q = 2, 11

    if p % 2 == 0:
        M = q
    else:
        M = 2 * q

    ####################################################################################################################

    count, nk = 0, 30

    eigval_bands_2 = np.zeros((2*M, 3 * nk))

    for i in range(nk):
        k = K1 - K1 * float(i) / float(nk - 1)
        k = np.matmul(k, bvec)
        eigvals = np.linalg.eigvals(hamiltonian2(k, M))
        idx = np.argsort(eigvals)
        for band in range(M):
            eigval_bands_2[band, count] = np.real(+np.sqrt(3 + eigvals[idx[band]]))
            eigval_bands_2[M + band, count] = np.real(-np.sqrt(3 + eigvals[idx[band]]))
        count += 1

    for i in range(nk):
        k = GA + (MM - GA) * float(i) / float(nk - 1)
        k = np.matmul(k, bvec)
        eigvals = np.linalg.eigvals(hamiltonian2(k, M))
        idx = np.argsort(eigvals)
        for band in range(M):
            eigval_bands_2[band, count] = np.real(+np.sqrt(3 + eigvals[idx[band]]))
            eigval_bands_2[M + band, count] = np.real(-np.sqrt(3 + eigvals[idx[band]]))
        count += 1

    for i in range(nk):
        k = MM + (K2 - MM) * float(i) / float(nk - 1)
        k = np.matmul(k, bvec)
        eigvals = np.linalg.eigvals(hamiltonian2(k, M))
        idx = np.argsort(eigvals)
        for band in range(M):
            eigval_bands_2[band, count] = np.real(+np.sqrt(3 + eigvals[idx[band]]))
            eigval_bands_2[M + band, count] = np.real(-np.sqrt(3 + eigvals[idx[band]]))
        count += 1

    ####################################################################################################################

    samples_x, samples_y = 101, 101
    max_idx_x, max_idx_y = samples_x - 1, samples_y - 1

    energy_matrix_2 = np.zeros((2*M, samples_x, samples_y))
    u_matrix_2 = np.zeros((M, M, samples_x, samples_y), dtype=np.complex128)

    for idx_x in range(samples_x):
        frac_kx = idx_x / max_idx_x
        for idx_y in range(samples_y):
            frac_ky = idx_y / max_idx_y

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

            eigvals, eigvecs = np.linalg.eig(hamiltonian2(k_u, M))
            idx = np.argsort(eigvals)
            for band in range(M):
                energy_matrix_2[band, idx_x, idx_y] = np.real(+np.sqrt(3 + eigvals[idx[band]]))
                energy_matrix_2[M + band, idx_x, idx_y] = np.real(-np.sqrt(3 + eigvals[idx[band]]))
                u_matrix_2[:, band, idx_x, idx_y] = eigvecs[:, idx[band]]

    berry_flux_matrix_2 = np.zeros((M, samples_x - 1, samples_y - 1))

    for band in range(M):
        for idx_x in range(max_idx_x):
            for idx_y in range(max_idx_y):
                berry_flux_matrix_2[band, idx_x, idx_y] = berry_curv(u_matrix_2[:, band, idx_x, idx_y],
                                                                     u_matrix_2[:, band, idx_x + 1, idx_y],
                                                                     u_matrix_2[:, band, idx_x, idx_y + 1],
                                                                     u_matrix_2[:, band, idx_x + 1, idx_y + 1])

    chern_numbers_2 = np.zeros(M)
    for band in range(M):
        chern_numbers_2[band] = np.sum(berry_flux_matrix_2[band, :, :]) / (2 * np.pi)

    ##############
    # Final Plot #
    ##############

    fig = plt.figure()

    gs = gridspec.GridSpec(2, 2, width_ratios=[2, 1], height_ratios=[1, 1])

    ax0 = plt.subplot(gs[0])

    # color index
    # colors = plt.cm.RdBu_r(np.linspace(0, 1, 2))
    cidx = ['b', 'b', 'b', 'b', 'r', 'r', 'r', 'r']

    for nb in range(num_bands_1):
        plt.scatter(np.linspace(0, 89, 90), eigval_bands_1[nb, :], c=cidx[nb])
    ax0.axvline(30, color='k', linewidth=0.5)
    ax0.axvline(60, color='k', linewidth=0.5)
    ax0.axhline(0, color='k', linewidth=0.5, ls='--')
    plt.xlim((0, 89))
    plt.setp(ax0.get_xticklabels(), visible=False)
    ax0.text(10, 0.1, '(a) $B_z=0$')

########################################################################################################################

    ax1 = plt.subplot(gs[1], sharey=ax0)
    n, bins, patches = ax1.hist(np.ndarray.flatten(energy_matrix_1), 100,
                                density=False, facecolor='k', alpha=0.75, orientation='horizontal')
    plt.setp(ax1.get_yticklabels(), visible=False)
    ax1.axhline(0, color='k', linewidth=0.5, ls='--')
    plt.xticks([1000, 2000], ["1000", "2000"])
    plt.setp(ax1.get_xticklabels(), visible=False)

    gs.update(wspace=0)

########################################################################################################################

    ax2 = plt.subplot(gs[2])

    # color index
    colors = plt.cm.RdBu_r(np.linspace(0, 1, 2 * M))
    cidx = np.zeros(2 * M, dtype=int)
    for i in range(M):
        cidx[i] = M + i
        cidx[M + i] = (M - 1) - i

    for nb in range(2*M):
        plt.scatter(np.linspace(0, 89, 90), eigval_bands_2[nb, :], c=colors[cidx[nb]])
    # ax2.set_xlabel('Path')
    ax2.axvline(30, color='k', linewidth=0.5)
    ax2.axvline(60, color='k', linewidth=0.5)
    ax2.axhline(0, color='k', linewidth=0.5, ls='--')
    plt.xlim((0, 89))
    plt.xticks([0, 30, 60, 89], ["K", "Γ", "M", "K'"])
    ax2.text(10, 0.2, '(b) $B_z \\neq 0$')

########################################################################################################################

    ax3 = plt.subplot(gs[3], sharey=ax2)
    n, bins, patches = ax3.hist(np.ndarray.flatten(energy_matrix_2), 100,
                                density=False, facecolor='k', alpha=0.75, orientation='horizontal')
    plt.setp(ax3.get_yticklabels(), visible=False)
    ax3.set_xlabel('# States')
    ax3.axhline(0, color='k', linewidth=0.5, ls='--')
    plt.xticks([1000, 2000], ["1000", "2000"])
    plt.setp(ax3.get_xticklabels(), visible=False)
    plt.tick_params(
        axis='x',  # changes apply to the x-axis
        which='both',  # both major and minor ticks are affected
        bottom=False,  # ticks along the bottom edge are off
        labelbottom=False)

    gs.update(wspace=0)
    gs.update(hspace=0)

    fig.text(0.02, 0.5, '$E$ / meV', va='center', rotation='vertical')

    plt.show()