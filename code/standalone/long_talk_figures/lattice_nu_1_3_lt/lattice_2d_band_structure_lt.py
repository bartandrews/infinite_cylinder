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


def min_hamiltonian(k, num_bands_1):

    t1, t2, t2dash = 1, -0.0075, 0.03

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


def hamiltonian(k, num_bands_1):

    t1, t2, t2dash = 1, -0.0075, 0.03
    # t1, t2, t2dash = 2, 0.05, 0.2

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


def hamiltonian2(k, M, p, q):

    t1, t2, t2dash = 1, -0.0075, 0.03
    a = 1
    c = np.sqrt(3) * a / 6  # ... / 6
    eta = 1 * k[0] * M * a / 2  # 3 * ...
    alpha = float(p / q)

    matrix = np.zeros(shape=(M, M), dtype=np.complex128)

    def A(phi, m):
        return t2 * (2 * np.cos(4 * np.pi * phi * m + 12 * k[1] * c)
                     + 2 * np.cos(2 * np.pi * phi * (m - 3 / 2) + 6 * k[1] * c)
                     + 2 * np.cos(2 * np.pi * phi * (m + 3 / 2) + 6 * k[1] * c) + 3) \
               + t2dash * (-2 * np.cos(4 * np.pi * phi * m + 12 * k[1] * c)
                           - 2 * np.cos(2 * np.pi * phi * (m - 3 / 2) + 6 * k[1] * c)
                           - 2 * np.cos(2 * np.pi * phi * (m + 3 / 2) + 6 * k[1] * c) - 9)

    def B(phi, m):
        return t1 * (2 * np.exp(1j * np.pi * phi / 3) * np.cos(np.pi * phi * (m + 1 / 2) + 3 * k[1] * c))

    def C(phi):
        return t1 * np.exp(1j * np.pi * phi / 3)

    def D(phi, m):
        return t2 * (2 * np.cos(np.pi * phi * (m - 5 / 2) + 3 * k[1] * c)
                     + 2 * np.cos(np.pi * phi * (m + 7 / 2) + 3 * k[1] * c)
                     + 2 * np.cos(3 * np.pi * phi * (m - 1 / 2) + 9 * k[1] * c)
                     + 2 * np.cos(3 * np.pi * phi * (m + 3 / 2) + 9 * k[1] * c)) \
               + t2dash * (-2 * np.cos(np.pi * phi * (m - 5 / 2) + 3 * k[1] * c)
                           - 2 * np.cos(np.pi * phi * (m + 7 / 2) + 3 * k[1] * c)
                           - 2 * np.cos(3 * np.pi * phi * (m - 1 / 2) + 9 * k[1] * c)
                           - 2 * np.cos(3 * np.pi * phi * (m + 3 / 2) + 9 * k[1] * c))

    def G(phi, m):
        return t2 * (2 * np.cos(2 * np.pi * phi * (m + 1) + 6 * k[1] * c)
                     + 2 * np.cos(3 * np.pi * phi)) \
               + t2dash * (-2 * np.cos(2 * np.pi * phi * (m + 1) + 6 * k[1] * c)
                           - 2 * np.cos(3 * np.pi * phi))

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
    matrix[0][M - 1] = np.conj(B(alpha, M)) * np.exp(-1j * eta)

    matrix[0][M - 2] = C(alpha) * np.exp(-1j * eta)
    matrix[1][M - 1] = C(alpha) * np.exp(-1j * eta)

    matrix[0][M - 3] = D(alpha, M - 1) * np.exp(-1j * eta)
    matrix[1][M - 2] = D(alpha, M) * np.exp(-1j * eta)
    matrix[2][M - 1] = D(alpha, 1) * np.exp(-1j * eta)

    matrix[0][M - 6] = G(alpha, M - 3) * np.exp(-1j * eta)
    matrix[1][M - 5] = G(alpha, M - 2) * np.exp(-1j * eta)
    matrix[2][M - 4] = G(alpha, M - 1) * np.exp(-1j * eta)
    matrix[3][M - 3] = G(alpha, M) * np.exp(-1j * eta)
    matrix[4][M - 2] = G(alpha, 1) * np.exp(-1j * eta)
    matrix[5][M - 1] = G(alpha, 2) * np.exp(-1j * eta)

    # bottom-left
    matrix[M - 1][0] = B(alpha, M) * np.exp(1j * eta)

    matrix[M - 2][0] = np.conj(C(alpha)) * np.exp(1j * eta)
    matrix[M - 1][1] = np.conj(C(alpha)) * np.exp(1j * eta)

    matrix[M - 3][0] = D(alpha, M - 1) * np.exp(1j * eta)
    matrix[M - 2][1] = D(alpha, M) * np.exp(1j * eta)
    matrix[M - 1][2] = D(alpha, 1) * np.exp(1j * eta)

    matrix[M - 6][0] = G(alpha, M - 3) * np.exp(1j * eta)
    matrix[M - 5][1] = G(alpha, M - 2) * np.exp(1j * eta)
    matrix[M - 4][2] = G(alpha, M - 1) * np.exp(1j * eta)
    matrix[M - 3][3] = G(alpha, M) * np.exp(1j * eta)
    matrix[M - 2][4] = G(alpha, 1) * np.exp(1j * eta)
    matrix[M - 1][5] = G(alpha, 2) * np.exp(1j * eta)

    return matrix


def berry_curv(ev, ev_alpha, ev_beta, ev_alpha_beta):

    bc = - np.imag(np.log(np.conj(ev).dot(ev_alpha) * np.conj(ev_alpha).dot(ev_alpha_beta)
                          * np.conj(ev_alpha_beta).dot(ev_beta) * np.conj(ev_beta).dot(ev)))

    return bc


def multi_berry_curv(ev1, ev1_alpha, ev1_beta, ev1_alpha_beta, ev2, ev2_alpha, ev2_beta, ev2_alpha_beta):

    matrix1 = np.zeros((2, 2), dtype=np.complex128)
    matrix1[0][0] = np.conj(ev1).dot(ev1_alpha)
    matrix1[0][1] = np.conj(ev1).dot(ev2_alpha)
    matrix1[1][0] = np.conj(ev2).dot(ev1_alpha)
    matrix1[1][1] = np.conj(ev2).dot(ev2_alpha)

    matrix2 = np.zeros((2, 2), dtype=np.complex128)
    matrix2[0][0] = np.conj(ev1_alpha).dot(ev1_alpha_beta)
    matrix2[0][1] = np.conj(ev1_alpha).dot(ev2_alpha_beta)
    matrix2[1][0] = np.conj(ev2_alpha).dot(ev1_alpha_beta)
    matrix2[1][1] = np.conj(ev2_alpha).dot(ev2_alpha_beta)

    matrix3 = np.zeros((2, 2), dtype=np.complex128)
    matrix3[0][0] = np.conj(ev1_alpha_beta).dot(ev1_beta)
    matrix3[0][1] = np.conj(ev1_alpha_beta).dot(ev2_beta)
    matrix3[1][0] = np.conj(ev2_alpha_beta).dot(ev1_beta)
    matrix3[1][1] = np.conj(ev2_alpha_beta).dot(ev2_beta)

    matrix4 = np.zeros((2, 2), dtype=np.complex128)
    matrix4[0][0] = np.conj(ev1_beta).dot(ev1)
    matrix4[0][1] = np.conj(ev1_beta).dot(ev2)
    matrix4[1][0] = np.conj(ev2_beta).dot(ev1)
    matrix4[1][1] = np.conj(ev2_beta).dot(ev2)

    multi_bc = - np.imag(np.log(np.linalg.det(matrix1) * np.linalg.det(matrix2) * np.linalg.det(matrix3) * np.linalg.det(matrix4)))

    return multi_bc


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

    p, q = 1, 3

    if p % 2 == 0:
        M = q
    else:
        M = 2 * q

    # reciprocal lattice vectors
    b1_2 = (2. * np.pi / q) * np.array([1, -1 / np.sqrt(3)])
    b2_2 = (2. * np.pi) * np.array([0, 2 / np.sqrt(3)])
    bvec_2 = np.vstack((b1_2, b2_2))

    ####################################################################################################################

    count, nk = 0, 30

    eigval_bands_2 = np.zeros((2*M, 3 * nk))

    for i in range(nk):
        k = K1 - K1 * float(i) / float(nk - 1)
        k = np.matmul(k, bvec)
        eigvals = np.linalg.eigvals(hamiltonian2(k, M, p, q))
        idx = np.argsort(eigvals)
        for band in range(M):
            eigval_bands_2[M+band, count] = np.real(+np.sqrt(3 + eigvals[idx[band]]))
            eigval_bands_2[(M-1) - band, count] = np.real(-np.sqrt(3 + eigvals[idx[band]]))
        count += 1

    for i in range(nk):
        k = GA + (MM - GA) * float(i) / float(nk - 1)
        k = np.matmul(k, bvec)
        eigvals = np.linalg.eigvals(hamiltonian2(k, M, p, q))
        idx = np.argsort(eigvals)
        for band in range(M):
            eigval_bands_2[M+band, count] = np.real(+np.sqrt(3 + eigvals[idx[band]]))
            eigval_bands_2[(M-1) - band, count] = np.real(-np.sqrt(3 + eigvals[idx[band]]))
        count += 1

    for i in range(nk):
        k = MM + (K2 - MM) * float(i) / float(nk - 1)
        k = np.matmul(k, bvec)
        eigvals = np.linalg.eigvals(hamiltonian2(k, M, p, q))
        idx = np.argsort(eigvals)
        for band in range(M):
            eigval_bands_2[M+band, count] = np.real(+np.sqrt(3 + eigvals[idx[band]]))
            eigval_bands_2[(M-1) - band, count] = np.real(-np.sqrt(3 + eigvals[idx[band]]))
        count += 1

    ####################################################################################################################

    samples_x, samples_y = 101, 101
    max_idx_x, max_idx_y = samples_x - 1, samples_y - 1

    energy_matrix_2 = np.zeros((2*M, samples_x, samples_y))
    u_matrix_2 = np.zeros((M, 2*M, samples_x, samples_y), dtype=np.complex128)

    for idx_x in range(samples_x):
        frac_kx = idx_x / max_idx_x
        for idx_y in range(samples_y):
            frac_ky = idx_y / max_idx_y

            k_u = np.matmul(np.array([frac_kx, frac_ky]), bvec_2)

            eigvals, eigvecs = np.linalg.eig(hamiltonian2(k_u, M, p, q))
            idx = np.argsort(eigvals)
            for band in range(M):
                energy_matrix_2[M + band, idx_x, idx_y] = np.real(+np.sqrt(3 + eigvals[idx[band]]))
                u_matrix_2[:, M + band, idx_x, idx_y] = eigvecs[:, idx[band]]
                energy_matrix_2[(M-1) - band, idx_x, idx_y] = np.real(-np.sqrt(3 + eigvals[idx[band]]))
                u_matrix_2[:, (M-1) - band, idx_x, idx_y] = eigvecs[:, idx[band]]

    berry_flux_matrix_2 = np.zeros((M, samples_x - 1, samples_y - 1))

    for band in range(0, M, 2):
        for idx_x in range(max_idx_x):
            for idx_y in range(max_idx_y):
                berry_flux_matrix_2[band, idx_x, idx_y] = multi_berry_curv(u_matrix_2[:, band, idx_x, idx_y],
                                                                           u_matrix_2[:, band, idx_x + 1, idx_y],
                                                                           u_matrix_2[:, band, idx_x, idx_y + 1],
                                                                           u_matrix_2[:, band, idx_x + 1, idx_y + 1],
                                                                           u_matrix_2[:, band+1, idx_x, idx_y],
                                                                           u_matrix_2[:, band+1, idx_x + 1, idx_y],
                                                                           u_matrix_2[:, band+1, idx_x, idx_y + 1],
                                                                           u_matrix_2[:, band+1, idx_x + 1, idx_y + 1])
                berry_flux_matrix_2[band+1, idx_x, idx_y] = berry_flux_matrix_2[band, idx_x, idx_y]

    chern_numbers_2 = np.zeros(M)
    for band in range(0, M, 2):
        chern_numbers_2[band] = np.sum(berry_flux_matrix_2[band, :, :]) / (2 * np.pi)
        print("Chern number ( band", band, ") = ", chern_numbers_2[band])

    ##############
    # Final Plot #
    ##############

    fig = plt.figure(figsize=(6, 3))
    ax2 = plt.subplot(111)

    for nb in range(2*M):
        if nb < M:
            plt.scatter(np.linspace(0, 89, 90), eigval_bands_2[nb, :], c='b', s=0.5)
        else:
            plt.scatter(np.linspace(0, 89, 90), eigval_bands_2[nb, :], c='r', s=0.5)
    ax2.set_ylabel('$E$ / meV', fontsize=11)
    ax2.tick_params(axis='both', which='major', labelsize=10)
    ax2.axvline(30, color='k', linewidth=0.5)
    ax2.axvline(60, color='k', linewidth=0.5)
    ax2.axhline(0, color='k', linewidth=0.5)

    for band in range(0, M-2, 2):
        ax2.text(14+band*15, eigval_bands_2[band, band*int(89/(M-1))], "$\mathbf{{{}}}$".format(int(round(chern_numbers_2[band]))), fontsize=11)

    plt.xlim((0, 89))
    plt.xticks([0, 30, 60, 89], ["K", "Γ", "M", "K'"], fontsize=11)

    plt.savefig("/home/bart/Documents/presentations/2020_02_21/figures/lattice_2d_band_structure_lt.png".format(p, q), bbox_inches='tight', dpi=300)
    plt.show()
