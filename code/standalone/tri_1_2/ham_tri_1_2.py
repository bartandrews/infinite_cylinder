"""Based on: Energy Spectrum of a Triangular Lattice in a Uniform Magnetic Field:
Effect of Next-Nearest-Neighbor Hopping"""

from mpl_toolkits.mplot3d import axes3d, Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np

p = 2
q = 9

# reciprocal lattice vectors

b1 = (2.*np.pi/q) * np.array([1, -1/np.sqrt(3)])
b2 = (2.*np.pi) * np.array([0, 2/np.sqrt(3)])
bvec = np.vstack((b1, b2))


def hamiltonian(k, M):

    t, tdash = 0, 1
    a = 1
    c = np.sqrt(3) * a / 2  # ... / 6
    delta = 1 * k[0] * M * a / 2  # 3 * ...
    alpha = float(p / q)

    matrix = np.zeros(shape=(M, M), dtype=np.complex128)

    def A(phi, m):
        return 2 * tdash * np.cos(2 * np.pi * phi * m - 2 * k[1] * c)

    def B(phi, m):
        return 2 * t * np.cos(np.pi * phi * (m + 1 / 2) - k[1] * c)  # correction: (m + 1/2)

    def D(phi, m):
        return 2 * tdash * np.cos(np.pi * phi * (m + 1 / 2) - k[1] * c)  # correction: (m + 1/2)

    for i in range(M):
        matrix[i][i] = A(alpha, i + 1)
    for i in range(M - 1):
        matrix[i, i + 1] = B(alpha, i + 1)
        matrix[i + 1, i] = B(alpha, i + 1)
    for i in range(M - 2):
        matrix[i, i + 2] = t
        matrix[i + 2, i] = t
    for i in range(M - 3):
        matrix[i, i + 3] = D(alpha, i + 2)
        matrix[i + 3, i] = D(alpha, i + 2)

    # top-right
    matrix[0][M - 3] = D(alpha, M - 1) * np.exp(-1j * delta)
    matrix[1][M - 2] = D(alpha, M) * np.exp(-1j * delta)
    matrix[2][M - 1] = D(alpha, 1) * np.exp(-1j * delta)
    matrix[0][M - 2] = t * np.exp(-1j * delta)
    matrix[1][M - 1] = t * np.exp(-1j * delta)
    matrix[0][M - 1] = B(alpha, M) * np.exp(-1j * delta)

    # bottom-left
    matrix[M - 3][0] = D(alpha, M - 1) * np.exp(1j * delta)
    matrix[M - 2][1] = D(alpha, M) * np.exp(1j * delta)
    matrix[M - 1][2] = D(alpha, 1) * np.exp(1j * delta)
    matrix[M - 2][0] = t * np.exp(1j * delta)
    matrix[M - 1][1] = t * np.exp(1j * delta)
    matrix[M - 1][0] = B(alpha, M) * np.exp(1j * delta)

    return matrix


def berry_curv(ev, ev_alpha, ev_beta, ev_alpha_beta):

    bc = - np.imag(np.log(np.conj(ev).dot(ev_alpha) * np.conj(ev_alpha).dot(ev_alpha_beta)
                          * np.conj(ev_alpha_beta).dot(ev_beta) * np.conj(ev_beta).dot(ev)))

    return bc


def triple_berry_curv(ev1, ev1_alpha, ev1_beta, ev1_alpha_beta, ev2, ev2_alpha, ev2_beta, ev2_alpha_beta,
                      ev3, ev3_alpha, ev3_beta, ev3_alpha_beta):

    matrix1 = np.zeros((3, 3), dtype=np.complex128)
    matrix1[0][0] = np.conj(ev1).dot(ev1_alpha)
    matrix1[0][1] = np.conj(ev1).dot(ev2_alpha)
    matrix1[0][2] = np.conj(ev1).dot(ev3_alpha)
    matrix1[1][0] = np.conj(ev2).dot(ev1_alpha)
    matrix1[1][1] = np.conj(ev2).dot(ev2_alpha)
    matrix1[1][2] = np.conj(ev2).dot(ev3_alpha)
    matrix1[2][0] = np.conj(ev3).dot(ev1_alpha)
    matrix1[2][1] = np.conj(ev3).dot(ev2_alpha)
    matrix1[2][2] = np.conj(ev3).dot(ev3_alpha)

    matrix2 = np.zeros((3, 3), dtype=np.complex128)
    matrix2[0][0] = np.conj(ev1_alpha).dot(ev1_alpha_beta)
    matrix2[0][1] = np.conj(ev1_alpha).dot(ev2_alpha_beta)
    matrix2[0][2] = np.conj(ev1_alpha).dot(ev3_alpha_beta)
    matrix2[1][0] = np.conj(ev2_alpha).dot(ev1_alpha_beta)
    matrix2[1][1] = np.conj(ev2_alpha).dot(ev2_alpha_beta)
    matrix2[1][2] = np.conj(ev2_alpha).dot(ev3_alpha_beta)
    matrix2[2][0] = np.conj(ev3_alpha).dot(ev1_alpha_beta)
    matrix2[2][1] = np.conj(ev3_alpha).dot(ev2_alpha_beta)
    matrix2[2][2] = np.conj(ev3_alpha).dot(ev3_alpha_beta)

    matrix3 = np.zeros((3, 3), dtype=np.complex128)
    matrix3[0][0] = np.conj(ev1_alpha_beta).dot(ev1_beta)
    matrix3[0][1] = np.conj(ev1_alpha_beta).dot(ev2_beta)
    matrix3[0][2] = np.conj(ev1_alpha_beta).dot(ev3_beta)
    matrix3[1][0] = np.conj(ev2_alpha_beta).dot(ev1_beta)
    matrix3[1][1] = np.conj(ev2_alpha_beta).dot(ev2_beta)
    matrix3[1][2] = np.conj(ev2_alpha_beta).dot(ev3_beta)
    matrix3[2][0] = np.conj(ev3_alpha_beta).dot(ev1_beta)
    matrix3[2][1] = np.conj(ev3_alpha_beta).dot(ev2_beta)
    matrix3[2][2] = np.conj(ev3_alpha_beta).dot(ev3_beta)

    matrix4 = np.zeros((3, 3), dtype=np.complex128)
    matrix4[0][0] = np.conj(ev1_beta).dot(ev1)
    matrix4[0][1] = np.conj(ev1_beta).dot(ev2)
    matrix4[0][2] = np.conj(ev1_beta).dot(ev3)
    matrix4[1][0] = np.conj(ev2_beta).dot(ev1)
    matrix4[1][1] = np.conj(ev2_beta).dot(ev2)
    matrix4[1][2] = np.conj(ev2_beta).dot(ev3)
    matrix4[2][0] = np.conj(ev3_beta).dot(ev1)
    matrix4[2][1] = np.conj(ev3_beta).dot(ev2)
    matrix4[2][2] = np.conj(ev3_beta).dot(ev3)

    triple_bc = - np.imag(np.log(np.linalg.det(matrix1) * np.linalg.det(matrix2) * np.linalg.det(matrix3) * np.linalg.det(matrix4)))

    return triple_bc


if __name__ == '__main__':

    numb_samples = 101

    if p % 2 == 0:
        M = q
    else:
        M = 2 * q

    eigenvalues = np.zeros((M, numb_samples, numb_samples))
    eigenvectors = np.zeros((M, M, numb_samples, numb_samples), dtype=np.complex128)

    for band in range(M):
        for idx_x in range(numb_samples):
            frac_kx = idx_x / (numb_samples-1)
            for idx_y in range(numb_samples):
                frac_ky = idx_y / (numb_samples-1)

                # # wavevector used for u (depends on boundary conditions)
                # if idx_x == (numb_samples-1) and idx_y == (numb_samples-1):
                #     k_u = np.matmul(np.array([0, 0]), bvec)
                # elif idx_x == (numb_samples-1):
                #     k_u = np.matmul(np.array([0, frac_ky]), bvec)
                # elif idx_y == (numb_samples-1):
                #     k_u = np.matmul(np.array([frac_kx, 0]), bvec)
                # else:
                #     k_u = np.matmul(np.array([frac_kx, frac_ky]), bvec)

                k_u = np.matmul(np.array([frac_kx, frac_ky]), bvec)

                eigvals, eigvecs = np.linalg.eig(hamiltonian(k_u, M))
                idx = np.argsort(eigvals)
                eigenvalues[band][idx_x][idx_y] = np.real(eigvals[idx[band]])
                eigenvectors[:, band, idx_x, idx_y] = eigvecs[:, idx[band]]

    berry_flux_matrix = np.zeros((M, numb_samples - 1, numb_samples - 1))

    # for band in range(M):
    #     for idx_x in range(numb_samples-1):
    #         for idx_y in range(numb_samples-1):
    #             berry_flux_matrix[band, idx_x, idx_y] = berry_curv(eigenvectors[:, band, idx_x, idx_y],
    #                                                                eigenvectors[:, band, idx_x + 1, idx_y],
    #                                                                eigenvectors[:, band, idx_x, idx_y + 1],
    #                                                                eigenvectors[:, band, idx_x + 1, idx_y + 1])

    band1 = 0
    band2 = 1
    band3 = 2

    for idx_x in range(numb_samples - 1):
        for idx_y in range(numb_samples - 1):
            berry_flux_matrix[band1, idx_x, idx_y] = triple_berry_curv(eigenvectors[:, band1, idx_x, idx_y],
                                                                   eigenvectors[:, band1, idx_x + 1, idx_y],
                                                                   eigenvectors[:, band1, idx_x, idx_y + 1],
                                                                   eigenvectors[:, band1, idx_x + 1, idx_y + 1],
                                                                   eigenvectors[:, band2, idx_x, idx_y],
                                                                   eigenvectors[:, band2, idx_x + 1, idx_y],
                                                                   eigenvectors[:, band2, idx_x, idx_y + 1],
                                                                   eigenvectors[:, band2, idx_x + 1, idx_y + 1],
                                                                   eigenvectors[:, band3, idx_x, idx_y],
                                                                   eigenvectors[:, band3, idx_x + 1, idx_y],
                                                                   eigenvectors[:, band3, idx_x, idx_y + 1],
                                                                   eigenvectors[:, band3, idx_x + 1, idx_y + 1])
            berry_flux_matrix[band2, idx_x, idx_y] = berry_flux_matrix[band1, idx_x, idx_y]
            berry_flux_matrix[band3, idx_x, idx_y] = berry_flux_matrix[band1, idx_x, idx_y]

    band1 = 3
    band2 = 4
    band3 = 5

    for idx_x in range(numb_samples - 1):
        for idx_y in range(numb_samples - 1):
            berry_flux_matrix[band1, idx_x, idx_y] = triple_berry_curv(eigenvectors[:, band1, idx_x, idx_y],
                                                                       eigenvectors[:, band1, idx_x + 1, idx_y],
                                                                       eigenvectors[:, band1, idx_x, idx_y + 1],
                                                                       eigenvectors[:, band1, idx_x + 1, idx_y + 1],
                                                                       eigenvectors[:, band2, idx_x, idx_y],
                                                                       eigenvectors[:, band2, idx_x + 1, idx_y],
                                                                       eigenvectors[:, band2, idx_x, idx_y + 1],
                                                                       eigenvectors[:, band2, idx_x + 1, idx_y + 1],
                                                                       eigenvectors[:, band3, idx_x, idx_y],
                                                                       eigenvectors[:, band3, idx_x + 1, idx_y],
                                                                       eigenvectors[:, band3, idx_x, idx_y + 1],
                                                                       eigenvectors[:, band3, idx_x + 1, idx_y + 1])
            berry_flux_matrix[band2, idx_x, idx_y] = berry_flux_matrix[band1, idx_x, idx_y]
            berry_flux_matrix[band3, idx_x, idx_y] = berry_flux_matrix[band1, idx_x, idx_y]

    band1 = 6
    band2 = 7
    band3 = 8

    for idx_x in range(numb_samples - 1):
        for idx_y in range(numb_samples - 1):
            berry_flux_matrix[band1, idx_x, idx_y] = triple_berry_curv(eigenvectors[:, band1, idx_x, idx_y],
                                                                       eigenvectors[:, band1, idx_x + 1, idx_y],
                                                                       eigenvectors[:, band1, idx_x, idx_y + 1],
                                                                       eigenvectors[:, band1, idx_x + 1, idx_y + 1],
                                                                       eigenvectors[:, band2, idx_x, idx_y],
                                                                       eigenvectors[:, band2, idx_x + 1, idx_y],
                                                                       eigenvectors[:, band2, idx_x, idx_y + 1],
                                                                       eigenvectors[:, band2, idx_x + 1, idx_y + 1],
                                                                       eigenvectors[:, band3, idx_x, idx_y],
                                                                       eigenvectors[:, band3, idx_x + 1, idx_y],
                                                                       eigenvectors[:, band3, idx_x, idx_y + 1],
                                                                       eigenvectors[:, band3, idx_x + 1, idx_y + 1])
            berry_flux_matrix[band2, idx_x, idx_y] = berry_flux_matrix[band1, idx_x, idx_y]
            berry_flux_matrix[band3, idx_x, idx_y] = berry_flux_matrix[band1, idx_x, idx_y]

    chern_numbers = np.zeros(M)

    for band in range(M):
        chern_numbers[band] = np.sum(berry_flux_matrix[band, :, :]) / (2 * np.pi)
        print("Chern number ( band", band, ") = ", chern_numbers[band])

    ##########
    # Figure #
    ##########

    fig = plt.figure()
    ax = Axes3D(fig)

    idx_x = np.linspace(0, numb_samples-1, numb_samples, dtype=int)
    idx_y = np.linspace(0, numb_samples-1, numb_samples, dtype=int)

    KX, KY = np.meshgrid(idx_x, idx_y)

    E = np.zeros(M, dtype=object)

    for i in range(M):
        E[i] = eigenvalues[i, KX, KY]
        ax.plot_surface(KX, KY, E[i])

    ax.set_aspect('equal', adjustable='box')

    ax.set_xlabel('$k_x / |\mathbf{B}_1|$', fontsize=14)
    ax.set_ylabel('$k_y / |\mathbf{B}_2|$', fontsize=14)
    ax.set_zlabel('$E$ / meV', fontsize=14)

    def custom(value, tick_number):

        if value == 0:
            return "0"
        elif value == 100:
            return "1"
        else:
            return "{}".format(value / 100)

    ax.xaxis.set_major_formatter(plt.FuncFormatter(custom))
    ax.yaxis.set_major_formatter(plt.FuncFormatter(custom))

    # for i in range(M):
    #     start = -0.05
    #     stop = 0.03
    #     interval = (stop-start)/(M-1)
    #     ax.text2D(-0.105, start+i*interval, "$C_{{{:2d}}}={:2d}$".format(i+1, int(round(chern_numbers[i]))), color="k", fontsize=14)

    for i in [0, 3, 6]:  # range(M)
        start = -0.05
        stop = 0.03
        interval = (stop-start)/(3-1)
        ax.text2D(-0.12, start+(i/3)*interval, "$C_{{{:2d}+{:2d}+{:2d}}}={:2d}$".format(i+1, i+2, i+3, int(round(chern_numbers[i]))), color="k", fontsize=14)

    ax.text2D(-0.107, 0.07, "(b) $\phi={}/{}$".format(p, q), color="k", fontsize=18)

    ax.tick_params(axis="x", labelsize=14)
    ax.tick_params(axis="y", labelsize=14)
    ax.tick_params(axis="z", labelsize=14)

    # plt.rc('xtick', labelsize=20)
    # plt.rc('ytick', labelsize=20)

    plt.show()
