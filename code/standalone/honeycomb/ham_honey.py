"""Based on: Energy Spectrum of a Honeycomb Lattice under Nonuniform Magnetic Fields"""

from mpl_toolkits.mplot3d import axes3d, Axes3D
import matplotlib.pyplot as plt
import numpy as np


p = 1
q = 2


def hamiltonian(k, M):

    t = 1
    tdash = 0

    matrix = np.zeros(shape=(M, M), dtype=np.complex128)

    alpha = float(p / q)

    def A(k, phi, m):
        return 2 * tdash * np.cos(2 * np.pi * phi * m - 2 * k[1] * (np.sqrt(3)/2))

    def B(k, phi, m):
        return 2 * t * np.cos(np.pi * phi * (m + 1 / 2) - k[1] * (np.sqrt(3)/2))

    def C(k, phi, m):
        return 2 * tdash * np.cos(np.pi * phi * (m + 1 / 2) - k[1] * (np.sqrt(3)/2))

    def delta(k, M):
        return k[0] * M / 2

    for i in range(M):
        matrix[i][i] = A(k, alpha, i + 1)
    for i in range(M - 1):
        matrix[i, i + 1] = B(k, alpha, i + 1)
        matrix[i + 1, i] = B(k, alpha, i + 1)
    for i in range(M - 2):
        matrix[i, i + 2] = t
        matrix[i + 2, i] = t
    for i in range(M - 3):
        matrix[i, i + 3] = C(k, alpha, i + 2)
        matrix[i + 3, i] = C(k, alpha, i + 2)

    # top-right
    matrix[0][M - 3] = C(k, alpha, M - 1) * np.exp(-1j * delta(k, M))
    matrix[1][M - 2] = C(k, alpha, M) * np.exp(-1j * delta(k, M))
    matrix[2][M - 1] = C(k, alpha, 1) * np.exp(-1j * delta(k, M))
    matrix[0][M - 2] = t * np.exp(-1j * delta(k, M))
    matrix[1][M - 1] = t * np.exp(-1j * delta(k, M))
    matrix[0][M - 1] = B(k, alpha, M) * np.exp(-1j * delta(k, M))

    # bottom-left
    matrix[M - 3][0] = C(k, alpha, M - 1) * np.exp(1j * delta(k, M))
    matrix[M - 2][1] = C(k, alpha, M) * np.exp(1j * delta(k, M))
    matrix[M - 1][2] = C(k, alpha, 1) * np.exp(1j * delta(k, M))
    matrix[M - 2][0] = t * np.exp(1j * delta(k, M))
    matrix[M - 1][1] = t * np.exp(1j * delta(k, M))
    matrix[M - 1][0] = B(k, alpha, M) * np.exp(1j * delta(k, M))

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

    numb_samples = 100

    if p % 2 == 0:
        M = q
    else:
        M = 2 * q

    eigenvalues = np.zeros((M, numb_samples, numb_samples))
    eigenvectors = np.zeros((M, M, numb_samples, numb_samples), dtype=np.complex128)

    for band in range(M):
        for i in range(numb_samples):
            kx = (np.pi/(M*0.5))*(-1 + i * 2/(numb_samples-1))
            for j in range(numb_samples):
                ky = (np.pi)*(-1 + j * 2/(numb_samples-1))
                eigvals, eigvecs = np.linalg.eig(hamiltonian([kx, ky], M))
                idx = np.argsort(eigvals)
                eigenvalues[band][i][j] = np.real(eigvals[idx[band]])
                eigenvectors[:, band, i, j] = eigvecs[:, idx[band]]

    berry_flux_matrix = np.zeros((M, numb_samples - 1, numb_samples - 1))

    # for band in range(M):
    #     for idx_x in range(numb_samples-1):
    #         for idx_y in range(numb_samples-1):
    #             berry_flux_matrix[band, idx_x, idx_y] = berry_curv(eigenvectors[:, band, idx_x, idx_y],
    #                                                                eigenvectors[:, band, idx_x + 1, idx_y],
    #                                                                eigenvectors[:, band, idx_x, idx_y + 1],
    #                                                                eigenvectors[:, band, idx_x + 1, idx_y + 1])

    band0 = 0
    band1 = 1

    for idx_x in range(numb_samples - 1):
        for idx_y in range(numb_samples - 1):
            berry_flux_matrix[band0, idx_x, idx_y] = multi_berry_curv(eigenvectors[:, band0, idx_x, idx_y],
                                                                  eigenvectors[:, band0, idx_x + 1, idx_y],
                                                                  eigenvectors[:, band0, idx_x, idx_y + 1],
                                                                  eigenvectors[:, band0, idx_x + 1, idx_y + 1],
                                                                  eigenvectors[:, band1, idx_x, idx_y],
                                                                  eigenvectors[:, band1, idx_x + 1, idx_y],
                                                                  eigenvectors[:, band1, idx_x, idx_y + 1],
                                                                  eigenvectors[:, band1, idx_x + 1, idx_y + 1])
            berry_flux_matrix[band1, idx_x, idx_y] = berry_flux_matrix[band0, idx_x, idx_y]

    band2 = 2
    band3 = 3

    for idx_x in range(numb_samples - 1):
        for idx_y in range(numb_samples - 1):
            berry_flux_matrix[band2, idx_x, idx_y] = multi_berry_curv(eigenvectors[:, band2, idx_x, idx_y],
                                                                      eigenvectors[:, band2, idx_x + 1, idx_y],
                                                                      eigenvectors[:, band2, idx_x, idx_y + 1],
                                                                      eigenvectors[:, band2, idx_x + 1, idx_y + 1],
                                                                      eigenvectors[:, band3, idx_x, idx_y],
                                                                      eigenvectors[:, band3, idx_x + 1, idx_y],
                                                                      eigenvectors[:, band3, idx_x, idx_y + 1],
                                                                      eigenvectors[:, band3, idx_x + 1, idx_y + 1])
            berry_flux_matrix[band3, idx_x, idx_y] = berry_flux_matrix[band2, idx_x, idx_y]

    chern_numbers = np.zeros(M)

    for band in range(M):
        chern_numbers[band] = np.sum(berry_flux_matrix[band, :, :]) / (2 * np.pi)
        print("Chern number ( band", band, ") = ", chern_numbers[band])

    ##########
    # Figure #
    ##########

    fig = plt.figure()
    ax = Axes3D(fig)

    kx = np.linspace(0, numb_samples-1, numb_samples, dtype=int)
    ky = np.linspace(0, numb_samples-1, numb_samples, dtype=int)

    KX, KY = np.meshgrid(kx, ky)

    E = np.zeros(M, dtype=object)

    for i in range(M):
        E[i] = eigenvalues[i, KX, KY]
        ax.plot_surface(KX, KY, E[i])

    ax.set_xlabel('kx')
    ax.set_ylabel('ky')
    ax.set_zlabel('E')

    plt.show()
