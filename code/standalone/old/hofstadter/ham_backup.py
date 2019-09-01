from mpl_toolkits.mplot3d import axes3d, Axes3D
import matplotlib.pyplot as plt
import numpy as np

def hamiltonian(k, alpha, q):

    J = 1

    def h(k, alpha, q):
        return 2 * np.cos(k[0] - q * 2 * np.pi * alpha)

    ham = np.zeros((q, q), dtype=np.complex128)

    for i in range(q):
        ham[i][i] = -J * h(k, alpha, i)

    for i in range(q-1):
        ham[i][i+1] = -J * np.exp(-1j * k[1])
        ham[i+1][i] = -J * np.exp(1j * k[1])

    ham[0][q-1] = -J * np.exp(1j * k[1])
    ham[q-1][0] = -J * np.exp(-1j * k[1])

    return ham


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
    q = 6

    eigenvalues = np.zeros((q, numb_samples, numb_samples))
    eigenvectors = np.zeros((q, q, numb_samples, numb_samples), dtype=np.complex128)

    for band in range(q):
        for i in range(numb_samples):
            #kx = np.pi*(-1 + i * 2/(numb_samples-1))
            kx = (2 * np.pi) * (i / (numb_samples - 1))
            for j in range(numb_samples):
                #ky = (np.pi/q)*(-1 + j * 2/(numb_samples-1))
                ky = (2 * np.pi / q) * (j / (numb_samples - 1))
                eigvals, eigvecs = np.linalg.eig(hamiltonian([kx, ky], 1/q, q))
                idx = np.argsort(eigvals)
                eigenvalues[band][i][j] = np.real(eigvals[idx[band]])
                eigenvectors[:, band, i, j] = eigvecs[:, idx[band]]

    berry_flux_matrix = np.zeros((q, numb_samples - 1, numb_samples - 1))

    # for band in range(q):
    #     for idx_x in range(numb_samples-1):
    #         for idx_y in range(numb_samples-1):
    #             berry_flux_matrix[band, idx_x, idx_y] = berry_curv(eigenvectors[:, band, idx_x, idx_y],
    #                                                                eigenvectors[:, band, idx_x + 1, idx_y],
    #                                                                eigenvectors[:, band, idx_x, idx_y + 1],
    #                                                                eigenvectors[:, band, idx_x + 1, idx_y + 1])

    for band in [0, 1, 4, 5]:
        for idx_x in range(numb_samples-1):
            for idx_y in range(numb_samples-1):
                berry_flux_matrix[band, idx_x, idx_y] = berry_curv(eigenvectors[:, band, idx_x, idx_y],
                                                                   eigenvectors[:, band, idx_x + 1, idx_y],
                                                                   eigenvectors[:, band, idx_x, idx_y + 1],
                                                                   eigenvectors[:, band, idx_x + 1, idx_y + 1])
    band1 = 2
    band2 = 3

    for idx_x in range(numb_samples-1):
        for idx_y in range(numb_samples-1):
            berry_flux_matrix[2, idx_x, idx_y] = multi_berry_curv(eigenvectors[:, band1, idx_x, idx_y],
                                                                  eigenvectors[:, band1, idx_x + 1, idx_y],
                                                                  eigenvectors[:, band1, idx_x, idx_y + 1],
                                                                  eigenvectors[:, band1, idx_x + 1, idx_y + 1],
                                                                  eigenvectors[:, band2, idx_x, idx_y],
                                                                  eigenvectors[:, band2, idx_x + 1, idx_y],
                                                                  eigenvectors[:, band2, idx_x, idx_y + 1],
                                                                  eigenvectors[:, band2, idx_x + 1, idx_y + 1])
            berry_flux_matrix[3, idx_x, idx_y] = berry_flux_matrix[2, idx_x, idx_y]

    chern_numbers = np.zeros(q)

    for band in range(q):
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

    E = np.zeros(q, dtype=object)

    for i in range(q):
        E[i] = eigenvalues[i, KX, KY]
        ax.plot_surface(KX, KY, E[i])

    ax.set_xlabel('kx')
    ax.set_ylabel('ky')
    ax.set_zlabel('E')

    plt.show()

    ######################################################
    # Berry fluxes along Wilson loops in the y direction #
    ######################################################

    # wl_berry_flux = np.ones((q, numb_samples), dtype=np.complex128)
    # hwcc = np.zeros((q, numb_samples))
    #
    # for band in [0, 1, 2, 3, 4]:
    #     for idx_x in range(numb_samples):
    #         for idx_y in range(numb_samples - 1):
    #             wl_berry_flux[band, idx_x] *= \
    #                 np.conj(eigenvectors[:, band, idx_x, idx_y]).dot(eigenvectors[:, band, idx_x, idx_y + 1])
    #             hwcc[band, idx_x] = -(1 / (2 * np.pi)) * np.imag(np.log(wl_berry_flux[band, idx_x]))
    #             print(band, hwcc[band, idx_x])
    #
    # fig, ax = plt.subplots()
    #
    # t = np.arange(0, numb_samples, 1)
    # ax.plot(t, hwcc[0, :])
    #
    # plt.show()
