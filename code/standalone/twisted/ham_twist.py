from mpl_toolkits.mplot3d import axes3d, Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np


p = 2
q = 7

# reciprocal lattice vectors

b1 = (2.*np.pi/q) * np.array([1, -1/np.sqrt(3)])
b2 = (2.*np.pi) * np.array([0, 2/np.sqrt(3)])
bvec = np.vstack((b1, b2))


def hamiltonian(k, M):

    # t1, t2 = 1, 1
    # a = 1
    # c = np.sqrt(3) * a / 6
    # delta = k[0] * M * a / 2
    #
    # matrix = np.zeros(shape=(M, M), dtype=np.complex128)
    #
    # alpha = float(p / q)
    #
    # def A(phi, m):
    #     return 2 * t1 * np.cos(6 * k[1] * c) + 2 * t2 * np.cos(2 * np.pi * phi * m + 18 * k[1] * c)  # 18
    #
    # def B():
    #     return 2 * t1 * np.cos(3 * k[1] * c)
    #
    # def C(phi, m):
    #     return 2 * t2 * np.cos(np.pi * phi * (m + 1 / 2) + 9 * k[1] * c)  # 9
    #
    # for i in range(M):
    #     matrix[i][i] = A(alpha, i + 1)
    # for i in range(M - 1):
    #     matrix[i, i + 1] = B()
    #     matrix[i + 1, i] = B()
    # for i in range(M - 3):
    #     matrix[i, i + 3] = C(alpha, i + 2)
    #     matrix[i + 3, i] = C(alpha, i + 2)
    #
    # # top-right
    # matrix[0][M - 3] = C(alpha, M - 1) * np.exp(-1j * delta)
    # matrix[1][M - 2] = C(alpha, M) * np.exp(-1j * delta)
    # matrix[2][M - 1] = C(alpha, 1) * np.exp(-1j * delta)
    # matrix[0][M - 2] = 0
    # matrix[1][M - 1] = 0
    # matrix[0][M - 1] = B() * np.exp(-1j * delta)
    #
    # # bottom-left
    # matrix[M - 3][0] = C(alpha, M - 1) * np.exp(1j * delta)
    # matrix[M - 2][1] = C(alpha, M) * np.exp(1j * delta)
    # matrix[M - 1][2] = C(alpha, 1) * np.exp(1j * delta)
    # matrix[M - 2][0] = 0
    # matrix[M - 1][1] = 0
    # matrix[M - 1][0] = B() * np.exp(1j * delta)

    t = 1
    tdash = 0

    matrix = np.zeros(shape=(M, M), dtype=np.complex128)

    alpha = float(p / q)

    def A(k, phi, m):
        return 2 * tdash * np.cos(2 * np.pi * phi * m + 2 * k[1] * (np.sqrt(3)/2))

    def B(k, phi, m):
        return 2 * t * np.cos(np.pi * phi * (m + 1 / 2) + k[1] * (np.sqrt(3)/2))

    def C(k, phi, m):
        return 2 * tdash * np.cos(np.pi * phi * (m + 1 / 2) + k[1] * (np.sqrt(3)/2))

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

    for band in range(M):
        for idx_x in range(numb_samples-1):
            for idx_y in range(numb_samples-1):
                berry_flux_matrix[band, idx_x, idx_y] = berry_curv(eigenvectors[:, band, idx_x, idx_y],
                                                                   eigenvectors[:, band, idx_x + 1, idx_y],
                                                                   eigenvectors[:, band, idx_x, idx_y + 1],
                                                                   eigenvectors[:, band, idx_x + 1, idx_y + 1])

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

    for i in range(M):
        start = -0.05
        stop = 0.03
        interval = (stop-start)/(M-1)
        ax.text2D(-0.105, start+i*interval, "$C_{:2d}={:2d}$".format(i+1, int(round(chern_numbers[i]))), color="k", fontsize=14)

    ax.text2D(-0.107, 0.07, "(a) $\phi={}/{}$".format(p, q), color="k", fontsize=18)

    ax.tick_params(axis="x", labelsize=14)
    ax.tick_params(axis="y", labelsize=14)
    ax.tick_params(axis="z", labelsize=14)

    # plt.rc('xtick', labelsize=20)
    # plt.rc('ytick', labelsize=20)

    plt.show()
