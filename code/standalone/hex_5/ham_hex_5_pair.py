"""Based on: Energy Spectrum of a Honeycomb Lattice under Nonuniform Magnetic Fields"""

from mpl_toolkits.mplot3d import axes3d, Axes3D
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec


def hamiltonian(k, M, p, q):

    a = 1
    c = np.sqrt(3) * a / 6
    eta = 1 * k[0] * M * a / 2
    alpha = float(p / q)

    matrix = np.zeros(shape=(M, M), dtype=np.complex128)

    def A(phi, m):
        return 2 * np.cos(4 * np.pi * phi * m + 12 * k[1] * c) \
               + 2 * np.cos(2 * np.pi * phi * (m - 3 / 2) + 6 * k[1] * c) \
               + 2 * np.cos(2 * np.pi * phi * (m + 3 / 2) + 6 * k[1] * c) + 3

    def D(phi, m):
        return 2 * np.cos(np.pi * phi * (m - 5 / 2) + 3 * k[1] * c) \
               + 2 * np.cos(np.pi * phi * (m + 7 / 2) + 3 * k[1] * c) \
               + 2 * np.cos(3 * np.pi * phi * (m - 1 / 2) + 9 * k[1] * c) \
               + 2 * np.cos(3 * np.pi * phi * (m + 3 / 2) + 9 * k[1] * c)

    def G(phi, m):
        return 2 * np.cos(2 * np.pi * phi * (m + 1) + 6 * k[1] * c) \
               + 2 * np.cos(3 * np.pi * phi)

    for i in range(M):
        matrix[i, i] = A(alpha, i + 1)
    for i in range(M - 3):
        matrix[i, i + 3] = D(alpha, i + 2)
        matrix[i + 3, i] = D(alpha, i + 2)
    for i in range(M - 6):
        matrix[i, i + 6] = G(alpha, i + 3)
        matrix[i + 6, i] = G(alpha, i + 3)

    # top-right
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


def double_berry_curv(ev1, ev1_alpha, ev1_beta, ev1_alpha_beta, ev2, ev2_alpha, ev2_beta, ev2_alpha_beta):

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

    double_bc = - np.imag(np.log(np.linalg.det(matrix1) * np.linalg.det(matrix2) * np.linalg.det(matrix3) * np.linalg.det(matrix4)))

    return double_bc


numb_samples = 101

if __name__ == '__main__':

    print("***Chern numbers***\n")

    ############
    # System 1 #########################################################################################################
    ############

    p1 = 1
    q1 = 4

    # reciprocal lattice vectors
    b1 = (2. * np.pi / q1) * np.array([1, -1 / np.sqrt(3)])
    b2 = (2. * np.pi) * np.array([0, 2 / np.sqrt(3)])
    bvec = np.vstack((b1, b2))

    if p1 % 2 == 0:
        M1 = q1
    else:
        M1 = 2 * q1

    eigenvalues1 = np.zeros((2*M1, numb_samples, numb_samples))
    eigenvectors1 = np.zeros((M1, M1, numb_samples, numb_samples), dtype=np.complex128)

    for band in range(M1):
        for idx_x in range(numb_samples):
            frac_kx = idx_x / (numb_samples-1)
            for idx_y in range(numb_samples):
                frac_ky = idx_y / (numb_samples-1)

                k_u = np.matmul(np.array([frac_kx, frac_ky]), bvec)

                eigvals, eigvecs = np.linalg.eig(hamiltonian(k_u, M1, p1, q1))
                idx = np.argsort(eigvals)
                eigenvalues1[band][idx_x][idx_y] = np.real(+np.sqrt(3 + eigvals[idx[band]]))
                eigenvalues1[M1 + band][idx_x][idx_y] = np.real(-np.sqrt(3 + eigvals[idx[band]]))
                eigenvectors1[:, band, idx_x, idx_y] = eigvecs[:, idx[band]]

    berry_flux_matrix1 = np.zeros((M1, numb_samples - 1, numb_samples - 1))

    band1 = 0
    band2 = 1

    for idx_x in range(numb_samples - 1):
        for idx_y in range(numb_samples - 1):
            berry_flux_matrix1[band1, idx_x, idx_y] = double_berry_curv(eigenvectors1[:, band1, idx_x, idx_y],
                                                                        eigenvectors1[:, band1, idx_x + 1, idx_y],
                                                                        eigenvectors1[:, band1, idx_x, idx_y + 1],
                                                                        eigenvectors1[:, band1, idx_x + 1, idx_y + 1],
                                                                        eigenvectors1[:, band2, idx_x, idx_y],
                                                                        eigenvectors1[:, band2, idx_x + 1, idx_y],
                                                                        eigenvectors1[:, band2, idx_x, idx_y + 1],
                                                                        eigenvectors1[:, band2, idx_x + 1, idx_y + 1])
            berry_flux_matrix1[band2, idx_x, idx_y] = berry_flux_matrix1[band1, idx_x, idx_y]

    band1 = 2
    band2 = 3

    for idx_x in range(numb_samples - 1):
        for idx_y in range(numb_samples - 1):
            berry_flux_matrix1[band1, idx_x, idx_y] = double_berry_curv(eigenvectors1[:, band1, idx_x, idx_y],
                                                                        eigenvectors1[:, band1, idx_x + 1, idx_y],
                                                                        eigenvectors1[:, band1, idx_x, idx_y + 1],
                                                                        eigenvectors1[:, band1, idx_x + 1, idx_y + 1],
                                                                        eigenvectors1[:, band2, idx_x, idx_y],
                                                                        eigenvectors1[:, band2, idx_x + 1, idx_y],
                                                                        eigenvectors1[:, band2, idx_x, idx_y + 1],
                                                                        eigenvectors1[:, band2, idx_x + 1, idx_y + 1])
            berry_flux_matrix1[band2, idx_x, idx_y] = berry_flux_matrix1[band1, idx_x, idx_y]

    band1 = 4
    band2 = 5

    for idx_x in range(numb_samples - 1):
        for idx_y in range(numb_samples - 1):
            berry_flux_matrix1[band1, idx_x, idx_y] = double_berry_curv(eigenvectors1[:, band1, idx_x, idx_y],
                                                                        eigenvectors1[:, band1, idx_x + 1, idx_y],
                                                                        eigenvectors1[:, band1, idx_x, idx_y + 1],
                                                                        eigenvectors1[:, band1, idx_x + 1, idx_y + 1],
                                                                        eigenvectors1[:, band2, idx_x, idx_y],
                                                                        eigenvectors1[:, band2, idx_x + 1, idx_y],
                                                                        eigenvectors1[:, band2, idx_x, idx_y + 1],
                                                                        eigenvectors1[:, band2, idx_x + 1, idx_y + 1])
            berry_flux_matrix1[band2, idx_x, idx_y] = berry_flux_matrix1[band1, idx_x, idx_y]

    band1 = 6
    band2 = 7

    for idx_x in range(numb_samples - 1):
        for idx_y in range(numb_samples - 1):
            berry_flux_matrix1[band1, idx_x, idx_y] = double_berry_curv(eigenvectors1[:, band1, idx_x, idx_y],
                                                                        eigenvectors1[:, band1, idx_x + 1, idx_y],
                                                                        eigenvectors1[:, band1, idx_x, idx_y + 1],
                                                                        eigenvectors1[:, band1, idx_x + 1, idx_y + 1],
                                                                        eigenvectors1[:, band2, idx_x, idx_y],
                                                                        eigenvectors1[:, band2, idx_x + 1, idx_y],
                                                                        eigenvectors1[:, band2, idx_x, idx_y + 1],
                                                                        eigenvectors1[:, band2, idx_x + 1, idx_y + 1])
            berry_flux_matrix1[band2, idx_x, idx_y] = berry_flux_matrix1[band1, idx_x, idx_y]

    chern_numbers1 = np.zeros(M1)

    for band in range(M1):
        chern_numbers1[band] = np.sum(berry_flux_matrix1[band, :, :]) / (2 * np.pi)
        print("I) Chern number ( band", band, ") = ", chern_numbers1[band])

    ########################################################################
    # Berry fluxes along Wilson loops in the (non-transformed) y direction #############################################
    ########################################################################

    wl_berry_flux = np.ones((8, numb_samples), dtype=np.complex128)
    hwcc = np.zeros((8, numb_samples))

    for band in [0, 1]:
        for idx_x in range(numb_samples):
            for idx_y in range(numb_samples - 1):

                # boundary conditions
                if idx_y == (numb_samples - 2) and idx_x == (numb_samples - 1):
                    idx_y_f, idx_x = 0, 0
                elif idx_y == (numb_samples - 2):
                    idx_y_f = 0
                elif idx_x == (numb_samples - 1):
                    idx_x = 0
                else:
                    idx_y_f = idx_y + 1

                wl_berry_flux[band, idx_x] *= \
                    np.conj(eigenvectors1[:, band, idx_x, idx_y]).dot(eigenvectors1[:, band, idx_x, idx_y_f])
                hwcc[band, idx_x] = -(1 / (2 * np.pi)) * np.imag(np.log(wl_berry_flux[band, idx_x]))

    for band in [2, 3]:
        for idx_x in range(numb_samples):
            for idx_y in range(numb_samples - 1):

                # boundary conditions
                if idx_y == (numb_samples - 2) and idx_x == (numb_samples - 1):
                    idx_y_f, idx_x = 0, 0
                elif idx_y == (numb_samples - 2):
                    idx_y_f = 0
                elif idx_x == (numb_samples - 1):
                    idx_x = 0
                else:
                    idx_y_f = idx_y + 1

                band1 = 2
                band2 = 3

                matrix = np.zeros((2, 2), dtype=np.complex128)
                matrix[0, 0] = np.conj(eigenvectors1[:, band1, idx_x, idx_y]).dot(
                    eigenvectors1[:, band1, idx_x, idx_y_f])
                matrix[0, 1] = -np.conj(eigenvectors1[:, band1, idx_x, idx_y]).dot(
                    eigenvectors1[:, band2, idx_x, idx_y_f])
                matrix[1, 0] = -np.conj(eigenvectors1[:, band2, idx_x, idx_y]).dot(
                    eigenvectors1[:, band1, idx_x, idx_y_f])
                matrix[1, 1] = np.conj(eigenvectors1[:, band2, idx_x, idx_y]).dot(
                    eigenvectors1[:, band2, idx_x, idx_y_f])

                wl_berry_flux[band, idx_x] *= np.linalg.det(matrix)
                hwcc[band, idx_x] = -(1 / (2 * np.pi)) * np.imag(np.log(wl_berry_flux[band, idx_y]))

    for band in [4, 5]:
        for idx_x in range(numb_samples):
            for idx_y in range(numb_samples - 1):

                # boundary conditions
                if idx_y == (numb_samples - 2) and idx_x == (numb_samples - 1):
                    idx_y_f, idx_x = 0, 0
                elif idx_y == (numb_samples - 2):
                    idx_y_f = 0
                elif idx_x == (numb_samples - 1):
                    idx_x = 0
                else:
                    idx_y_f = idx_y + 1

                band1 = 4
                band2 = 5

                matrix = np.zeros((2, 2), dtype=np.complex128)
                matrix[0, 0] = np.conj(eigenvectors1[:, band1, idx_x, idx_y]).dot(
                    eigenvectors1[:, band1, idx_x, idx_y_f])
                matrix[0, 1] = -np.conj(eigenvectors1[:, band1, idx_x, idx_y]).dot(
                    eigenvectors1[:, band2, idx_x, idx_y_f])
                matrix[1, 0] = -np.conj(eigenvectors1[:, band2, idx_x, idx_y]).dot(
                    eigenvectors1[:, band1, idx_x, idx_y_f])
                matrix[1, 1] = np.conj(eigenvectors1[:, band2, idx_x, idx_y]).dot(
                    eigenvectors1[:, band2, idx_x, idx_y_f])

                wl_berry_flux[band, idx_x] *= np.linalg.det(matrix)
                hwcc[band, idx_x] = -(1 / (2 * np.pi)) * np.imag(np.log(wl_berry_flux[band, idx_y]))

    for band in [6, 7]:
        for idx_x in range(numb_samples):
            for idx_y in range(numb_samples - 1):

                # boundary conditions
                if idx_y == (numb_samples - 2) and idx_x == (numb_samples - 1):
                    idx_y_f, idx_x = 0, 0
                elif idx_y == (numb_samples - 2):
                    idx_y_f = 0
                elif idx_x == (numb_samples - 1):
                    idx_x = 0
                else:
                    idx_y_f = idx_y + 1

                band1 = 6
                band2 = 7

                matrix = np.zeros((2, 2), dtype=np.complex128)
                matrix[0, 0] = np.conj(eigenvectors1[:, band1, idx_x, idx_y]).dot(
                    eigenvectors1[:, band1, idx_x, idx_y_f])
                matrix[0, 1] = -np.conj(eigenvectors1[:, band1, idx_x, idx_y]).dot(
                    eigenvectors1[:, band2, idx_x, idx_y_f])
                matrix[1, 0] = -np.conj(eigenvectors1[:, band2, idx_x, idx_y]).dot(
                    eigenvectors1[:, band1, idx_x, idx_y_f])
                matrix[1, 1] = np.conj(eigenvectors1[:, band2, idx_x, idx_y]).dot(
                    eigenvectors1[:, band2, idx_x, idx_y_f])

                wl_berry_flux[band, idx_x] *= np.linalg.det(matrix)
                hwcc[band, idx_x] = -(1 / (2 * np.pi)) * np.imag(np.log(wl_berry_flux[band, idx_y]))

    ##########
    # Figure ###########################################################################################################
    ##########

    fig = plt.figure(figsize=(12, 6))
    gs = gridspec.GridSpec(1, 2, width_ratios=[1, 1])

    ###############
    # Subfigure 1 ######################################################################################################
    ###############

    idx_x = np.linspace(0, numb_samples - 1, numb_samples, dtype=int)
    idx_y = np.linspace(0, numb_samples - 1, numb_samples, dtype=int)

    KX, KY = np.meshgrid(idx_x, idx_y)

    ax = plt.subplot(gs[0], projection='3d')

    E1 = np.zeros(2*M1, dtype=object)

    # # color index
    # colors = plt.cm.RdBu_r(np.linspace(0, 1, 2 * M1))
    # cidx = np.zeros(2 * M1, dtype=int)
    # for i in range(M1):
    #     cidx[i] = M1 + i
    #     cidx[M1 + i] = (M1 - 1) - i

    for i in range(2*M1):
        E1[i] = eigenvalues1[i, KX, KY]
        ax.plot_surface(KX, KY, E1[i])

    ax.set_aspect('equal', adjustable='box')

    ax.tick_params(axis='x', which='major', pad=0.5)

    ax.set_xlabel('\n$k_x / |\mathbf{B}_1|$', fontsize=14, linespacing=0.3)
    ax.set_ylabel('\n$k_y / |\mathbf{B}_2|$', fontsize=14, linespacing=1.5)
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

    for i in [0, 2, 4, 6]:
        start = -0.05
        stop = 0.03
        interval = (stop - start) / (3 - 1)
        if i == 0:
            ax.text2D(-0.12, start+i*interval, "$C_{{\pm{:2d}, \pm{:2d}}}={:2d}$".format(i+1, i+2, 0*int(round(chern_numbers1[i]))),
                      color='C{}'.format(i), fontsize=14)
        else:
            ax.text2D(-0.12, start + (i / 3) * interval,
                      "$C_{{{:2d},{:2d}}}={:2d}$".format(i + 1, i + 2, int(round(chern_numbers1[i]))), color='C{}'.format(i),
                      fontsize=14)

    ax.text2D(-0.12, 0.09, "(a)                   $\phi={}/{}$".format(p1, q1), color="k", fontsize=18)
    ax.text2D(0.15, 0.09, "(b)", color="k", fontsize=18)

    ax.tick_params(axis="x", labelsize=14)
    ax.tick_params(axis="y", labelsize=14)
    ax.tick_params(axis="z", labelsize=14)

    ###############
    # Subfigure 2 ######################################################################################################
    ###############

    ax1 = plt.subplot(gs[1])

    ax1.set_aspect('equal', adjustable='box')

    ax1.set_xlim([0, 1])
    ax1.set_ylim([-0.5, 0.5])
    ax1.axhline(0, color='k', linewidth=0.25, ls='--')
    ax1.axvline(0.5, color='k', linewidth=0.25, ls='--')
    ax1.set_xlabel('$k_x / |\mathbf{B}_1|$', fontsize=14)
    ax1.set_ylabel('$\Sigma \mathrm{HWCC} / 2 \pi$', fontsize=14)

    ax1.tick_params(axis="x", labelsize=14)
    ax1.tick_params(axis="y", labelsize=14)

    for band in [0, 2, 4, 6]:
        ax1.scatter(np.arange(numb_samples) / 100, hwcc[band, :], s=5)

    gs.update(wspace=0.5)

    print("\n***Band touching***\n")

    ########################
    # Band touching report #############################################################################################
    ########################

    threshold = 1e-5

    ######

    maxima1 = np.zeros(2 * M1)
    minima1 = np.zeros(2 * M1)
    for i in range(2 * M1):
        maxima1[i] = np.max(E1[i])
        minima1[i] = np.min(E1[i])

    # band index
    bidx1 = np.zeros(2 * M1, dtype=int)
    for i in range(M1):
        bidx1[i] = (2 * M1 - 1) - i
        bidx1[M1 + i] = i

    gap1 = np.zeros(2 * M1 - 1)
    for i in range(2 * M1 - 1):
        gap1[i] = minima1[bidx1[i + 1]] - maxima1[bidx1[i]]

    for i in range(len(gap1)):
        if gap1[i] < threshold:
            print("I) Bands", i + 1, "and", i + 2, "have overlapping extrema within a threshold of", threshold)
    print(gap1)

    ######

    for i in range(2 * M1 - 1):
        for kx in idx_x:
            for ky in idx_y:
                if abs(eigenvalues1[bidx1[i], kx, ky] - eigenvalues1[bidx1[i + 1], kx, ky]) < threshold:
                    print("I) Bands", i + 1, "and", i + 2, "touch within a threshold of", threshold)
                    break
            if abs(eigenvalues1[bidx1[i], kx, ky] - eigenvalues1[bidx1[i + 1], kx, ky]) < threshold:
                break

    ############

    plt.savefig("/home/bart/Documents/papers/TBG/figures/hex_5_bands_pair.png", bbox_inches='tight', dpi=300)
    plt.show()
