from mpl_toolkits.mplot3d import axes3d, Axes3D
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec


def hamiltonian(k, M, p, q):

    t1, t2, t2dash = 1, -0.025, 1/18.8
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
                     + 2 * np.cos(3 * np.pi * phi * (m + 3 / 2) + 9 * k[1] * c))\
               + t2dash * (-2 * np.cos(np.pi * phi * (m - 5 / 2) + 3 * k[1] * c)
                     - 2 * np.cos(np.pi * phi * (m + 7 / 2) + 3 * k[1] * c)
                     - 2 * np.cos(3 * np.pi * phi * (m - 1 / 2) + 9 * k[1] * c)
                     - 2 * np.cos(3 * np.pi * phi * (m + 3 / 2) + 9 * k[1] * c))

    def G(phi, m):
        return t2 * (2 * np.cos(2 * np.pi * phi * (m + 1) + 6 * k[1] * c)
                     + 2 * np.cos(3 * np.pi * phi))\
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


numb_samples = 101

if __name__ == '__main__':

    print("***Chern numbers***\n")

    ############
    # System 1 #########################################################################################################
    ############

    p1 = 1
    q1 = 3

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
                # eigenvalues1[band][idx_x][idx_y] = np.real(+np.sqrt(3 + eigvals[idx[band]]))
                eigenvalues1[(M1-1) - band][idx_x][idx_y] = np.real(-np.sqrt(3 + eigvals[idx[band]]))
                eigenvectors1[:, (M1-1) - band, idx_x, idx_y] = eigvecs[:, idx[band]]

    berry_flux_matrix1 = np.zeros((M1, numb_samples - 1, numb_samples - 1))

    for band in range(0, M1, 2):
        for idx_x in range(numb_samples-1):
            for idx_y in range(numb_samples-1):
                berry_flux_matrix1[band, idx_x, idx_y] = multi_berry_curv(eigenvectors1[:, band, idx_x, idx_y],
                                                                          eigenvectors1[:, band, idx_x + 1, idx_y],
                                                                          eigenvectors1[:, band, idx_x, idx_y + 1],
                                                                          eigenvectors1[:, band, idx_x + 1, idx_y + 1],
                                                                          eigenvectors1[:, band+1, idx_x, idx_y],
                                                                          eigenvectors1[:, band+1, idx_x + 1, idx_y],
                                                                          eigenvectors1[:, band+1, idx_x, idx_y + 1],
                                                                          eigenvectors1[:, band+1, idx_x + 1, idx_y + 1])
                berry_flux_matrix1[band+1, idx_x, idx_y] = berry_flux_matrix1[band, idx_x, idx_y]

    chern_numbers1 = np.zeros(M1)

    for band in range(0, M1, 2):
        chern_numbers1[band] = np.sum(berry_flux_matrix1[band, :, :]) / (2 * np.pi)
        print("I) Chern number ( band", band, ") = ", chern_numbers1[band])

    ########################################################################
    # Berry fluxes along Wilson loops in the (non-transformed) y direction #############################################
    ########################################################################

    wl_berry_flux = np.ones((M1, numb_samples), dtype=np.complex128)
    hwcc = np.zeros((M1, numb_samples))

    for band in range(0, M1, 2):
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

                matrix1 = np.zeros((2, 2), dtype=np.complex128)
                matrix1[0][0] = np.conj(eigenvectors1[:, band, idx_x, idx_y]).dot(eigenvectors1[:, band, idx_x, idx_y_f])
                matrix1[0][1] = np.conj(eigenvectors1[:, band, idx_x, idx_y]).dot(eigenvectors1[:, band+1, idx_x, idx_y_f])
                matrix1[1][0] = np.conj(eigenvectors1[:, band+1, idx_x, idx_y]).dot(eigenvectors1[:, band, idx_x, idx_y_f])
                matrix1[1][1] = np.conj(eigenvectors1[:, band+1, idx_x, idx_y]).dot(eigenvectors1[:, band+1, idx_x, idx_y_f])

                wl_berry_flux[band, idx_x] *= \
                    np.linalg.det(matrix1)
                hwcc[band, idx_x] = -(1 / (2 * np.pi)) * np.imag(np.log(wl_berry_flux[band, idx_x]))

    ##########
    # Figure ###########################################################################################################
    ##########

    fig = plt.figure()
    gs = gridspec.GridSpec(1, 2, width_ratios=[1, 1])

    bands_to_study = 6

    ###############
    # Subfigure 1 ######################################################################################################
    ###############

    idx_x = np.linspace(0, numb_samples - 1, numb_samples, dtype=int)
    idx_y = np.linspace(0, numb_samples - 1, numb_samples, dtype=int)

    KX, KY = np.meshgrid(idx_x, idx_y)

    ax = plt.subplot(gs[0], projection='3d')

    E1 = np.zeros(M1, dtype=object)

    # # color index
    # colors = plt.cm.RdBu_r(np.linspace(0, 1, 2 * M1))
    # cidx = np.zeros(2 * M1, dtype=int)
    # for i in range(M1):
    #     cidx[i] = M1 + i
    #     cidx[M1 + i] = (M1 - 1) - i

    for i in range(bands_to_study):
        E1[i] = eigenvalues1[i, KX, KY]
        ax.plot_surface(KX, KY, E1[i])

    ax.set_aspect('equal', adjustable='box')

    ax.tick_params(axis='x', which='major', pad=-4)
    ax.tick_params(axis='y', which='major', pad=-3)
    ax.tick_params(axis='z', which='major', pad=3)

    ax.set_xlabel('$k_x / |\mathbf{B}_1|$', fontsize=11, linespacing=-1)
    ax.set_ylabel('\n$k_y / |\mathbf{B}_2|$', fontsize=11, linespacing=1)
    ax.set_zlabel('\n$E$ / meV', fontsize=11, linespacing=1)

    ax.xaxis.labelpad = -2
    ax.yaxis.labelpad = -2

    def custom(value, tick_number):

        if value == 0:
            return "0"
        elif value == numb_samples-1:
            return "1"
        else:
            return "{}".format(round(value / (numb_samples-1), 2))

    ax.xaxis.set_major_formatter(plt.FuncFormatter(custom))
    ax.yaxis.set_major_formatter(plt.FuncFormatter(custom))

    for i in range(0, bands_to_study, 2):
        start = -0.05
        stop = 0.03
        interval = (stop - start) / (bands_to_study - 1)
        ax.text2D(-0.18, start + i * interval, "$C_{{{:2d},{:2d}}}={:2d}$".format(-M1 + i+1, -M1 + i, int(round(chern_numbers1[i]))),
                  color='C{}'.format(i+1 % 10), fontsize=11)

    ax.text2D(-0.16, 0.09, "(c)                 $n_\phi={}/{}$".format(p1, q1), color="k", fontsize=12)
    ax.text2D(0.16, 0.09, "(d)", color="k", fontsize=12)

    ax.tick_params(axis="x", labelsize=10)
    ax.tick_params(axis="y", labelsize=10)
    ax.tick_params(axis="z", labelsize=10)

    ###############
    # Subfigure 2 ######################################################################################################
    ###############

    ax1 = plt.subplot(gs[1])

    ax1.set_aspect('equal', adjustable='box')

    ax1.set_xticks(np.arange(0, 1.05, 0.2))
    ax1.set_xlim([0, 1])
    ax1.set_yticks(np.arange(-0.5, 0.55, 0.2))
    ax1.set_ylim([-0.5, 0.5])
    # ax1.axhline(-0.17, color='k', linewidth=0.5, ls='--')
    # ax1.axvline(0.33, color='k', linewidth=0.5, ls='--')
    ax1.set_xlabel('$k_x / |\mathbf{B}_1|$', fontsize=11)
    ax1.set_ylabel('$\prod \\theta_\mathrm{B} / 2 \pi$', fontsize=11)

    ax1.tick_params(axis="x", labelsize=10)
    ax1.tick_params(axis="y", labelsize=10)

    for band in range(0, bands_to_study, 2):
        print(np.arange(numb_samples) / (numb_samples-1), hwcc[band, :])
        hwcc[band, 0] = None
        hwcc[band, numb_samples-1] = None
        ax1.scatter(np.arange(numb_samples) / (numb_samples-1), hwcc[band, :], s=5, zorder=bands_to_study-band, color='C{}'.format(band+1 % 10))

    # ax1.text(0.55, -0.35, "$C=-1$", fontsize=11)

    gs.update(wspace=0.75)

    print("\n***Band touching***\n")

    ########################
    # Band touching report #############################################################################################
    ########################

    threshold = 1e-5

    ######

    maxima1 = np.zeros(bands_to_study)
    minima1 = np.zeros(bands_to_study)
    width1 = np.zeros(bands_to_study)
    for i in range(bands_to_study):
        maxima1[i] = np.max(E1[i])
        minima1[i] = np.min(E1[i])
        width1[i] = maxima1[i] - minima1[i]

    print("widths = ", width1)

    # # band index
    # bidx1 = np.zeros(M1, dtype=int)
    # for i in range(M1):
    #     bidx1[i] = (2 * M1 - 1) - i
    #     bidx1[M1 + i] = i

    gap1 = np.zeros(bands_to_study - 1)
    for i in range(bands_to_study - 1):
        gap1[i] = minima1[i + 1] - maxima1[i]

    for i in range(len(gap1)):
        if gap1[i] < threshold:
            print("I) Bands", i + 1, "and", i + 2, "have overlapping extrema within a threshold of", threshold)

    print("gaps = ", gap1)

    ######

    for i in range(M1 - 1):
        for kx in idx_x:
            for ky in idx_y:
                if abs(eigenvalues1[i, kx, ky] - eigenvalues1[i + 1, kx, ky]) < threshold:
                    print("I) Bands", i + 1, "and", i + 2, "touch within a threshold of", threshold)
                    break
            if abs(eigenvalues1[i, kx, ky] - eigenvalues1[i + 1, kx, ky]) < threshold:
                break

    ############

    # plt.savefig("/home/bart/Documents/papers/TBG/figures/3d_band_structure_phi_{}_{}.png".format(p1, q1), bbox_inches='tight', dpi=300)
    plt.show()
