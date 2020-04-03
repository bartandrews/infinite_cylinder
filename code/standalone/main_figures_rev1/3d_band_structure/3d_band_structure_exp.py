# --- python imports
from mpl_toolkits.mplot3d import axes3d, Axes3D
import matplotlib.pyplot as plt
import numpy as np
import sys
import matplotlib.gridspec as gridspec


def hamiltonian(k, M, p_val, q_val):

    t1 = 1  # hopping parameter
    a = 1  # from Fig. 1
    c = np.sqrt(3) * a / 6  # from Fig. 1
    eta = 1 * k[0] * M * a / 2  # from hex1 paper
    nphi = float(p_val / q_val)  # flux density

    matrix = np.zeros(shape=(M, M), dtype=np.complex128)  # initialize the matrix

    def B(nphi_val, m):
        return t1 * (2 * np.exp(1j * np.pi * nphi_val / 3) * np.cos(np.pi * nphi_val * (m + 1 / 2) + 3 * k[1] * c))

    def C(nphi_val):
        return t1 * np.exp(1j * np.pi * nphi_val / 3)

    for j in range(M - 1):
        matrix[j, j + 1] = B(nphi, j + 1)
        matrix[j + 1, j] = np.conj(B(nphi, j + 1))
    for j in range(M - 2):
        matrix[j, j + 2] = np.conj(C(nphi))
        matrix[j + 2, j] = C(nphi)

    # top-right
    matrix[0][M - 1] = np.conj(B(nphi, M)) * np.exp(-1j * eta)

    matrix[0][M - 2] = C(nphi) * np.exp(-1j * eta)
    matrix[1][M - 1] = C(nphi) * np.exp(-1j * eta)

    # bottom-left
    matrix[M - 1][0] = B(nphi, M) * np.exp(1j * eta)

    matrix[M - 2][0] = np.conj(C(nphi)) * np.exp(1j * eta)
    matrix[M - 1][1] = np.conj(C(nphi)) * np.exp(1j * eta)

    return matrix


def berry_curv(ev, ev_alpha, ev_beta, ev_alpha_beta):

    bc = - np.imag(np.log(np.conj(ev).dot(ev_alpha) * np.conj(ev_alpha).dot(ev_alpha_beta)
                          * np.conj(ev_alpha_beta).dot(ev_beta) * np.conj(ev_beta).dot(ev)))

    return bc


numb_samples = 101

if __name__ == '__main__':

    p, q = 1, 3

    # reciprocal lattice vectors
    b1 = (2. * np.pi / q) * np.array([1, -1 / np.sqrt(3)])
    b2 = (2. * np.pi) * np.array([0, 2 / np.sqrt(3)])
    bvec = np.vstack((b1, b2))

    if p % 2 == 0:
        M = q
    else:
        M = 2 * q

    eigenvalues = np.zeros((2*M, numb_samples, numb_samples))
    eigenvectors = np.zeros((M, M, numb_samples, numb_samples), dtype=np.complex128)

    for band in range(M):
        for idx_x in range(numb_samples):
            frac_kx = idx_x / (numb_samples-1)
            for idx_y in range(numb_samples):
                frac_ky = idx_y / (numb_samples-1)

                k_u = np.matmul(np.array([frac_kx, frac_ky]), bvec)

                eigvals, eigvecs = np.linalg.eig(hamiltonian(k_u, M, p, q))
                idx = np.argsort(eigvals)
                # print(eigvals[idx[0]], eigvals[idx[1]], eigvals[idx[2]], eigvals[idx[3]], eigvals[idx[4]], eigvals[idx[5]])
                if band % 2 != 0:
                    eigenvalues[(q-1) - int(band/2)][idx_x][idx_y] = np.real(-np.sqrt(3 + eigvals[idx[band]]))
                    eigenvectors[:, (q-1) - int(band/2), idx_x, idx_y] = eigvecs[:, idx[band]]

    print(eigenvalues[0][0][0], eigenvalues[1][0][0], eigenvalues[2][0][0])

    berry_flux_matrix = np.zeros((q, numb_samples - 1, numb_samples - 1))

    for band in range(0, q):
        for idx_x in range(numb_samples-1):
            for idx_y in range(numb_samples-1):
                berry_flux_matrix[band, idx_x, idx_y] = berry_curv(eigenvectors[:, band, idx_x, idx_y],
                                                                   eigenvectors[:, band, idx_x + 1, idx_y],
                                                                   eigenvectors[:, band, idx_x, idx_y + 1],
                                                                   eigenvectors[:, band, idx_x + 1, idx_y + 1])

    chern_numbers = np.zeros(q)

    for band in range(0, q):
        chern_numbers[band] = np.sum(berry_flux_matrix[band, :, :]) / (2 * np.pi)
        print("Chern number ( band", band, ") = ", chern_numbers[band])

    # ########################################################################
    # # Berry fluxes along Wilson loops in the (non-transformed) y direction #############################################
    # ########################################################################
    #
    # wl_berry_flux = np.ones((M, numb_samples), dtype=np.complex128)
    # hwcc = np.zeros((M, numb_samples))
    #
    # band = 0
    #
    # # for band in range(0, M):
    # for idx_x in range(numb_samples):
    #     for idx_y in range(numb_samples - 1):
    #
    #         # boundary conditions
    #         if idx_y == (numb_samples - 2) and idx_x == (numb_samples - 1):
    #             idx_y_f, idx_x = 0, 0
    #         elif idx_y == (numb_samples - 2):
    #             idx_y_f = 0
    #         elif idx_x == (numb_samples - 1):
    #             idx_x = 0
    #         else:
    #             idx_y_f = idx_y + 1
    #
    #         wl_berry_flux[band, idx_x] *= \
    #             np.conj(eigenvectors[:, band, idx_x, idx_y]).dot(eigenvectors[:, band, idx_x, idx_y_f])
    #
    #     hwcc[band, idx_x] = -(1 / (2 * np.pi)) * np.imag(np.log(wl_berry_flux[band, idx_x]))

    # for band in range(M, 2*M):
    #     for idx_x in range(numb_samples):
    #         for idx_y in range(numb_samples - 1):
    #             hwcc[band, idx_x] = -(1 / (2 * np.pi)) * np.imag(np.log(wl_berry_flux[(M - 1) - (band - M), idx_x]))

    ##########
    # Figure ###########################################################################################################
    ##########

    fig = plt.figure()
    # gs = gridspec.GridSpec(1, 2, width_ratios=[1, 1])

    bands_to_study = 3

    ###############
    # Subfigure 1 ######################################################################################################
    ###############

    idx_x = np.linspace(0, numb_samples - 1, numb_samples, dtype=int)
    idx_y = np.linspace(0, numb_samples - 1, numb_samples, dtype=int)

    KX, KY = np.meshgrid(idx_x, idx_y)

    ax = plt.subplot(111, projection='3d')

    E = np.zeros(M, dtype=object)

    # # color index
    # colors = plt.cm.RdBu_r(np.linspace(0, 1, 2 * M))
    # cidx = np.zeros(2 * M, dtype=int)
    # for i in range(M):
    #     cidx[i] = M + i
    #     cidx[M + i] = (M - 1) - i

    for i in range(bands_to_study):
        E[i] = eigenvalues[i, KX, KY]
        ax.plot_surface(KX, KY, E[i])

    # ax.set_aspect('equal', adjustable='box')

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

    for i in range(bands_to_study):
        start = -0.05
        stop = 0.03
        interval = (stop - start) / (bands_to_study - 1)
        ax.text2D(-0.15, start + i * interval, "$C_{{{:2d}}}={:2d}$".format(i, int(round(chern_numbers[i]))),
              color='C{}'.format(i % 10), fontsize=11)
    ax.text2D(-0.16, 0.09, "(c)                 $n_\phi={}/{}$".format(p, q), color="k", fontsize=12)
    ax.text2D(0.16, 0.09, "(d)", color="k", fontsize=12)

    ax.tick_params(axis="x", labelsize=10)
    ax.tick_params(axis="y", labelsize=10)
    ax.tick_params(axis="z", labelsize=10)

    # ###############
    # # Subfigure 2 ######################################################################################################
    # ###############
    #
    # ax1 = plt.subplot(gs[1])
    #
    # ax1.set_aspect('equal', adjustable='box')
    #
    # ax1.set_xticks(np.arange(0, 1.05, 0.2))
    # ax1.set_xlim([0, 1])
    # ax1.set_yticks(np.arange(-0.5, 0.55, 0.2))
    # ax1.set_ylim([-0.5, 0.5])
    # # ax1.axhline(-0.17, color='k', linewidth=0.5, ls='--')
    # # ax1.axvline(0.33, color='k', linewidth=0.5, ls='--')
    # ax1.set_xlabel('$k_x / |\mathbf{B}_1|$', fontsize=11)
    # ax1.set_ylabel('$\prod \\theta_\mathrm{B} / 2 \pi$', fontsize=11)
    #
    # ax1.tick_params(axis="x", labelsize=10)
    # ax1.tick_params(axis="y", labelsize=10)
    #
    # band = 0
    #
    # # for band in range(0, bands_to_study):
    #
    # print(np.arange(numb_samples) / (numb_samples-1), hwcc[band, :])
    # hwcc[band, 0] = None
    # hwcc[band, numb_samples-1] = None
    # ax1.scatter(np.arange(numb_samples) / (numb_samples-1), hwcc[band, :], s=5, zorder=bands_to_study-band, color='C{}'.format(band % 10))
    #
    # # ax1.text(0.55, -0.35, "$C=-1$", fontsize=11)
    #
    # gs.update(wspace=0.75)

    # print("\n***Band touching***\n")

    # ########################
    # # Band touching report #############################################################################################
    # ########################
    #
    # threshold = 1e-5
    #
    # ######
    #
    # maxima = np.zeros(bands_to_study)
    # minima = np.zeros(bands_to_study)
    # width = np.zeros(bands_to_study)
    # for i in range(bands_to_study):
    #     maxima[i] = np.max(E[i])
    #     minima[i] = np.min(E[i])
    #     width[i] = maxima[i] - minima[i]
    #
    # print("widths = ", width)
    #
    # # # band index
    # # bidx1 = np.zeros(M, dtype=int)
    # # for i in range(M):
    # #     bidx1[i] = (2 * M - 1) - i
    # #     bidx1[M + i] = i
    #
    # gap = np.zeros(bands_to_study - 1)
    # for i in range(bands_to_study - 1):
    #     gap[i] = minima[i + 1] - maxima[i]
    #
    # for i in range(len(gap)):
    #     if gap[i] < threshold:
    #         print("Bands", i + 1, "and", i + 2, "have overlapping extrema within a threshold of", threshold)
    #
    # print("gaps = ", gap)
    #
    # ######
    #
    # for i in range(M - 1):
    #     for kx in idx_x:
    #         for ky in idx_y:
    #             if abs(eigenvalues[i, kx, ky] - eigenvalues[i + 1, kx, ky]) < threshold:
    #                 print("Bands", i + 1, "and", i + 2, "touch within a threshold of", threshold)
    #                 break
    #         if abs(eigenvalues[i, kx, ky] - eigenvalues[i + 1, kx, ky]) < threshold:
    #             break
    #
    # ############

    # plt.savefig("/home/bart/Documents/papers/TBG_rev1/figures/3d_band_structure_phi_{}_{}.png".format(p, q), bbox_inches='tight', dpi=300)
    plt.show()
