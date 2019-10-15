from mpl_toolkits.mplot3d import axes3d, Axes3D
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec

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
    q = 4

    # reciprocal lattice vectors
    b1 = (2. * np.pi) * np.array([1, 0])
    b2 = (2. * np.pi / q) * np.array([0, 1])
    bvec = np.vstack((b1, b2))

    eigenvalues = np.zeros((q, numb_samples, numb_samples))
    eigenvectors = np.zeros((q, q, numb_samples, numb_samples), dtype=np.complex128)

    for band in range(q):
        for idx_x in range(numb_samples):
            frac_kx = idx_x / (numb_samples - 1)
            for idx_y in range(numb_samples):
                frac_ky = idx_y / (numb_samples - 1)

                # # wavevector used for u (depends on boundary conditions)
                # if idx_x == (numb_samples - 1) and idx_y == (numb_samples - 1):
                #     k_u = np.matmul(np.array([0, 0]), bvec)
                # elif idx_x == (numb_samples - 1):
                #     k_u = np.matmul(np.array([0, frac_ky]), bvec)
                # elif idx_y == (numb_samples - 1):
                #     k_u = np.matmul(np.array([frac_kx, 0]), bvec)
                # else:
                #     k_u = np.matmul(np.array([frac_kx, frac_ky]), bvec)

                k_u = np.matmul(np.array([frac_kx, frac_ky]), bvec)

                eigvals, eigvecs = np.linalg.eig(hamiltonian(k_u, 1/q, q))
                idx = np.argsort(eigvals)
                eigenvalues[band][idx_x][idx_y] = np.real(eigvals[idx[band]])
                eigenvectors[:, band, idx_x, idx_y] = eigvecs[:, idx[band]]

    berry_flux_matrix = np.zeros((q, numb_samples - 1, numb_samples - 1))

    for band in range(q):
        for idx_x in range(numb_samples-1):
            for idx_y in range(numb_samples-1):
                berry_flux_matrix[band, idx_x, idx_y] = berry_curv(eigenvectors[:, band, idx_x, idx_y],
                                                                   eigenvectors[:, band, idx_x + 1, idx_y],
                                                                   eigenvectors[:, band, idx_x, idx_y + 1],
                                                                   eigenvectors[:, band, idx_x + 1, idx_y + 1])

    # for band in [0, 1, 4, 5]:
    #     for idx_x in range(numb_samples-1):
    #         for idx_y in range(numb_samples-1):
    #             berry_flux_matrix[band, idx_x, idx_y] = berry_curv(eigenvectors[:, band, idx_x, idx_y],
    #                                                                eigenvectors[:, band, idx_x + 1, idx_y],
    #                                                                eigenvectors[:, band, idx_x, idx_y + 1],
    #                                                                eigenvectors[:, band, idx_x + 1, idx_y + 1])
    # band1 = 2
    # band2 = 3
    #
    # for idx_x in range(numb_samples-1):
    #     for idx_y in range(numb_samples-1):
    #         berry_flux_matrix[2, idx_x, idx_y] = multi_berry_curv(eigenvectors[:, band1, idx_x, idx_y],
    #                                                               eigenvectors[:, band1, idx_x + 1, idx_y],
    #                                                               eigenvectors[:, band1, idx_x, idx_y + 1],
    #                                                               eigenvectors[:, band1, idx_x + 1, idx_y + 1],
    #                                                               eigenvectors[:, band2, idx_x, idx_y],
    #                                                               eigenvectors[:, band2, idx_x + 1, idx_y],
    #                                                               eigenvectors[:, band2, idx_x, idx_y + 1],
    #                                                               eigenvectors[:, band2, idx_x + 1, idx_y + 1])
    #         berry_flux_matrix[3, idx_x, idx_y] = berry_flux_matrix[2, idx_x, idx_y]

    chern_numbers = np.zeros(q)

    for band in range(q):
        chern_numbers[band] = np.sum(berry_flux_matrix[band, :, :]) / (2 * np.pi)
        print("Chern number ( band", band, ") = ", chern_numbers[band])

    ##########
    # Figure #
    ##########

    fig = plt.figure(figsize=(12, 6))
    gs = gridspec.GridSpec(1, 2, width_ratios=[1, 1])

    ax = plt.subplot(gs[0], projection='3d')

    idx_x = np.linspace(0, numb_samples-1, numb_samples, dtype=int)
    idx_y = np.linspace(0, numb_samples-1, numb_samples, dtype=int)

    KX, KY = np.meshgrid(idx_x, idx_y)

    E = np.zeros(q, dtype=object)

    for i in range(q):
        E[i] = eigenvalues[i, KX, KY]
        ax.plot_surface(KX, KY, E[i])

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

    for i in range(q):
        start = -0.05
        stop = 0.03
        interval = (stop - start) / (q - 1)
        ax.text2D(-0.12, start + i * interval, "$C_{{{:2d}}}={:2d}$".format(i + 1, int(round(chern_numbers[i]))),
                  color="k", fontsize=14)

    ax.text2D(-0.12, 0.09, "(a)                   $\phi=1/{}$".format(q), color="k", fontsize=18)
    ax.text2D(0.15, 0.09, "(b)", color="k", fontsize=18)

    ax.tick_params(axis="x", labelsize=14)
    ax.tick_params(axis="y", labelsize=14)
    ax.tick_params(axis="z", labelsize=14)

    ########################################################################
    # Berry fluxes along Wilson loops in the (non-transformed) x direction #
    ########################################################################

    ax1 = plt.subplot(gs[1])

    wl_berry_flux = np.ones((q, numb_samples), dtype=np.complex128)
    hwcc = np.zeros((q, numb_samples))

    for band in range(q):
        for idx_y in range(numb_samples):
            for idx_x in range(numb_samples - 1):

                # boundary conditions
                if idx_x == (numb_samples - 2) and idx_y == (numb_samples - 1):
                    idx_x_f, idx_y = 0, 0
                elif idx_x == (numb_samples - 2):
                    idx_x_f = 0
                elif idx_y == (numb_samples - 1):
                    idx_y = 0
                else:
                    idx_x_f = idx_x + 1

                wl_berry_flux[band, idx_y] *= \
                    np.conj(eigenvectors[:, band, idx_x, idx_y]).dot(eigenvectors[:, band, idx_x_f, idx_y])
                hwcc[band, idx_y] = -(1 / (2 * np.pi)) * np.imag(np.log(wl_berry_flux[band, idx_y]))
                #print(band, hwcc[band, idx_y])

    ax1.set_aspect('equal', adjustable='box')

    ax1.set_xlim([0, 1])
    ax1.set_ylim([-0.5, 0.5])
    ax1.axhline(0, color='k', linewidth=0.25, ls='--')
    ax1.axvline(0.5, color='k', linewidth=0.25, ls='--')
    ax1.set_xlabel('$k_y / |\mathbf{B}_2|$', fontsize=14)
    ax1.set_ylabel('$\Sigma \mathrm{HWCC} / 2 \pi$', fontsize=14)

    ax1.tick_params(axis="x", labelsize=14)
    ax1.tick_params(axis="y", labelsize=14)

    for band in range(q):
        ax1.scatter(np.arange(numb_samples)/100, hwcc[band, :], s=1)

    gs.update(wspace=0.5)

    # plt.savefig("/home/bart/Documents/papers/TBG/figures/hofstadter_bands_1_3.png", bbox_inches='tight', dpi=300)
    plt.show()
