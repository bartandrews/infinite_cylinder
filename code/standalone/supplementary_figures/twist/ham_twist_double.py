from mpl_toolkits.mplot3d import axes3d, Axes3D
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec


def hamiltonian(k, M, p, q):

    t1, t2 = 0.5, 0.5
    a = 1
    c = np.sqrt(3) * a / 6  # ... / 6
    eta = 1 * k[0] * M * a / 2  # 3 * ...
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


numb_samples = 101

if __name__ == '__main__':

    print("***Chern numbers***\n")

    ############
    # System 1 #########################################################################################################
    ############

    p1 = 2
    q1 = 7

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

    for band in range(M1):
        for idx_x in range(numb_samples-1):
            for idx_y in range(numb_samples-1):
                berry_flux_matrix1[band, idx_x, idx_y] = berry_curv(eigenvectors1[:, band, idx_x, idx_y],
                                                                    eigenvectors1[:, band, idx_x + 1, idx_y],
                                                                    eigenvectors1[:, band, idx_x, idx_y + 1],
                                                                    eigenvectors1[:, band, idx_x + 1, idx_y + 1])

    chern_numbers1 = np.zeros(M1)

    for band in range(M1):
        chern_numbers1[band] = np.sum(berry_flux_matrix1[band, :, :]) / (2 * np.pi)
        print("I) Chern number ( band", band, ") = ", chern_numbers1[band])

    ############
    # System 2 #########################################################################################################
    ############

    p2 = 2
    q2 = 9

    # reciprocal lattice vectors
    b1 = (2. * np.pi / q2) * np.array([1, -1 / np.sqrt(3)])
    b2 = (2. * np.pi) * np.array([0, 2 / np.sqrt(3)])
    bvec = np.vstack((b1, b2))

    if p2 % 2 == 0:
        M2 = q2
    else:
        M2 = 2 * q2

    eigenvalues2 = np.zeros((2 * M2, numb_samples, numb_samples))
    eigenvectors2 = np.zeros((M2, M2, numb_samples, numb_samples), dtype=np.complex128)

    for band in range(M2):
        for idx_x in range(numb_samples):
            frac_kx = idx_x / (numb_samples - 1)
            for idx_y in range(numb_samples):
                frac_ky = idx_y / (numb_samples - 1)

                k_u = np.matmul(np.array([frac_kx, frac_ky]), bvec)

                eigvals, eigvecs = np.linalg.eig(hamiltonian(k_u, M2, p2, q2))
                idx = np.argsort(eigvals)
                eigenvalues2[band][idx_x][idx_y] = np.real(+np.sqrt(3 + eigvals[idx[band]]))
                eigenvalues2[M2 + band][idx_x][idx_y] = np.real(-np.sqrt(3 + eigvals[idx[band]]))
                eigenvectors2[:, band, idx_x, idx_y] = eigvecs[:, idx[band]]

    berry_flux_matrix2 = np.zeros((M2, numb_samples - 1, numb_samples - 1))

    for band in range(M2):
        for idx_x in range(numb_samples - 1):
            for idx_y in range(numb_samples - 1):
                berry_flux_matrix2[band, idx_x, idx_y] = berry_curv(eigenvectors2[:, band, idx_x, idx_y],
                                                                    eigenvectors2[:, band, idx_x + 1, idx_y],
                                                                    eigenvectors2[:, band, idx_x, idx_y + 1],
                                                                    eigenvectors2[:, band, idx_x + 1, idx_y + 1])

    chern_numbers2 = np.zeros(M2)

    for band in range(M2):
        chern_numbers2[band] = np.sum(berry_flux_matrix2[band, :, :]) / (2 * np.pi)
        print("II) Chern number ( band", band, ") = ", chern_numbers2[band])

    ##########
    # Figure ###########################################################################################################
    ##########

    fig = plt.figure(figsize=(12, 6))
    gs = gridspec.GridSpec(1, 2, width_ratios=[1, 1])

    idx_x = np.linspace(0, numb_samples - 1, numb_samples, dtype=int)
    idx_y = np.linspace(0, numb_samples - 1, numb_samples, dtype=int)

    KX, KY = np.meshgrid(idx_x, idx_y)

    ###############
    # Subfigure 1 ######################################################################################################
    ###############

    ax = plt.subplot(gs[0], projection='3d')

    E1 = np.zeros(2*M1, dtype=object)

    # color index
    colors = plt.cm.RdBu_r(np.linspace(0, 1, 2 * M1))
    cidx = np.zeros(2 * M1, dtype=int)
    for i in range(M1):
        cidx[i] = M1 + i
        cidx[M1 + i] = (M1 - 1) - i

    for i in range(2*M1):
        E1[i] = eigenvalues1[i, KX, KY]
        ax.plot_surface(KX, KY, E1[i], color=colors[cidx[i]])

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

    for i in range(M1):
        start = -0.05
        stop = 0.03
        interval = (stop - start) / (M1 - 1)
        ax.text2D(-0.105, start + i * interval, "$C_{{{:2d}}}={:2d}$".format(i + 1, int(round(chern_numbers1[i]))),
                  color="k", fontsize=14)

    ax.text2D(-0.107, 0.07, "(a) $\phi={}/{}$".format(p1, q1), color="k", fontsize=18)

    ax.tick_params(axis="x", labelsize=14)
    ax.tick_params(axis="y", labelsize=14)
    ax.tick_params(axis="z", labelsize=14)

    ###############
    # Subfigure 2 ######################################################################################################
    ###############

    ax1 = plt.subplot(gs[1], projection='3d')

    E2 = np.zeros(2 * M2, dtype=object)

    # color index
    colors = plt.cm.RdBu_r(np.linspace(0, 1, 2 * M2))
    cidx = np.zeros(2 * M2, dtype=int)
    for i in range(M2):
        cidx[i] = M2 + i
        cidx[M2 + i] = (M2 - 1) - i

    for i in range(2 * M2):
        E2[i] = eigenvalues2[i, KX, KY]
        ax1.plot_surface(KX, KY, E2[i], color=colors[cidx[i]])

    ax1.set_aspect('equal', adjustable='box')

    ax1.tick_params(axis='x', which='major', pad=0.5)

    ax1.set_xlabel('\n$k_x / |\mathbf{B}_1|$', fontsize=14, linespacing=0.3)
    ax1.set_ylabel('\n$k_y / |\mathbf{B}_2|$', fontsize=14, linespacing=1.5)
    ax1.set_zlabel('$E$ / meV', fontsize=14)

    def custom(value, tick_number):

        if value == 0:
            return "0"
        elif value == 100:
            return "1"
        else:
            return "{}".format(value / 100)

    ax1.xaxis.set_major_formatter(plt.FuncFormatter(custom))
    ax1.yaxis.set_major_formatter(plt.FuncFormatter(custom))

    for i in range(M2):
        start = -0.05
        stop = 0.03
        interval = (stop - start) / (M2 - 1)
        ax1.text2D(-0.105, start + i * interval, "$C_{{{:2d}}}={:2d}$".format(i + 1, int(round(chern_numbers2[i]))),
                  color="k", fontsize=14)

    ax1.text2D(-0.107, 0.07, "(b) $\phi={}/{}$".format(p2, q2), color="k", fontsize=18)

    ax1.tick_params(axis="x", labelsize=14)
    ax1.tick_params(axis="y", labelsize=14)
    ax1.tick_params(axis="z", labelsize=14)

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

    maxima2 = np.zeros(2 * M2)
    minima2 = np.zeros(2 * M2)
    for i in range(2 * M2):
        maxima2[i] = np.max(E2[i])
        minima2[i] = np.min(E2[i])

    # band index
    bidx2 = np.zeros(2 * M2, dtype=int)
    for i in range(M2):
        bidx2[i] = (2 * M2 - 1) - i
        bidx2[M2 + i] = i

    gap2 = np.zeros(2 * M2 - 1)
    for i in range(2 * M2 - 1):
        gap2[i] = minima2[bidx2[i + 1]] - maxima2[bidx2[i]]

    for i in range(len(gap2)):
        if gap2[i] < threshold:
            print("II) Bands", i + 1, "and", i + 2, "have overlapping extrema within a threshold of", threshold)
    print(gap2)

    ######

    for i in range(2 * M2 - 1):
        for kx in idx_x:
            for ky in idx_y:
                if abs(eigenvalues2[bidx2[i], kx, ky] - eigenvalues2[bidx2[i + 1], kx, ky]) < threshold:
                    print("II) Bands", i + 1, "and", i + 2, "touch within a threshold of", threshold)
                    break
            if abs(eigenvalues2[bidx2[i], kx, ky] - eigenvalues2[bidx2[i + 1], kx, ky]) < threshold:
                break

    ############

    plt.savefig("/home/bart/Documents/papers/TBG/figures/twist_bands.png", bbox_inches='tight', dpi=300)
    plt.show()
