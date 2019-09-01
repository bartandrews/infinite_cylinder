"""Based on: Energy Spectrum of a Honeycomb Lattice under Nonuniform Magnetic Fields"""

from mpl_toolkits.mplot3d import axes3d, Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import sys


p = 2
q = 5

# reciprocal lattice vectors

b1 = (2.*np.pi/q) * np.array([1, -1/np.sqrt(3)])
b2 = (2.*np.pi) * np.array([0, 2/np.sqrt(3)])
bvec = np.vstack((b1, b2))


def hamiltonian(k, M):

    a = 1
    c = np.sqrt(3) * a / 6  # ... / 6
    eta = 1 * k[0] * M * a / 2  # 3 * ...
    alpha = float(p / q)

    matrix = np.zeros(shape=(M, M), dtype=np.complex128)

    def B(phi, m):
        return np.exp(1j * np.pi * phi / 3) * 2 * np.cos(np.pi * phi * (m + 1 / 2) + 3 * k[1] * c)

    def C(phi):
        return np.exp(1j * np.pi * phi / 3)

    for i in range(M - 1):
        matrix[i, i + 1] = B(alpha, i + 1)
        matrix[i + 1, i] = np.conj(B(alpha, i + 1))
    for i in range(M - 2):
        matrix[i, i + 2] = np.conj(C(alpha))
        matrix[i + 2, i] = C(alpha)

    # top-right
    matrix[0][M - 1] = np.conj(B(alpha, M)) * np.exp(-1j * eta)
    matrix[0][M - 2] = C(alpha) * np.exp(-1j * eta)
    matrix[1][M - 1] = C(alpha) * np.exp(-1j * eta)

    # bottom-left
    matrix[M - 1][0] = B(alpha, M) * np.exp(1j * eta)
    matrix[M - 2][0] = np.conj(C(alpha)) * np.exp(1j * eta)
    matrix[M - 1][1] = np.conj(C(alpha)) * np.exp(1j * eta)

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

    eigenvalues = np.zeros((2*M, numb_samples, numb_samples))
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
                eigenvalues[band][idx_x][idx_y] = np.real(+np.sqrt(3 + eigvals[idx[band]]))
                eigenvalues[M + band][idx_x][idx_y] = np.real(-np.sqrt(3 + eigvals[idx[band]]))
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

    E = np.zeros(2*M, dtype=object)

    # color index
    colors = plt.cm.RdBu_r(np.linspace(0, 1, 2 * M))
    cidx = np.zeros(2 * M, dtype=int)
    for i in range(M):
        cidx[i] = M + i
        cidx[M + i] = (M - 1) - i

    for i in range(2*M):
        E[i] = eigenvalues[i, KX, KY]
        ax.plot_surface(KX, KY, E[i], color=colors[cidx[i]])

    ########################
    # Band touching report #
    ########################

    # Using maxima and minima

    maxima = np.zeros(2 * M)
    minima = np.zeros(2 * M)
    for i in range(2 * M):
        maxima[i] = np.max(E[i])
        minima[i] = np.min(E[i])

    bidx = np.zeros(2 * M, dtype=int)
    for i in range(M):
        bidx[i] = (2 * M - 1) - i
        bidx[M + i] = i

    gap = np.zeros(2 * M - 1)
    for i in range(2 * M - 1):
        gap[i] = minima[bidx[i + 1]] - maxima[bidx[i]]

    for i in range(len(gap)):
        if gap[i] < 1e-5:
            print("Bands", i + 1, "and", i + 2, "have overlapping extrema.")
    print(gap)

    # Using robust method

    threshold = 1e-5

    for i in range(2 * M - 1):
        for kx in idx_x:
            for ky in idx_y:
                if abs(eigenvalues[bidx[i], kx, ky] - eigenvalues[bidx[i + 1], kx, ky]) < threshold:
                    print("Bands", i + 1, "and", i + 2, "touch within a threshold of", threshold)
                    break
            if abs(eigenvalues[bidx[i], kx, ky] - eigenvalues[bidx[i + 1], kx, ky]) < threshold:
                break

    ########################

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

    for i in range(M):
        start = -0.05
        stop = 0.03
        interval = (stop-start)/(M-1)
        if i == 0:
            ax.text2D(-0.105, start+i*interval, "$C_{{\pm{:2d}}}={:2d}$".format(i+1, 2*int(round(chern_numbers[i]))), color="k", fontsize=14)
        else:
            ax.text2D(-0.105, start + i * interval, "$C_{{{:2d}}}={:2d}$".format(i + 1, int(round(chern_numbers[i]))),
                      color="k", fontsize=14)

    ax.text2D(-0.107, 0.07, "(a) $\phi={}/{}$".format(p, q), color="k", fontsize=18)

    ax.tick_params(axis="x", labelsize=14)
    ax.tick_params(axis="y", labelsize=14)
    ax.tick_params(axis="z", labelsize=14)

    # plt.rc('xtick', labelsize=20)
    # plt.rc('ytick', labelsize=20)

    # plt.savefig("/home/bart/Documents/papers/TBG/figures/hex_1_bands_{}_{}.png".format(p, q), bbox_inches='tight', dpi=300)
    plt.show()
