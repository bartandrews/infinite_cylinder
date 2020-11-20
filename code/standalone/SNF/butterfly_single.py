import numpy as np
import matplotlib.pyplot as plt
import sys
import matplotlib.gridspec as gridspec
from matplotlib.ticker import FormatStrFormatter

q = 299  # prime > 7


def matrix_eigenvalues(p, M, t1, t2, t2dash):

    k = [0, 0]
    a = 1
    c = np.sqrt(3) * a / 6
    eta = 3 * k[0] * M * a / 2

    matrix = np.zeros(shape=(M, M), dtype=np.complex128)

    alpha = float(p / q)

    def A(phi, m):
        return t2 * (2 * np.cos(4 * np.pi * phi * m + 12 * k[1] * c)
                     + 2 * np.cos(2 * np.pi * phi * (m - 3 / 2) + 6 * k[1] * c)
                     + 2 * np.cos(2 * np.pi * phi * (m + 3 / 2) + 6 * k[1] * c) + 3) \
               + t2dash * (-2 * np.cos(4 * np.pi * phi * m + 12 * k[1] * c)
                     - 2 * np.cos(2 * np.pi * phi * (m - 3 / 2) + 6 * k[1] * c)
                     - 2 * np.cos(2 * np.pi * phi * (m + 3 / 2) + 6 * k[1] * c) - 9)

    def B(phi, m):
        return t1 * (2*np.exp(1j * np.pi*phi/3)*np.cos(np.pi*phi*(m+1/2)+3*k[1]*c))

    def C(phi):
        return t1 * np.exp(1j*np.pi*phi/3)

    def D(phi, m):
        return t2 * (2 * np.cos(np.pi * phi * (m - 5 / 2) + 3 * k[1] * c)
                     + 2 * np.cos(np.pi * phi * (m + 7 / 2) + 3 * k[1] * c)
                     + 2 * np.cos(3 * np.pi * phi * (m - 1 / 2) + 9 * k[1] * c)
                     + 2 * np.cos(3 * np.pi * phi * (m + 3 / 2) + 9 * k[1] * c)) \
               + t2dash * (-2 * np.cos(np.pi * phi * (m - 5 / 2) + 3 * k[1] * c)
                     - 2 * np.cos(np.pi * phi * (m + 7 / 2) + 3 * k[1] * c)
                     - 2 * np.cos(3 * np.pi * phi * (m - 1 / 2) + 9 * k[1] * c)
                     - 2 * np.cos(3 * np.pi * phi * (m + 3 / 2) + 9 * k[1] * c))

    def G(phi, m):
        return t2 * (2 * np.cos(2 * np.pi * phi * (m + 1) + 6 * k[1] * c)
                     + 2 * np.cos(3 * np.pi * phi)) \
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
    matrix[0][M - 1] = np.conj(B(alpha, M))

    matrix[0][M - 2] = C(alpha)
    matrix[1][M - 1] = C(alpha)

    matrix[0][M - 3] = D(alpha, M - 1) * np.exp(-1j * eta)
    matrix[1][M - 2] = D(alpha, M) * np.exp(-1j * eta)
    matrix[2][M - 1] = D(alpha, 1) * np.exp(-1j * eta)

    matrix[0][M - 6] = G(alpha, M - 3)
    matrix[1][M - 5] = G(alpha, M - 2)
    matrix[2][M - 4] = G(alpha, M - 1)
    matrix[3][M - 3] = G(alpha, M)
    matrix[4][M - 2] = G(alpha, 1)
    matrix[5][M - 1] = G(alpha, 2)

    # bottom-left
    matrix[M - 1][0] = B(alpha, M) * np.exp(1j * eta)

    matrix[M - 2][0] = np.conj(C(alpha))
    matrix[M - 1][1] = np.conj(C(alpha))

    matrix[M - 3][0] = D(alpha, M - 1) * np.exp(1j * eta)
    matrix[M - 2][1] = D(alpha, M) * np.exp(1j * eta)
    matrix[M - 1][2] = D(alpha, 1) * np.exp(1j * eta)

    matrix[M - 6][0] = G(alpha, M - 3) * np.exp(1j * eta)
    matrix[M - 5][1] = G(alpha, M - 2) * np.exp(1j * eta)
    matrix[M - 4][2] = G(alpha, M - 1) * np.exp(1j * eta)
    matrix[M - 3][3] = G(alpha, M) * np.exp(1j * eta)
    matrix[M - 2][4] = G(alpha, 1)
    matrix[M - 1][5] = G(alpha, 2)

    lmbda = np.real(np.linalg.eigvals(matrix))

    eenergies = np.zeros(2 * len(lmbda))
    for i in range(len(lmbda)):
        eenergies[i] = +np.sqrt(3 + lmbda[i])
        eenergies[len(lmbda) + i] = -np.sqrt(3 + lmbda[i])

    return eenergies


if __name__ == '__main__':

    values0 = []

    for p in range(q):
        if p % 2 == 0:
            M = q
        else:
            M = 2 * q

        alpha = float(p/q)
        alpha_list = [alpha for i in range(2*M)]
        eigenvalues0 = matrix_eigenvalues(p, M, 1, -0.0075, 0.03)
        values0.append((eigenvalues0, alpha_list))

    ##########
    # Figure #
    ##########

    fig = plt.figure(figsize=(4, 3.9))

    ax0 = plt.subplot(111)
    ax0.tick_params('x', direction='out', bottom=True)

    for eigenvalues, alphas in values0:
        ax0.plot(alphas, eigenvalues, '.', color='red', markersize=0.2)

    ax0.tick_params(axis='both', which='major', labelsize=15)

    ax0.axvline(1/3, color='k', linewidth=0.5, ls="--")
    ax0.axvline(2/3, color='k', linewidth=0.5, ls="--")

    ax0.set_xlim([0, 1])
    ax0.set_xlabel('$n_{\phi}$', fontsize=17, labelpad=9)
    ax0.xaxis.set_major_formatter(FormatStrFormatter('%g'))
    ax0.set_ylabel('$E$ / meV', fontsize=17)
    #ax0.set_title("hello", fontsize=14)

    plt.savefig("/home/bart/Documents/Berkeley/research_statement/figures/butterfly.png", bbox_inches='tight', dpi=1000)
    plt.show()
