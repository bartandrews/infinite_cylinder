import numpy as np
import matplotlib.pyplot as plt
import sys
import matplotlib.gridspec as gridspec


q = 199  # prime > 7


def matrix_eigenvalues(p, M, t1):

    t2 = 1
    kx, ky = 0, 0
    a = 1
    c = np.sqrt(3) * a / 6
    delta = kx * M * a / 2

    matrix = np.zeros(shape=(M, M), dtype=np.complex128)

    alpha = float(p/q)

    def A(phi, m):
        return 2 * t1 * np.cos(6 * ky * c) + 2 * t2 * np.cos(2 * np.pi * phi * m + 18 * ky * c)

    def B():
        return 2 * t1 * np.cos(3 * ky * c)

    def C(phi, m):
        return 2 * t2 * np.cos(np.pi * phi * (m + 1/2) + 9 * ky * c)

    for i in range(M):
        matrix[i][i] = A(alpha, i+1)
    for i in range(M - 1):
        matrix[i, i + 1] = B()
        matrix[i + 1, i] = B()
    for i in range(M - 3):
        matrix[i, i + 3] = C(alpha, i+2)
        matrix[i + 3, i] = C(alpha, i+2)

    # top-right
    matrix[0][M - 3] = C(alpha, M - 1) * np.exp(-1j * delta)
    matrix[1][M - 2] = C(alpha, M) * np.exp(-1j * delta)
    matrix[2][M - 1] = C(alpha, 1) * np.exp(-1j * delta)
    matrix[0][M - 2] = 0
    matrix[1][M - 1] = 0
    matrix[0][M - 1] = B() * np.exp(-1j * delta)

    # bottom-left
    matrix[M - 3][0] = C(alpha, M - 1) * np.exp(1j * delta)
    matrix[M - 2][1] = C(alpha, M) * np.exp(1j * delta)
    matrix[M - 1][2] = C(alpha, 1) * np.exp(1j * delta)
    matrix[M - 2][0] = 0
    matrix[M - 1][1] = 0
    matrix[M - 1][0] = B() * np.exp(1j * delta)

    return np.real(np.linalg.eigvals(matrix))


if __name__ == '__main__':

    values0 = []
    values1 = []
    values2 = []

    for p in range(q):
        if p % 2 == 0:
            M = q
        else:
            M = 2 * q

        alpha = float(p/q)
        alpha_list = [alpha for i in range(M)]
        eigenvalues0 = matrix_eigenvalues(p, M, 0)
        values0.append((eigenvalues0, alpha_list))
        eigenvalues1 = matrix_eigenvalues(p, M, 0.5)
        values1.append((eigenvalues1, alpha_list))
        eigenvalues2 = matrix_eigenvalues(p, M, 1)
        values2.append((eigenvalues2, alpha_list))

    ##########
    # Figure #
    ##########

    fig = plt.figure()

    gs = gridspec.GridSpec(3, 1, height_ratios=[1, 1, 1])

    ax0 = plt.subplot(gs[0])
    ax0.tick_params('x', direction='in', bottom=True)

    for eigenvalues, alphas in values0:
        ax0.plot(alphas, eigenvalues, '.', color='r', markersize=1.8)

    ax1 = plt.subplot(gs[1], sharex=ax0)
    ax1.tick_params('x', direction='in', bottom=True)

    for eigenvalues, alphas in values1:
        ax1.plot(alphas, eigenvalues, '.', color='r', markersize=1.8)

    plt.setp(ax0.get_xticklabels(), visible=False)

    ax2 = plt.subplot(gs[2], sharex=ax0)

    for eigenvalues, alphas in values2:
        ax2.plot(alphas, eigenvalues, '.', color='r', markersize=1.8)

    plt.setp(ax1.get_xticklabels(), visible=False)

    ax0.axvline(1/3, color='k', linewidth=0.5, ls="--")
    ax0.axvline(2/3, color='k', linewidth=0.5, ls="--")
    ax1.axvline(1 / 3, color='k', linewidth=0.5, ls="--")
    ax1.axvline(2 / 3, color='k', linewidth=0.5, ls="--")
    ax2.axvline(1 / 3, color='k', linewidth=0.5, ls="--")
    ax2.axvline(2 / 3, color='k', linewidth=0.5, ls="--")

    gs.update(hspace=0)

    ax0.set_xlim([0, 1])
    ax2.set_xlabel('$\phi$')
    ax1.set_ylabel('$E$ / meV')

    ax0.text(0.4, 4.5, '(a) $t_1 = 0$', fontsize=10)
    ax1.text(0.4, 7, '(b) $t_1 = 0.5$', fontsize=10)
    ax2.text(0.4, 9.5, '(c) $t_1 = 1$', fontsize=10)

    plt.show()
