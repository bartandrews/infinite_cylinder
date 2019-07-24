import numpy as np
import matplotlib.pyplot as plt
import sys


q = 199  # prime > 7


def matrix_eigenvalues(p, M):

    k = [0, 0]
    a = 1
    c = np.sqrt(3) * a / 6
    eta = 3 * k[0] * M * a / 2
    alpha = float(p / q)
    t2, t2dash = 0, 0.097

    matrix = np.zeros(shape=(M, M), dtype=np.complex128)

    def A(phi, m):
        return t2 * (2*np.cos(4*np.pi*phi*m + 12*k[1]*c) + 2*np.cos(2*np.pi*phi*(m-3/2)+6*k[1]*c) + 2*np.cos(2*np.pi*phi*(m+3/2)+6*k[1]*c) + 3) \
               + t2dash * (-2*np.cos(4*np.pi*phi*m + 12*k[1]*c) - 2*np.cos(2*np.pi*phi*(m-3/2)+6*k[1]*c) - 2*np.cos(2*np.pi*phi*(m+3/2)+6*k[1]*c) - 9)

    def D(phi, m):
        return t2 * (2*np.cos(np.pi*phi*(m-5/2)+3*k[1]*c) + 2*np.cos(np.pi*phi*(m+7/2)+3*k[1]*c) + 2*np.cos(3*np.pi*phi*(m-1/2)+9*k[1]*c) + 2*np.cos(3*np.pi*phi*(m+3/2)+9*k[1]*c)) \
               + t2dash * (-2*np.cos(np.pi*phi*(m-5/2)+3*k[1]*c) - 2*np.cos(np.pi*phi*(m+7/2)+3*k[1]*c) - 2*np.cos(3*np.pi*phi*(m-1/2)+9*k[1]*c) - 2*np.cos(3*np.pi*phi*(m+3/2)+9*k[1]*c))

    def G(phi, m):
        return t2 * (2*np.cos(2*np.pi*phi*(m+1)+6*k[1]*c) + 2*np.cos(3*np.pi*phi)) \
               + t2dash * (-2*np.cos(2*np.pi*phi*(m+1)+6*k[1]*c) - 2*np.cos(3*np.pi*phi))

    for i in range(M):
        matrix[i, i] = A(alpha, i + 1)
    for i in range(M - 3):
        matrix[i, i + 3] = D(alpha, i + 2)
        matrix[i + 3, i] = D(alpha, i + 2)
    for i in range(M - 6):
        matrix[i, i + 6] = G(alpha, i + 3)
        matrix[i + 6, i] = G(alpha, i + 3)

    # top-right
    matrix[0][M - 3] = D(alpha, M-1) * np.exp(-1j * eta)
    matrix[1][M - 2] = D(alpha, M) * np.exp(-1j * eta)
    matrix[2][M - 1] = D(alpha, 1) * np.exp(-1j * eta)

    matrix[0][M - 6] = G(alpha, M-3)
    matrix[1][M - 5] = G(alpha, M-2)
    matrix[2][M - 4] = G(alpha, M-1)
    matrix[3][M - 3] = G(alpha, M)
    matrix[4][M - 2] = G(alpha, 1)
    matrix[5][M - 1] = G(alpha, 2)

    # bottom-left
    matrix[M - 3][0] = D(alpha, M-1) * np.exp(1j * eta)
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

    values = []

    for p in range(q):
        if p % 2 == 0:
            M = q
        else:
            M = 2 * q

        alpha = float(p/q)
        alpha_list = [alpha for i in range(2*M)]
        eigenvalues = matrix_eigenvalues(p, M)
        values.append((eigenvalues, alpha_list))

    ##########
    # Figure #
    ##########

    fig = plt.figure()
    ax = plt.subplot(1, 1, 1)
    for eigenvalues, alphas in values:
        ax.plot(alphas, eigenvalues, '.', color='r', markersize=0.5)

    ax.set_xlim([0, 1])
    ax.set_xlabel('$\phi$')
    ax.set_ylabel('$E$ / meV')

    # plt.savefig("/home/bart/Documents/papers/TBG/figures/hex_5_dash.png", bbox_inches='tight', dpi=300)
    plt.show()
