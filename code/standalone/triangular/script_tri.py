"""Based on: Energy Spectrum of a Triangular Lattice in a Uniform Magnetic Field:
Effect of Next-Nearest-Neighbor Hopping"""

import numpy as np
import matplotlib.pyplot as plt
import sys


q = 199  # prime > 7


def matrix_eigenvalues(p, M):

    t = 0
    tdash = 1

    matrix = np.zeros(shape=(M, M))

    alpha = float(p/q)

    def A(phi, m):
        return 2 * tdash * np.cos(2 * np.pi * phi * m)

    def B(phi, m):
        return 2 * t * np.cos(np.pi * phi * (m + 1/2))  # correction: (m + 1/2)

    def D(phi, m):
        return 2 * tdash * np.cos(np.pi * phi * (m + 1/2))  # correction: (m + 1/2)

    for i in range(M):
        matrix[i][i] = A(alpha, i+1)
    for i in range(M - 1):
        matrix[i, i + 1] = B(alpha, i+1)
        matrix[i + 1, i] = B(alpha, i+1)
    for i in range(M - 2):
        matrix[i, i + 2] = t
        matrix[i + 2, i] = t
    for i in range(M - 3):
        matrix[i, i + 3] = D(alpha, i+2)
        matrix[i + 3, i] = D(alpha, i+2)

    # top-right
    matrix[0][M - 3] = D(alpha, M - 1)
    matrix[1][M - 2] = D(alpha, M)
    matrix[2][M - 1] = D(alpha, 1)
    matrix[0][M - 2] = t
    matrix[1][M - 1] = t
    matrix[0][M - 1] = B(alpha, M)

    # bottom-left
    matrix[M - 3][0] = D(alpha, M - 1)
    matrix[M - 2][1] = D(alpha, M)
    matrix[M - 1][2] = D(alpha, 1)
    matrix[M - 2][0] = t
    matrix[M - 1][1] = t
    matrix[M - 1][0] = B(alpha, M)

    return np.real(np.linalg.eigvals(matrix))

    # lmbda = np.real(np.linalg.eigvals(matrix))
    #
    # eenergies = np.zeros(2 * len(lmbda))
    # for i in range(len(lmbda)):
    #     eenergies[i] = + lmbda[i]
    #     eenergies[len(lmbda) + i] = -lmbda[i]
    #
    # return eenergies


if __name__ == '__main__':

    values = []

    for p in range(q):
        if p % 2 == 0:
            M = q
        else:
            M = 2 * q

        alpha = float(p/q)
        alpha_list = [alpha for i in range(M)]
        # alpha_list = [alpha for i in range(2*M)]
        eigenvalues = matrix_eigenvalues(p, M)
        values.append((eigenvalues, alpha_list))

    ##########
    # Figure #
    ##########

    fig = plt.figure()
    ax = plt.subplot(1, 1, 1)
    for eigenvalues, alphas in values:
        ax.plot(alphas, eigenvalues, '.', color='r', markersize=1.8)

    ax.set_xlim([0, 1])
    ax.set_xlabel('$\phi$')
    ax.set_ylabel('$E$ / meV')

    plt.show()
