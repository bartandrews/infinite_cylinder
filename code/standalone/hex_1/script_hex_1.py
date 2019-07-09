"""Based on: Energy Spectrum of a Honeycomb Lattice under Nonuniform Magnetic Fields"""

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

    matrix = np.zeros(shape=(M, M), dtype=np.complex128)

    def B(phi, m):
        return np.exp(1j*np.pi*phi/3)*2*np.cos(np.pi*phi*(m+1/2) + 3*k[1]*c)

    def C(phi):
        return np.exp(1j * np.pi * phi / 3)

    for i in range(M - 1):
        matrix[i, i + 1] = B(alpha, i + 1)
        matrix[i + 1, i] = np.conj(B(alpha, i + 1))
    for i in range(M - 2):
        matrix[i, i + 2] = np.conj(C(alpha))
        matrix[i + 2, i] = C(alpha)

    # top-right
    matrix[0][M - 1] = np.conj(B(alpha, M))
    matrix[0][M - 2] = C(alpha)
    matrix[1][M - 1] = C(alpha)

    # bottom-left
    matrix[M - 1][0] = B(alpha, M) * np.exp(1j * eta)
    matrix[M - 2][0] = np.conj(C(alpha))
    matrix[M - 1][1] = np.conj(C(alpha))

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
        ax.plot(alphas, eigenvalues, '.', color='r', markersize=1.8)

    ax.set_xlim([0, 1])
    ax.set_xlabel('$\phi$')
    ax.set_ylabel('$E$ / meV')

    plt.show()
