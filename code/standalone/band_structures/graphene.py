import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d

# 0) define the lattice constant
a0 = 1

# 1) define the lattice vectors
a1 = (a0/2) * np.array([3, np.sqrt(3)])
a2 = (a0/2) * np.array([3, -np.sqrt(3)])
avec = np.vstack((a1, a2))

# 2) define the reciprocal lattice vectors
b1 = (2.*np.pi)/(3*a0) * np.array([1, np.sqrt(3)])
b2 = (2.*np.pi)/(3*a0) * np.array([1, -np.sqrt(3)])
bvec = np.vstack((b1, b2))


# 3-5) write down the Hamiltonian matrix
def hamiltonian(k_val, num_bands_val):

    # hopping amplitude
    t = -1

    # nearest neighbors
    delta = np.zeros((3, 2))
    delta[0, :] = (1 / 3) * avec[0, :] + (1 / 3) * avec[1, :]
    delta[1, :] = (-2 / 3) * avec[0, :] + (1 / 3) * avec[1, :]
    delta[2, :] = (1 / 3) * avec[0, :] + (-2 / 3) * avec[1, :]

    # initialize the Hamiltonian
    Hamiltonian = np.zeros((num_bands_val, num_bands_val), dtype=np.complex128)

    # construct the f function
    f = 0
    for i in range(3):
        f += t * np.exp(1j * k_val.dot(delta[i, :]))

    # define the Hamiltonian
    Hamiltonian[0][0] = 0
    Hamiltonian[0][1] = f
    Hamiltonian[1][0] = np.conj(f)
    Hamiltonian[1][1] = 0

    return Hamiltonian


if __name__ == '__main__':

    num_bands = 2
    num_samples = 101

    eigenvalues = np.zeros((num_bands, num_samples, num_samples))
    eigenvectors = np.zeros((num_bands, num_bands, num_samples, num_samples), dtype=np.complex128)

    for band in range(num_bands):
        for idx_x in range(num_samples):
            frac_kx = idx_x / (num_samples-1)
            for idx_y in range(num_samples):
                frac_ky = idx_y / (num_samples-1)

                k = np.matmul(np.array([frac_kx, frac_ky]), bvec)

                eigvals, eigvecs = np.linalg.eig(hamiltonian(k, num_bands))
                idx = np.argsort(eigvals)
                eigenvalues[band][idx_x][idx_y] = eigvals[idx[band]]
                eigenvectors[:, band, idx_x, idx_y] = eigvecs[:, idx[band]]

    # plot figure
    fig = plt.figure()
    ax = fig.gca(projection='3d')

    idx_x = np.linspace(0, num_samples - 1, num_samples, dtype=int)
    idx_y = np.linspace(0, num_samples - 1, num_samples, dtype=int)
    KX, KY = np.meshgrid(idx_x, idx_y)

    for i in range(num_bands):
        ax.plot_surface(KX, KY, eigenvalues[i, KX, KY])

    # label figure
    ax.set_xlabel('$k_x$')
    ax.set_ylabel('$k_y$')
    ax.set_zlabel('$E$ / meV')

    plt.show()
