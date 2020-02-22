from mpl_toolkits.mplot3d import axes3d, Axes3D
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec


# lattice vectors (normalized)

a1 = (1/2) * np.array([np.sqrt(3), 1])
a2 = (1/2) * np.array([np.sqrt(3), -1])
avec = np.vstack((a1, a2))

# reciprocal lattice vectors

b1 = (2.*np.pi) * np.array([1/np.sqrt(3), 1])
b2 = (2.*np.pi) * np.array([1/np.sqrt(3), -1])
bvec = np.vstack((b1, b2))

# normed euclidean unit vectors

e1 = (2.*np.pi) * np.array([1, 0])
e2 = (2.*np.pi) * np.array([0, 1])
evec = np.vstack((e1, e2))

# high symmetry points

K1 = np.array([2/3, 1/3])
K2 = np.array([1/3, 2/3])
GA = np.array([0., 0.])
MM = np.array([0.5, 0.5])


def hamiltonian(k, num_bands_1):

    t, t0, t3 = 3, 0.3, 0.3

    delta = np.zeros((3, 2))
    delta[0, :] = (1 / 3) * avec[0, :] + (1 / 3) * avec[1, :]
    delta[1, :] = (-2 / 3) * avec[0, :] + (1 / 3) * avec[1, :]
    delta[2, :] = (1 / 3) * avec[0, :] + (-2 / 3) * avec[1, :]

    Hamiltonian = np.zeros((num_bands_1, num_bands_1), dtype=np.complex128)

    f = 0
    for i in range(3):
        f += np.exp(1j * k.dot(delta[i, :]))

    Hamiltonian[0][0] = 0
    Hamiltonian[0][1] = 0
    Hamiltonian[0][2] = -t*f
    Hamiltonian[0][3] = t3*np.conj(f)

    Hamiltonian[1][1] = 0
    Hamiltonian[1][2] = t0
    Hamiltonian[1][3] = -t*f

    Hamiltonian[2][2] = 0
    Hamiltonian[2][3] = 0

    Hamiltonian[3][3] = 0

    # H.c.

    Hamiltonian[1][0] = np.conj(Hamiltonian[0][1])

    Hamiltonian[2][0] = np.conj(Hamiltonian[0][2])
    Hamiltonian[2][1] = np.conj(Hamiltonian[1][2])

    Hamiltonian[3][0] = np.conj(Hamiltonian[0][3])
    Hamiltonian[3][1] = np.conj(Hamiltonian[1][3])
    Hamiltonian[3][2] = np.conj(Hamiltonian[2][3])

    return Hamiltonian


if __name__ == '__main__':

    num_bands = 4
    num_samples = 101

    ############
    # System 1 #########################################################################################################
    ############

    eigenvalues1 = np.zeros((num_bands, num_samples, num_samples))
    eigenvectors1 = np.zeros((num_bands, num_bands, num_samples, num_samples), dtype=np.complex128)

    for band in range(num_bands):
        for idx_x in range(num_samples):
            frac_kx = -1 + 2*idx_x / (num_samples-1)
            for idx_y in range(num_samples):
                frac_ky = -1 + 2*idx_y / (num_samples-1)

                k_u = np.matmul(np.array([frac_kx, frac_ky]), evec)

                eigvals, eigvecs = np.linalg.eig(hamiltonian(k_u, num_bands))
                idx = np.argsort(eigvals)
                # eigenvalues1[band][idx_x][idx_y] = np.real(+np.sqrt(3 + eigvals[idx[band]]))
                eigenvalues1[(num_bands-1) - band][idx_x][idx_y] = np.real(-np.sqrt(3 + eigvals[idx[band]]))
                eigenvectors1[:, (num_bands-1) - band, idx_x, idx_y] = eigvecs[:, idx[band]]


    ##########
    # Figure ###########################################################################################################
    ##########

    fig = plt.figure()

    ###############
    # Subfigure 1 ######################################################################################################
    ###############

    idx_x = np.linspace(-1, num_samples - 1, num_samples, dtype=int)
    idx_y = np.linspace(-1, num_samples - 1, num_samples, dtype=int)

    KX, KY = np.meshgrid(idx_x, idx_y)

    ax = plt.subplot(111, projection='3d')

    E1 = np.zeros(num_bands, dtype=object)

    # # color index
    # colors = plt.cm.RdBu_r(np.linspace(0, 1, 2 * M1))
    # cidx = np.zeros(2 * M1, dtype=int)
    # for i in range(M1):
    #     cidx[i] = M1 + i
    #     cidx[M1 + i] = (M1 - 1) - i

    for i in range(num_bands):
        E1[i] = eigenvalues1[i, KX, KY]
        ax.plot_surface(KX, KY, E1[i])

    ax.set_aspect('equal', adjustable='box')

    ax.tick_params(axis='x', which='major', pad=-4)
    ax.tick_params(axis='y', which='major', pad=-3)
    ax.tick_params(axis='z', which='major', pad=3)

    ax.set_xlabel('$k_x$', fontsize=11, linespacing=-1)
    ax.set_ylabel('\n$k_y$', fontsize=11, linespacing=1)
    ax.set_zlabel('\n$E$ / meV', fontsize=11, linespacing=1)

    ax.xaxis.labelpad = -2
    ax.yaxis.labelpad = -2

    def custom(value, tick_number):

        if value == 0:
            return "0"
        elif value == num_samples-1:
            return "1"
        else:
            return "{}".format(round(value / (num_samples-1), 2))

    ax.xaxis.set_major_formatter(plt.FuncFormatter(custom))
    ax.yaxis.set_major_formatter(plt.FuncFormatter(custom))

    ax.tick_params(axis="x", labelsize=10)
    ax.tick_params(axis="y", labelsize=10)
    ax.tick_params(axis="z", labelsize=10)

    plt.savefig("/home/bart/Documents/FSK2020/proposal/figures/bilayer_graphene/3d_band_structure.png", bbox_inches='tight', dpi=300)
    plt.show()
