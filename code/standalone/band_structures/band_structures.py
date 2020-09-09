import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d

a0 = 1  # lattice constant
t = 1  # hopping amplitude


def define_unit_cell(lattice_val, model_val):

    if lattice_val is 'square':
        if model_val is 'HalSquC1':
            num_bands_val = 2
            # lattice vectors
            a1 = a0 * np.array([2, 0])
            a2 = a0 * np.array([0, 1])
            avec_val = np.vstack((a1, a2))
            # reciprocal lattice vectors
            b1 = (2. * np.pi) / a0 * np.array([1 / 2, 0])
            b2 = (2. * np.pi) / a0 * np.array([0, 1])
            bvec_val = np.vstack((b1, b2))
        elif model_val is 'HalSquC2':
            num_bands_val = 2
            # lattice vectors
            a1 = a0 * np.array([1, 0])
            a2 = a0 * np.array([0, 1])
            avec_val = np.vstack((a1, a2))
            # reciprocal lattice vectors
            b1 = (2. * np.pi) / a0 * np.array([1, 0])
            b2 = (2. * np.pi) / a0 * np.array([0, 1])
            bvec_val = np.vstack((b1, b2))
        else:
            return ValueError("Requested model is not implemented in define_unit_cell function.")
    elif lattice_val is 'triangular':
        if model_val is 'HalTriC3':
            num_bands_val = 2
            # lattice vectors
            a1 = a0 * np.array([1, 0])
            a2 = a0 * np.array([1/2, np.sqrt(3)/2])
            avec_val = np.vstack((a1, a2))
            # reciprocal lattice vectors
            b1 = (2. * np.pi) / a0 * np.array([1, -1/np.sqrt(3)])
            b2 = (2. * np.pi) / a0 * np.array([0, 2/np.sqrt(3)])
            bvec_val = np.vstack((b1, b2))
        else:
            return ValueError("Requested model is not implemented in define_unit_cell function.")
    elif lattice_val is 'honeycomb':
        num_bands_val = 2
        # lattice vectors
        a1 = (a0 / 2) * np.array([3, np.sqrt(3)])
        a2 = (a0 / 2) * np.array([3, -np.sqrt(3)])
        avec_val = np.vstack((a1, a2))
        # reciprocal lattice vectors
        b1 = (2. * np.pi) / (3 * a0) * np.array([1, np.sqrt(3)])
        b2 = (2. * np.pi) / (3 * a0) * np.array([1, -np.sqrt(3)])
        bvec_val = np.vstack((b1, b2))
    else:
        return ValueError("Requested lattice is not implemented in define_unit_cell function.")

    return num_bands_val, avec_val, bvec_val


def hamiltonian(lattice_val, model_val, k_val, num_bands_val, avec_val):

    # initialize the Hamiltonian
    Hamiltonian = np.zeros((num_bands_val, num_bands_val), dtype=np.complex128)

    if lattice_val is 'square':

        # nearest neighbors
        delta = np.zeros((2, 2))
        delta[0, :] = a0 * np.array([1, 0])
        delta[1, :] = a0 * np.array([0, 1])

        if model_val in ['HalSquC1', 'HalSquC2']:
            tphi = t * np.exp(1j * np.pi/4)
            tdash = t/(2+np.sqrt(2))
            tddash = t/(2+2*np.sqrt(2))

            # first-nearest neighbors
            firstNN = np.zeros((4, 2))
            # positive direction for AdB
            firstNN[0, :] = delta[1, :]  # up
            firstNN[1, :] = -delta[1, :]  # down
            # positive direction for BdA
            firstNN[2, :] = -delta[0, :]  # left
            firstNN[3, :] = delta[0, :]  # right

            # second-nearest neighbors
            secondNN = np.zeros((2, 2))
            secondNN[0, :] = (delta[0, :] + delta[1, :])
            secondNN[1, :] = (delta[0, :] - delta[1, :])

            # third-nearest neighbors
            thirdNN = np.zeros((2, 2))
            thirdNN[0, :] = 2*delta[0, :]
            thirdNN[1, :] = 2*delta[1, :]

            f_AdB = 0
            for m in range(0, 2):
                f_AdB += tphi * np.exp(1j * k_val.dot(firstNN[m, :]))
            f_BdA = 0
            for m in range(2, 4):
                f_BdA += tphi * np.exp(1j * k_val.dot(firstNN[m, :]))

            fdash_A = tdash * np.exp(1j * k_val.dot(secondNN[0, :])) - tdash * np.exp(1j * k_val.dot(secondNN[1, :]))
            fdash_B = -tdash * np.exp(1j * k_val.dot(secondNN[0, :])) + tdash * np.exp(1j * k_val.dot(secondNN[1, :]))

            fddash = 0
            for m in range(0, 2):
                fddash += tddash * np.exp(1j * k_val.dot(thirdNN[m, :]))

            Hamiltonian[0][0] = fdash_A + np.conj(fdash_A) + fddash + np.conj(fddash)
            Hamiltonian[0][1] = f_AdB + np.conj(f_BdA)
            Hamiltonian[1][0] = f_BdA + np.conj(f_AdB)
            Hamiltonian[1][1] = fdash_B + np.conj(fdash_B) + fddash + np.conj(fddash)
        else:
            return ValueError("Requested model is not implemented in hamiltonian function.")

    elif lattice_val is 'triangular':

        # nearest neighbors
        delta = np.zeros((3, 2))
        delta[0, :] = a0 * np.array([1, 0])
        delta[1, :] = a0 * np.array([-1/2, np.sqrt(3)/2])
        delta[2, :] = a0 * np.array([-1/2, -np.sqrt(3)/2])

        if model_val is 'HalTriC3':
            t1 = t
            t2 = 0.39 * t * np.exp(1j * np.pi/2)
            t3 = -0.34 * t

            # first-nearest neighbors
            firstNN = np.zeros((3, 2))
            firstNN[0, :] = delta[0, :]
            firstNN[1, :] = delta[1, :]
            firstNN[2, :] = delta[2, :]

            # second-nearest neighbors
            secondNN = np.zeros((6, 2))
            # positive direction for A
            secondNN[0, :] = a0 * np.array([0, np.sqrt(3)])
            secondNN[1, :] = a0 * np.array([1.5, -np.sqrt(3)/2])
            secondNN[2, :] = a0 * np.array([-1.5, -np.sqrt(3)/2])
            # positive direction for B
            secondNN[3, :] = a0 * np.array([1.5, np.sqrt(3)/2])
            secondNN[4, :] = a0 * np.array([-1.5, np.sqrt(3)/2])
            secondNN[5, :] = a0 * np.array([0, -np.sqrt(3)])

            # third-nearest neighbors
            thirdNN = np.zeros((3, 2))
            thirdNN[0, :] = 2 * delta[0, :]
            thirdNN[1, :] = 2 * delta[1, :]
            thirdNN[2, :] = 2 * delta[2, :]

            f1 = 0
            for m in range(0, 3):
                f1 += -t1 * np.exp(1j * k_val.dot(firstNN[m, :]))

            f2_A = 0
            for m in range(0, 3):
                f2_A += -t2 * np.exp(1j * k_val.dot(secondNN[m, :]))
            f2_B = 0
            for m in range(3, 6):
                f2_B += -t2 * np.exp(1j * k_val.dot(secondNN[m, :]))

            f3 = 0
            for m in range(0, 3):
                f3 += -t3 * np.exp(1j * k_val.dot(thirdNN[m, :]))

            Hamiltonian[0][0] = f2_A + np.conj(f2_A)
            Hamiltonian[0][1] = f1 + np.conj(f3)
            Hamiltonian[1][0] = np.conj(f1) + f3
            Hamiltonian[1][1] = f2_B + np.conj(f2_B)
        else:
            return ValueError("Requested model is not implemented in hamiltonian function.")

    elif lattice_val is 'honeycomb':

        # nearest neighbors
        delta = np.zeros((3, 2))
        delta[0, :] = (1 / 3) * avec[0, :] + (1 / 3) * avec[1, :]
        delta[1, :] = (-2 / 3) * avec[0, :] + (1 / 3) * avec[1, :]
        delta[2, :] = (1 / 3) * avec[0, :] + (-2 / 3) * avec[1, :]

        if model_val is 'graphene':
            f = 0
            for m in range(3):
                f += -t * np.exp(1j * k_val.dot(delta[m, :]))

            Hamiltonian[0][0] = 0
            Hamiltonian[0][1] = f
            Hamiltonian[1][0] = np.conj(f)
            Hamiltonian[1][1] = 0
        elif model_val is 'HalHexC1':
            t2 = t * np.sqrt(129) / 36
            phi = np.arccos(3*np.sqrt(3/43))

            # next-nearest neighbors
            secondNN = np.zeros((6, 2))
            # positive direction for A sites / negative direction for B sites
            secondNN[0, :] = avec_val[0, :]
            secondNN[1, :] = -avec_val[0, :] + avec_val[1, :]
            secondNN[2, :] = -avec_val[1, :]
            # negative direction for A sites / positive direction for B sites
            secondNN[3, :] = avec_val[1, :]
            secondNN[4, :] = -avec_val[0, :]
            secondNN[5, :] = avec_val[0, :] - avec_val[1, :]

            f = 0
            for m in range(0, 3):
                f += t * np.exp(1j * k_val.dot(delta[m, :]))
            f1 = 0
            for m in range(0, 3):
                f1 += t2 * np.exp(1j * k_val.dot(secondNN[m, :]))
            f2 = 0
            for m in range(3, 6):
                f2 += t2 * np.exp(1j * k_val.dot(secondNN[m, :]))

            Hamiltonian[0][0] = np.exp(1j * phi) * f1 + np.exp(-1j * phi) * f2
            Hamiltonian[0][1] = f
            Hamiltonian[1][0] = np.conj(f)
            Hamiltonian[1][1] = np.exp(1j * phi) * f2 + np.exp(-1j * phi) * f1
        else:
            return ValueError("Requested model is not implemented in hamiltonian function.")
    else:
        return ValueError("Requested lattice is not implemented in hamiltonian function.")

    return Hamiltonian


def berry_curv(eigvec, eigvec_alpha, eigvec_beta, eigvec_alpha_beta):

    Berry_curv = - np.imag(np.log(np.conj(eigvec).dot(eigvec_alpha) * np.conj(eigvec_alpha).dot(eigvec_alpha_beta)
                           * np.conj(eigvec_alpha_beta).dot(eigvec_beta) * np.conj(eigvec_beta).dot(eigvec)))

    return Berry_curv


if __name__ == '__main__':

    # initialization
    num_samples = 101
    lattice = 'square'  # square, triangular, honeycomb
    model = 'HalSquC2'  # HalSquC1, HalSquC2, HalTriC3, graphene, HalHexC1

    # define unit cell
    num_bands, avec, bvec = define_unit_cell(lattice, model)

    # construct bands
    eigenvalues = np.zeros((num_bands, num_samples, num_samples))  # real
    eigenvectors = np.zeros((num_bands, num_bands, num_samples, num_samples), dtype=np.complex128)  # complex
    for band in range(num_bands):
        for idx_x in range(num_samples):
            frac_kx = idx_x / (num_samples-1)
            for idx_y in range(num_samples):
                frac_ky = idx_y / (num_samples-1)
                k = np.matmul(np.array([frac_kx, frac_ky]), bvec)
                eigvals, eigvecs = np.linalg.eig(hamiltonian(lattice, model, k, num_bands, avec))
                idx = np.argsort(eigvals)
                eigenvalues[band][idx_x][idx_y] = np.real(eigvals[idx[band]])
                eigenvectors[:, band, idx_x, idx_y] = eigvecs[:, idx[band]]

    # compute Chern numbers
    berry_fluxes = np.zeros((num_bands, num_samples-1, num_samples-1))  # real
    for band in range(num_bands):
        for idx_x in range(num_samples-1):
            for idx_y in range(num_samples-1):
                berry_fluxes[band, idx_x, idx_y] = berry_curv(eigenvectors[:, band, idx_x, idx_y],
                                                              eigenvectors[:, band, idx_x+1, idx_y],
                                                              eigenvectors[:, band, idx_x, idx_y+1],
                                                              eigenvectors[:, band, idx_x+1, idx_y+1])
    chern_numbers = np.zeros(num_bands)
    for band in range(num_bands):
        chern_numbers[band] = np.sum(berry_fluxes[band, :, :]) / (2*np.pi)
        print(f"Chern number ({band}) = {chern_numbers[band]}")

    # band analysis
    band_gap = np.min(eigenvalues[1]) - np.max(eigenvalues[0])
    band_width = np.max(eigenvalues[0]) - np.min(eigenvalues[0])
    print(f"band width (0) = {band_width}")
    print(f"band gap (0-1) = {band_gap}")
    print(f"gap-to-width ratio (0-1) = {band_gap/band_width}")

    # plot figure
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    idx_x = np.linspace(0, num_samples-1, num_samples, dtype=int)
    idx_y = np.linspace(0, num_samples-1, num_samples, dtype=int)
    kx, ky = np.meshgrid(idx_x, idx_y)
    for i in range(num_bands):
        ax.plot_surface(kx, ky, eigenvalues[i, kx, ky])
    ax.set_xlabel('$k_x$')
    ax.set_ylabel('$k_y$')
    ax.set_zlabel('$E$')
    plt.show()
