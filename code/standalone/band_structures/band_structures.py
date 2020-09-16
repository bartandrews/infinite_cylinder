# --- python imports
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d

a0 = 1  # lattice constant
t = 1  # hopping amplitude


def define_unit_cell(model_val, q_val=5):

    if 'Squ' in model_val:
        if model_val is 'HofSqu1':
            num_bands_val = q_val
            # lattice vectors
            a1 = a0 * np.array([num_bands_val, 0])
            a2 = a0 * np.array([0, 1])
            avec_val = np.vstack((a1, a2))
            # reciprocal lattice vectors
            b1 = (2. * np.pi) / a0 * np.array([1 / num_bands_val, 0])
            b2 = (2. * np.pi) / a0 * np.array([0, 1])
            bvec_val = np.vstack((b1, b2))
            # symmetry points
            GA = np.array([0, 0])
            Y = np.array([0, 0.5])
            S = np.array([0.5, 0.5])
            X = np.array([0.5, 0])
            sym_points_val = [GA, Y, S, X]
        elif model_val is 'HalSquC1':
            num_bands_val = 2
            # lattice vectors
            a1 = a0 * np.array([2, 0])
            a2 = a0 * np.array([0, 1])
            avec_val = np.vstack((a1, a2))
            # reciprocal lattice vectors
            b1 = (2. * np.pi) / a0 * np.array([1 / 2, 0])
            b2 = (2. * np.pi) / a0 * np.array([0, 1])
            bvec_val = np.vstack((b1, b2))
            # symmetry points
            GA = np.array([0, 0])
            M = np.array([0.5, 0.5])
            X = np.array([0.5, 0])
            sym_points_val = [GA, M, X]
        elif model_val in ['HalSquC2', 'HalSquC3', 'HalSquC4', 'HalSquC5']:
            if model_val is 'HalSquC2':
                num_bands_val = 2
            elif model_val is 'HalSquC3':
                num_bands_val = 3
            elif model_val is 'HalSquC4':
                num_bands_val = 4
            else:
                num_bands_val = 5
            # lattice vectors
            a1 = a0 * np.array([1, 0])
            a2 = a0 * np.array([0, 1])
            avec_val = np.vstack((a1, a2))
            # reciprocal lattice vectors
            b1 = (2. * np.pi) / a0 * np.array([1, 0])
            b2 = (2. * np.pi) / a0 * np.array([0, 1])
            bvec_val = np.vstack((b1, b2))
            # symmetry points
            GA = np.array([0, 0])
            M = np.array([0.5, 0.5])
            X = np.array([0.5, 0])
            sym_points_val = [GA, M, X]
        else:
            return ValueError("Requested model is not implemented in define_unit_cell function.")
    elif 'Tri' in model_val:
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
            # symmetry points
            sym_points_val = []  # to implement
        else:
            return ValueError("Requested model is not implemented in define_unit_cell function.")
    elif 'Hex' in model_val or model_val is 'graphene':
        num_bands_val = 2
        # lattice vectors
        a1 = (a0 / 2) * np.array([3, np.sqrt(3)])
        a2 = (a0 / 2) * np.array([3, -np.sqrt(3)])
        avec_val = np.vstack((a1, a2))
        # reciprocal lattice vectors
        b1 = (2. * np.pi) / (3 * a0) * np.array([1, np.sqrt(3)])
        b2 = (2. * np.pi) / (3 * a0) * np.array([1, -np.sqrt(3)])
        bvec_val = np.vstack((b1, b2))
        # symmetry points
        K1 = np.array([2/3, 1/3])
        GA = np.array([0., 0.])
        MM = np.array([0.5, 0.5])
        K2 = np.array([1/3, 2/3])
        sym_points_val = [K1, GA, MM, K2]
    else:
        return ValueError("Requested lattice cannot be read from model name.")

    return num_bands_val, avec_val, bvec_val, sym_points_val


def hamiltonian(model_val, k_val, num_bands_val, avec_val, p_val=1, tx_factor_val=1):

    # initialize the Hamiltonian
    Hamiltonian = np.zeros((num_bands_val, num_bands_val), dtype=np.complex128)

    if 'Squ' in model_val:

        # nearest neighbors
        delta = np.zeros((2, 2))
        delta[0, :] = a0 * np.array([1, 0])
        delta[1, :] = a0 * np.array([0, 1])

        if model_val is 'HofSqu1':

            q_val = num_bands_val
            nphi = p_val/q_val

            def h(k_val_val, m_val):
                return 2 * np.cos(2 * np.pi * nphi * m_val + k_val_val[1] * a0)

            for n in range(q_val):
                Hamiltonian[n][n] = t * h(k_val, n)

            for n in range(q_val-1):
                Hamiltonian[n][n+1] = tx_factor_val*t * np.exp(+1j*k_val[0]*a0)
                Hamiltonian[n+1][n] = tx_factor_val*t * np.exp(-1j*k_val[0]*a0)

            Hamiltonian[0][q_val-1] = tx_factor_val*t * np.exp(-1j*k_val[0]*a0)
            Hamiltonian[q_val-1][0] = tx_factor_val*t * np.exp(+1j*k_val[0]*a0)

        elif model_val in ['HalSquC1', 'HalSquC2']:
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
        # elif model_val is 'HalSquC3':
        #     t1 = t
        #     t2 = - t / np.sqrt(3)
        #     phi = np.pi / 3
        #
        #     # first-nearest neighbors
        #     firstNN = np.zeros((2, 2))
        #     # positive phase direction for inter-orbital creation
        #     firstNN[0, :] = -delta[1, :]  # down
        #     # neutral direction
        #     firstNN[1, :] = delta[0, :]  # right
        #
        #     # second-nearest neighbors
        #     secondNN = np.zeros((3, 2))
        #     # positive phase direction for inter-orbital creation
        #     secondNN[0, :] = (delta[0, :] - delta[1, :])  # bottom right
        #     # positive phase direction for same-orbital creation
        #     secondNN[1, :] = (-delta[0, :] - delta[1, :])  # bottom left
        #     # negative phase direction for same orbital creation
        #     secondNN[2, :] = (delta[0, :] + delta[1, :])  # top right
        #
        #     f_BdA = t1 * np.exp(1j * 2 * phi) * np.exp(1j * k_val.dot(firstNN[0, :]))
        #     f_CdB = t1 * np.exp(1j * 4 * phi) * np.exp(1j * k_val.dot(firstNN[0, :]))
        #     f_AdC = t1 * np.exp(1j * 6 * phi) * np.exp(1j * k_val.dot(firstNN[0, :]))
        #     f_horizontal = t1 * np.exp(1j * k_val.dot(firstNN[1, :]))
        #
        #     fdash_AdA = t2 * np.exp(1j * phi) * np.exp(1j * k_val.dot(secondNN[1, :]))
        #     fdash_BdB = t2 * np.exp(1j * 3 * phi) * np.exp(1j * k_val.dot(secondNN[1, :]))
        #     fdash_CdC = t2 * np.exp(1j * 5 * phi) * np.exp(1j * k_val.dot(secondNN[1, :]))
        #
        #     fdash_CdA = t2 * np.exp(1j * 3 * phi) * np.exp(1j * k_val.dot(secondNN[0, :]))
        #     fdash_AdB = t2 * np.exp(1j * 5 * phi) * np.exp(1j * k_val.dot(secondNN[0, :]))
        #     fdash_BdC = t2 * np.exp(1j * 7 * phi) * np.exp(1j * k_val.dot(secondNN[0, :]))
        #
        #     Hamiltonian[0][0] = fdash_AdA
        #     Hamiltonian[0][1] = fdash_AdB
        #     Hamiltonian[0][2] = f_AdC + f_horizontal
        #
        #     Hamiltonian[1][0] = f_BdA + f_horizontal
        #     Hamiltonian[1][1] = fdash_BdB
        #     Hamiltonian[1][2] = fdash_BdC
        #
        #     Hamiltonian[2][0] = fdash_CdA
        #     Hamiltonian[2][1] = f_CdB + f_horizontal
        #     Hamiltonian[2][2] = fdash_CdC
        #
        #     Hamiltonian += Hamiltonian.conj().transpose()
        elif model_val in ['HalSquC3', 'HalSquC4', 'HalSquC5']:
            if model_val is 'HalSquC3':
                N = 3
            elif model_val is 'HalSquC4':
                N = 4
            else:
                N = 5

            t1 = t
            t2 = - t / np.sqrt(N)
            phi = np.pi / 3

            # first-nearest neighbors
            firstNN = np.zeros((2, 2))
            # positive phase direction for inter-orbital creation
            firstNN[0, :] = -delta[1, :]  # down
            # neutral direction
            firstNN[1, :] = delta[0, :]  # right

            # second-nearest neighbors
            secondNN = np.zeros((3, 2))
            # positive phase direction for inter-orbital creation
            secondNN[0, :] = (delta[0, :] - delta[1, :])  # bottom right
            # positive phase direction for same-orbital creation
            secondNN[1, :] = (-delta[0, :] - delta[1, :])  # bottom left
            # negative phase direction for same orbital creation
            secondNN[2, :] = (delta[0, :] + delta[1, :])  # top right

            for a in range(N):
                Hamiltonian[(a+1) % N][a] = t1 * (np.exp(1j * k_val.dot(firstNN[1, :])) + np.exp(1j * 2 * (a+1) * phi) * np.exp(1j * k_val.dot(firstNN[0, :])))
                Hamiltonian[a][a] = t2 * np.exp(1j * (2 * (a+1) - 1) * phi) * np.exp(1j * k_val.dot(secondNN[1, :]))
                Hamiltonian[(a+2) % N][a] = t2 * np.exp(1j * (2 * (a+1) + 1) * phi) * np.exp(1j * k_val.dot(secondNN[0, :]))

            # f_BdA = t1 * np.exp(1j * 2 * phi) * np.exp(1j * k_val.dot(firstNN[0, :]))
            # f_CdB = t1 * np.exp(1j * 4 * phi) * np.exp(1j * k_val.dot(firstNN[0, :]))
            # f_AdC = t1 * np.exp(1j * 6 * phi) * np.exp(1j * k_val.dot(firstNN[0, :]))
            # f_horizontal = t1 * np.exp(1j * k_val.dot(firstNN[1, :]))
            #
            # fdash_AdA = t2 * np.exp(1j * phi) * np.exp(1j * k_val.dot(secondNN[1, :]))
            # fdash_BdB = t2 * np.exp(1j * 3 * phi) * np.exp(1j * k_val.dot(secondNN[1, :]))
            # fdash_CdC = t2 * np.exp(1j * 5 * phi) * np.exp(1j * k_val.dot(secondNN[1, :]))
            #
            # fdash_CdA = t2 * np.exp(1j * 3 * phi) * np.exp(1j * k_val.dot(secondNN[0, :]))
            # fdash_AdB = t2 * np.exp(1j * 5 * phi) * np.exp(1j * k_val.dot(secondNN[0, :]))
            # fdash_BdC = t2 * np.exp(1j * 7 * phi) * np.exp(1j * k_val.dot(secondNN[0, :]))
            #
            # Hamiltonian[0][0] = fdash_AdA
            # Hamiltonian[0][1] = fdash_AdB
            # Hamiltonian[0][2] = f_AdC + f_horizontal
            #
            # Hamiltonian[1][0] = f_BdA + f_horizontal
            # Hamiltonian[1][1] = fdash_BdB
            # Hamiltonian[1][2] = fdash_BdC
            #
            # Hamiltonian[2][0] = fdash_CdA
            # Hamiltonian[2][1] = f_CdB + f_horizontal
            # Hamiltonian[2][2] = fdash_CdC

            Hamiltonian += Hamiltonian.conj().transpose()
        else:
            return ValueError("Requested model is not implemented in hamiltonian function.")

    elif 'Tri' in model_val:

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

    elif 'Hex' in model_val or model_val is 'graphene':

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
        return ValueError("Requested lattice cannot be read from model name.")

    return Hamiltonian


def berry_curv(eigvec, eigvec_alpha, eigvec_beta, eigvec_alpha_beta):

    Berry_curv = - np.imag(np.log(np.conj(eigvec).dot(eigvec_alpha) * np.conj(eigvec_alpha).dot(eigvec_alpha_beta)
                           * np.conj(eigvec_alpha_beta).dot(eigvec_beta) * np.conj(eigvec_beta).dot(eigvec)))

    return Berry_curv


if __name__ == '__main__':

    # initialization
    num_samples = 101
    model = 'HalSquC3'  # (HofSqu1, HalSquC1, HalSquC2), (HalTriC3), (graphene, HalHexC1)
    mining = False  # data mining mode for 2D band structures
    if mining:
        flag_3D = False  # only works in 2D mode
        p = [1, 4, 4, 4, 4]
        q = [4, 7, 11, 15, 19]  # lists need to be the same length
        if 'Hof' in model:
            path_to_file = "code/models/hofstadter/squ_1_hoppings_ratio.dat"
        elif 'Hal' in model:
            path_to_file = "code/models/haldane/squ_hoppings_ratio.dat"
        else:
            raise ValueError("Chosen model is not implemented for mining.")
    else:
        flag_3D = False  # choose between 3D or 2D band structure
        p, q = 1, 4  # for Hofstadter model only

    if flag_3D:
        # define unit cell
        num_bands, avec, bvec, sym_points = define_unit_cell(model, q_val=q)

        # construct bands
        eigenvalues = np.zeros((num_bands, num_samples, num_samples))  # real
        eigenvectors = np.zeros((num_bands, num_bands, num_samples, num_samples), dtype=np.complex128)  # complex
        for band in range(num_bands):
            for idx_x in range(num_samples):
                frac_kx = idx_x / (num_samples-1)
                for idx_y in range(num_samples):
                    frac_ky = idx_y / (num_samples-1)
                    k = np.matmul(np.array([frac_kx, frac_ky]), bvec)
                    eigvals, eigvecs = np.linalg.eig(hamiltonian(model, k, num_bands, avec, p_val=p, tx_factor_val=1))
                    idx = np.argsort(eigvals)
                    eigenvalues[band][idx_x][idx_y] = np.real(eigvals[idx[band]])
                    eigenvectors[:, band, idx_x, idx_y] = eigvecs[:, idx[band]]

        # compute Chern numbers
        berry_fluxes = np.zeros((num_bands, num_samples - 1, num_samples - 1))  # real
        for band in range(num_bands):
            for idx_x in range(num_samples - 1):
                for idx_y in range(num_samples - 1):
                    berry_fluxes[band, idx_x, idx_y] = berry_curv(eigenvectors[:, band, idx_x, idx_y],
                                                                  eigenvectors[:, band, idx_x + 1, idx_y],
                                                                  eigenvectors[:, band, idx_x, idx_y + 1],
                                                                  eigenvectors[:, band, idx_x + 1, idx_y + 1])
        chern_numbers = np.zeros(num_bands)
        for band in range(num_bands):
            chern_numbers[band] = np.sum(berry_fluxes[band, :, :]) / (2 * np.pi)
            print(f"Chern number ({band}) = {chern_numbers[band]}")

        # construct figure
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        idx_x = np.linspace(0, num_samples - 1, num_samples, dtype=int)
        idx_y = np.linspace(0, num_samples - 1, num_samples, dtype=int)
        kx, ky = np.meshgrid(idx_x, idx_y)
        for i in range(num_bands):
            ax.plot_surface(kx, ky, eigenvalues[i, kx, ky])
        ax.set_xlabel('$k_x$')
        ax.set_ylabel('$k_y$')
        ax.set_zlabel('$E$')
    elif mining:
        total_num_chern = len(p)
        ratio = np.zeros((total_num_chern, num_samples))
        # loop through all chern numbers
        for num_chern in range(total_num_chern):
            for num_sample, tx_factor in enumerate(np.linspace(0, 10, num_samples)):
                # define unit cell
                num_bands, avec, bvec, sym_points = define_unit_cell(model, q_val=q[num_chern])
                # construct bands
                num_paths = len(sym_points)
                points_per_path = int(num_samples / num_paths)
                num_points = num_paths * points_per_path
                eigenvalues = np.zeros((num_bands, num_points))  # real
                count = 0
                for i in range(num_paths):
                    for j in range(points_per_path):
                        k = sym_points[i] + (sym_points[(i + 1) % num_paths] - sym_points[i]) * float(j) / float(
                            points_per_path - 1)
                        k = np.matmul(k, bvec)
                        eigvals = np.linalg.eigvals(hamiltonian(model, k, num_bands, avec, p_val=p[num_chern], tx_factor_val=tx_factor))
                        idx = np.argsort(eigvals)
                        for band in range(num_bands):
                            eigenvalues[band, count] = np.real(eigvals[idx[band]])
                        count += 1
                band_gap = np.min(eigenvalues[1]) - np.max(eigenvalues[0])
                band_width = np.max(eigenvalues[0]) - np.min(eigenvalues[0])
                ratio[num_chern][num_sample] = band_gap/band_width

        with open(path_to_file, 'w') as file:
            for num_sample, tx_factor in enumerate(np.linspace(0, 10, num_samples)):
                file.write(f"{tx_factor}\t{ratio[0][num_sample]}\t{ratio[1][num_sample]}\t{ratio[2][num_sample]}\t{ratio[3][num_sample]}\t{ratio[4][num_sample]}\n")
    else:
        # define unit cell
        num_bands, avec, bvec, sym_points = define_unit_cell(model, q_val=q)

        # construct bands
        num_paths = len(sym_points)
        points_per_path = int(num_samples/num_paths)
        num_points = num_paths*points_per_path
        eigenvalues = np.zeros((num_bands, num_points))  # real
        count = 0
        for i in range(num_paths):
            for j in range(points_per_path):
                k = sym_points[i] + (sym_points[(i+1) % num_paths] - sym_points[i]) * float(j) / float(points_per_path - 1)
                k = np.matmul(k, bvec)
                eigvals = np.linalg.eigvals(hamiltonian(model, k, num_bands, avec, p_val=p, tx_factor_val=1))
                idx = np.argsort(eigvals)
                for band in range(num_bands):
                    eigenvalues[band, count] = np.real(eigvals[idx[band]])
                count += 1

        # construct figure
        fig = plt.figure()
        ax = plt.subplot(111)
        for i in range(num_bands):
            ax.plot(eigenvalues[i])
        for i in range(1, num_paths):
            ax.axvline(i*points_per_path, color='k', linewidth=0.5, ls='--')
        ax.set_xlim([0, num_points])
        ax.set_xlabel('symmetry path')
        ax.set_ylabel('$E$')

    # band analysis
    band_gap = np.min(eigenvalues[1]) - np.max(eigenvalues[0])
    band_width = np.max(eigenvalues[0]) - np.min(eigenvalues[0])
    print(f"band width (0) = {band_width}")
    print(f"band gap (0-1) = {band_gap}")
    print(f"gap-to-width ratio (0-1) = {band_gap / band_width}")

    plt.show()
