#######################################################
# CHERN NUMBER FOR THE FOUR BAND MODEL (CONVENTION I) #
#######################################################

import numpy as np
import matplotlib.pyplot as plt
import sys

# lattice vectors (normalized)

a1 = (1/2) * np.array([np.sqrt(3), 1])
a2 = (1/2) * np.array([np.sqrt(3), -1])
avec = np.vstack((a1, a2))

# reciprocal lattice vectors

b1 = (2.*np.pi) * np.array([1/np.sqrt(3), 1])
b2 = (2.*np.pi) * np.array([1/np.sqrt(3), -1])
bvec = np.vstack((b1, b2))

# high symmetry points

K1 = np.array([2/3, 1/3])
K2 = np.array([1/3, 2/3])
GA = np.array([0., 0.])
MM = np.array([0.5, 0.5])

# print the lattice

print("avec = ", avec)
print("bvec = ", bvec)
print("K1 = ", np.matmul(K1, bvec))
print("K2 = ", np.matmul(K2, bvec))
print("GA = ", np.matmul(GA, bvec))
print("MM = ", np.matmul(MM, bvec))

# sanity check the lattice

print("|avec| = ", np.sqrt(avec[0, 0]**2+avec[0, 1]**2))
print("a0 . b0 = ", avec[0, :].dot(bvec[0, :]))
print("a1 . b1 = ", avec[1, :].dot(bvec[1, :]))
print("a1 . b0 = ", avec[1, :].dot(bvec[0, :]))
print("a0 . b1 = ", avec[0, :].dot(bvec[1, :]))

# tau vectors for positions of orbitals within the home unit cell

tau = np.zeros((2, 2))
tau[0, :] = np.array([0, 0])
tau[1, :] = (1/3)*np.array([1, 1])

# plot the lattice

# plt.scatter(np.matmul(K1, bvec)[0], np.matmul(K1, bvec)[1], color='blue')
# plt.scatter(np.matmul(K2, bvec)[0], np.matmul(K2, bvec)[1], color='green')
# plt.scatter(np.matmul(GA, bvec)[0], np.matmul(GA, bvec)[1], color='orange')
# plt.scatter(np.matmul(MM, bvec)[0], np.matmul(MM, bvec)[1], color='red')
# plt.plot([0, avec[0, 0]], [0, avec[0, 1]], color='black')
# plt.plot([0, avec[1, 0]], [0, avec[1, 1]], color='black')
# plt.plot([0, bvec[0, 0]], [0, bvec[0, 1]], color='purple')
# plt.plot([0, bvec[1, 0]], [0, bvec[1, 1]], color='purple')
# plt.axis('equal')  # plt.gca().set_aspect('equal', adjustable='box')
# plt.show()
# sys.exit()


def hamiltonian(k):

    t1 = 0.331
    t2 = 0.01
    t2dash = 0.097

    delta = np.zeros((3, 2))
    # delta[0, :] = (1 / 2) * np.array([1, np.sqrt(3)])
    # delta[1, :] = (1 / 2) * np.array([1, -np.sqrt(3)])
    # delta[2, :] = np.array([-1, 0])
    delta[0, :] = (1 / 3) * avec[0, :] + (1 / 3) * avec[1, :]
    delta[1, :] = (-2 / 3) * avec[0, :] + (1 / 3) * avec[1, :]
    delta[2, :] = (1 / 3) * avec[0, :] + (-2 / 3) * avec[1, :]
    
    fifthNN = np.zeros((6, 2))
    fifthNN[0, :] = -avec[0, :] - avec[1, :]
    fifthNN[1, :] = -avec[0, :] + 2*avec[1, :]
    fifthNN[2, :] = 2*avec[0, :] - avec[1, :]
    fifthNN[3, :] = avec[0, :] + avec[1, :]
    fifthNN[4, :] = -2*avec[0, :] + avec[1, :]
    fifthNN[5, :] = avec[0, :] - 2*avec[1, :]

    # plot the neighbors

    # plt.scatter(delta[:, 0], delta[:, 1], color='blue')
    # plt.scatter(fifthNN[0:3, 0], fifthNN[0:3, 1], color='red')
    # plt.scatter(fifthNN[3:6, 0], fifthNN[3:6, 1], color='red', marker='x')
    #
    # plt.plot([0, avec[0, 0]], [0, avec[0, 1]], color='black')
    # plt.plot([0, avec[1, 0]], [0, avec[1, 1]], color='black')
    # plt.axis('equal')  # plt.gca().set_aspect('equal', adjustable='box')
    # plt.show()
    # sys.exit()

    Hamiltonian = np.zeros((4, 4), dtype=complex)

    f = 0
    for i in range(0, 3):
        f += t1 * np.exp(1j * k.dot(delta[i, :]))
    f2 = 0
    for i in range(0, 3):
        f2 += t2 * np.exp(1j * k.dot(fifthNN[i, :]))
    xi = 0
    for i in range(0, 3):
        xi += t2dash * np.exp(1j * k.dot(fifthNN[i, :]))

    # Zeeman coefficients

    Bx = 0.5
    By = 0.5
    Bz = 0.5

    Hamiltonian[0][0] = f2 + np.conj(f2) + Bz
    Hamiltonian[0][1] = f
    Hamiltonian[1][0] = np.conj(f)
    Hamiltonian[1][1] = f2 + np.conj(f2) + Bz

    Hamiltonian[2][2] = f2 + np.conj(f2) - Bz
    Hamiltonian[2][3] = f
    Hamiltonian[3][2] = np.conj(f)
    Hamiltonian[3][3] = f2 + np.conj(f2) - Bz

    Hamiltonian[0][2] = xi - np.conj(xi) + Bx - By*1j
    Hamiltonian[2][0] = -xi + np.conj(xi) + Bx - By*1j
    Hamiltonian[1][3] = xi - np.conj(xi) + Bx + By*1j
    Hamiltonian[3][1] = -xi + np.conj(xi) + Bx + By*1j

    return Hamiltonian


def get_eigen(k):

    eigval, eigvec = np.linalg.eigh(hamiltonian(k))
    idx = np.argsort(eigval)[::-1]
    eigvals = np.real(eigval[idx])
    eigvecs = eigvec[:, idx]

    return eigvals[0], eigvecs[0], eigvals[1], eigvecs[1], eigvals[2], eigvecs[2], eigvals[3], eigvecs[3]


def berry_curv(ev, ev_alpha, ev_beta, ev_alpha_beta):

    bc = - np.imag(np.log(np.conj(ev).dot(ev_alpha) * np.conj(ev_alpha).dot(ev_alpha_beta)
                              * np.conj(ev_alpha_beta).dot(ev_beta) * np.conj(ev_beta).dot(ev)))

    return bc


def bloch_factor(k):

    bloch_factor = np.zeros(2, dtype=complex)

    bloch_factor[0] = np.exp(-1j * k.dot(np.matmul(tau[0, :], avec)))
    bloch_factor[1] = np.exp(-1j * k.dot(np.matmul(tau[1, :], avec)))

    return bloch_factor


if __name__ == '__main__':

    ##################################
    # Plot 2D Haldane Band Structure #
    ##################################

    open("2D_four_band_structure.txt", "w")
    band_structure_2d = open("2D_four_band_structure.txt", "a", buffering=1)

    count = 0
    nk = 100

    for i in range(nk):
        k = K1 - K1 * float(i) / float(nk - 1)
        k = np.matmul(k, bvec)
        a, a_vec, b, b_vec, c, c_vec, d, d_vec = get_eigen(k)
        band_structure_2d.write("{} {} {} {} {} {} {} {} {}\n"
                                .format(count, a, b, c, d,
                                        abs(a_vec[0]) ** 2, abs(b_vec[0]) ** 2, abs(c_vec[0]) ** 2, abs(d_vec[0]) ** 2))
        count += 1

    for i in range(nk):
        k = GA + (MM - GA) * float(i) / float(nk - 1)
        k = np.matmul(k, bvec)
        a, a_vec, b, b_vec, c, c_vec, d, d_vec = get_eigen(k)
        band_structure_2d.write("{} {} {} {} {} {} {} {} {}\n"
                                .format(count, a, b, c, d,
                                        abs(a_vec[0]) ** 2, abs(b_vec[0]) ** 2, abs(c_vec[0]) ** 2, abs(d_vec[0]) ** 2))
        count += 1

    for i in range(nk):
        k = MM + (K2 - MM) * float(i) / float(nk - 1)
        k = np.matmul(k, bvec)
        a, a_vec, b, b_vec, c, c_vec, d, d_vec = get_eigen(k)
        band_structure_2d.write("{} {} {} {} {} {} {} {} {}\n"
                                .format(count, a, b, c, d,
                                        abs(a_vec[0]) ** 2, abs(b_vec[0]) ** 2, abs(c_vec[0]) ** 2, abs(d_vec[0]) ** 2))
        count += 1

    ###############################
    # Calculate the Chern numbers #
    ###############################

    num_bands = len(hamiltonian(np.array([0, 0])))  # get the dimension of the square hamiltonian

    samples_x = 101
    samples_y = 101
    max_idx_x = samples_x - 1
    max_idx_y = samples_y - 1

    u_matrix = np.zeros((num_bands, num_bands, samples_x, samples_y), dtype=complex)

    for idx_x in range(samples_x):
        for idx_y in range(samples_y):

            frac_kx = idx_x / max_idx_x
            frac_ky = idx_y / max_idx_y

            # wavevector used for the Bloch factor (always the same)
            k_bloch = np.matmul(np.array([frac_kx, frac_ky]), bvec)

            # wavevector used for u (depends on boundary conditions)
            if idx_x == max_idx_x and idx_y == max_idx_y:
                k_u = np.matmul(np.array([0, 0]), bvec)
                factor = bloch_factor(k_bloch)
            elif idx_x == max_idx_x:
                k_u = np.matmul(np.array([0, frac_ky]), bvec)
                factor = bloch_factor(k_bloch)
            elif idx_y == max_idx_y:
                k_u = np.matmul(np.array([frac_kx, 0]), bvec)
                factor = bloch_factor(k_bloch)
            else:
                k_u = np.matmul(np.array([frac_kx, frac_ky]), bvec)
                factor = 1

            eigval, u_eigvec = np.linalg.eig(hamiltonian(k_bloch))
            idx0 = np.argsort(eigval)[::-1]
            for band in range(num_bands):
                u_matrix[:, band, idx_x, idx_y] = u_eigvec[:, idx0[band]]

    berry_flux_matrix = np.zeros((num_bands, samples_x-1, samples_y-1))

    for band in range(num_bands):
        for idx_x in range(samples_x-1):
            for idx_y in range(samples_y-1):
                berry_flux_matrix[band, idx_x, idx_y] = berry_curv(u_matrix[:, band, idx_x, idx_y],
                                                                   u_matrix[:, band, idx_x + 1, idx_y],
                                                                   u_matrix[:, band, idx_x, idx_y + 1],
                                                                   u_matrix[:, band, idx_x + 1, idx_y + 1])

    print("")
    for band in range(num_bands):
        print("Chern number (", band, ") = ", np.sum(berry_flux_matrix[band, :, :])/(2*np.pi))

