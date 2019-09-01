#####################################################
# CHERN NUMBER FOR THE HALDANE MODEL (CONVENTION I) #
#####################################################

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
#
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


def hamiltonian(k, phi_input, M_input):

    t1 = 1
    t2 = 1
    phi = phi_input
    M = M_input

    delta = np.zeros((3, 2))
    #* BUG is here (I think)
    #* this makes crystal structure to be non-equivalent to the in convetion_II script
    # delta[0, :] = (2/3)*avec[0, :] - (1/3)*avec[1, :]
    # delta[1, :] = (-1/3)*avec[0, :] + (2/3)*avec[1, :]
    # delta[2, :] = (-1/3)*(avec[0, :] + avec[1, :])
    #* I change it:
    delta[0, :] = (1/3) *avec[0, :] + (1/3 )*avec[1, :]
    delta[1, :] = (-2/3)*avec[0, :] + (1/3 )*avec[1, :]
    delta[2, :] = (1/3) *avec[0, :] + (-2/3)*avec[1, :]
 

    secondNN = np.zeros((6, 2))
    # positive direction for A sites / negative direction for B sites
    secondNN[0, :] = avec[0, :]
    secondNN[1, :] = -avec[0, :] + avec[1, :]
    secondNN[2, :] = -avec[1, :]
    # negative direction for A sites / positive direction for B sites
    secondNN[3, :] = avec[1, :]
    secondNN[4, :] = -avec[0, :]
    secondNN[5, :] = avec[0, :] - avec[1, :]

    # plot the neighbors
    #
    # plt.scatter(delta[:, 0], delta[:, 1], color='blue')
    # plt.scatter(secondNN[0:3, 0], secondNN[0:3, 1], color='red')
    # plt.scatter(secondNN[3:6, 0], secondNN[3:6, 1], color='red', marker='x')
    #
    # plt.plot([0, avec[0, 0]], [0, avec[0, 1]], color='black')
    # plt.plot([0, avec[1, 0]], [0, avec[1, 1]], color='black')
    # plt.axis('equal')  # plt.gca().set_aspect('equal', adjustable='box')
    # plt.show()
    # sys.exit()

    Hamiltonian = np.zeros((2, 2), dtype=np.complex_)

    f = 0
    for i in range(0, 3):
        f += -t1 * np.exp(1j * k.dot(delta[i, :]))
    f1 = 0
    for i in range(0, 3):
        f1 += -t2 * np.exp(1j * k.dot(secondNN[i, :]))
    f2 = 0
    for i in range(3, 6):
        f2 += -t2 * np.exp(1j * k.dot(secondNN[i, :]))

    Hamiltonian[0][0] = np.exp(1j * phi) * f1 + np.exp(-1j * phi) * f2 + M
    Hamiltonian[0][1] = f
    Hamiltonian[1][0] = np.conj(f)
    Hamiltonian[1][1] = np.exp(1j * phi) * f2 + np.exp(-1j * phi) * f1 - M

    return Hamiltonian


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

    open("haldane_I_phase_diagram_with_bloch_original_neighbours.dat", "w")
    phase_diagram = open("haldane_I_phase_diagram_with_bloch_original_neighbours.dat", "a", buffering=1)

    samples_x = 101
    samples_y = 101
    max_idx_x = samples_x - 1
    max_idx_y = samples_y - 1

    for phi in np.linspace(-np.pi, np.pi, 20):
        for M in np.linspace(-10, 10, 20):

            ####################################################################################
            # Matrix of u eigenvectors for each site, size = (samples_x)*(samples_y) = 101*101 #
            ####################################################################################

            u_matrix = np.zeros((2, 2, samples_x, samples_y), dtype=np.complex_)
            eig_matrix = np.zeros((2, samples_x, samples_y), dtype=np.float_)

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

                    eigval, u_eigvec = np.linalg.eig(hamiltonian(k_bloch, phi, M))
                    idx0 = np.argsort(eigval)[::-1]
                    for i in range(2):
                        eig_matrix[i, idx_x, idx_y] = np.real(eigval[idx0[i]])
                        u_matrix[:, i, idx_x, idx_y] = factor * u_eigvec[:, idx0[i]]

            ###############################################################################################
            # Matrix of Berry fluxes for each plaquette, size = (samples_x - 1)*(samples_y - 1) = 100*100 #
            ###############################################################################################

            berry_flux_matrix = np.zeros((samples_x-1, samples_y-1), dtype=np.ndarray)

            for idx_x in range(samples_x-1):
                for idx_y in range(samples_y-1):

                    berry_flux_0 = berry_curv(u_matrix[:, 0, idx_x, idx_y], u_matrix[:, 0, idx_x + 1, idx_y],
                                              u_matrix[:, 0, idx_x, idx_y + 1], u_matrix[:, 0, idx_x + 1, idx_y + 1])

                    berry_flux_1 = berry_curv(u_matrix[:, 1, idx_x, idx_y], u_matrix[:, 1, idx_x + 1, idx_y],
                                              u_matrix[:, 1, idx_x, idx_y + 1], u_matrix[:, 1, idx_x + 1, idx_y + 1])

                    berry_flux_matrix[idx_x, idx_y] = np.array([berry_flux_0, berry_flux_1])

            berry_flux_matrix_0 = \
                [[berry_flux_matrix[idx_x, idx_y][0] for idx_x in range(samples_x-1)] for idx_y in range(samples_y-1)]
            berry_flux_matrix_1 = \
                [[berry_flux_matrix[idx_x, idx_y][1] for idx_x in range(samples_x-1)] for idx_y in range(samples_y-1)]

            #* write eigenvalues in all available k-points *#
            f = open('eigvals_I.dat', 'w')
            for i in range(len(eig_matrix[:, 0, 0])):
                count = 0
                for j in range(len(eig_matrix[0, :, 0])):
                    for k in range(len(eig_matrix[0, 0, :])):
                        count = count+1
                        e = eig_matrix[i, j, k]
                        f.write(str(count)+' '+str(e)+'\n')
                f.write('\n')
            f.close()

            print(phi, M, np.sum(berry_flux_matrix_0) / (2 * np.pi))
            phase_diagram.write("{} {} {}\n".format(phi, M, np.sum(berry_flux_matrix_0) / (2 * np.pi)))

        phase_diagram.write("\n")
