"""Calculate the band structure in Fig. 5.a) of
"Maximally Localized Wannier Orbitals and the Extended Hubbard Model for Twisted Bilayer Graphene"""

import numpy as np
import matplotlib.pyplot as plt
import sys

# lat_const = 1.

# lattice vectors (unnormalized)

a1 = (1/2) * np.array([3, np.sqrt(3)])
a2 = (1/2) * np.array([3, -np.sqrt(3)])
avec = np.vstack((a1, a2))

# reciprocal lattice vectors

b1 = (2.*np.pi/3) * np.array([1, np.sqrt(3)])
b2 = (2.*np.pi/3) * np.array([1, -np.sqrt(3)])
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
    delta[0, :] = (1 / 2) * np.array([1, np.sqrt(3)])
    delta[1, :] = (1 / 2) * np.array([1, -np.sqrt(3)])
    delta[2, :] = np.array([-1, 0])
    
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

    Bx = 0.2
    By = 0.2
    Bz = 0.2

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

    bc = - 2 * np.imag(np.log(np.conj(ev).dot(ev_alpha) * np.conj(ev_alpha).dot(ev_alpha_beta)
                              * np.conj(ev_alpha_beta).dot(ev_beta) * np.conj(ev_beta).dot(ev)))

    return bc


def berry_curv3(ev, ev_alpha, ev_beta, ev_alpha_beta):

    gammaAB = -np.angle((np.conj(ev).dot(ev_alpha))/(np.abs(np.conj(ev).dot(ev_alpha))))
    gammaBC = -np.angle((np.conj(ev_alpha).dot(ev_alpha_beta))/(np.abs(np.conj(ev_alpha).dot(ev_alpha_beta))))
    gammaCD = -np.angle((np.conj(ev_alpha_beta).dot(ev_beta))/(np.abs(np.conj(ev_alpha_beta).dot(ev_beta))))
    gammaDA = -np.angle((np.conj(ev_beta).dot(ev))/(np.abs(np.conj(ev_beta).dot(ev))))

    bc = gammaAB + gammaBC + gammaCD + gammaDA

    return bc


def bottom_line(kx_value):

    if kx < 0:
        ky_value = np.matmul(K2, bvec)[1] + (np.matmul(K2, bvec)[1]/np.matmul(K2, bvec)[0]) * \
                   (kx_value - (-np.matmul(K2, bvec)[0]))
    else:
        ky_value = 2*np.matmul(K2, bvec)[1] - \
                   (np.matmul(K2, bvec)[1] / np.matmul(K2, bvec)[0]) * kx_value

    return ky_value


def top_line(kx_value):

    if kx < 0:
        ky_value = np.matmul(K1, bvec)[1] + (np.matmul(K1, bvec)[1] / np.matmul(K1, bvec)[0]) * \
                   (kx_value - (-np.matmul(K1, bvec)[0]))
    else:
        ky_value = 2 * np.matmul(K1, bvec)[1] \
                   - (np.matmul(K1, bvec)[1] / np.matmul(K1, bvec)[0]) * kx_value

    return ky_value


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

    ############################################################
    # Plot 3D Haldane Band Structure & Haldane Berry Curvature #
    ############################################################

    open("3D_four_band_structure.txt", "w")
    band_structure_3d = open("3D_four_band_structure.txt", "a", buffering=1)
    open("four_berry_curvature.txt", "w")
    berry_curvature = open("four_berry_curvature.txt", "a", buffering=1)

    # Need ~1000 samples to get accurate Chern number!
    min_x = -np.pi
    max_x = np.pi
    samples_x = 200
    min_y = -np.pi
    max_y = np.pi
    samples_y = 200
    delta_alpha = (max_x - min_x) / samples_x
    delta_beta = (max_y - min_y) / samples_y

    tot_ber_curv_0, tot_ber_curv_1, tot_ber_curv_2, tot_ber_curv_3 = 0, 0, 0, 0

    for kx in np.linspace(min_x, max_x, samples_x):
        if abs(kx) <= abs(np.matmul(K1, bvec)[0]):
            # if -abs(np.matmul(K1, bvec)[0]) <= kx < abs(np.matmul(K1, bvec)[0]):
            for ky in np.linspace(min_y, max_y, samples_y):
                if abs(ky) <= abs(top_line(kx)):
                    # if bottom_line(kx) <= ky < top_line(kx):

                    eigval, eigvec = np.linalg.eig(hamiltonian(np.array([kx, ky])))
                    idx0 = np.argsort(eigval)[::-1]
                    eigvals = np.real(eigval[idx0])
                    eigvecs = eigvec[:, idx0]

                    eigval_alpha, eigvec_alpha = np.linalg.eig(hamiltonian(np.array([kx + delta_alpha, ky])))
                    idx1 = np.argsort(eigval_alpha)[::-1]
                    eigvals_alpha = np.real(eigval_alpha[idx1])
                    eigvecs_alpha = eigvec_alpha[:, idx1]

                    eigval_beta, eigvec_beta = np.linalg.eig(hamiltonian(np.array([kx, ky + delta_beta])))
                    idx2 = np.argsort(eigval_beta)[::-1]
                    eigvals_beta = np.real(eigval_beta[idx2])
                    eigvecs_beta = eigvec_beta[:, idx2]

                    eigval_alpha_beta, eigvec_alpha_beta = np.linalg.eig(
                        hamiltonian(np.array([kx + delta_alpha, ky + delta_beta])))
                    idx3 = np.argsort(eigval_alpha_beta)[::-1]
                    eigvals_alpha_beta = np.real(eigval_alpha_beta[idx3])
                    eigvecs_alpha_beta = eigvec_alpha_beta[:, idx3]

                    ber_curv_0 = berry_curv3(eigvecs[0], eigvecs_alpha[0], eigvecs_beta[0], eigvecs_alpha_beta[0])
                    ber_curv_1 = berry_curv3(eigvecs[1], eigvecs_alpha[1], eigvecs_beta[1], eigvecs_alpha_beta[1])
                    ber_curv_2 = berry_curv3(eigvecs[2], eigvecs_alpha[2], eigvecs_beta[2], eigvecs_alpha_beta[2])
                    ber_curv_3 = berry_curv3(eigvecs[3], eigvecs_alpha[3], eigvecs_beta[3], eigvecs_alpha_beta[3])

                    tot_ber_curv_0 += ber_curv_0
                    tot_ber_curv_1 += ber_curv_1
                    tot_ber_curv_2 += ber_curv_2
                    tot_ber_curv_3 += ber_curv_3

                    band_structure_3d.write("{} {} {} {} {} {} {} {} {} {}\n"
                                            .format(kx, ky, eigvals[0], eigvals[1], eigvals[2], eigvals[3],
                                                    abs(eigvecs[0][0]) ** 2, abs(eigvecs[1][0]) ** 2,
                                                    abs(eigvecs[2][0]) ** 2, abs(eigvecs[3][0]) ** 2))
                    berry_curvature.write("{} {} {} {} {} {}\n".format(kx, ky, ber_curv_0, ber_curv_1,
                                                                       ber_curv_2, ber_curv_3))

            berry_curvature.write("\n")

    print("")
    print("Chern number (0) = ", tot_ber_curv_0 / (2 * np.pi))
    print("Chern number (1) = ", tot_ber_curv_1 / (2 * np.pi))
    print("Chern number (2) = ", tot_ber_curv_2 / (2 * np.pi))
    print("Chern number (3) = ", tot_ber_curv_3 / (2 * np.pi))

