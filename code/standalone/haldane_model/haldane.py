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


# def hamiltonian(k):
#
#     ##################
#     # Initialization #
#     ##################
#
#     t = 1
#     t2 = 0
#     # t2 = np.sqrt(129)/36
#     M = 0
#     phi = 0
#     # phi = np.arccos(3*np.sqrt(3/43))
#
#     # Pauli spin
#     s0 = np.array([[1, 0], [0, 1]])
#     sx = np.array([[0, 1], [1, 0]])
#     sy = np.array([[0, -1j], [1j, 0]])
#     sz = np.array([[1, 0], [0, -1]])
#
#     a = np.array([[(np.sqrt(3) / 2), (1 / 2)],
#                   [(-np.sqrt(3) / 2), (1 / 2)],
#                   [0, -1]])
#
#     b = np.array([[(-np.sqrt(3) / 2), (3 / 2)],
#                   [(-np.sqrt(3) / 2), (-3 / 2)],
#                   [-np.sqrt(3), 0]])
#
#     def h0(k):
#
#         h0 = 0
#         for i in range(3):
#             h0 += 2*t2*np.cos(phi)*np.cos(k.dot(b[i]))
#
#         return h0
#
#     def hx(k):
#
#         hx = 0
#         for i in range(3):
#             hx += t * np.cos(k.dot(a[i]))
#
#         return hx
#
#     def hy(k):
#
#         hy = 0
#         for i in range(3):
#             hy += -t * np.sin(k.dot(a[i]))
#
#         return hy
#
#     def hz(k):
#
#         hz = M
#         for i in range(3):
#             hz += M + 2 * t2 * np.sin(phi) * np.sin(k.dot(b[i]))
#
#         return hz
#
#     def hamiltonian_internal(k):
#
#         return h0(k)*s0 + hx(k)*sx + hy(k)*sy + hz(k)*sz
#
#     return hamiltonian_internal(k)


def hamiltonian(k):

    t1 = 1
    t2 = np.sqrt(129)/36
    phi = np.arccos(3*np.sqrt(3/43))

    delta = np.zeros((3, 2))
    delta[0, :] = (1/2)*np.array([1, np.sqrt(3)])
    delta[1, :] = (1/2)*np.array([1, -np.sqrt(3)])
    delta[2, :] = np.array([-1, 0])

    secondNN = np.zeros((6, 2))
    # positive direction for A sites / negative direction for B sites
    secondNN[0, :] = avec[0, :]
    secondNN[1, :] = -avec[0, :] + avec[1, :]
    secondNN[2, :] = -avec[1, :]
    # negative direction for A sites / positive direction for B sites
    secondNN[3, :] = avec[1, :]
    secondNN[4, :] = -avec[0, :]
    secondNN[5, :] = avec[0, :] - avec[1, :]

    # plt.scatter(delta[:, 0], delta[:, 1], color='blue')
    # plt.scatter(secondNN[0:3, 0], secondNN[0:3, 1], color='red')
    # plt.scatter(secondNN[3:6, 0], secondNN[3:6, 1], color='red', marker='x')
    # #
    # plt.plot([0, avec[0, 0]], [0, avec[0, 1]], color='black')
    # plt.plot([0, avec[1, 0]], [0, avec[1, 1]], color='black')
    # plt.axis('equal')  # plt.gca().set_aspect('equal', adjustable='box')
    # plt.show()
    # sys.exit()

    Hamiltonian = np.zeros((2, 2), dtype=complex)

    f = 0
    for i in range(0, 3):
        f += -t1 * np.exp(1j * k.dot(delta[i, :]))
    f1 = 0
    for i in range(0, 3):
        f1 += t2 * np.exp(1j * k.dot(secondNN[i, :]))
    f2 = 0
    for i in range(3, 6):
        f2 += t2 * np.exp(1j * k.dot(secondNN[i, :]))

    Hamiltonian[0][0] = np.exp(1j * phi) * f1 + np.exp(-1j * phi) * f2
    Hamiltonian[0][1] = f
    Hamiltonian[1][0] = np.conj(f)
    Hamiltonian[1][1] = np.exp(1j * phi) * f2 + np.exp(-1j * phi) * f1

    return Hamiltonian


def get_eigen(k):

    eigval, eigvec = np.linalg.eigh(hamiltonian(k))
    idx = np.argsort(eigval)[::-1]
    eigvals = np.real(eigval[idx])
    eigvecs = eigvec[:, idx]

    return eigvals[0], eigvecs[0], eigvals[1], eigvecs[1]


def berry_curv(ev, ev_alpha, ev_beta, ev_alpha_beta):

    bc = - 2 * np.imag(np.log(np.conj(ev).dot(ev_alpha) * np.conj(ev_alpha).dot(ev_alpha_beta)
                              * np.conj(ev_alpha_beta).dot(ev_beta) * np.conj(ev_beta).dot(ev)))

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

    # print(get_eigen(np.matmul(K1, bvec)))
    # sys.exit()

    ##################################
    # Plot 2D Haldane Band Structure #
    ##################################

    open("2D_haldane_band_structure.txt", "w")
    band_structure_2d = open("2D_haldane_band_structure.txt", "a", buffering=1)

    count = 0
    nk = 100

    for i in range(nk):

        k = K1-K1*float(i)/float(nk-1)
        k = np.matmul(k, bvec)
        a, a_vec, b, b_vec = get_eigen(k)
        band_structure_2d.write("{} {} {} {} {}\n".format(count, a, b, abs(a_vec[0])**2, abs(b_vec[0])**2))
        count += 1

    for i in range(nk):

        k = GA+(MM-GA)*float(i)/float(nk-1)
        k = np.matmul(k, bvec)
        a, a_vec, b, b_vec = get_eigen(k)
        band_structure_2d.write("{} {} {} {} {}\n".format(count, a, b, abs(a_vec[0])**2, abs(b_vec[0])**2))
        count += 1

    for i in range(nk):

        k = MM+(K2-MM)*float(i)/float(nk-1)
        k = np.matmul(k, bvec)
        a, a_vec, b, b_vec = get_eigen(k)
        band_structure_2d.write("{} {} {} {} {}\n".format(count, a, b, abs(a_vec[0])**2, abs(b_vec[0])**2))
        count += 1

    ############################################################
    # Plot 3D Haldane Band Structure & Haldane Berry Curvature #
    ############################################################

    open("3D_haldane_band_structure.txt", "w")
    band_structure_3d = open("3D_haldane_band_structure.txt", "a", buffering=1)
    open("haldane_berry_curvature.txt", "w")
    berry_curvature = open("haldane_berry_curvature.txt", "a", buffering=1)

    min_x = -np.pi
    max_x = np.pi
    samples_x = 100
    min_y = -np.pi
    max_y = np.pi
    samples_y = 100
    delta_alpha = (max_x - min_x) / samples_x
    delta_beta = (max_y - min_y) / samples_y

    tot_ber_curv_0, tot_ber_curv_1 = 0, 0

    for kx in np.linspace(min_x, max_x, samples_x):
        if abs(kx) <= abs(np.matmul(K1, bvec)[0]):
            for ky in np.linspace(min_y, max_y, samples_y):
                if abs(ky) <= abs(top_line(kx)):

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

                    eigval_alpha_beta, eigvec_alpha_beta = np.linalg.eig(hamiltonian(np.array([kx + delta_alpha, ky + delta_beta])))
                    idx3 = np.argsort(eigval_alpha_beta)[::-1]
                    eigvals_alpha_beta = np.real(eigval_alpha_beta[idx3])
                    eigvecs_alpha_beta = eigvec_alpha_beta[:, idx3]

                    ber_curv_0 = berry_curv(eigvecs[0], eigvecs_alpha[0], eigvecs_beta[0], eigvecs_alpha_beta[0])
                    ber_curv_1 = berry_curv(eigvecs[1], eigvecs_alpha[1], eigvecs_beta[1], eigvecs_alpha_beta[1])

                    tot_ber_curv_0 += ber_curv_0
                    tot_ber_curv_1 += ber_curv_1

                    #print(kx, ky, ber_curv_0, ber_curv_1)
                    band_structure_3d.write("{} {} {} {} {} {}\n".format(kx, ky, eigvals[0], eigvals[1],
                                                                         abs(eigvecs[0][0])**2, abs(eigvecs[1][0])**2))
                    berry_curvature.write("{} {} {} {}\n".format(kx, ky, ber_curv_0, ber_curv_1))

            berry_curvature.write("\n")

    print("")
    print("Chern number (0) = ", tot_ber_curv_0 / (2 * np.pi))
    print("Chern number (1) = ", tot_ber_curv_1 / (2 * np.pi))
