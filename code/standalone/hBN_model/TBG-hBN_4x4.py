"""Calculate the band structure in Fig. 1 of
"Twisted Bilayer Graphene Aligned with Hexagonal Boron Nitride: Anomalous Hall Effect and a Lattice Model"""

import numpy as np
import matplotlib.pyplot as plt
import sys

lat_const = 1.
avec = lat_const*np.array([[0.5*np.sqrt(3), 0.5], [0, 1]])
b2 = 2.*np.pi/lat_const * np.array([-1./np.sqrt(3), 1.])
b1 = 2.*np.pi/lat_const * np.array([ 2./np.sqrt(3), 0.])
bvec = np.vstack((b1, b2))

# print(np.sqrt(avec[0, 0]**2+avec[0, 1]**2))
# print(avec[0, :].dot(bvec[0, :]))
# print(avec[1, :].dot(bvec[1, :]))
# print(avec[1, :].dot(bvec[0, :]))
# print(avec[0, :].dot(bvec[1, :]))

K1 = np.array([1/3, 2/3])
K2 = np.array([2/3, 1/3])
GA = np.array([0., 0.])
MM = np.array([0.5, 0.5])


def hamiltonian(k):

    # # --- twist of 1.08 degrees
    # m0 = 1.36
    # t = 1.22
    # tAd = 0.670 * np.exp(1j * 0.366 * np.pi)
    # tBd = 0.731 * np.exp(-1j * 0.657 * np.pi)
    # tdd = 0.801 * np.exp(-1j * 0.685 * np.pi)
    # t1ddd = 0.123 * np.exp(-1j * 0.48 * np.pi)
    # t2ddd = 0.355 * np.exp(-1j * 0.411 * np.pi)

    # --- twist of 1.20 degrees
    m0 = 0.076
    t = 3.056
    tAd = 0.837 * np.exp(1j * 0.56 * np.pi)
    tBd = 0.828 * np.exp(-1j * 0.469 * np.pi)
    tdd = 2.062 * np.exp(-1j * 0.54 * np.pi)
    t1ddd = 0.915 * np.exp(-1j * 0.337 * np.pi)
    t2ddd = 0.815 * np.exp(-1j * 0.434 * np.pi)

    delta = np.zeros((3, 2))
    delta[0, :] = np.matmul(np.array([-1/3, 2/3]), avec)
    delta[1, :] = np.matmul(np.array([-1/3, -1/3]), avec)
    delta[2, :] = np.matmul(np.array([2/3, -1/3]), avec)

    secondNN = np.zeros((6, 2))
    secondNN[0, :] = avec[0, :]
    secondNN[1, :] = avec[1, :]
    secondNN[2, :] = -avec[0, :] + avec[1, :]
    secondNN[3, :] = -avec[0, :]
    secondNN[4, :] = -avec[1, :]
    secondNN[5, :] = avec[0, :] - avec[1, :]

    thirdNN = np.zeros((3, 2))
    thirdNN[0, :] = avec[0, :] + delta[0, :]
    thirdNN[1, :] = -avec[0, :] + delta[0, :]
    thirdNN[2, :] = -avec[1, :] + delta[2, :]

    fourthNN = np.zeros((6, 2))
    # t1ddd
    fourthNN[0, :] = -avec[1, :] + delta[1, :]
    fourthNN[1, :] = -avec[0, :] + avec[1, :] + delta[0, :]
    fourthNN[2, :] = avec[0, :] + delta[2, :]
    # t2ddd
    fourthNN[3, :] = -avec[0, :] + delta[1, :]
    fourthNN[4, :] = avec[1, :] + delta[0, :]
    fourthNN[5, :] = -avec[1, :] + avec[0, :] + delta[2, :]

    # plt.scatter(delta[:, 0], delta[:, 1], color='blue')
    # plt.scatter(secondNN[:, 0], secondNN[:, 1], color='green')
    # plt.scatter(thirdNN[:, 0], thirdNN[:, 1], color='orange')
    # plt.scatter(fourthNN[0:3, 0], fourthNN[0:3, 1], color='red')
    # plt.scatter(fourthNN[3:6, 0], fourthNN[3:6, 1], color='red', marker='x')
    # #
    # plt.plot([0, avec[0, 0]], [0, avec[0, 1]], color='black')
    # plt.plot([0, avec[1, 0]], [0, avec[1, 1]], color='black')
    # plt.axis('equal')  # plt.gca().set_aspect('equal', adjustable='box')
    # plt.show()
    # sys.exit()

    Hamiltonian = np.zeros((2, 2), dtype=complex)

    f = 0
    for i in range(0, 3):
        f += np.exp(1j * k.dot(delta[i, :]))
    fd = 0
    for i in range(0, 6):
        fd += np.exp(1j * k.dot(secondNN[i, :]))
    fdd = 0
    for i in range(0, 3):
        fdd += np.exp(1j * k.dot(thirdNN[i, :]))
    f1ddd = 0
    for i in range(0, 3):
        f1ddd += np.exp(1j * k.dot(fourthNN[i, :]))
    f2ddd = 0
    for i in range(3, 6):
        f2ddd += np.exp(1j * k.dot(fourthNN[i, :]))

    Hamiltonian[0][0] = -tAd*fd + m0 + np.conj(-tAd*fd)
    Hamiltonian[0][1] = np.conj(-t1ddd*f1ddd - t2ddd*f2ddd - t*f - tdd*fdd)
    Hamiltonian[1][0] = -t1ddd*f1ddd - t2ddd*f2ddd - t*f - tdd*fdd
    Hamiltonian[1][1] = -tBd*fd - m0 + np.conj(-tBd*fd)

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


if __name__ == '__main__':

    ##############################
    # Plot the 2D Band Structure #
    ##############################

    open("band_structure.txt", "w")
    band_structure = open("band_structure.txt", "a", buffering=1)

    nk = 101
    count = 0

    for i in range(nk):

        count = count+1
        k = K1-K1*float(i)/float(nk-1)
        k = np.matmul(k, bvec)
        count = count+1
        a, a_vec, b, b_vec = get_eigen(k)
        band_structure.write("{} {} {} {} {}\n".format(count, a, b, abs(a_vec[0])**2, abs(b_vec[0])**2))

    for i in range(nk):

        count = count+1
        k = GA+(MM-GA)*float(i)/float(nk-1)
        k = np.matmul(k, bvec)
        count = count+1
        a, a_vec, b, b_vec = get_eigen(k)
        band_structure.write("{} {} {} {} {}\n".format(count, a, b, abs(a_vec[0])**2, abs(b_vec[0])**2))

    for i in range(nk):

        count = count+1
        k = MM+(K2-MM)*float(i)/float(nk-1)
        k = np.matmul(k, bvec)
        count = count+1
        a, a_vec, b, b_vec = get_eigen(k)
        band_structure.write("{} {} {} {} {}\n".format(count, a, b, abs(a_vec[0])**2, abs(b_vec[0])**2))

    ############################
    # Plot the Berry Curvature #
    ############################

    open("berry_curvature.txt", "w")
    berry_curvature = open("berry_curvature.txt", "a", buffering=1)

    min_x = -np.pi
    max_x = np.pi
    samples_x = 300
    min_y = -np.pi
    max_y = np.pi
    samples_y = 300
    delta_alpha = (max_x - min_x) / samples_x
    delta_beta = (max_y - min_y) / samples_y

    tot_ber_curv_0, tot_ber_curv_1 = 0, 0

    for kx in np.linspace(min_x, max_x, samples_x):
        for ky in np.linspace(min_y, max_y, samples_y):

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

            print(kx, ky, ber_curv_0, ber_curv_1)
            berry_curvature.write("{} {} {} {}\n".format(kx, ky, ber_curv_0, ber_curv_1))

        berry_curvature.write("\n")

    print("")
    print("Chern number (0) = ", tot_ber_curv_0 / (2 * np.pi))
    print("Chern number (1) = ", tot_ber_curv_1 / (2 * np.pi))

