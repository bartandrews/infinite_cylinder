"""Calculate the band structure in Fig. 5.a) of
"Maximally Localized Wannier Orbitals and the Extended Hubbard Model for Twisted Bilayer Graphene"""

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

    t1 = 0.331
    t2 = 0.01
    t2dash = 0.097

    delta = np.zeros((3, 2))
    delta[0, :] = np.matmul(np.array([1/3, 1/3]), avec)
    delta[1, :] = np.matmul(np.array([-2/3, 1/3]), avec)
    delta[2, :] = np.matmul(np.array([1/3, -2/3]), avec)
    
    fifthNN = np.zeros((6, 2))
    fifthNN[0, :] = -avec[0, :] - avec[1, :]
    fifthNN[1, :] = -avec[0, :] + 2*avec[1, :]
    fifthNN[2, :] = 2*avec[0, :] - avec[1, :]
    fifthNN[3, :] = avec[0, :] + avec[1, :]
    fifthNN[4, :] = -2*avec[0, :] + avec[1, :]
    fifthNN[5, :] = avec[0, :] - 2*avec[1, :]

    # plt.scatter(delta[:, 0], delta[:, 1], color='blue')
    # plt.scatter(fifthNN[3, 0], fifthNN[3, 1], color='red')
    # plt.scatter(fifthNN[4, 0], fifthNN[4, 1], color='red')
    # plt.scatter(fifthNN[5, 0], fifthNN[5, 1], color='red')
    #
    # plt.plot([0, avec[0, 0]], [0, avec[0, 1]], color='k')
    # plt.plot([0, avec[1, 0]], [0, avec[1, 1]], color='k')
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

    Hamiltonian[0][0] = f2 + np.conj(f2) #+ 3
    Hamiltonian[0][1] = f
    Hamiltonian[1][0] = np.conj(f)
    Hamiltonian[1][1] = f2 + np.conj(f2) #- 1

    Hamiltonian[2][2] = f2 + np.conj(f2) #+ 1
    Hamiltonian[2][3] = f
    Hamiltonian[3][2] = np.conj(f)
    Hamiltonian[3][3] = f2 + np.conj(f2) #- 3

    Hamiltonian[0][2] = xi - np.conj(xi)
    Hamiltonian[2][0] = -xi + np.conj(xi)
    Hamiltonian[1][3] = xi - np.conj(xi)
    Hamiltonian[3][1] = -xi + np.conj(xi)

    return Hamiltonian


def hamiltonian_Haldane(k):

    ##################
    # Initialization #
    ##################

    t = 1
    t2 = np.sqrt(129)/36
    # M=0
    M = 3*np.sqrt(3)*t2 - 1
    # phi = np.arccos(3*np.sqrt(3/43))
    phi = np.pi/2

    # pauli spin
    s0 = np.array([[1, 0], [0, 1]])
    sx = np.array([[0, 1], [1, 0]])
    sy = np.array([[0, -1j], [1j, 0]])
    sz = np.array([[1, 0], [0, -1]])

    a = np.array([[(np.sqrt(3) / 2), (1 / 2)],
                  [(-np.sqrt(3) / 2), (1 / 2)],
                  [0, -1]])

    b = np.array([[(-np.sqrt(3) / 2), (3 / 2)],
                  [(-np.sqrt(3) / 2), (-3 / 2)],
                  [-np.sqrt(3), 0]])

    def h0(k):

        h0 = 0
        for i in range(3):
            h0 += 2*t2*np.cos(phi)*np.cos(k.dot(b[i]))

        return h0

    def hx(k):

        hx = 0
        for i in range(3):
            hx += t * np.cos(k.dot(a[i]))

        return hx

    def hy(k):

        hy = 0
        for i in range(3):
            hy += -t * np.sin(k.dot(a[i]))

        return hy

    def hz(k):

        hz = M
        for i in range(3):
            hz += M + 2 * t2 * np.sin(phi) * np.sin(k.dot(b[i]))

        return hz

    def hamiltonian_internal(k):

        return h0(k)*s0 + hx(k)*sx + hy(k)*sy + hz(k)*sz

    return hamiltonian_internal(k)


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
        a, a_vec, b, b_vec, c, c_vec, d, d_vec = get_eigen(k)
        band_structure.write("{} {} {} {} {} {} {} {} {}\n"
                             .format(count, a, b, c, d,
                                     abs(a_vec[0])**2, abs(b_vec[0])**2, abs(c_vec[0])**2, abs(d_vec[0])**2))

    for i in range(nk):

        count = count+1
        k = GA+(MM-GA)*float(i)/float(nk-1)
        k = np.matmul(k, bvec)
        count = count+1
        a, a_vec, b, b_vec, c, c_vec, d, d_vec = get_eigen(k)
        band_structure.write("{} {} {} {} {} {} {} {} {}\n"
                             .format(count, a, b, c, d,
                                     abs(a_vec[0])**2, abs(b_vec[0])**2, abs(c_vec[0])**2, abs(d_vec[0])**2))

    for i in range(nk):

        count = count+1
        k = MM+(K2-MM)*float(i)/float(nk-1)
        k = np.matmul(k, bvec)
        count = count+1
        a, a_vec, b, b_vec, c, c_vec, d, d_vec = get_eigen(k)
        band_structure.write("{} {} {} {} {} {} {} {} {}\n"
                             .format(count, a, b, c, d,
                                     abs(a_vec[0])**2, abs(b_vec[0])**2, abs(c_vec[0])**2, abs(d_vec[0])**2))

    ############################
    # Plot the Berry Curvature #
    ############################

    open("berry_curvature.txt", "w")
    berry_curvature = open("berry_curvature.txt", "a", buffering=1)

    min_x = -np.pi
    max_x = np.pi
    samples_x = 100
    min_y = -np.pi
    max_y = np.pi
    samples_y = 100
    delta_alpha = (max_x - min_x) / samples_x
    delta_beta = (max_y - min_y) / samples_y

    tot_ber_curv_0, tot_ber_curv_1, tot_ber_curv_2, tot_ber_curv_3 = 0, 0, 0, 0

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
            ber_curv_2 = berry_curv(eigvecs[2], eigvecs_alpha[2], eigvecs_beta[2], eigvecs_alpha_beta[2])
            ber_curv_3 = berry_curv(eigvecs[3], eigvecs_alpha[3], eigvecs_beta[3], eigvecs_alpha_beta[3])

            tot_ber_curv_0 += ber_curv_0
            tot_ber_curv_1 += ber_curv_1
            tot_ber_curv_2 += ber_curv_2
            tot_ber_curv_3 += ber_curv_3

            print(kx, ky, ber_curv_0, ber_curv_1, ber_curv_2, ber_curv_3)
            berry_curvature.write("{} {} {} {}\n".format(kx, ky, ber_curv_0, ber_curv_1, ber_curv_2, ber_curv_3))

        berry_curvature.write("\n")

    print("")
    print("Chern number (0) = ", tot_ber_curv_0 / (2 * np.pi))
    print("Chern number (1) = ", tot_ber_curv_1 / (2 * np.pi))
    print("Chern number (2) = ", tot_ber_curv_2 / (2 * np.pi))
    print("Chern number (3) = ", tot_ber_curv_3 / (2 * np.pi))

