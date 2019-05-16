"""Calculate the Chern number of the bands in Fig. 1.c) of
"Chareterization and stability of a fermionic nu=1/3 fractional Chern insulator"""

import numpy as np


def hamiltonian(kx, ky):

    ##################
    # Initialization #
    ##################

    t = 1
    t2 = np.sqrt(129)/36
    M = 0
    phi = np.arccos(3*np.sqrt(3/43))

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

    def h0(kx, ky):

        k = np.array([kx, ky])

        h0 = 0
        for i in range(3):
            h0 += 2*t2*np.cos(phi)*np.cos(k.dot(b[i]))

        return h0

    def hx(kx, ky):

        k = np.array([kx, ky])

        hx = 0
        for i in range(3):
            hx += t * np.cos(k.dot(a[i]))

        return hx

    def hy(kx, ky):

        k = np.array([kx, ky])

        hy = 0
        for i in range(3):
            hy += -t * np.sin(k.dot(a[i]))

        return hy

    def hz(kx, ky):

        k = np.array([kx, ky])

        hz = M
        for i in range(3):
            hz += 2 * t2 * np.sin(phi) * np.sin(k.dot(b[i]))

        return hz

    def hamiltonian(kx, ky):

        return h0(kx, ky)*s0 + hx(kx, ky)*sx + hy(kx, ky)*sy + hz(kx, ky)*sz

    return hamiltonian(kx, ky)


def berry_curv(ev, ev_alpha, ev_beta, ev_alpha_beta):

    bc = - 2 * np.imag(np.log(np.conj(ev).dot(ev_alpha) * np.conj(ev_alpha).dot(ev_alpha_beta)
                              * np.conj(ev_alpha_beta).dot(ev_beta) * np.conj(ev_beta).dot(ev)))

    return bc


if __name__ == '__main__':

    #######################
    # File initialization #
    #######################

    open("band_structure.txt", "w")
    band_structure = open("band_structure.txt", "a", buffering=1)
    open("berry_curvature.txt", "w")
    berry_curvature = open("berry_curvature.txt", "a", buffering=1)

    #######################
    # Grid discretization #
    #######################

    min_x = -np.pi
    max_x = np.pi
    samples_x = 100
    min_y = -np.pi
    max_y = np.pi
    samples_y = 100
    delta_alpha = (max_x - min_x) / samples_x
    delta_beta = (max_y - min_y) / samples_y

    ##################
    # Initialization #
    ##################

    tot_ber_curv_p = 0
    tot_ber_curv_m = 0

    for kx in np.linspace(min_x, max_x, samples_x):
        for ky in np.linspace(min_y, max_y, samples_y):

            eigval, eigvec = np.linalg.eig(hamiltonian(kx, ky))
            idx1 = np.argsort(eigval)[::-1]
            eigval_p = np.real(eigval[idx1])[0]
            eigval_m = np.real(eigval[idx1])[1]
            eigvec_p = eigvec[:, idx1][0]
            eigvec_m = eigvec[:, idx1][1]
            # alpha
            eigval_alpha, eigvec_alpha = np.linalg.eig(hamiltonian(kx+delta_alpha, ky))
            idx2 = np.argsort(eigval_alpha)[::-1]
            eigval_p_alpha = np.real(eigval_alpha[idx2])[0]
            eigval_m_alpha = np.real(eigval_alpha[idx2])[1]
            eigvec_p_alpha = eigvec_alpha[:, idx2][0]
            eigvec_m_alpha = eigvec_alpha[:, idx2][1]
            # beta
            eigval_beta, eigvec_beta = np.linalg.eig(hamiltonian(kx, ky+delta_beta))
            idx3 = np.argsort(eigval_beta)[::-1]
            eigval_p_beta = np.real(eigval_beta[idx3])[0]
            eigval_m_beta = np.real(eigval_beta[idx3])[1]
            eigvec_p_beta = eigvec_beta[:, idx3][0]
            eigvec_m_beta = eigvec_beta[:, idx3][1]
            # alpha_beta
            eigval_alpha_beta, eigvec_alpha_beta = np.linalg.eig(hamiltonian(kx+delta_alpha, ky+delta_beta))
            idx4 = np.argsort(eigval_alpha_beta)[::-1]
            eigval_p_alpha_beta = np.real(eigval_alpha_beta[idx4])[0]
            eigval_m_alpha_beta = np.real(eigval_alpha_beta[idx4])[1]
            eigvec_p_alpha_beta = eigvec_alpha_beta[:, idx4][0]
            eigvec_m_alpha_beta = eigvec_alpha_beta[:, idx4][1]

            ber_curv_p = berry_curv(eigvec_p, eigvec_p_alpha, eigvec_p_beta, eigvec_p_alpha_beta)
            ber_curv_m = berry_curv(eigvec_m, eigvec_m_alpha, eigvec_m_beta, eigvec_m_alpha_beta)

            tot_ber_curv_p += ber_curv_p
            tot_ber_curv_m += ber_curv_m

            print(kx, ky, ber_curv_p, ber_curv_m)
            band_structure.write("%.15f\t%.15f\t%.15f\t%.15f\n" % (kx, ky, eigval_p, eigval_m))
            berry_curvature.write("%.15f\t%.15f\t%.15f\t%.15f\n" % (kx, ky, ber_curv_p, ber_curv_m))
        print(" ")
        band_structure.write("\n")
        berry_curvature.write("\n")

    print("Chern number (+) = ", tot_ber_curv_p / (2 * np.pi))
    print("Chern number (-) = ", tot_ber_curv_m / (2 * np.pi))
