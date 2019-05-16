import numpy as np


def hamiltonian(kx, ky, a, b, c, d):

    ##################
    # Initialization #
    ##################

    k = np.array([kx, ky])
    lat_const = 1
    t1 = 1
    t2 = 0.1
    t2dash = 0.4

    delta = np.array([[(lat_const / 2), (lat_const / 2) * np.sqrt(3)],
                      [(lat_const / 2), -(lat_const / 2) * np.sqrt(3)],
                      [-lat_const, 0]])

    fifthNN = np.array([[3 * lat_const, 0],
                        [(3 / 2) * lat_const, (3 / 2) * lat_const * np.sqrt(3)],
                        [-(3 / 2) * lat_const, (3 / 2) * lat_const * np.sqrt(3)],
                        [-(3 / 2) * lat_const, -(3 / 2) * lat_const * np.sqrt(3)],
                        [(3 / 2) * lat_const, -(3 / 2) * lat_const * np.sqrt(3)],
                        [-3 * lat_const, 0]])

    ######################
    # Define Hamiltonian #
    ######################

    Hamiltonian = np.zeros((8, 8), dtype=complex)

    f = 0
    for i in range(0, 3):
        f += t1 * np.exp(1j * k.dot(delta[i]))
    for i in range(0, 6):
        f += t2 * np.exp(1j * k.dot(fifthNN[i]))

    Hamiltonian[0][0] = a
    Hamiltonian[0][1] = f
    Hamiltonian[1][0] = np.conj(f)
    Hamiltonian[1][1] = -a

    Hamiltonian[2][2] = b
    Hamiltonian[2][3] = f
    Hamiltonian[3][2] = np.conj(f)
    Hamiltonian[3][3] = -b

    Hamiltonian[4][4] = c
    Hamiltonian[4][5] = f
    Hamiltonian[5][4] = np.conj(f)
    Hamiltonian[5][5] = -c

    Hamiltonian[6][6] = d
    Hamiltonian[6][7] = f
    Hamiltonian[7][6] = np.conj(f)
    Hamiltonian[7][7] = -d

    # xi = 0
    # for i in range(0, 6):
    #     xi += t2dash * np.exp(1j * k.dot(fifthNN[i]))
    #
    # Hamiltonian[0][3] = xi
    # Hamiltonian[1][2] = -np.conj(xi)
    # Hamiltonian[2][1] = -xi
    # Hamiltonian[3][0] = np.conj(xi)
    #
    # Hamiltonian[4][7] = xi
    # Hamiltonian[5][6] = -np.conj(xi)
    # Hamiltonian[6][5] = -xi
    # Hamiltonian[7][4] = np.conj(xi)

    return Hamiltonian


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

    a, b, c, d = 10, 20, 30, 40
    tot_ber_curv_a_p, tot_ber_curv_a_m, tot_ber_curv_b_p, tot_ber_curv_b_m,\
        tot_ber_curv_c_p, tot_ber_curv_c_m, tot_ber_curv_d_p, tot_ber_curv_d_m = 0, 0, 0, 0, 0, 0, 0, 0

    for kx in np.linspace(min_x, max_x, samples_x):
        for ky in np.linspace(min_y, max_y, samples_y):

            eigval, eigvec = np.linalg.eig(hamiltonian(kx, ky, a, b, c, d))
            idx1 = np.argsort(eigval)[::-1]
            eigval_a_p = np.real(eigval[idx1])[3]
            eigval_a_m = np.real(eigval[idx1])[4]
            eigvec_a_p = eigvec[:, idx1][3]
            eigvec_a_m = eigvec[:, idx1][4]
            eigval_b_p = np.real(eigval[idx1])[2]
            eigval_b_m = np.real(eigval[idx1])[5]
            eigvec_b_p = eigvec[:, idx1][2]
            eigvec_b_m = eigvec[:, idx1][5]
            eigval_c_p = np.real(eigval[idx1])[1]
            eigval_c_m = np.real(eigval[idx1])[6]
            eigvec_c_p = eigvec[:, idx1][1]
            eigvec_c_m = eigvec[:, idx1][6]
            eigval_d_p = np.real(eigval[idx1])[0]
            eigval_d_m = np.real(eigval[idx1])[7]
            eigvec_d_p = eigvec[:, idx1][0]
            eigvec_d_m = eigvec[:, idx1][7]
            # alpha
            eigval_alpha, eigvec_alpha = np.linalg.eig(hamiltonian(kx + delta_alpha, ky, a, b, c, d))
            idx2 = np.argsort(eigval_alpha)[::-1]
            eigval_a_p_alpha = np.real(eigval_alpha[idx2])[3]
            eigval_a_m_alpha = np.real(eigval_alpha[idx2])[4]
            eigvec_a_p_alpha = eigvec_alpha[:, idx2][3]
            eigvec_a_m_alpha = eigvec_alpha[:, idx2][4]
            eigval_b_p_alpha = np.real(eigval_alpha[idx2])[2]
            eigval_b_m_alpha = np.real(eigval_alpha[idx2])[5]
            eigvec_b_p_alpha = eigvec_alpha[:, idx2][2]
            eigvec_b_m_alpha = eigvec_alpha[:, idx2][5]
            eigval_c_p_alpha = np.real(eigval_alpha[idx2])[1]
            eigval_c_m_alpha = np.real(eigval_alpha[idx2])[6]
            eigvec_c_p_alpha = eigvec_alpha[:, idx2][1]
            eigvec_c_m_alpha = eigvec_alpha[:, idx2][6]
            eigval_d_p_alpha = np.real(eigval_alpha[idx2])[0]
            eigval_d_m_alpha = np.real(eigval_alpha[idx2])[7]
            eigvec_d_p_alpha = eigvec_alpha[:, idx2][0]
            eigvec_d_m_alpha = eigvec_alpha[:, idx2][7]
            # beta
            eigval_beta, eigvec_beta = np.linalg.eig(hamiltonian(kx, ky + delta_beta, a, b, c, d))
            idx3 = np.argsort(eigval_beta)[::-1]
            eigval_a_p_beta = np.real(eigval_beta[idx3])[3]
            eigval_a_m_beta = np.real(eigval_beta[idx3])[4]
            eigvec_a_p_beta = eigvec_beta[:, idx3][3]
            eigvec_a_m_beta = eigvec_beta[:, idx3][4]
            eigval_b_p_beta = np.real(eigval_beta[idx3])[2]
            eigval_b_m_beta = np.real(eigval_beta[idx3])[5]
            eigvec_b_p_beta = eigvec_beta[:, idx3][2]
            eigvec_b_m_beta = eigvec_beta[:, idx3][5]
            eigval_c_p_beta = np.real(eigval_beta[idx3])[1]
            eigval_c_m_beta = np.real(eigval_beta[idx3])[6]
            eigvec_c_p_beta = eigvec_beta[:, idx3][1]
            eigvec_c_m_beta = eigvec_beta[:, idx3][6]
            eigval_d_p_beta = np.real(eigval_beta[idx3])[0]
            eigval_d_m_beta = np.real(eigval_beta[idx3])[7]
            eigvec_d_p_beta = eigvec_beta[:, idx3][0]
            eigvec_d_m_beta = eigvec_beta[:, idx3][7]
            # alpha_beta
            eigval_alpha_beta, eigvec_alpha_beta = np.linalg.eig(hamiltonian(kx + delta_alpha, ky + delta_beta, a, b, c, d))
            idx4 = np.argsort(eigval_alpha_beta)[::-1]
            eigval_a_p_alpha_beta = np.real(eigval_alpha_beta[idx4])[3]
            eigval_a_m_alpha_beta = np.real(eigval_alpha_beta[idx4])[4]
            eigvec_a_p_alpha_beta = eigvec_alpha_beta[:, idx4][3]
            eigvec_a_m_alpha_beta = eigvec_alpha_beta[:, idx4][4]
            eigval_b_p_alpha_beta = np.real(eigval_alpha_beta[idx4])[2]
            eigval_b_m_alpha_beta = np.real(eigval_alpha_beta[idx4])[5]
            eigvec_b_p_alpha_beta = eigvec_alpha_beta[:, idx4][2]
            eigvec_b_m_alpha_beta = eigvec_alpha_beta[:, idx4][5]
            eigval_c_p_alpha_beta = np.real(eigval_alpha_beta[idx4])[1]
            eigval_c_m_alpha_beta = np.real(eigval_alpha_beta[idx4])[6]
            eigvec_c_p_alpha_beta = eigvec_alpha_beta[:, idx4][1]
            eigvec_c_m_alpha_beta = eigvec_alpha_beta[:, idx4][6]
            eigval_d_p_alpha_beta = np.real(eigval_alpha_beta[idx4])[0]
            eigval_d_m_alpha_beta = np.real(eigval_alpha_beta[idx4])[7]
            eigvec_d_p_alpha_beta = eigvec_alpha_beta[:, idx4][0]
            eigvec_d_m_alpha_beta = eigvec_alpha_beta[:, idx4][7]

            ber_curv_a_p = berry_curv(eigvec_a_p, eigvec_a_p_alpha, eigvec_a_p_beta, eigvec_a_p_alpha_beta)
            ber_curv_a_m = berry_curv(eigvec_a_m, eigvec_a_m_alpha, eigvec_a_m_beta, eigvec_a_m_alpha_beta)
            ber_curv_b_p = berry_curv(eigvec_b_p, eigvec_b_p_alpha, eigvec_b_p_beta, eigvec_b_p_alpha_beta)
            ber_curv_b_m = berry_curv(eigvec_b_m, eigvec_b_m_alpha, eigvec_b_m_beta, eigvec_b_m_alpha_beta)
            ber_curv_c_p = berry_curv(eigvec_c_p, eigvec_c_p_alpha, eigvec_c_p_beta, eigvec_c_p_alpha_beta)
            ber_curv_c_m = berry_curv(eigvec_c_m, eigvec_c_m_alpha, eigvec_c_m_beta, eigvec_c_m_alpha_beta)
            ber_curv_d_p = berry_curv(eigvec_d_p, eigvec_d_p_alpha, eigvec_d_p_beta, eigvec_d_p_alpha_beta)
            ber_curv_d_m = berry_curv(eigvec_d_m, eigvec_d_m_alpha, eigvec_d_m_beta, eigvec_d_m_alpha_beta)

            tot_ber_curv_a_p += ber_curv_a_p
            tot_ber_curv_a_m += ber_curv_a_m
            tot_ber_curv_b_p += ber_curv_b_p
            tot_ber_curv_b_m += ber_curv_b_m
            tot_ber_curv_c_p += ber_curv_c_p
            tot_ber_curv_c_m += ber_curv_c_m
            tot_ber_curv_d_p += ber_curv_d_p
            tot_ber_curv_d_m += ber_curv_d_m

            print(kx, ky, ber_curv_a_p, ber_curv_a_m, ber_curv_b_p, ber_curv_b_m,
                  ber_curv_c_p, ber_curv_c_m, ber_curv_d_p, ber_curv_d_m)
            band_structure.write("%.15f\t%.15f\t%.15f\t%.15f\t%.15f\t%.15f\t"
                                 "%.15f\t%.15f\t%.15f\t%.15f\n" %
                                 (kx, ky, eigval_a_p, eigval_a_m, eigval_b_p, eigval_b_m,
                                  eigval_c_p, eigval_c_m, eigval_d_p, eigval_d_m))
            berry_curvature.write("%.15f\t%.15f\t%.15f\t%.15f\t%.15f\t%.15f\t"
                                  "%.15f\t%.15f\t%.15f\t%.15f\n" %
                                  (kx, ky, ber_curv_a_p, ber_curv_a_m, ber_curv_b_p, ber_curv_b_m,
                                   ber_curv_c_p, ber_curv_c_m, ber_curv_d_p, ber_curv_d_m))
        print(" ")
        band_structure.write("\n")
        berry_curvature.write("\n")

    print("Chern number (a+) = ", tot_ber_curv_a_p / (2 * np.pi))
    print("Chern number (a-) = ", tot_ber_curv_a_m / (2 * np.pi))
    print("Chern number (b+) = ", tot_ber_curv_a_p / (2 * np.pi))
    print("Chern number (b-) = ", tot_ber_curv_a_m / (2 * np.pi))
    print("Chern number (c+) = ", tot_ber_curv_a_p / (2 * np.pi))
    print("Chern number (c-) = ", tot_ber_curv_a_m / (2 * np.pi))
    print("Chern number (d+) = ", tot_ber_curv_a_p / (2 * np.pi))
    print("Chern number (d-) = ", tot_ber_curv_a_m / (2 * np.pi))
