################################
# GRAPHENE IN A MAGNETIC FIELD #
################################

import numpy as np
import matplotlib.pyplot as plt
import sys
from itertools import product

# variables to set (when changing the Hamiltonian)

num_bands = 2
band_structure_file_name = "2D_graphene_band_structure.dat"
A_sites = [0]
wilson_loop_file_name = "graphene_wilson_loop.dat"
energies_file_name = "graphene_energies.dat"

# lattice vectors (normalized)

a1 = np.array([np.sqrt(3), 0])
a2 = (1/2) * np.array([np.sqrt(3), 3])
avec = np.vstack((a1, a2))

# reciprocal lattice vectors

b1 = (2.*np.pi) * np.array([1/np.sqrt(3), -1/3])
b2 = (2.*np.pi) * np.array([0, 2/3])
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
tau[1, :] = (1/3)*np.array([0, 0])

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
    t = 1

    delta = np.zeros((3, 2))
    # * BUG is here (I think)
    # * this makes crystal structure to be non-equivalent to the in convetion_II script
    # delta[0, :] = (2/3)*avec[0, :] - (1/3)*avec[1, :]
    # delta[1, :] = (-1/3)*avec[0, :] + (2/3)*avec[1, :]
    # delta[2, :] = (-1/3)*(avec[0, :] + avec[1, :])
    # * I change it:
    delta[0, :] = - (1 / 3) * avec[0, :] + (2 / 3) * avec[1, :]
    delta[1, :] = -(1 / 3) * avec[0, :] - (1 / 3) * avec[1, :]
    delta[2, :] = (2 / 3) * avec[0, :] - (1 / 3) * avec[1, :]

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
    for i in range(1, 3):
        f += - t * np.exp(1j * k.dot(delta[i, :]))

    phi = -0.3 * np.pi

    Hamiltonian[0][0] = 0
    Hamiltonian[0][1] = - t * np.exp(1j * k.dot(delta[0, :])) - t * np.exp(1j * k.dot(delta[1, :])) * np.exp(1j * phi) - t * np.exp(1j * k.dot(delta[2, :])) * np.exp(-1j * phi)
    Hamiltonian[1][0] = np.conj(- t * np.exp(1j * k.dot(delta[0, :])) - t * np.exp(1j * k.dot(delta[1, :])) * np.exp(1j * phi) - t * np.exp(1j * k.dot(delta[2, :])) * np.exp(-1j * phi))
    Hamiltonian[1][1] = 0

    # Hamiltonian[0][3] = np.conj(-t * np.exp(1j * k.dot(delta[0, :])))
    # Hamiltonian[1][2] = np.conj(-t * np.exp(1j * k.dot(delta[0, :])))
    # Hamiltonian[2][1] = -t * np.exp(1j * k.dot(delta[0, :]))
    # Hamiltonian[3][0] = -t * np.exp(1j * k.dot(delta[0, :]))
    #
    # Hamiltonian[2][2] = 0
    # Hamiltonian[2][3] = np.conj(- t * np.exp(1j * k.dot(delta[1, :])) + - t * np.exp(1j * k.dot(delta[2, :])))
    # Hamiltonian[3][2] = - t * np.exp(1j * k.dot(delta[1, :])) + - t * np.exp(1j * k.dot(delta[2, :]))
    # Hamiltonian[3][3] = 0

    return Hamiltonian


def berry_curv(ev, ev_alpha, ev_beta, ev_alpha_beta):

    bc = - np.imag(np.log(np.conj(ev).dot(ev_alpha) * np.conj(ev_alpha).dot(ev_alpha_beta)
                          * np.conj(ev_alpha_beta).dot(ev_beta) * np.conj(ev_beta).dot(ev)))

    return bc


def bloch_factor(k):

    bloch_factor = np.zeros(num_bands, dtype=np.complex128)

    # using the A, B, A, B,... Hamiltonian convention
    for component in range(0, num_bands-1, 2):
        bloch_factor[component] = np.exp(-1j * k.dot(np.matmul(tau[0, :], avec)))
        bloch_factor[component+1] = np.exp(-1j * k.dot(np.matmul(tau[1, :], avec)))

    return bloch_factor


def amp(vec, character):

    amplitude = 0

    if character == 'site':
        for component in A_sites:
            amplitude += abs(vec[component]) ** 2

    return amplitude


def write_eigensystem(file, count, k, eigval):

    eigvec = np.zeros((num_bands, num_bands), dtype=np.complex128)

    eigvals, eigvecs = np.linalg.eig(hamiltonian(k))

    idx = np.argsort(eigvals)[::-1]

    for band in range(num_bands):
        eigval[band, count] = np.real(eigvals[idx[band]])
        eigvec[:, band] = eigvecs[:, idx[band]]

    file.write("{} ".format(count))
    for band in range(num_bands):
        file.write("{} ".format(eigval[band, count]))
    for band in range(num_bands):
        if band == num_bands - 1:
            file.write("{}\n".format(amp(eigvec[:, band], 'site')))
        else:
            file.write("{} ".format(amp(eigvec[:, band], 'site')))

    return eigval[:, count]


def eval_eigensystem(count, k, eigval, B1, B2, B3, B4, B5, B6, B7, B8):

    eigvec = np.zeros((num_bands, num_bands), dtype=np.complex128)

    eigvals, eigvecs = np.linalg.eig(hamiltonian(k, B1, B2, B3, B4, B5, B6, B7, B8))

    idx = np.argsort(eigvals)[::-1]

    for band in range(num_bands):
        eigval[band, count] = np.real(eigvals[idx[band]])
        eigvec[:, band] = eigvecs[:, idx[band]]

    return eigval[:, count]


if __name__ == '__main__':

    parameter_search = False

    if parameter_search == False:

        ##########################
        # Plot 2D band structure #
        ##########################

        open(band_structure_file_name, "w")
        band_structure_2d = open(band_structure_file_name, "a", buffering=1)

        count, nk = 0, 30

        eigval_bands = np.zeros((num_bands, 3*nk))

        for i in range(nk):
            k = K1 - K1 * float(i) / float(nk - 1)
            k = np.matmul(k, bvec)
            eigval_bands[:, count] = write_eigensystem(band_structure_2d, count, k, eigval_bands)
            count += 1

        for i in range(nk):
            k = GA + (MM - GA) * float(i) / float(nk - 1)
            k = np.matmul(k, bvec)
            eigval_bands[:, count] = write_eigensystem(band_structure_2d, count, k, eigval_bands)
            count += 1

        for i in range(nk):
            k = MM + (K2 - MM) * float(i) / float(nk - 1)
            k = np.matmul(k, bvec)
            eigval_bands[:, count] = write_eigensystem(band_structure_2d, count, k, eigval_bands)
            count += 1

        ###############################
        # Calculate the Chern numbers #
        ###############################

        samples_x, samples_y = 101, 101
        max_idx_x, max_idx_y = samples_x - 1, samples_y - 1

        energy_matrix = np.zeros((num_bands, samples_x, samples_y))
        u_matrix = np.zeros((num_bands, num_bands, samples_x, samples_y), dtype=np.complex128)
        chern_numbers = np.zeros(num_bands)

        for idx_x in range(samples_x):
            for idx_y in range(samples_y):

                frac_kx, frac_ky = idx_x/max_idx_x, idx_y/max_idx_y

                # wavevector used for u (depends on boundary conditions)
                if idx_x == max_idx_x and idx_y == max_idx_y:
                    k_u = np.matmul(np.array([0, 0]), bvec)
                elif idx_x == max_idx_x:
                    k_u = np.matmul(np.array([0, frac_ky]), bvec)
                elif idx_y == max_idx_y:
                    k_u = np.matmul(np.array([frac_kx, 0]), bvec)
                else:
                    k_u = np.matmul(np.array([frac_kx, frac_ky]), bvec)

                eigval, u_eigvec = np.linalg.eig(hamiltonian(k_u))
                idx0 = np.argsort(eigval)[::-1]
                for band in range(num_bands):
                    energy_matrix[band, idx_x, idx_y] = np.real(eigval[idx0[band]])
                    u_matrix[:, band, idx_x, idx_y] = u_eigvec[:, idx0[band]]

        berry_flux_matrix = np.zeros((num_bands, samples_x-1, samples_y-1))

        for band in range(num_bands):
            for idx_x in range(max_idx_x):
                for idx_y in range(max_idx_y):
                    berry_flux_matrix[band, idx_x, idx_y] = berry_curv(u_matrix[:, band, idx_x, idx_y],
                                                                       u_matrix[:, band, idx_x + 1, idx_y],
                                                                       u_matrix[:, band, idx_x, idx_y + 1],
                                                                       u_matrix[:, band, idx_x + 1, idx_y + 1])

        print("")
        for band in range(num_bands):
            chern_numbers[band] = np.sum(berry_flux_matrix[band, :, :]) / (2 * np.pi)
            print("Chern number ( band", band, ") = ", chern_numbers[band])

        ######################################################
        # Berry fluxes along Wilson loops in the y direction #
        ######################################################

        open(wilson_loop_file_name, "w")
        wilson_loop = open(wilson_loop_file_name, "a", buffering=1)

        wl_berry_flux = np.ones((num_bands, samples_x), dtype=np.complex128)

        for band in range(num_bands):
            for idx_x in range(samples_x):
                for idx_y in range(samples_y - 1):
                    wl_berry_flux[band, idx_x] *= \
                        np.conj(u_matrix[:, band, idx_x, idx_y]).dot(u_matrix[:, band, idx_x, idx_y+1])
                wilson_loop.write("{} {}\n"
                                  .format(idx_x/max_idx_x, -(1/(2*np.pi)) * np.imag(np.log(wl_berry_flux[band, idx_x]))))
            wilson_loop.write("\n\n")

        #####################
        # 2D electronic DOS #
        #####################

        open(energies_file_name, "w")
        eigenenergies = open(energies_file_name, "a", buffering=1)

        for band in range(num_bands):
            for idx_x in range(samples_x):
                for idx_y in range(samples_y):
                    eigenenergies.write("{}\n".format(energy_matrix[band, idx_x, idx_y]))
            wilson_loop.write("\n\n")

    elif parameter_search == True:

        ####################################
        # Report non-trivial Chern numbers #
        ####################################

        open("non_trivial_chern.dat", "w")
        non_trivial_chern = open("non_trivial_chern.dat", "a", buffering=1)

        for Bx1, Bx2, Bx3, Bx4, By1, By2, By3, By4 in product([-1, 0, 1], repeat=8):

            # band_sep = 1.5
            # Bz1, Bz2, Bz3, Bz4, Bz5, Bz6, Bz7, Bz8 = (7/2)*band_sep, (5/2)*band_sep, (3/2)*band_sep, (1/2)*band_sep, \
            #                                          (-1/2)*band_sep, (-3/2)*band_sep, (-5/2)*band_sep, (-7/2)*band_sep

            print(Bx1, Bx2, Bx3, Bx4, By1, By2, By3, By4)

            ###################################################

            count, nk = 0, 30

            eigval_bands = np.zeros((num_bands, 3 * nk))

            for i in range(nk):
                k = K1 - K1 * float(i) / float(nk - 1)
                k = np.matmul(k, bvec)
                eigval_bands[:, count] = eval_eigensystem(count, k, eigval_bands, Bx1, Bx2, Bx3, Bx4, By1, By2, By3, By4)
                count += 1

            for i in range(nk):
                k = GA + (MM - GA) * float(i) / float(nk - 1)
                k = np.matmul(k, bvec)
                eigval_bands[:, count] = eval_eigensystem(count, k, eigval_bands, Bx1, Bx2, Bx3, Bx4, By1, By2, By3, By4)
                count += 1

            for i in range(nk):
                k = MM + (K2 - MM) * float(i) / float(nk - 1)
                k = np.matmul(k, bvec)
                eigval_bands[:, count] = eval_eigensystem(count, k, eigval_bands, Bx1, Bx2, Bx3, Bx4, By1, By2, By3, By4)
                count += 1

            num_degen = 0

            minima = np.zeros(num_bands)
            maxima = np.zeros(num_bands)

            for band in range(num_bands):
                minima[band] = np.min(eigval_bands[band, :])
                maxima[band] = np.max(eigval_bands[band, :])

            for ref_band in range(num_bands):
                for band in range(num_bands):
                    if ref_band != band:
                        if minima[ref_band] <= minima[band] <= maxima[ref_band] or \
                                minima[ref_band] <= maxima[band] <= maxima[ref_band]:
                            # print("Band", ref_band, "is degenerate with band", band)
                            num_degen += 1

            if num_degen > 0:
                continue

            ###################################################

            samples_x, samples_y = 101, 101
            max_idx_x, max_idx_y = samples_x - 1, samples_y - 1

            u_matrix = np.zeros((num_bands, num_bands, samples_x, samples_y), dtype=np.complex128)
            chern_numbers = np.zeros(num_bands)

            for idx_x in range(samples_x):
                for idx_y in range(samples_y):

                    frac_kx, frac_ky = idx_x / max_idx_x, idx_y / max_idx_y

                    # wavevector used for u (depends on boundary conditions)
                    if idx_x == max_idx_x and idx_y == max_idx_y:
                        k_u = np.matmul(np.array([0, 0]), bvec)
                    elif idx_x == max_idx_x:
                        k_u = np.matmul(np.array([0, frac_ky]), bvec)
                    elif idx_y == max_idx_y:
                        k_u = np.matmul(np.array([frac_kx, 0]), bvec)
                    else:
                        k_u = np.matmul(np.array([frac_kx, frac_ky]), bvec)

                    eigval, u_eigvec = np.linalg.eig(hamiltonian(k_u, Bx1, Bx2, Bx3, Bx4, By1, By2, By3, By4))
                    idx0 = np.argsort(eigval)[::-1]
                    for band in range(num_bands):
                        u_matrix[:, band, idx_x, idx_y] = u_eigvec[:, idx0[band]]

            berry_flux_matrix = np.zeros((num_bands, samples_x - 1, samples_y - 1))

            for band in range(num_bands):
                for idx_x in range(max_idx_x):
                    for idx_y in range(max_idx_y):
                        berry_flux_matrix[band, idx_x, idx_y] = berry_curv(u_matrix[:, band, idx_x, idx_y],
                                                                           u_matrix[:, band, idx_x + 1, idx_y],
                                                                           u_matrix[:, band, idx_x, idx_y + 1],
                                                                           u_matrix[:, band, idx_x + 1, idx_y + 1])

            print("")
            for band in range(num_bands):
                chern_numbers[band] = np.sum(berry_flux_matrix[band, :, :]) / (2 * np.pi)

            ###################################################

            non_trivial_chern = 0

            for band in range(num_bands):
                if abs(chern_numbers[band]) > 0.5:
                    #print("Band", band, "has a non-trivial Chern number")
                    non_trivial_chern += 1

            if non_trivial_chern > 0:
                print("Non-trivial Chern number found!")
                non_trivial_chern.write("{} {} {} {} {} {} {} {}".format(Bx1, Bx2, Bx3, Bx4, By1, By2, By3, By4))
