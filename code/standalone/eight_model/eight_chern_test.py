########################################################
# CHERN NUMBER FOR THE EIGHT BAND MODEL (CONVENTION I) #
########################################################

import numpy as np
import matplotlib.pyplot as plt
import sys
from itertools import product
import matplotlib.gridspec as gridspec

# variables to set (when changing the Hamiltonian)

num_bands = 8
band_structure_file_name = "2D_eight_band_structure.dat"
A_sites = [0, 2, 4, 6]
x_sites = [0, 1, 4, 5]
up_sites = [0, 1, 2, 3]
wilson_loop_file_name = "eight_wilson_loop.dat"
energies_file_name = "eight_energies.dat"

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

    # values from arXiv:1805.06819
    t1 = 0.331
    t2 = np.sqrt(129)/36 * t1  # 0.01
    t2dash = 0.097
    # Haldane hoppings in each valley
    M = 0
    phi = np.arccos(3*np.sqrt(3/43))
    tsec2 = -0.097

    delta = np.zeros((3, 2))
    delta[0, :] = (1 / 3) * avec[0, :] + (1 / 3) * avec[1, :]
    delta[1, :] = (-2 / 3) * avec[0, :] + (1 / 3) * avec[1, :]
    delta[2, :] = (1 / 3) * avec[0, :] + (-2 / 3) * avec[1, :]

    secondNN = np.zeros((6, 2))
    # positive direction for A sites / negative direction for B sites
    secondNN[0, :] = avec[0, :]
    secondNN[1, :] = -avec[0, :] + avec[1, :]
    secondNN[2, :] = -avec[1, :]
    # negative direction for A sites / positive direction for B sites
    secondNN[3, :] = avec[1, :]
    secondNN[4, :] = -avec[0, :]
    secondNN[5, :] = avec[0, :] - avec[1, :]

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

    Hamiltonian = np.zeros((num_bands, num_bands), dtype=np.complex128)

    f = 0
    for i in range(3):
        f += t1 * np.exp(1j * k.dot(delta[i, :]))
    f2 = 0
    for i in range(0, 3):
        f2 += t2 * np.exp(1j * k.dot(secondNN[i, :]))
    # f2_2 = 0
    # for i in range(3, 6):
    #     f2_2 += t2 * np.exp(1j * k.dot(secondNN[i, :]))
    xi = 0
    for i in range(3):
        xi += t2dash * np.exp(1j * k.dot(fifthNN[i, :]))
    ###
    fsec1 = 0
    for i in range(0, 3):
        fsec1 += - tsec2 * np.exp(1j * k.dot(secondNN[i, :]))
    fsec2 = 0
    for i in range(3, 6):
        fsec2 += - tsec2 * np.exp(1j * k.dot(secondNN[i, :]))

    # complex hoppings

    t2_phase = np.exp(1j * phi)
    t2dash_phase = 1  # np.exp(1j * phi)

    # spin up block

    Hamiltonian[0][0] = t2_phase * f2 + np.conj(t2_phase) * np.conj(f2)
    Hamiltonian[0][1] = f
    Hamiltonian[1][0] = np.conj(f)
    Hamiltonian[1][1] = t2_phase * np.conj(f2) + np.conj(t2_phase) * f2

    Hamiltonian[0][2] = t2dash_phase * xi - np.conj(t2dash_phase) * np.conj(xi)
    Hamiltonian[1][3] = np.conj(t2dash_phase) * xi - t2dash_phase * np.conj(xi)

    Hamiltonian[2][0] = np.conj(Hamiltonian[0][2])
    Hamiltonian[3][1] = np.conj(Hamiltonian[1][3])

    Hamiltonian[2][2] = Hamiltonian[0][0]
    Hamiltonian[2][3] = Hamiltonian[0][1]
    Hamiltonian[3][2] = Hamiltonian[1][0]
    Hamiltonian[3][3] = Hamiltonian[1][1]

    # spin down block (the minimal model is spin degenerate)

    Hamiltonian[4][4] = Hamiltonian[0][0]
    Hamiltonian[4][5] = Hamiltonian[0][1]
    Hamiltonian[5][4] = Hamiltonian[1][0]
    Hamiltonian[5][5] = Hamiltonian[1][1]

    Hamiltonian[4][6] = Hamiltonian[0][2]
    Hamiltonian[5][7] = Hamiltonian[1][3]

    Hamiltonian[6][4] = Hamiltonian[2][0]
    Hamiltonian[7][5] = Hamiltonian[3][1]

    Hamiltonian[6][6] = Hamiltonian[2][2]
    Hamiltonian[6][7] = Hamiltonian[2][3]
    Hamiltonian[7][6] = Hamiltonian[3][2]
    Hamiltonian[7][7] = Hamiltonian[3][3]

    # orbital Zeeman contribution (Btau)

    Btau_x, Btau_y, Btau_z = 0, 0, 3  # 3

    for i in [0, 1, 4, 5]:
        Hamiltonian[i][i] += Btau_z
        Hamiltonian[i + 2][i + 2] -= Btau_z
        Hamiltonian[i][i + 2] += (Btau_x - 1j * Btau_y)
        Hamiltonian[i + 2][i] += (Btau_x + 1j * Btau_y)

    # spin Zeeman contribution (Bsigma)

    Bsigma_x, Bsigma_y, Bsigma_z = 0, 0, 6  # 6

    for i in range(4):

        Hamiltonian[i][i] += Bsigma_z
        Hamiltonian[i+4][i+4] -= Bsigma_z
        Hamiltonian[i+4][i] += (Bsigma_x - 1j*Bsigma_y)
        Hamiltonian[i][i+4] += (Bsigma_x + 1j*Bsigma_y)

    return Hamiltonian


def berry_curv(ev, ev_alpha, ev_beta, ev_alpha_beta):

    bc = - np.imag(np.log(np.conj(ev).dot(ev_alpha) * np.conj(ev_alpha).dot(ev_alpha_beta)
                          * np.conj(ev_alpha_beta).dot(ev_beta) * np.conj(ev_beta).dot(ev)))

    return bc


def amp(vec, character):

    amplitude = 0

    if character == 'site':
        for component in A_sites:
            amplitude += abs(vec[component]) ** 2
    elif character == 'orbital':
        for component in x_sites:
            amplitude += abs(vec[component]) ** 2
    elif character == 'spin':
        for component in up_sites:
            amplitude += abs(vec[component]) ** 2

    return amplitude


def write_eigensystem(file, count, k, eigval, amp_bands):

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
        file.write("{} ".format(amp(eigvec[:, band], 'site')))
        amp_bands[band, count] = amp(eigvec[:, band], 'site')
    for band in range(num_bands):
        file.write("{} ".format(amp(eigvec[:, band], 'orbital')))
    for band in range(num_bands):
        if band == num_bands - 1:
            file.write("{}\n".format(amp(eigvec[:, band], 'spin')))
        else:
            file.write("{} ".format(amp(eigvec[:, band], 'spin')))

    return eigval[:, count], amp_bands[:, count]


if __name__ == '__main__':

    ##########################
    # Plot 2D band structure #
    ##########################

    open(band_structure_file_name, "w")
    band_structure_2d = open(band_structure_file_name, "a", buffering=1)

    count, nk = 0, 30

    eigval_bands = np.zeros((num_bands, 3*nk))
    amp_bands = np.zeros((num_bands, 3*nk))

    for i in range(nk):
        k = K1 - K1 * float(i) / float(nk - 1)
        k = np.matmul(k, bvec)
        eigval_bands[:, count], amp_bands[:, count] = write_eigensystem(band_structure_2d, count, k, eigval_bands, amp_bands)
        count += 1

    for i in range(nk):
        k = GA + (MM - GA) * float(i) / float(nk - 1)
        k = np.matmul(k, bvec)
        eigval_bands[:, count], amp_bands[:, count] = write_eigensystem(band_structure_2d, count, k, eigval_bands, amp_bands)
        count += 1

    for i in range(nk):
        k = MM + (K2 - MM) * float(i) / float(nk - 1)
        k = np.matmul(k, bvec)
        eigval_bands[:, count], amp_bands[:, count] = write_eigensystem(band_structure_2d, count, k, eigval_bands, amp_bands)
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

    ##################################
    # Calculate the gap/width ratios #
    ##################################

    width = np.zeros(num_bands)

    for i in range(num_bands):
        width[i] = max(eigval_bands[i, :]) - min(eigval_bands[i, :])

    gap = np.zeros(num_bands-1)

    for i in range(num_bands-1):
        gap[i] = min(eigval_bands[i, :]) - max(eigval_bands[i+1, :])

    ratio = np.zeros(num_bands-1)

    for i in range(num_bands-1):
        ratio[i] = gap[i]/width[i+1]

    print("")
    for band in range(num_bands):
        chern_numbers[band] = np.sum(berry_flux_matrix[band, :, :]) / (2 * np.pi)
        print("Chern number ( band", band, ") = ", chern_numbers[band])
        if (band + 1) % 2 == 0:
            print("Ratio of", band - 1, "-", band, "gap to", band, "width = ", ratio[band - 1])
            print("")

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

    ##############
    # Final Plot #
    ##############

    fig = plt.figure()

    gs = gridspec.GridSpec(1, 3, width_ratios=[2, 1, 1])

    ax0 = plt.subplot(gs[0])
    # colors = ['red', 'blue', 'green', 'yellow', 'purple', 'orange', 'black', 'brown']
    markers = ["o", "^", "s", "D", "p", "P", "*", "X"]
    # fills = ['red', 'blue', 'green', 'yellow', 'none', 'none', 'none', 'none']

    for nb in range(num_bands):
        plt.scatter(np.linspace(0, 89, 90), eigval_bands[nb, :], c=amp_bands[nb, :], cmap=plt.cm.coolwarm,
                    marker=markers[nb])
    ax0.legend(['Ax↑ (C={})'.format(int(chern_numbers[0])), 'Bx↑ (C={})'.format(int(chern_numbers[1])),
                'Ay↑ (C={})'.format(int(chern_numbers[2])), 'By↑ (C={})'.format(int(chern_numbers[3])),
                'Ax↓ (C={})'.format(int(chern_numbers[4])), 'Bx↓ (C={})'.format(int(chern_numbers[5])),
                'Ay↓ (C={})'.format(int(chern_numbers[6])), 'By↓ (C={})'.format(int(chern_numbers[7]))],
               loc='upper center', bbox_to_anchor=(0.5, -0.05), ncol=4)
    ax0.set_title('Band Structure')
    ax0.set_ylabel('Energy / meV')
    ax0.axvline(30, color='k', linewidth=0.5)
    ax0.axvline(60, color='k', linewidth=0.5)
    ax0.axhline(0, color='k', linewidth=0.5, ls='--')
    plt.xlim((0, 89))
    plt.xticks([0, 30, 60, 89], ["K", "Γ", "M", "K"])
    plt.colorbar(label='$|c_\mathrm{A}|^2$')

########################################################################################################################

    ax1 = plt.subplot(gs[1], sharey=ax0)
    n, bins, patches = ax1.hist(np.ndarray.flatten(energy_matrix), 100,
                                density=False, facecolor='k', alpha=0.75, orientation='horizontal')
    plt.setp(ax1.get_yticklabels(), visible=False)
    ax1.set_title('Density of States')
    ax1.set_xlabel('# States')
    ax1.axhline(0, color='k', linewidth=0.5, ls='--')

    # plt.subplots_adjust(wspace=.0)
    gs.update(wspace=0)

########################################################################################################################

    left, bottom, width, height = [0.69, 0.58, 0.3, 0.3]
    ax2 = fig.add_axes([left, bottom, width, height])
    colors = ['red', 'blue', 'green', 'yellow', 'purple', 'orange', 'black', 'brown']
    markers = ["o", "^", "s", "D", "p", "P", "*", "X"]
    for nb in range(num_bands):
        ax2.scatter(np.linspace(0, 1, 101), -(1 / (2 * np.pi)) * np.imag(np.log(wl_berry_flux[nb, :])),
                    marker=markers[nb], facecolors='none', edgecolors=colors[nb])
    ax2.legend(['Ax↑ (C={})'.format(int(chern_numbers[0])), 'Bx↑ (C={})'.format(int(chern_numbers[1])),
                'Ay↑ (C={})'.format(int(chern_numbers[2])), 'By↑ (C={})'.format(int(chern_numbers[3])),
                'Ax↓ (C={})'.format(int(chern_numbers[4])), 'Bx↓ (C={})'.format(int(chern_numbers[5])),
                'Ay↓ (C={})'.format(int(chern_numbers[6])), 'By↓ (C={})'.format(int(chern_numbers[7]))],
               loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=2)
    plt.xlim((0, 1))
    plt.ylim((-0.5, 0.5))
    ax2.set_aspect('equal', adjustable='box')
    ax2.set_title('Berry flux around Wilson Loop')
    ax2.set_xlabel('$k_\mathrm{x} / \mathrm{b}_1$')
    ax2.set_ylabel('$\Sigma \mathrm{HWCC} / 2 \pi$')
    ax2.axhline(0, color='k', linewidth=0.5, ls='--')
    ax2.axvline(0.5, color='k', linewidth=0.5, ls='--')

########################################################################################################################

    left, bottom, width, height = [0.715, 0.11, 0.25, 0.25]
    ax3 = fig.add_axes([left, bottom, width, height])

    delta = np.zeros((3, 2))
    delta[0, :] = (1 / 3) * avec[0, :] + (1 / 3) * avec[1, :]
    delta[1, :] = (-2 / 3) * avec[0, :] + (1 / 3) * avec[1, :]
    delta[2, :] = (1 / 3) * avec[0, :] + (-2 / 3) * avec[1, :]

    secondNN = np.zeros((6, 2))
    # positive direction for A sites / negative direction for B sites
    secondNN[0, :] = avec[0, :]
    secondNN[1, :] = -avec[0, :] + avec[1, :]
    secondNN[2, :] = -avec[1, :]
    # negative direction for A sites / positive direction for B sites
    secondNN[3, :] = avec[1, :]
    secondNN[4, :] = -avec[0, :]
    secondNN[5, :] = avec[0, :] - avec[1, :]

    fifthNN = np.zeros((6, 2))
    fifthNN[0, :] = -avec[0, :] - avec[1, :]
    fifthNN[1, :] = -avec[0, :] + 2 * avec[1, :]
    fifthNN[2, :] = 2 * avec[0, :] - avec[1, :]
    fifthNN[3, :] = avec[0, :] + avec[1, :]
    fifthNN[4, :] = -2 * avec[0, :] + avec[1, :]
    fifthNN[5, :] = avec[0, :] - 2 * avec[1, :]

    for m in [-3, -2, -1, 0, 1, 2, 3]:
        for n in [-3, -2, -1, 0, 1, 2, 3]:
            xcoord = 0 + n*avec[0, 0]
            ycoord = 0 + m + n*avec[0, 1]
            for i in range(3):
                ax3.arrow(xcoord, ycoord, delta[i, 0], delta[i, 1], color='gray', width=0.0001, ls='--', alpha=0.5)

    ax3.scatter(delta[:, 0], delta[:, 1], color='blue')
    ax3.scatter(secondNN[0:3, 0], secondNN[0:3, 1], color='green')
    ax3.scatter(secondNN[3:6, 0], secondNN[3:6, 1], color='green', marker='x')
    ax3.scatter(fifthNN[0:3, 0], fifthNN[0:3, 1], color='red')
    ax3.scatter(fifthNN[3:6, 0], fifthNN[3:6, 1], color='red', marker='x')
    ax3.set_title('Real Space Lattice')
    ax3.set_xlabel('$x$')
    ax3.set_ylabel('$y$')

    ax3.arrow(0, 0, avec[0, 0], avec[0, 1], color='black', head_width=0.1, length_includes_head=True)
    ax3.arrow(0, 0, avec[1, 0], avec[1, 1], color='black', head_width=0.1, length_includes_head=True)
    ax3.text(avec[0, 0]/3, 0.4, "$\mathbf{a}_1$")
    ax3.text(avec[0, 0]/3, -0.5, "$\mathbf{a}_2$")

    ax3.set_aspect('equal', adjustable='box')

    plt.show()
