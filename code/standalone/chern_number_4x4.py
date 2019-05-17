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
    t2 = 0.01  # 1
    t2dash = 0.097

    delta = np.zeros((3,2))
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

    Hamiltonian[0][0] = f2+np.conj(f2)
    Hamiltonian[0][1] = f
    Hamiltonian[1][0] = np.conj(f)
    Hamiltonian[1][1] = f2+np.conj(f2)

    Hamiltonian[2][2] = f2+np.conj(f2)
    Hamiltonian[2][3] = f
    Hamiltonian[3][2] = np.conj(f)
    Hamiltonian[3][3] = f2+np.conj(f2)

    Hamiltonian[0][2] = xi-np.conj(xi)
    Hamiltonian[2][0] = -xi+np.conj(xi)
    Hamiltonian[1][3] = xi-np.conj(xi)
    Hamiltonian[3][1] = -xi+np.conj(xi)

    # print('02,20:')
    # print(Hamiltonian[0,2])
    # print(np.conj(Hamiltonian[2,0]))

    return Hamiltonian


def get_eigen(k):

    v, w = np.linalg.eigh(hamiltonian(k))
    eig = np.sort(v.real)

    return eig[0], eig[1], eig[2], eig[3]


if __name__ == '__main__':

    open("band_structure.txt", "w")
    band_structure = open("band_structure.txt", "a", buffering=1)

    nk = 101
    count = 0

    for i in range(nk):

        count = count+1
        k = K1-K1*float(i)/float(nk-1)
        k = np.matmul(k, bvec)
        count = count+1
        a, b, c, d = get_eigen(k)
        band_structure.write("{} {} {} {} {}\n".format(count, a, b, c, d))

    for i in range(nk):

        count = count+1
        k = GA+(MM-GA)*float(i)/float(nk-1)
        k = np.matmul(k, bvec)
        count = count+1
        a, b, c, d = get_eigen(k)
        band_structure.write("{} {} {} {} {}\n".format(count, a, b, c, d))

    for i in range(nk):

        count = count+1
        k = MM+(K2-MM)*float(i)/float(nk-1)
        k = np.matmul(k, bvec)
        count = count+1
        a, b, c, d = get_eigen(k)
        band_structure.write("{} {} {} {} {}\n".format(count, a, b, c, d))

