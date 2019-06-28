from matrix import *
import matplotlib.pyplot as plt
import sys

values = []

# numb. energy bands 
sz = 100

# alpha samples
maximum = 100


def _make_alphas(a, size):
    return [a for i in range(0, size)]


if __name__ == '__main__':

    for index in range(0, maximum + 1):
        
        alpha = float(index)/float(maximum)
        
        alpha_list = _make_alphas(alpha, sz)

        eigenvalues = matrix_eigenvalues(0, alpha, sz)

        #print(eigenvalues)

        values.append((eigenvalues, alpha_list))

    fig = plt.figure()
    ax = plt.subplot(1,1,1)
    for eigenvalues, alphas in values:
        ax.plot(alphas, eigenvalues, '.', color='r', markersize=1.8)

    ax.set_title('Square Butterfly')
    ax.set_xlabel("Flux per Plaquette, Phi / 2*pi = alpha")
    ax.set_ylabel('Energy / J')

    #plt.axis([-4.5, 4.5, 0, 1])
    plt.show()