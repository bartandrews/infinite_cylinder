import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

if __name__ == '__main__':

    # initialization
    dr, dx = 0.01, 0.01
    max_it = 10000
    R = np.arange(0, 4, dr)
    X = np.arange(0, 1, dx)
    Z = np.zeros((len(X), len(R)))

    # evaluation
    for ir, r in enumerate(R):
        if ir % 100 == 0:
            print(ir, "of", len(R))
        z = 0.25
        for i in range(1000):  # stabilization
            z = r * z * (1 - z)
        for i in range(max_it):  # iteration
            z = r * z * (1 - z)
            Z[int(z / dx), ir] += 1
        Z[:, ir] *= np.count_nonzero(Z[:, ir])  # even out intensities
    Z = np.where(Z > 0, Z, np.NaN)

    print(Z)

    # plotting
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.imshow(Z, interpolation='none', origin='lower', norm=LogNorm(), extent=(min(R), max(R), min(X), max(X)))
    ax.set_xlabel("r")
    ax.set_ylabel("x")
    plt.show()
