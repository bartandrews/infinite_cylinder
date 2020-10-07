import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import time
from sys import getsizeof


if __name__ == '__main__':

    dr, dx = 0.001/4, 0.001/4
    max_it = 10000
    it_cyc = 10
    R = np.arange(2, 4, dr)
    X = np.arange(0, 1, dx)
    Z = np.zeros((len(X), len(R)))
    z = np.zeros((len(R)))
    out = np.zeros((max_it, len(R)))

    for cy in range(it_cyc):
        z[:] = (cy + .5) / it_cyc  # use different starting values, not just 0.25
        for i in range(1000):  # stabilization
            z[:] = R[:] * z[:] * (1 - z[:])
        for i in range(max_it):  # iteration
            z[:] = R[:] * z[:] * (1 - z[:])
            out[i, :] = z[:]

        for ir, r in enumerate(R):
            h = np.histogram(out[:, ir], bins=list(X))[0]
            Z[1:, ir] += h[::-1]

    for ir, r in enumerate(R):  # even out intensities
        Z[:, ir] *= np.count_nonzero(Z[:, ir])

    Z = np.where(Z > 0, Z, np.NaN)

    # pick color bar range:
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.imshow(Z, interpolation='none', origin='lower', norm=LogNorm(), extent=(min(R), max(R), min(X), max(X)))
    ax.set_xlabel("r")
    ax.set_ylabel("x")
    plt.show()
