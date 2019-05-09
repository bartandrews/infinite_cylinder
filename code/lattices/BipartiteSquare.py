"""Bipartite square lattice for the chiral-pi-flux model.
Lattice based on: "Fractional Quantum Hall States at Zero Magnetic Field"."""

import numpy as np
from tenpy.models import lattice
from tenpy.networks import site


class BipartiteSquare(lattice.Lattice):
    def __init__(self, Lx, Ly, siteA, **kwargs):

        basis = np.array(([2, 0.], [0, 2]))

        pos = np.array(([0, 0], [1, 0], [0, 1], [1, 1]))

        kwargs.setdefault('order', 'default')
        kwargs.setdefault('bc', 'periodic')
        kwargs.setdefault('bc_MPS', 'infinite')
        kwargs.setdefault('basis', basis)
        kwargs.setdefault('positions', pos)

        super().__init__([Lx, Ly], [siteA, siteA, siteA, siteA], **kwargs)

        self.NN = [(0, 1, np.array([0, 0])), (1, 3, np.array([0, 0])),
                   (3, 2, np.array([0, 0])), (2, 0, np.array([0, 0])),
                   (2, 0, np.array([0, 1])), (1, 3, np.array([0, -1])),
                   (0, 1, np.array([-1, 0])), (3, 2, np.array([1, 0]))]
        self.nNNdashed = [(0, 3, np.array([0, 0])), (2, 1, np.array([0, 0])),
                          (3, 0, np.array([1, 1])), (1, 2, np.array([1, -1]))]
        self.nNNdotted = [(1, 2, np.array([1, 0])), (3, 0, np.array([1, 0])),
                          (2, 1, np.array([0, 1])), (3, 0, np.array([0, 1]))]


def plot_lattice():
    import matplotlib.pyplot as plt
    ax = plt.gca()
    fs = site.FermionSite()
    lat = BipartiteSquare(3, 3, fs, basis=[[2, 0], [0, 2]])
    lat.plot_sites(ax)
    lat.plot_coupling(ax, lat.NN, linestyle='-', color='green')
    lat.plot_coupling(ax, lat.nNNdashed, linestyle='--', color='black')
    lat.plot_coupling(ax, lat.nNNdotted, linestyle='--', color='red')
    ax.set_aspect('equal')
    plt.show()


if __name__ == "__main__":
    plot_lattice()
