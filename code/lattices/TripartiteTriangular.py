"""Tripartite triangular lattice for the C=3 Haldane model.
Lattice based on: "Topological Flat Band models with arbitrary Chern numbers"."""

import numpy as np
from tenpy.models import lattice
from tenpy.networks import site


class TripartiteTriangular(lattice.Lattice):
    def __init__(self, Lx, Ly, siteA, **kwargs):

        basis = np.array(([3., 0.], [0.5, 0.5*np.sqrt(3)]))

        pos = np.array(([0., 0.], [1., 0.], [2., 0.]))

        kwargs.setdefault('order', 'default')
        kwargs.setdefault('bc', 'periodic')
        kwargs.setdefault('bc_MPS', 'infinite')
        kwargs.setdefault('basis', basis)
        kwargs.setdefault('positions', pos)

        super().__init__([Lx, Ly], [siteA, siteA, siteA], **kwargs)

        self.NN = [(0, 2, np.array([-1, 1])), (0, 1, np.array([0, 0])), (0, 0, np.array([0, -1])),
                   (1, 0, np.array([0, 1])), (1, 2, np.array([0, 0])), (1, 1, np.array([0, -1])),
                   (2, 1, np.array([0, 1])), (2, 0, np.array([1, 0])), (2, 2, np.array([0, -1]))]

        self.nNNA = [(0, 2, np.array([-1, 2])), (0, 2, np.array([0, -1])), (0, 2, np.array([-1, -1])),
                     (1, 0, np.array([0, 2])), (1, 0, np.array([1, -1])), (1, 0, np.array([0, -1])),
                     (2, 1, np.array([0, 2])), (2, 1, np.array([1, -1])), (2, 1, np.array([0, -1]))]

        self.nNNB = [(0, 1, np.array([0, 1])), (0, 1, np.array([-1, 1])), (0, 1, np.array([0, -2])),
                     (1, 2, np.array([0, 1])), (1, 2, np.array([-1, 1])), (1, 2, np.array([0, -2])),
                     (2, 0, np.array([1, 1])), (2, 0, np.array([0, 1])), (2, 0, np.array([1, -2]))]

        self.nnNN = [(0, 1, np.array([-1, 2])), (0, 2, np.array([0, 0])), (0, 0, np.array([0, -2])),
                     (1, 2, np.array([-1, 2])), (1, 0, np.array([1, 0])), (1, 1, np.array([0, -2])),
                     (2, 0, np.array([0, 2])), (2, 1, np.array([1, 0])), (2, 2, np.array([0, -2]))]


def plot_lattice():
    import matplotlib.pyplot as plt
    ax = plt.gca()
    fs = site.FermionSite()
    lat = TripartiteTriangular(3, 3, fs)
    lat.plot_sites(ax)
    lat.plot_coupling(ax, lat.NN, linestyle='--', color='green')
    lat.plot_coupling(ax, lat.nNNA, linestyle='--', color='red')
    lat.plot_coupling(ax, lat.nNNB, linestyle='--', color='blue')
    lat.plot_coupling(ax, lat.nnNN, linestyle='--', color='black')
    ax.set_aspect('equal')
    plt.show()


if __name__ == "__main__":
    plot_lattice()
