"""Square lattice in a perpendicular magnetic field (Landau gauge in y direction with alpha=1/5).
Lattice based on: "Square Lattice with Magnetic Field", Aidelsburger PhD thesis, Fig.2.4.(a)."""

import numpy as np
from tenpy.models import lattice
from tenpy.networks import site


class MagneticSquare(lattice.Lattice):

    def __init__(self, Lx, Ly, siteA, **kwargs):

        basis = np.array(([1, 0.], [0, 3]))

        pos = np.array(([0, 0], [0, 1], [0, 2]))

        kwargs.setdefault('order', 'default')
        kwargs.setdefault('bc', 'periodic')
        kwargs.setdefault('bc_MPS', 'infinite')
        kwargs.setdefault('basis', basis)
        kwargs.setdefault('positions', pos)

        super().__init__([Lx, Ly], [siteA, siteA, siteA], **kwargs)

        self.NN_vert = [(0, 1, np.array([0, 0])), (1, 2, np.array([0, 0])), (2, 0, np.array([0, 1]))]

        self.NN_h0 = [(0, 0, np.array([1, 0]))]
        self.NN_h1 = [(1, 1, np.array([1, 0]))]
        self.NN_h2 = [(2, 2, np.array([1, 0]))]


def plot_lattice():

    import matplotlib.pyplot as plt

    ax = plt.gca()
    fs = site.FermionSite()
    lat = MagneticSquare(2, 1, fs, basis=[[1, 0], [0, 3]])
    lat.plot_sites(ax)
    lat.plot_coupling(ax, lat.NN_vert, linestyle='-', color='black')
    lat.plot_coupling(ax, lat.NN_h2, linestyle='-', color='green')
    lat.plot_coupling(ax, lat.NN_h1, linestyle='-', color='blue')
    lat.plot_coupling(ax, lat.NN_h0, linestyle='-', color='red')
    ax.set_aspect('equal')
    plt.show()


if __name__ == "__main__":
    plot_lattice()
