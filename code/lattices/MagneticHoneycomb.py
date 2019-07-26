"""Honeycomb lattice in a perpendicular magnetic field (Landau gauge in x direction with alpha=2/5)."""

import numpy as np
from tenpy.models import lattice
from tenpy.networks import site


class MagneticHoneycomb(lattice.Lattice):

    def __init__(self, Lx, Ly, siteA, **kwargs):

        basis = np.array(([2.5 * np.sqrt(3), 2.5], [0., 1]))
        delta = np.array([1 / (2. * np.sqrt(3.)), 0.5])
        pos = (-delta / 2., delta / 2, [0.5 * np.sqrt(3), 0.5] - delta / 2., [0.5 * np.sqrt(3), 0.5] + delta / 2,
               [1 * np.sqrt(3), 1] - delta / 2., [1 * np.sqrt(3), 1] + delta / 2, [1.5 * np.sqrt(3), 1.5] - delta / 2., [1.5 * np.sqrt(3), 1.5] + delta / 2,
               [2 * np.sqrt(3), 2] - delta / 2., [2 * np.sqrt(3), 2] + delta / 2)

        kwargs.setdefault('order', 'default')
        kwargs.setdefault('bc', 'periodic')
        kwargs.setdefault('bc_MPS', 'infinite')
        kwargs.setdefault('basis', basis)
        kwargs.setdefault('positions', pos)

        super().__init__([Lx, Ly], [siteA, siteA, siteA, siteA, siteA, siteA, siteA, siteA, siteA, siteA], **kwargs)

        self.NN0d = [(0, 1, np.array([0, -1]))]
        self.NN2d = [(2, 3, np.array([0, -1]))]
        self.NN4d = [(4, 5, np.array([0, -1]))]
        self.NN6d = [(6, 7, np.array([0, -1]))]
        self.NN8d = [(8, 9, np.array([0, -1]))]

        self.NN0ul = [(0, 9, np.array([-1, 0]))]
        self.NN2ul = [(2, 1, np.array([0, 0]))]
        self.NN4ul = [(4, 3, np.array([0, 0]))]
        self.NN6ul = [(6, 5, np.array([0, 0]))]
        self.NN8ul = [(8, 7, np.array([0, 0]))]

        self.NN0ur = [(0, 1, np.array([0, 0]))]
        self.NN2ur = [(2, 3, np.array([0, 0]))]
        self.NN4ur = [(4, 5, np.array([0, 0]))]
        self.NN6ur = [(6, 7, np.array([0, 0]))]
        self.NN8ur = [(8, 9, np.array([0, 0]))]


def plot_lattice():

    import matplotlib.pyplot as plt

    ax = plt.gca()
    fs = site.FermionSite()
    lat = MagneticHoneycomb(5, 5, fs, basis=[[2.5 *np.sqrt(3), 2.5], [0, 1]])
    lat.plot_sites(ax)
    # down
    lat.plot_coupling(ax, lat.NN0d, linestyle='-', color='black')
    lat.plot_coupling(ax, lat.NN2d, linestyle='-', color='red')
    lat.plot_coupling(ax, lat.NN4d, linestyle='-', color='blue')
    lat.plot_coupling(ax, lat.NN6d, linestyle='-', color='green')
    lat.plot_coupling(ax, lat.NN8d, linestyle='-', color='orange')
    # upper left
    lat.plot_coupling(ax, lat.NN0ul, linestyle='-', color='black')
    lat.plot_coupling(ax, lat.NN2ul, linestyle='-', color='red')
    lat.plot_coupling(ax, lat.NN4ul, linestyle='-', color='blue')
    lat.plot_coupling(ax, lat.NN6ul, linestyle='-', color='green')
    lat.plot_coupling(ax, lat.NN8ul, linestyle='-', color='orange')
    # upper right
    lat.plot_coupling(ax, lat.NN0ur, linestyle='-', color='black')
    lat.plot_coupling(ax, lat.NN2ur, linestyle='-', color='red')
    lat.plot_coupling(ax, lat.NN4ur, linestyle='-', color='blue')
    lat.plot_coupling(ax, lat.NN6ur, linestyle='-', color='green')
    lat.plot_coupling(ax, lat.NN8ur, linestyle='-', color='orange')
    ax.set_aspect('equal')
    plt.show()


if __name__ == "__main__":
    plot_lattice()
