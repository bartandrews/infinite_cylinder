"""Honeycomb lattice in a perpendicular magnetic field (Landau gauge in x direction with alpha=2/5)."""

import numpy as np
from tenpy.models import lattice
from tenpy.networks import site


class MagneticHoneycomb3(lattice.Lattice):

    def __init__(self, Lx, Ly, siteA, **kwargs):

        numb_sites = 8

        basis = np.array(([(numb_sites/4) * np.sqrt(3), numb_sites/4], [0., 1]))
        delta = np.array([1 / (2. * np.sqrt(3.)), 0.5])

        pos_list = []
        for i in range(0, int(numb_sites / 2)):
            pos_list.append([(i / 2) * np.sqrt(3), i / 2] - delta / 2.)  # A site
            pos_list.append([(i / 2) * np.sqrt(3), i / 2] + delta / 2.)  # B site
        pos = tuple(pos_list)

        kwargs.setdefault('order', 'default')
        kwargs.setdefault('bc', 'periodic')
        kwargs.setdefault('bc_MPS', 'infinite')
        kwargs.setdefault('basis', basis)
        kwargs.setdefault('positions', pos)

        super().__init__([Lx, Ly], [siteA] * numb_sites, **kwargs)

        for i in range(numb_sites):
            if i % 2 != 0:
                setattr(self, "NN{}ul".format(i), [(i, i-1, np.array([0, 1]))])

        for i in range(numb_sites):
            if i % 2 != 0:
                setattr(self, "NN{}bl".format(i), [(i, i-1, np.array([0, 0]))])

        for i in range(numb_sites):
            if i % 2 != 0:
                if i == numb_sites-1:
                    setattr(self, "NN{}r".format(i), [(i, 0, np.array([1, 0]))])
                else:
                    setattr(self, "NN{}r".format(i), [(i, i+1, np.array([0, 0]))])


def plot_lattice():

    numb_sites = 8

    import matplotlib.pyplot as plt

    ax = plt.gca()
    fs = site.FermionSite()
    lat = MagneticHoneycomb3(1, 4, fs, basis=[[(numb_sites/4) *np.sqrt(3), numb_sites/4], [0, 1]])
    lat.plot_sites(ax)

    for i in range(1, numb_sites, 2):
        lat.plot_coupling(ax, getattr(lat, "NN{}ul".format(i)), linestyle='-', color='red')
        lat.plot_coupling(ax, getattr(lat, "NN{}bl".format(i)), linestyle='-', color='blue')
        lat.plot_coupling(ax, getattr(lat, "NN{}r".format(i)), linestyle='-', color='green')

    ax.set_aspect('equal')
    plt.show()


if __name__ == "__main__":
    plot_lattice()
