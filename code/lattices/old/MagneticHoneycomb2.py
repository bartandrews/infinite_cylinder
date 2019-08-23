"""Honeycomb lattice in a perpendicular magnetic field (Landau gauge in x direction with alpha=2/5)."""

import numpy as np
from tenpy.models import lattice
from tenpy.networks import site


class MagneticHoneycomb2(lattice.Lattice):

    def __init__(self, Lx, Ly, siteA, **kwargs):

        numb_sites = 4*8

        basis = np.array(([(8/4) * np.sqrt(3), 8/4], [0., 4]))
        delta = np.array([1 / (2. * np.sqrt(3.)), 0.5])

        pos_list = []
        for j in range(4):
            for i in range(0, int(8 / 2)):
                pos_list.append([(i / 2) * np.sqrt(3), i / 2+j] - delta / 2.)  # A site
                pos_list.append([(i / 2) * np.sqrt(3), i / 2+j] + delta / 2.)  # B site
        pos = tuple(pos_list)

        # print(len(pos))

        kwargs.setdefault('order', 'default')
        kwargs.setdefault('bc', 'periodic')
        kwargs.setdefault('bc_MPS', 'infinite')
        kwargs.setdefault('basis', basis)
        kwargs.setdefault('positions', pos)

        super().__init__([Lx, Ly], [siteA] * numb_sites, **kwargs)

        for i in range(numb_sites):
            if i % 2 == 0:
                if i < 8:
                    setattr(self, "NN{}d".format(i), [(i, 24+1+i, np.array([0, -1]))])
                elif 8 <= i < 16:
                    setattr(self, "NN{}d".format(i), [(i, i-8+1, np.array([0, 0]))])
                elif 16 <= i < 24:
                    setattr(self, "NN{}d".format(i), [(i, i-8+1, np.array([0, 0]))])
                elif 24 <= i < 32:
                    setattr(self, "NN{}d".format(i), [(i, i-8+1, np.array([0, 0]))])
            else:
                if i < 8:
                    setattr(self, "NN{}u".format(i), [(i, i + 8 - 1, np.array([0, 0]))])
                elif 8 <= i < 16:
                    setattr(self, "NN{}u".format(i), [(i, i + 8 - 1, np.array([0, 0]))])
                elif 16 <= i < 24:
                    setattr(self, "NN{}u".format(i), [(i, i + 8 - 1, np.array([0, 0]))])
                elif 24 <= i < 32:
                    setattr(self, "NN{}u".format(i), [(i, (i+7)%8, np.array([0, 1]))])

        for i in range(numb_sites):
            if i % 2 == 0:
                if i % 8 == 0:
                    setattr(self, "NN{}ul".format(i), [(i, i+8-1, np.array([-1, 0]))])
                else:
                    setattr(self, "NN{}ul".format(i), [(i, i - 1, np.array([0, 0]))])
            else:
                setattr(self, "NN{}bl".format(i), [(i, i-1, np.array([0, 0]))])

        for i in range(numb_sites):
            if i % 2 == 0:
                setattr(self, "NN{}ur".format(i), [(i, i + 1, np.array([0, 0]))])
            else:
                if (i+1) % 8 == 0:
                    setattr(self, "NN{}br".format(i), [(i, i + 1 - 8, np.array([1, 0]))])
                else:
                    setattr(self, "NN{}br".format(i), [(i, i + 1, np.array([0, 0]))])


def plot_lattice():

    numb_sites = 4*8

    import matplotlib.pyplot as plt

    ax = plt.gca()
    fs = site.FermionSite()
    lat = MagneticHoneycomb2(1, 1, fs, basis=[[(8/4) *np.sqrt(3), 8/4], [0, 4]])
    lat.plot_sites(ax)

    for i in range(0, numb_sites, 2):
        lat.plot_coupling(ax, getattr(lat, "NN{}d".format(i)), linestyle='-', color='red')
        lat.plot_coupling(ax, getattr(lat, "NN{}ul".format(i)), linestyle='-', color='blue')
        lat.plot_coupling(ax, getattr(lat, "NN{}ur".format(i)), linestyle='-', color='green')

    # for i in range(1, numb_sites, 2):
        # lat.plot_coupling(ax, getattr(lat, "NN{}u".format(i)), linestyle='-', color='red')
        # lat.plot_coupling(ax, getattr(lat, "NN{}bl".format(i)), linestyle='-', color='blue')
        # lat.plot_coupling(ax, getattr(lat, "NN{}br".format(i)), linestyle='-', color='green')

    ax.set_aspect('equal')
    plt.show()


if __name__ == "__main__":
    plot_lattice()
