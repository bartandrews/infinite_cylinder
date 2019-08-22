"""Honeycomb lattice in a perpendicular magnetic field (Landau gauge in x direction with alpha=2/5)."""

import numpy as np
from tenpy.models import lattice
from tenpy.networks import site


class MagneticHoneycomb(lattice.Lattice):

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

        # nearest neighbors

        for i in range(numb_sites):
            if i % 2 == 0:
                setattr(self, "NN{}d".format(i), [(i, i+1, np.array([0, -1]))])
            else:
                setattr(self, "NN{}u".format(i), [(i, i-1, np.array([0, 1]))])

        for i in range(numb_sites):
            if i % 2 == 0:
                if i == 0:
                    setattr(self, "NN0ul", [(i, numb_sites-1, np.array([-1, 0]))])
                else:
                    setattr(self, "NN{}ul".format(i), [(i, i - 1, np.array([0, 0]))])
            else:
                setattr(self, "NN{}bl".format(i), [(i, i-1, np.array([0, 0]))])

        for i in range(numb_sites):
            if i % 2 == 0:
                setattr(self, "NN{}ur".format(i), [(i, i + 1, np.array([0, 0]))])
            else:
                if i == numb_sites-1:
                    setattr(self, "NN{}br".format(i), [(i, 0, np.array([1, 0]))])
                else:
                    setattr(self, "NN{}br".format(i), [(i, i+1, np.array([0, 0]))])

        # fifth nearest neighbors

        for i in range(numb_sites):
            if i % 2 == 0:
                if i == 0:
                    setattr(self, "fifthNN{}u".format(i), [(i, numb_sites-2, np.array([-1, 2]))])
                else:
                    setattr(self, "fifthNN{}u".format(i), [(i, i-2, np.array([0, 2]))])
            else:
                if i == 1:
                    setattr(self, "fifthNN{}u".format(i), [(i, numb_sites - 1, np.array([-1, 2]))])
                else:
                    setattr(self, "fifthNN{}u".format(i), [(i, i - 2, np.array([0, 2]))])

        for i in range(numb_sites):
            if i % 2 == 0:
                if i == 0:
                    setattr(self, "fifthNN{}ul".format(i), [(i, numb_sites-4, np.array([-1, 1]))])
                elif i == 2:
                    setattr(self, "fifthNN{}ul".format(i), [(i, numb_sites-2, np.array([-1, 1]))])
                else:
                    setattr(self, "fifthNN{}ul".format(i), [(i, i-4, np.array([0, 1]))])
            else:
                if i == 1:
                    setattr(self, "fifthNN{}ul".format(i), [(i, numb_sites - 3, np.array([-1, 1]))])
                elif i == 3:
                    setattr(self, "fifthNN{}ul".format(i), [(i, numb_sites - 1, np.array([-1, 1]))])
                else:
                    setattr(self, "fifthNN{}ul".format(i), [(i, i - 4, np.array([0, 1]))])

        for i in range(numb_sites):
            if i % 2 == 0:
                if i == numb_sites-2:
                    setattr(self, "fifthNN{}ur".format(i), [(i, 0, np.array([1, 1]))])
                else:
                    setattr(self, "fifthNN{}ur".format(i), [(i, i+2, np.array([0, 1]))])
            else:
                if i == numb_sites-1:
                    setattr(self, "fifthNN{}ur".format(i), [(i, 1, np.array([1, 1]))])
                else:
                    setattr(self, "fifthNN{}ur".format(i), [(i, i + 2, np.array([0, 1]))])


def plot_lattice():

    numb_sites = 8

    import matplotlib.pyplot as plt

    ax = plt.gca()
    fs = site.FermionSite()
    lat = MagneticHoneycomb(1, 1, fs, basis=[[(numb_sites/4) *np.sqrt(3), numb_sites/4], [0, 1]])
    lat.plot_sites(ax)

    # for i in range(0, numb_sites, 2):
    #     lat.plot_coupling(ax, getattr(lat, "NN{}d".format(i)), linestyle='-', color='red')
    #     lat.plot_coupling(ax, getattr(lat, "NN{}ul".format(i)), linestyle='-', color='blue')
    #     lat.plot_coupling(ax, getattr(lat, "NN{}ur".format(i)), linestyle='-', color='green')
    # for i in range(1, numb_sites, 2):
    #     lat.plot_coupling(ax, getattr(lat, "NN{}u".format(i)), linestyle='-', color='red')
    #     lat.plot_coupling(ax, getattr(lat, "NN{}bl".format(i)), linestyle='-', color='blue')
    #     lat.plot_coupling(ax, getattr(lat, "NN{}br".format(i)), linestyle='-', color='green')

    for i in range(numb_sites):
        lat.plot_coupling(ax, getattr(lat, "fifthNN{}u".format(i)), linestyle='-', color='red')
        lat.plot_coupling(ax, getattr(lat, "fifthNN{}ul".format(i)), linestyle='-', color='blue')
        lat.plot_coupling(ax, getattr(lat, "fifthNN{}ur".format(i)), linestyle='-', color='green')

    ax.set_aspect('equal')
    plt.show()


if __name__ == "__main__":
    plot_lattice()
