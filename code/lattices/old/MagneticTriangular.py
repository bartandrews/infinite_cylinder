"""Triangular lattice in a perpendicular magnetic field"""

import numpy as np
from tenpy.models import lattice
from tenpy.networks import site
import matplotlib.pyplot as plt


class MagneticTriangular(lattice.Lattice):

    def __init__(self, Lx, Ly, siteA, **kwargs):

        numb_sites = 3

        basis = np.array(([(numb_sites/2)*np.sqrt(3), numb_sites/2], [0, 1]))

        pos_list = []
        for i in range(numb_sites):
            pos_list.append([np.sqrt(3)/2 * i, i/2])
        pos = tuple(pos_list)

        kwargs.setdefault('order', 'default')
        kwargs.setdefault('bc', 'periodic')
        kwargs.setdefault('bc_MPS', 'infinite')
        kwargs.setdefault('basis', basis)
        kwargs.setdefault('positions', pos)

        super().__init__([Lx, Ly], [siteA] * numb_sites, **kwargs)

        # NN ###
        for i in range(numb_sites):
            if i == numb_sites - 1:
                setattr(self, "NN_ru{}".format(i), [(i, 0, np.array([1, 0]))])
            else:
                setattr(self, "NN_ru{}".format(i), [(i, i+1, np.array([0, 0]))])

            if i == numb_sites - 1:
                setattr(self, "NN_rd{}".format(i), [(i, 0, np.array([1, -1]))])
            else:
                setattr(self, "NN_rd{}".format(i), [(i, i + 1, np.array([0, -1]))])

            setattr(self, "NN_u{}".format(i), [(i, i, np.array([0, 1]))])

        # nNN ###
        for i in range(numb_sites):
            if i == 0:
                setattr(self, "nNN_u{}".format(i), [(i, numb_sites-1, np.array([-1, 2]))])
            else:
                setattr(self, "nNN_u{}".format(i), [(i, i-1, np.array([0, 2]))])

            # if i == 0:
            #     setattr(self, "nNN_ul{}".format(i), [(i, numb_sites-2, np.array([-1, 1]))])
            # elif i == 1:
            #     setattr(self, "nNN_ul{}".format(i), [(i, numb_sites-1, np.array([-1, 1]))])
            # else:
            #     setattr(self, "nNN_ul{}".format(i), [(i, i - 2, np.array([0, 1]))])

            if i == numb_sites-2:
                setattr(self, "nNN_br{}".format(i), [(i, 0, np.array([1, -1]))])
            elif i == numb_sites-1:
                setattr(self, "nNN_br{}".format(i), [(i, 1, np.array([1, -1]))])
            else:
                setattr(self, "nNN_br{}".format(i), [(i, i + 2, np.array([0, -1]))])

            if i == numb_sites-1:
                setattr(self, "nNN_ur{}".format(i), [(i, 0, np.array([1, 1]))])
            else:
                setattr(self, "nNN_ur{}".format(i), [(i, i+1, np.array([0, 1]))])


def plot_lattice():

    numb_sites = 3

    ax = plt.gca()
    fs = site.FermionSite()
    lat = MagneticTriangular(1, 6, fs, basis=[[(numb_sites/2)*np.sqrt(3), numb_sites/2], [0, 1]])
    lat.plot_sites(ax)

    for i in range(numb_sites):
        # lat.plot_coupling(ax, getattr(lat, "NN_u{}".format(i)), linestyle='-', color='C{}'.format(i))
        # lat.plot_coupling(ax, getattr(lat, "NN_ru{}".format(i)), linestyle='-', color='C{}'.format(i))
        # lat.plot_coupling(ax, getattr(lat, "NN_rd{}".format(i)), linestyle='-', color='C{}'.format(i))

        lat.plot_coupling(ax, getattr(lat, "nNN_u{}".format(i)), linestyle='-', color='C{}'.format(i))
        # lat.plot_coupling(ax, getattr(lat, "nNN_ul{}".format(i)), linestyle='-', color='C{}'.format(i))
        lat.plot_coupling(ax, getattr(lat, "nNN_br{}".format(i)), linestyle='-', color='C{}'.format(i))
        lat.plot_coupling(ax, getattr(lat, "nNN_ur{}".format(i)), linestyle='-', color='C{}'.format(i))

    ax.set_aspect('equal')
    plt.show()


if __name__ == "__main__":
    plot_lattice()
