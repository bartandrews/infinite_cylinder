"""Square lattice in a perpendicular magnetic field (Landau gauge in y direction with alpha=1/5).
Lattice based on: "Square Lattice with Magnetic Field", Aidelsburger PhD thesis, Fig.2.4.(a)."""

import numpy as np
from tenpy.models import lattice
from tenpy.networks import site


class MagneticSquare(lattice.Lattice):

    def __init__(self, Lx, Ly, siteA, **kwargs):

        numb_sites = 11

        basis = np.array(([1, 0.], [0, numb_sites]))

        pos_list = []
        for i in range(numb_sites):
            pos_list.append([0, i])
        pos = tuple(pos_list)

        kwargs.setdefault('order', 'default')
        kwargs.setdefault('bc', 'periodic')
        kwargs.setdefault('bc_MPS', 'infinite')
        kwargs.setdefault('basis', basis)
        kwargs.setdefault('positions', pos)

        super().__init__([Lx, Ly], [siteA] * numb_sites, **kwargs)

        NN_vert_list=[]
        for i in range(numb_sites):
            if i == numb_sites-1:
                NN_vert_list.append((i, 0, np.array([0, 1])))
            else:
                NN_vert_list.append((i, i+1, np.array([0, 0])))
        self.NN_vert = NN_vert_list

        for i in range(numb_sites):
            setattr(self, "NN_h{}".format(i), [(i, i, np.array([1, 0]))])


def plot_lattice():

    numb_sites = 11

    import matplotlib.pyplot as plt

    ax = plt.gca()
    fs = site.FermionSite()
    lat = MagneticSquare(2, 1, fs, basis=[[1, 0], [0, numb_sites]])
    lat.plot_sites(ax)
    lat.plot_coupling(ax, lat.NN_vert, linestyle='-', color='black')

    for i in range(numb_sites):
        lat.plot_coupling(ax, getattr(lat, "NN_h{}".format(i)), linestyle='-', color='C{}'.format(i))

    ax.set_aspect('equal')
    plt.show()


if __name__ == "__main__":
    plot_lattice()
