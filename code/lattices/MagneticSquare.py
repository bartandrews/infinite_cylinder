"""Square lattice in a perpendicular magnetic field (Landau gauge in y direction with alpha=1/5).
Lattice based on: "Square Lattice with Magnetic Field", Aidelsburger PhD thesis, Fig.2.4.(a)."""

import numpy as np
from tenpy.models import lattice
from tenpy.networks import site


class MagneticSquare(lattice.Lattice):

    def __init__(self, Lx, Ly, siteA, **kwargs):

        numb_sites = 4

        basis = np.array(([numb_sites, 0.], [0, 1]))

        pos_list = []
        for i in range(numb_sites):
            pos_list.append([i, 0])
        pos = tuple(pos_list)

        kwargs.setdefault('order', 'default')
        kwargs.setdefault('bc', 'periodic')
        kwargs.setdefault('bc_MPS', 'infinite')
        kwargs.setdefault('basis', basis)
        kwargs.setdefault('positions', pos)

        super().__init__([Lx, Ly], [siteA] * numb_sites, **kwargs)

        NN_horiz_list=[]
        for i in range(numb_sites):
            if i == numb_sites-1:
                NN_horiz_list.append((i, 0, np.array([1, 0])))
            else:
                NN_horiz_list.append((i, i+1, np.array([0, 0])))
        self.NN_horiz = NN_horiz_list

        for i in range(numb_sites):
            setattr(self, "NN_v{}".format(i), [(i, i, np.array([0, 1]))])


def plot_lattice():

    numb_sites = 4

    import matplotlib.pyplot as plt

    ax = plt.gca()
    fs = site.FermionSite()
    lat = MagneticSquare(1, 4, fs, basis=[[numb_sites, 0], [0, 1]])
    lat.plot_sites(ax)
    lat.plot_coupling(ax, lat.NN_horiz, linestyle='-', color='black')

    for i in range(numb_sites):
        lat.plot_coupling(ax, getattr(lat, "NN_v{}".format(i)), linestyle='-', color='C{}'.format(i))

    ax.set_aspect('equal')
    plt.show()


if __name__ == "__main__":
    plot_lattice()
