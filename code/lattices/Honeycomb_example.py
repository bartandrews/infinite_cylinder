# Example Honeycomb lattice class file, which is an executable python script used to plot a test lattice.

import numpy as np
from tenpy.models import lattice
from tenpy.networks import site


class MyLattice(lattice.Lattice):
    def __init__(self, Lx, Ly, siteA, siteB, **kwargs):
        basis = np.array(([0.5 * np.sqrt(3), 0.5], [0., 1]))
        delta = np.array([1 / (2. * np.sqrt(3.)), 0.5])
        pos = (-delta / 2., delta / 2)
        kwargs.setdefault('basis', basis)
        kwargs.setdefault('positions', pos)
        super().__init__([Lx, Ly], [siteA, siteB], **kwargs)
        self.NN = [(0, 1, np.array([0, 0])), (1, 0, np.array([1, 0])), (1, 0, np.array([0, 1]))]
        self.nNN = [(0, 0, np.array([1, 0])), (0, 0, np.array([0, 1])), (0, 0, np.array([1, -1])),
                    (1, 1, np.array([1, 0])), (1, 1, np.array([0, 1])), (1, 1, np.array([1, -1]))]
        self.nnNN = [(1, 0, np.array([1, 1])), (0, 1, np.array([-1, 1])), (0, 1, np.array([1, -1]))]
        self.fifthNN = [(0, 0, np.array([2, -1])), (0, 0, np.array([1, 1])), (0, 0, np.array([-1, 2])),
                        (0, 0, np.array([-2, 1])), (0, 0, np.array([-1, -1])), (0, 0, np.array([1, -2]))]


def plot_lattice():
    import matplotlib.pyplot as plt
    ax = plt.gca()
    site1 = site.FermionSite()
    lat = MyLattice(3, 3, site1, site1, basis=[[0.5 *np.sqrt(3), 0.5], [0, 1]])
    lat.plot_sites(ax)
    lat.plot_coupling(ax, lat.fifthNN, linestyle='-', color='green')
    # lat.plot_coupling(ax, lat.nNN, linestyle='-', color='blue')
    # lat.plot_coupling(ax, lat.nnNN, linestyle='-', color='red')
    ax.set_aspect('equal')
    plt.show()


if __name__ == "__main__":
    plot_lattice()
