# Example lattice class file, which is an executable python script used to plot a test lattice.

import numpy as np
from tenpy.models import lattice
from tenpy.networks import site


class MyLattice(lattice.Lattice):
    def __init__(self, Lx, Ly, siteA, siteB, **kwargs):
        pos = np.array([[0., 0.],  # center
                        [1., 0.],  # triangle right
                        [-0.5, 0.5*np.sqrt(3)],  # triangle top left
                        [-0.5, -0.5*np.sqrt(3)]])  # triangle bottom left
        pos *= 0.3 # scale triangle smaller
        site.multi_sites_combine_charges([siteA, siteB])  # independently conserved charges
        super().__init__([Lx, Ly], [siteA, siteB, siteB, siteB], positions=pos, **kwargs)
        self.BB_bonds = [(1, 2, np.array([0, 0])),
                         (2, 3, np.array([0, 0])),
                         (3, 1, np.array([0, 0]))]
        self.AB_bonds = [(1, 2, np.array([0, 0])),
                         (2, 3, np.array([0, 0])),
                         (3, 1, np.array([0, 0]))]
        self.AA_NN = [(0, 0, np.array([1, 0])),
                      (0, 0, np.array([0, 1]))]


def plot_lattice():
    import matplotlib.pyplot as plt
    ax = plt.gca()
    siteA = site.SpinHalfSite()
    siteB = site.FermionSite()
    lat = MyLattice(3, 3, siteA, siteB, basis=[[0.5 *np.sqrt(3), 0.5], [0, 1]])
    lat.plot_sites(ax)
    # lat.plot_coupling(ax, lat.AB_bonds, linestyle='-', color='green')
    lat.plot_coupling(ax, lat.AA_NN, linestyle='-', color='blue')
    # lat.plot_coupling(ax, lat.BB_bonds, linestyle='--', color='red')
    plt.show()


if __name__ == "__main__":
    plot_lattice()
