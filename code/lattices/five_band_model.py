"""Interpenetrating lattice for the five-band model of twisted bilayer graphene.
Lattice based on: "Faithful Tight-binding Models and Fragile Topology of Magic-angle Bilayer Graphene"."""

import numpy as np
from tenpy.models import lattice
from tenpy.networks import site


class FiveBandLattice(lattice.Lattice):

    def __init__(self, Lx, Ly, site_a1, site_a2, site_d, **kwargs):

        pos = np.array([[-0.75, -0.25*np.sqrt(3)],  # left
                        [-0.25, +0.25*np.sqrt(3)],  # up
                        [+0.25, -0.25*np.sqrt(3)],  # down
                        [+0.75, +0.25*np.sqrt(3)]])  # right

        site.multi_sites_combine_charges([site_a1, site_a2, site_d])  # independently conserved charges

        super().__init__([Lx, Ly], [site_a1, site_d, site_d, site_a2], positions=pos, **kwargs)

        self.a1_d = [(1, 0, np.array([0, 1])),
                     (1, 0, np.array([0, 0])),
                     (2, 0, np.array([0, 0])),
                     (2, 0, np.array([1, 0]))]
        self.a2_d = [(2, 3, np.array([0, 0])),
                     (1, 3, np.array([0, 0])),
                     (3, 2, np.array([0, 1])),
                     (1, 3, np.array([1, 0]))]
        self.d_d = [(1, 2, np.array([0, 1])),
                    (1, 2, np.array([0, 0]))]
        self.a1_a2 = [(3, 0, np.array([1, 1])),
                      (3, 0, np.array([1, 0]))]


def plot_lattice():
    import matplotlib.pyplot as plt
    ax = plt.gca()
    site_a1 = site.FermionSite()
    site_a2 = site.FermionSite()
    site_d = site.FermionSite()
    lat = FiveBandLattice(3, 3, site_a1, site_a2, site_d, basis=[[2, 0], [0, np.sqrt(3)]])
    lat.plot_sites(ax)
    lat.plot_coupling(ax, lat.a1_d, linestyle='-', color='red')
    lat.plot_coupling(ax, lat.a2_d, linestyle='-', color='blue')
    lat.plot_coupling(ax, lat.d_d, linestyle='-', color='green')
    lat.plot_coupling(ax, lat.a1_a2, linestyle='-', color='orange')
    plt.show()


if __name__ == "__main__":
    plot_lattice()
