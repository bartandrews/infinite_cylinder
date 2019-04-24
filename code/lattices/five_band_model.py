"""Interpenetrating lattice for the five-band model of twisted bilayer graphene.
Lattice based on: "Faithful Tight-binding Models and Fragile Topology of Magic-angle Bilayer Graphene"."""

import numpy as np
from tenpy.models import lattice
from tenpy.networks import site


class FiveBandLattice(lattice.Lattice):

    dim = 2

    def __init__(self, Lx, Ly, site_a, site_d, **kwargs):

        basis = [[3, 0], [0, np.sqrt(3)]]

        pos = np.array([[0, 0],
                        [0.5, 0.5*np.sqrt(3)],
                        [1, 0],
                        [1.5, 0.5*np.sqrt(3)],
                        [2, 0],
                        [2.5, 0.5*np.sqrt(3)]])

        kwargs.setdefault('basis', basis)
        kwargs.setdefault('positions', pos)

        site.multi_sites_combine_charges([site_a, site_d])  # independently conserved charges

        super().__init__([Lx, Ly], [site_a, site_d, site_d, site_a, site_d, site_d], **kwargs)

        self.a1_d = [(0, 1, np.array([0, 0])),
                     (0, 2, np.array([0, 0])),
                     (0, 1, np.array([0, -1])),
                     (0, 4, np.array([-1, 0])),
                     (0, 5, np.array([-1, 0])),
                     (0, 5, np.array([-1, -1]))]
        self.a2_d = [(3, 1, np.array([0, 0])),
                     (3, 2, np.array([0, 0])),
                     (3, 4, np.array([0, 0])),
                     (3, 5, np.array([0, 0])),
                     (3, 2, np.array([0, 1])),
                     (3, 4, np.array([0, 1]))]
        self.d_d = [(1, 2, np.array([0, 0])),
                    (1, 2, np.array([0, 1])),
                    (2, 4, np.array([0, 0])),
                    (4, 5, np.array([0, 0])),
                    (5, 4, np.array([0, 1])),
                    (5, 1, np.array([1, 0]))]


def plot_lattice():
    import matplotlib.pyplot as plt
    ax = plt.gca()
    site_a = site.FermionSite()
    site_d = site.FermionSite()
    lat = FiveBandLattice(2, 2, site_a, site_d)
    lat.plot_sites(ax)
    lat.plot_coupling(ax, lat.a1_d, linestyle='-', color='red')
    lat.plot_coupling(ax, lat.a2_d, linestyle='-', color='blue')
    lat.plot_coupling(ax, lat.d_d, linestyle='--', color='green')
    plt.show()


if __name__ == "__main__":
    plot_lattice()
