import numpy as np
from tenpy.models import lattice
from tenpy.networks import site


class BipartiteChain(lattice.Lattice):
    def __init__(self, Lx, Ly, siteA, **kwargs):
        basis = np.array(([2, 0.], [0, 1]))
        pos = np.array(([0, 0], [1, 0]))
        kwargs.setdefault('basis', basis)
        kwargs.setdefault('positions', pos)

        kwargs.setdefault('order', 'default')
        kwargs.setdefault('bc', ['periodic', 'open'])
        kwargs.setdefault('bc_MPS', 'infinite')

        super().__init__([Lx, Ly], [siteA, siteA], **kwargs)
        self.site = siteA

        self.intra_cell = [(0, 1, np.array([0, 0]))]
        self.inter_cell = [(1, 0, np.array([1, 0]))]

    def plot_lattice(self):
        import matplotlib.pyplot as plt
        ax = plt.gca()
        lat = BipartiteChain(6, 1, self.site)
        lat.plot_sites(ax)
        lat.plot_coupling(ax, lat.inter_cell, linestyle='-', color='green')
        lat.plot_coupling(ax, lat.intra_cell, linestyle='-', color='red')
        ax.set_aspect('equal')
        plt.show()
