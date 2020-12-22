# --- python imports
import numpy as np
import matplotlib.pyplot as plt
# --- tenpy imports
from tenpy.models.model import CouplingMPOModel
from tenpy.networks.site import FermionSite
from tenpy.models.lattice import Chain


class SSHModel(CouplingMPOModel):

    def __init__(self, params):
        CouplingMPOModel.__init__(self, params)

    def init_sites(self, params):
        conserve = params.get('conserve', 'N')
        n = params.get('n', (1, 2))
        n = n[0] / n[1]
        site = FermionSite(conserve=conserve, filling=n)
        return site

    def init_lattice(self, params):
        L = params.get('L', 8)
        site = self.init_sites(params)
        lat = Chain(L, site, bc='periodic', bc_MPS='infinite')
        return lat

    def init_terms(self, params):
        t1 = params.get('t1', 1)
        t2 = params.get('t2', 1)

        # u1, u2, dx = (0, 0, 1)
        self.add_coupling([np.conj(-t1), -t2], 0, 'Cd', 0, 'C', 1, plus_hc=True)
        # self.add_coupling([-t1, np.conj(-t2)], u1, creation, u2, annihilation, -dx)  # H.c.


if __name__ == "__main__":

    model_params = dict(t1=-1, t2=0, n=(int(1), int(2)), L=8)
    M = SSHModel(model_params)

    print("max MPO bond dimension = ", max(M.H_MPO.chi))

    ax = plt.gca()
    M.lat.plot_sites(ax)
    M.lat.plot_coupling(ax, linestyle='-', color='red')
    ax.set_aspect('equal')
    M.lat.plot_order(ax, linestyle='dotted', color='k')
    plt.show()
