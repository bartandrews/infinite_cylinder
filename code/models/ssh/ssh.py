# --- python imports
import numpy as np
# --- tenpy imports
from tenpy.models.model import CouplingMPOModel
from tenpy.networks.site import FermionSite
# --- infinite_cylinder imports
from lattices.BipartiteChain import BipartiteChain


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
        Lx = params.get('LxMUC', 6)
        Ly = params.get('Ly', 1)
        site = self.init_sites(params)
        order = params.get('order', 'default')
        bc_MPS = params.get('bc_MPS', 'infinite')
        bc_x = params.get('bc_x', 'periodic')
        bc_y = params.get('bc_y', 'open')
        bc = [bc_x, bc_y]
        lat = BipartiteChain(Lx, Ly, site, order=order, bc_MPS=bc_MPS, bc=bc)
        return lat

    def init_terms(self, params):
        creation, annihilation = 'Cd', 'C'
        t1 = params.get('t1', 1)
        t2 = params.get('t2', 1)

        # inter-cell hopping
        u1, u2, dx = (0, 1, np.array([0, 0]))
        self.add_coupling(t1, u1, creation, u2, annihilation, dx)
        self.add_coupling(np.conj(t1), u2, creation, u1, annihilation, -dx)  # H.c.
        # intra-cell hopping
        u1, u2, dx = (1, 0, np.array([1, 0]))
        self.add_coupling(t2, u1, creation, u2, annihilation, dx)
        self.add_coupling(np.conj(t2), u2, creation, u1, annihilation, -dx)  # H.c.


if __name__ == "__main__":

    model_params = dict(t1=0.5, t2=1, n=(int(1), int(2)), LxMUC=6, Ly=1)
    M = SSHModel(model_params)

    print("max MPO bond dimension = ", max(M.H_MPO.chi))
