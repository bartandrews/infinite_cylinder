# --- python imports
import numpy as np
# --- tenpy imports
from tenpy.models.model import CouplingMPOModel
from tenpy.networks.site import SpinSite
from tenpy.models.lattice import Chain


class HeisenbergModel(CouplingMPOModel):

    def __init__(self, params):
        CouplingMPOModel.__init__(self, params)

    def init_sites(self, params):
        site = SpinSite(S=1, conserve=None)
        return site

    def init_lattice(self, params):
        L = params.get('L', 6)
        site = self.init_sites(params)
        lat = Chain(L, site, bc='periodic', bc_MPS='infinite')
        return lat

    def init_terms(self, params):
        J = params.get('J', 1.)
        D = params.get('D', 1.)
        for u in range(len(self.lat.unit_cell)):
            self.add_onsite(0.3, u, 'Sx')
            self.add_onsite(D, u, 'Sz Sz')
        for u1, u2, dx in self.lat.pairs['nearest_neighbors']:
            self.add_coupling((J + J) / 4., u1, 'Sp', u2, 'Sm', dx, plus_hc=True)
            self.add_coupling((J - J) / 4., u1, 'Sp', u2, 'Sp', dx, plus_hc=True)
            self.add_coupling(J, u1, 'Sz', u2, 'Sz', dx)


if __name__ == "__main__":

    model_params = dict(J=1, D=1, L=6)
    M = HeisenbergModel(model_params)

    print("max MPO bond dimension = ", max(M.H_MPO.chi))
