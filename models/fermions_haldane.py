"""Fermionic Haldane model: (spinless fermions with NNN hoppings t and U=V=0 and some arbitrary mu)."""

from tenpy.models.model import CouplingMPOModel, NearestNeighborModel
from tenpy.tools.params import get_parameter
from tenpy.networks.site import FermionSite
from tenpy.tools.misc import to_array
import numpy as np


class FermionicHaldaneModel(CouplingMPOModel):

    def __init__(self, model_params):
        CouplingMPOModel.__init__(self, model_params)

    def init_sites(self, model_params):
        conserve = get_parameter(model_params, 'conserve', 'N', self.name)
        filling = get_parameter(model_params, 'filling', 1/3, self.name)
        site = FermionSite(conserve=conserve, filling=filling)
        return site

    def add_phi_ext(self, strength, dx, phase):
        phi_ext = np.exp(1.j * np.pi * phase)
        (_, c_shape) = self.lat._coupling_shape(dx)
        strength = to_array(strength, c_shape) * (1. + 0.j)
        dy = dx[1]
        if dy> 0:
            strength[:, -dy:] *= phi_ext
        elif dy < 0:
            strength[:, :dy] *= np.conj(phi_ext)

        return strength

    def init_terms(self, model_params):
        # 0) Read out/set default parameters.
        Lx = get_parameter(model_params, 'Lx', 2, self.name)
        Ly = get_parameter(model_params, 'Ly', 2, self.name)
        t = get_parameter(model_params, 't', 1., self.name)
        V = get_parameter(model_params, 'V', 0, self.name)
        mu = get_parameter(model_params, 'mu', 0., self.name)
        phi_ext = get_parameter(model_params, 'phi_ext', 0., self.name)

        phi = np.arccos(3*np.sqrt(3/43))
        t2 = (np.sqrt(129)/36)*t * np.exp(1j * phi)

        for u in range(len(self.lat.unit_cell)):
            self.add_onsite(mu, 0, 'N')
            self.add_onsite(-mu, 1, 'N')
        for u1, u2, dx in self.lat.nearest_neighbors:
            #print(self.lat.nearest_neighbors)

            # phi_ext = 1

            t_phi = self.add_phi_ext(t, dx, phi_ext)

            self.add_coupling(t_phi, u1, 'Cd', u2, 'C', dx, 'JW', True)
            self.add_coupling(np.conj(t_phi), u1, 'Cd', u2, 'C', -dx, 'JW', True)  # h.c.
            self.add_coupling(V, u1, 'N', u2, 'N', dx)
        for u1, u2, dx in [(0, 0, np.array([1, 0])), (0, 0, np.array([0, -1])), (0, 0, np.array([-1, 1])),
                            (1, 1, np.array([1, 0])), (1, 1, np.array([0, -1])), (1, 1, np.array([1, -1]))]:
            t2_phi = self.add_phi_ext(t2, dx, phi_ext)
            self.add_coupling(t2_phi, u1, 'Cd', u2, 'C', dx, 'JW', True)
            self.add_coupling(np.conj(t2_phi), u1, 'Cd', u2, 'C', -dx, 'JW', True)  # h.c.


class FermionicHaldaneChain(FermionicHaldaneModel, NearestNeighborModel):

    def __init__(self, model_params):
        model_params.setdefault('lattice', "Chain")
        CouplingMPOModel.__init__(self, model_params)
