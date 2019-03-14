"""Fermionic Haldane model: (spinless fermions with NNN hoppings t and U=V=0 and some arbitrary mu)."""

from tenpy.models.model import CouplingMPOModel, NearestNeighborModel
from tenpy.tools.params import get_parameter
from tenpy.networks.site import FermionSite
import numpy as np


class FermionicHaldaneModel(CouplingMPOModel):

    def __init__(self, model_params):
        model_params.setdefault('lattice', 'Honeycomb')
        CouplingMPOModel.__init__(self, model_params)

    def init_sites(self, model_params):

        conserve = get_parameter(model_params, 'conserve', 'N', self.name)
        filling = get_parameter(model_params, 'filling', 1/3, self.name)
        site = FermionSite(conserve=conserve, filling=filling)
        return site

    def init_terms(self, model_params):

        t = get_parameter(model_params, 't', 1., self.name, True)
        V = get_parameter(model_params, 'V', 0, self.name, True)
        mu = get_parameter(model_params, 'mu', 0., self.name, True)
        phi_ext = get_parameter(model_params, 'phi_ext', 0., self.name)

        phi = np.arccos(3*np.sqrt(3/43))
        t2 = (np.sqrt(129)/36)*t * np.exp(1j * phi)
        # t2 = 0

        for u in range(len(self.lat.unit_cell)):

            self.add_onsite(mu, 0, 'N')
            self.add_onsite(-mu, 1, 'N')

        for u1, u2, dx in self.lat.nearest_neighbors:

            t_phi = self.coupling_strength_add_ext_flux(t, dx, [0, phi_ext])

            self.add_coupling(t_phi, u1, 'Cd', u2, 'C', dx, 'JW', True)
            self.add_coupling(np.conj(t_phi), u2, 'Cd', u1, 'C', -dx, 'JW', True)  # h.c.
            self.add_coupling(V, u1, 'N', u2, 'N', dx)

        for u1, u2, dx in [(0, 0, np.array([1, 0])), (0, 0, np.array([0, -1])), (0, 0, np.array([-1, 1])),
                           (1, 1, np.array([-1, 0])), (1, 1, np.array([0, 1])), (1, 1, np.array([1, -1]))]:

            t2_phi = self.coupling_strength_add_ext_flux(t2, dx, [0, phi_ext])
            self.add_coupling(t2_phi, u1, 'Cd', u2, 'C', dx, 'JW', True)
            self.add_coupling(np.conj(t2_phi), u2, 'Cd', u1, 'C', -dx, 'JW', True)  # h.c.


class FermionicHaldaneChain(FermionicHaldaneModel, NearestNeighborModel):

    def __init__(self, model_params):
        model_params.setdefault('lattice', "Chain")
        CouplingMPOModel.__init__(self, model_params)
