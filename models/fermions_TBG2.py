"""Spinless fermions with two orbitals.
Hamiltonian based on: "Model for the metal-insulator transition in graphene superlattices and beyond" (Sec. III)."""

import numpy as np

from tenpy.models.model import CouplingMPOModel, NearestNeighborModel
from tenpy.tools.params import get_parameter
from tenpy.networks.site import FermionSite, GroupedSite


class FermionicTBG2Model(CouplingMPOModel):

    def __init__(self, model_params):
        CouplingMPOModel.__init__(self, model_params)

    def init_sites(self, model_params):

        conserve = get_parameter(model_params, 'conserve', 'N', self.name)
        filling = get_parameter(model_params, 'filling', 1/2, self.name)
        fs = FermionSite(conserve=conserve, filling=filling)

        gs = GroupedSite([fs, fs], labels=['px', 'py'], charges='same')
        gs.add_op('Ntot', gs.Npx + gs.Npy, False)

        return gs

    def init_terms(self, model_params):

        t = get_parameter(model_params, 't', 1., self.name)
        mu = get_parameter(model_params, 'mu', 0., self.name)

        t1 = 0.331*t  # real
        t2 = (-0.010 + 0.097 * 1j)*t  # complex

        for u in range(len(self.lat.unit_cell)):

            self.add_onsite(-mu, 0, 'Ntot')

        for u1, u2, dx in self.lat.nearest_neighbors:

            self.add_coupling(t1, u1, 'Cdpx', u2, 'Cpx', dx, 'JW', True)
            self.add_coupling(t1, u2, 'Cdpx', u1, 'Cpx', -dx, 'JW', True)  # h.c.
            self.add_coupling(t1, u1, 'Cdpy', u2, 'Cpy', dx, 'JW', True)
            self.add_coupling(t1, u2, 'Cdpy', u1, 'Cpy', -dx, 'JW', True)  # h.c.

        for u1, u2, dx in self.lat.next_nearest_neighbors:

            self.add_coupling(np.real(t2), u1, 'Cdpx', u2, 'Cpx', dx, 'JW', True)
            self.add_coupling(np.real(t2), u2, 'Cdpx', u1, 'Cpx', -dx, 'JW', True)  # h.c.
            self.add_coupling(np.real(t2), u1, 'Cdpy', u2, 'Cpy', dx, 'JW', True)
            self.add_coupling(np.real(t2), u2, 'Cdpy', u1, 'Cpy', -dx, 'JW', True)  # h.c.

            self.add_coupling(np.imag(t2), u1, 'Cdpx', u2, 'Cpy', dx, 'JW', True)
            self.add_coupling(np.imag(t2), u2, 'Cdpy', u1, 'Cpx', -dx, 'JW', True)  # h.c.
            self.add_coupling(-np.imag(t2), u1, 'Cdpy', u2, 'Cpx', dx, 'JW', True)
            self.add_coupling(-np.imag(t2), u2, 'Cdpx', u1, 'Cpy', -dx, 'JW', True)  # h.c.


class FermionicTBG1Chain(FermionicTBG2Model, NearestNeighborModel):

    def __init__(self, model_params):
        model_params.setdefault('lattice', "Chain")
        CouplingMPOModel.__init__(self, model_params)
