"""Spinless fermions with three orbitals.
Hamiltonian based on: "Model for the metal-insulator transition in graphene superlattices and beyond" (Sec. III)."""

import numpy as np

from tenpy.models.model import CouplingMPOModel, NearestNeighborModel
from tenpy.tools.params import get_parameter
from tenpy.networks.site import FermionSite, GroupedSite


class FermionicTBG6Model(CouplingMPOModel):

    def __init__(self, model_params):
        CouplingMPOModel.__init__(self, model_params)

    def init_sites(self, model_params):

        conserve = get_parameter(model_params, 'conserve', 'N', self.name)
        fs = FermionSite(conserve=conserve)

        gs = GroupedSite([fs, fs, fs], labels=['px', 'py', 'z'], charges='same')
        gs.add_op('Ntot', gs.Npx + gs.Npy + gs.Nz, False)

        return gs

    def init_terms(self, model_params):

        t = get_parameter(model_params, 't', 1., self.name)
        mu = get_parameter(model_params, 'mu', 0., self.name)
        U = get_parameter(model_params, 'U', 0., self.name)
        V = get_parameter(model_params, 'V', 1, self.name)
        phi_ext = 2 * np.pi * get_parameter(model_params, 'phi_ext', 0., self.name)

        t1 = t
        phi = np.arccos(3 * np.sqrt(3 / 43))
        t2 = (np.sqrt(129) / 36) * t * np.exp(1j * phi)
        t2dash = 0.4*t2
        t3 = 10*t
        mu3 = 50

        for u in range(len(self.lat.unit_cell)):

            self.add_onsite(mu3, 0, 'Nz')
            self.add_onsite(-mu3, 1, 'Nz')

        for u1, u2, dx in self.lat.nearest_neighbors:

            t1_phi = self.coupling_strength_add_ext_flux(t1, dx, [0, phi_ext])
            self.add_coupling(t1_phi, u1, 'Cdpx', u2, 'Cpx', dx, 'JW', True)
            self.add_coupling(np.conj(t1_phi), u2, 'Cdpx', u1, 'Cpx', -dx, 'JW', True)  # h.c.
            self.add_coupling(t1_phi, u1, 'Cdpy', u2, 'Cpy', dx, 'JW', True)
            self.add_coupling(np.conj(t1_phi), u2, 'Cdpy', u1, 'Cpy', -dx, 'JW', True)  # h.c.

            t3_phi = self.coupling_strength_add_ext_flux(t3, dx, [0, phi_ext])
            self.add_coupling(t3_phi, u1, 'Cdz', u2, 'Cz', dx, 'JW', True)
            self.add_coupling(np.conj(t3_phi), u2, 'Cdz', u1, 'Cz', -dx, 'JW', True)  # h.c.

            self.add_coupling(V, u1, 'Ntot', u2, 'Ntot', dx)

        for u1, u2, dx in [(0, 0, np.array([-1, 1])), (0, 0, np.array([1, 0])), (0, 0, np.array([0, -1])),
                           (1, 1, np.array([0, 1])), (1, 1, np.array([1, -1])), (1, 1, np.array([-1, 0]))]:

            t2_phi = self.coupling_strength_add_ext_flux(t2, dx, [0, phi_ext])
            self.add_coupling(t2_phi, u1, 'Cdpx', u2, 'Cpx', dx, 'JW', True)
            self.add_coupling(np.conj(t2_phi), u2, 'Cdpx', u1, 'Cpx', -dx, 'JW', True)  # h.c.
            self.add_coupling(t2_phi, u1, 'Cdpy', u2, 'Cpy', dx, 'JW', True)
            self.add_coupling(np.conj(t2_phi), u2, 'Cdpy', u1, 'Cpy', -dx, 'JW', True)  # h.c.

            t2dash_phi = self.coupling_strength_add_ext_flux(t2dash, dx, [0, phi_ext])
            self.add_coupling(t2dash_phi, u1, 'Cdpx', u2, 'Cpy', dx, 'JW', True)
            self.add_coupling(np.conj(t2dash_phi), u2, 'Cdpy', u1, 'Cpx', -dx, 'JW', True)  # h.c.
            self.add_coupling(-t2dash_phi, u1, 'Cdpy', u2, 'Cpx', dx, 'JW', True)
            self.add_coupling(-np.conj(t2dash_phi), u2, 'Cdpx', u1, 'Cpy', -dx, 'JW', True)  # h.c.


class FermionicTBG6Chain(FermionicTBG6Model, NearestNeighborModel):

    def __init__(self, model_params):
        model_params.setdefault('lattice', "Chain")
        CouplingMPOModel.__init__(self, model_params)
