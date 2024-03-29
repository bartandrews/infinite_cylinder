"""Spinful fermions with two orbitals.
Hamiltonian based on: "Model for the metal-insulator transition in graphene superlattices and beyond" (Sec. IV)."""

import numpy as np

from tenpy.models.model import CouplingMPOModel, NearestNeighborModel
from tenpy.tools.params import get_parameter
from tenpy.networks.site import SpinHalfFermionSite, GroupedSite


class FermionicTBG1Model(CouplingMPOModel):

    def __init__(self, model_params):
        CouplingMPOModel.__init__(self, model_params)

    def init_sites(self, model_params):

        cons_N = get_parameter(model_params, 'cons_N', 'N', self.name)
        cons_Sz = get_parameter(model_params, 'cons_Sz', 'Sz', self.name)
        shfs = SpinHalfFermionSite(cons_N=cons_N, cons_Sz=cons_Sz)

        gs = GroupedSite([shfs, shfs], labels=['px', 'py'], charges='same')
        gs.add_op('Ntot', gs.Ntotpx + gs.Ntotpy, False)

        return gs

    def init_terms(self, model_params):

        t = get_parameter(model_params, 't', 1., self.name, True)
        U = get_parameter(model_params, 'U', 0, self.name, True)
        V = get_parameter(model_params, 'V', 0, self.name, True)
        mu = get_parameter(model_params, 'mu', 0., self.name, True)

        t1 = 0.331*t  # real
        t2 = (-0.010 + 0.097 * 1j)*t  # complex

        for u in range(len(self.lat.unit_cell)):

            self.add_onsite(-mu, 0, 'Ntot')
            self.add_onsite(U/2, 0, 'Ntot Ntot')
            self.add_onsite(-U, 0, 'Ntot')

        for u1, u2, dx in self.lat.nearest_neighbors:

            # spin up
            self.add_coupling(t1, u1, 'Cdupx', u2, 'Cupx', dx, 'JW', True)
            self.add_coupling(t1, u2, 'Cdupx', u1, 'Cupx', -dx, 'JW', True)  # h.c.
            self.add_coupling(t1, u1, 'Cdupy', u2, 'Cupy', dx, 'JW', True)
            self.add_coupling(t1, u2, 'Cdupy', u1, 'Cupy', -dx, 'JW', True)  # h.c.
            # spin down
            self.add_coupling(t1, u1, 'Cddpx', u2, 'Cdpx', dx, 'JW', True)
            self.add_coupling(t1, u2, 'Cddpx', u1, 'Cdpx', -dx, 'JW', True)  # h.c.
            self.add_coupling(t1, u1, 'Cddpy', u2, 'Cdpy', dx, 'JW', True)
            self.add_coupling(t1, u2, 'Cddpy', u1, 'Cdpy', -dx, 'JW', True)  # h.c.

        for u1, u2, dx in self.lat.next_nearest_neighbors:

            # spin up
            self.add_coupling(np.real(t2), u1, 'Cdupx', u2, 'Cupx', dx, 'JW', True)
            self.add_coupling(np.real(t2), u2, 'Cdupx', u1, 'Cupx', -dx, 'JW', True)  # h.c.
            self.add_coupling(np.real(t2), u1, 'Cdupy', u2, 'Cupy', dx, 'JW', True)
            self.add_coupling(np.real(t2), u2, 'Cdupy', u1, 'Cupy', -dx, 'JW', True)  # h.c.
            # spin down
            self.add_coupling(np.real(t2), u1, 'Cddpx', u2, 'Cdpx', dx, 'JW', True)
            self.add_coupling(np.real(t2), u2, 'Cddpx', u1, 'Cdpx', -dx, 'JW', True)  # h.c.
            self.add_coupling(np.real(t2), u1, 'Cddpy', u2, 'Cdpy', dx, 'JW', True)
            self.add_coupling(np.real(t2), u2, 'Cddpy', u1, 'Cdpy', -dx, 'JW', True)  # h.c.

            # spin up
            self.add_coupling(np.imag(t2), u1, 'Cdupx', u2, 'Cupy', dx, 'JW', True)
            self.add_coupling(np.imag(t2), u2, 'Cdupy', u1, 'Cupx', -dx, 'JW', True)  # h.c.
            self.add_coupling(-np.imag(t2), u1, 'Cdupy', u2, 'Cupx', dx, 'JW', True)
            self.add_coupling(-np.imag(t2), u2, 'Cdupx', u1, 'Cupy', -dx, 'JW', True)  # h.c.
            # spin down
            self.add_coupling(np.imag(t2), u1, 'Cddpx', u2, 'Cdpy', dx, 'JW', True)
            self.add_coupling(np.imag(t2), u2, 'Cddpy', u1, 'Cdpx', -dx, 'JW', True)  # h.c.
            self.add_coupling(-np.imag(t2), u1, 'Cddpy', u2, 'Cdpx', dx, 'JW', True)
            self.add_coupling(-np.imag(t2), u2, 'Cddpx', u1, 'Cdpy', -dx, 'JW', True)  # h.c.

            self.add_coupling(V, u1, 'Ntot', u2, 'Ntot', dx, 'JW', True)


class FermionicTBG1Chain(FermionicTBG1Model, NearestNeighborModel):

    def __init__(self, model_params):
        model_params.setdefault('lattice', "Chain")
        CouplingMPOModel.__init__(self, model_params)
