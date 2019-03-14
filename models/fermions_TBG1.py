"""Spinful fermions with 2 orbitals, Hamiltonian based on the paper: "Model for the metal-insulator transition in
graphene superlattices and beyond" (Sec. IV)
"""

from tenpy.models.model import CouplingMPOModel, NearestNeighborModel
from tenpy.tools.params import get_parameter
from tenpy.networks.site import SpinHalfFermionSite, GroupedSite
import numpy as np


class FermionicTBG1Model(CouplingMPOModel):

    def __init__(self, model_params):
        CouplingMPOModel.__init__(self, model_params)

    def init_sites(self, model_params):

        cons_N = get_parameter(model_params, 'cons_N', 'N', self.name)
        cons_Sz = get_parameter(model_params, 'cons_Sz', 'Sz', self.name)
        shfs = SpinHalfFermionSite(cons_N=cons_N, cons_Sz=cons_Sz)

        gs = GroupedSite([shfs, shfs], labels=['px', 'py'], charges='same')
        gs.add_op('Ntot', gs.Ntotpx + gs.Ntotpy, False)

        print(sorted(gs.opnames))
        print(gs.state_labels)

        return gs

    def init_terms(self, model_params):

        t = get_parameter(model_params, 't', 1., self.name)
        U = get_parameter(model_params, 'U', 0, self.name)
        V = get_parameter(model_params, 'V', 0, self.name)
        mu = get_parameter(model_params, 'mu', 0., self.name)

        t1 = 0.331*t  # real
        t2 = (-0.010 + 0.097 * 1j)*t  # complex

        for u in range(len(self.lat.unit_cell)):

            self.add_onsite(-mu, 0, 'Ntot')
            self.add_onsite(U/2, 0, 'Ntot Ntot')
            self.add_onsite(-U, 0, 'Ntot')

        for u1, u2, dx in self.lat.nearest_neighbors:

            self.add_coupling(t1, u1, 'Cdupx', u2, 'Cupx', dx, 'JW', True)
            self.add_coupling(t1, u2, 'Cdupx', u1, 'Cupx', -dx, 'JW', True)  # h.c.
            self.add_coupling(t1, u1, 'Cdupy', u2, 'Cupy', dx, 'JW', True)
            self.add_coupling(t1, u2, 'Cdupy', u1, 'Cupy', -dx, 'JW', True)  # h.c.

        for u1, u2, dx in self.lat.next_nearest_neighbors:

            self.add_coupling(np.real(t2), u1, 'Cdupx', u2, 'Cupx', dx, 'JW', True)
            self.add_coupling(np.real(t2), u2, 'Cdupx', u1, 'Cupx', -dx, 'JW', True)  # h.c.
            self.add_coupling(np.real(t2), u1, 'Cdupy', u2, 'Cupy', dx, 'JW', True)
            self.add_coupling(np.real(t2), u2, 'Cdupy', u1, 'Cupy', -dx, 'JW', True)  # h.c.

            self.add_coupling(np.imag(t2), u1, 'Cdupx', u2, 'Cupy', dx, 'JW', True)
            self.add_coupling(np.imag(t2), u2, 'Cdupy', u1, 'Cupx', -dx, 'JW', True)  # h.c.
            self.add_coupling(-np.imag(t2), u1, 'Cdupy', u2, 'Cupx', dx, 'JW', True)
            self.add_coupling(-np.imag(t2), u2, 'Cdupx', u1, 'Cupy', -dx, 'JW', True)  # h.c.

            self.add_coupling(V, u1, 'Ntot', u2, 'Ntot', dx, 'JW', True)


class FermionicTBG1Chain(FermionicTBG1Model, NearestNeighborModel):

    def __init__(self, model_params):
        model_params.setdefault('lattice', "Chain")
        CouplingMPOModel.__init__(self, model_params)
