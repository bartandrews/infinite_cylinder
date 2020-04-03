"""Spinful fermions with two orbitals.
Hamiltonian based on: "Spin/orbital density wave and Mott insulator in two-orbital Hubbard model on honeycomb lattice"
"""

import numpy as np

from tenpy.models.model import CouplingMPOModel, NearestNeighborModel
from tenpy.tools.params import get_parameter
from tenpy.networks.site import SpinHalfFermionSite, GroupedSite


class FermionicTBG3Model(CouplingMPOModel):

    def __init__(self, model_params):
        CouplingMPOModel.__init__(self, model_params)

    def init_sites(self, model_params):

        cons_N = get_parameter(model_params, 'cons_N', 'N', self.name)
        cons_Sz = get_parameter(model_params, 'cons_Sz', 'Sz', self.name)
        shfs = SpinHalfFermionSite(cons_N=cons_N, cons_Sz=cons_Sz)

        gs = GroupedSite([shfs, shfs], labels=['px', 'py'], charges='same')
        gs.add_op('Ntot', gs.Ntotpx + gs.Ntotpy, False)

        # print(sorted(gs.opnames))
        # print(gs.state_labels)

        return gs

    def init_terms(self, model_params):

        t = get_parameter(model_params, 't', -1., self.name, True)
        U = get_parameter(model_params, 'U', 0, self.name, True)
        V = get_parameter(model_params, 'V', 0, self.name, True)
        mu = get_parameter(model_params, 'mu', 0., self.name, True)

        for u in range(len(self.lat.unit_cell)):

            self.add_onsite(U, 0, 'Ntot Ntot')
            self.add_onsite(-2*U, 0, 'Ntot')
            self.add_onsite(U, 0, 'Id')

        for u1, u2, dx in self.lat.nearest_neighbors:

            # spin up
            self.add_coupling(t, u1, 'Cdupx', u2, 'Cupx', dx, 'JW', True)
            self.add_coupling(t, u2, 'Cdupx', u1, 'Cupx', -dx, 'JW', True)  # h.c.
            self.add_coupling(t, u1, 'Cdupy', u2, 'Cupy', dx, 'JW', True)
            self.add_coupling(t, u2, 'Cdupy', u1, 'Cupy', -dx, 'JW', True)  # h.c.
            # spin down
            self.add_coupling(t, u1, 'Cddpx', u2, 'Cdpx', dx, 'JW', True)
            self.add_coupling(t, u2, 'Cddpx', u1, 'Cdpx', -dx, 'JW', True)  # h.c.
            self.add_coupling(t, u1, 'Cddpy', u2, 'Cdpy', dx, 'JW', True)
            self.add_coupling(t, u2, 'Cddpy', u1, 'Cdpy', -dx, 'JW', True)  # h.c.


class FermionicTBG3Chain(FermionicTBG3Model, NearestNeighborModel):

    def __init__(self, model_params):
        model_params.setdefault('lattice', "Chain")
        CouplingMPOModel.__init__(self, model_params)
