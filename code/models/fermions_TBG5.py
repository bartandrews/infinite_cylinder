"""Spin valley model based on two-orbital Hubbard model with Hund's coupling.
Hamiltonian based on: "Ferromagnetism and Spin-Valley states in Moire Correlated Insulators"
"""

import numpy as np

from tenpy.models.model import CouplingMPOModel, NearestNeighborModel
from tenpy.tools.params import get_parameter
from tenpy.networks.site import SpinSite, GroupedSite


class FermionicTBG5Model(CouplingMPOModel):

    def __init__(self, model_params):
        CouplingMPOModel.__init__(self, model_params)

    def init_sites(self, model_params):

        cons_parity = get_parameter(model_params, 'cons_parity', 'parity', self.name)
        cons_Sz = get_parameter(model_params, 'cons_Sz', 'Sz', self.name)
        ss = SpinSite(cons_parity=cons_parity, cons_Sz=cons_Sz)

        gs = GroupedSite([ss, ss], labels=['spin', 'valley'], charges='same')
        # gs.add_op('Ntot', gs.Ntotpx + gs.Ntotpy, False)

        print(sorted(gs.opnames))
        print(gs.state_labels)

        return gs

    def init_terms(self, model_params):

        J = get_parameter(model_params, 'J', 1., self.name, True)
        Js = get_parameter(model_params, 'Js', 1., self.name, True)
        Jv = get_parameter(model_params, 'Jv', 1., self.name, True)

        for u1, u2, dx in self.lat.nearest_neighbors:

            # term 1
            self.add_coupling(J, u1, 'Sxspin Sxvalley', u2, 'Sxspin Sxvalley', dx)
            self.add_coupling(J, u1, 'Sxspin Syvalley', u2, 'Sxspin Syvalley', dx)
            self.add_coupling(J, u1, 'Sxspin Szvalley', u2, 'Sxspin Szvalley', dx)

            self.add_coupling(J, u1, 'Syspin Sxvalley', u2, 'Syspin Sxvalley', dx)
            self.add_coupling(J, u1, 'Syspin Syvalley', u2, 'Syspin Syvalley', dx)
            self.add_coupling(J, u1, 'Syspin Szvalley', u2, 'Syspin Szvalley', dx)

            self.add_coupling(J, u1, 'Szspin Sxvalley', u2, 'Szspin Sxvalley', dx)
            self.add_coupling(J, u1, 'Szspin Syvalley', u2, 'Szspin Syvalley', dx)
            self.add_coupling(J, u1, 'Szspin Szvalley', u2, 'Szspin Szvalley', dx)

            # term 2
            self.add_coupling(3*Js, u1, 'Sxspin', u2, 'Sxspin', dx)
            self.add_coupling(3*Js, u1, 'Syspin', u2, 'Syspin', dx)
            self.add_coupling(3*Js, u1, 'Szspin', u2, 'Szspin', dx)

            # term 3
            self.add_coupling(3*Jv, u1, 'Sxvalley', u2, 'Sxvalley', dx)
            self.add_coupling(3*Jv, u1, 'Syvalley', u2, 'Syvalley', dx)
            self.add_coupling(3*Jv, u1, 'Szvalley', u2, 'Szvalley', dx)


class FermionicTBG5Chain(FermionicTBG5Model, NearestNeighborModel):

    def __init__(self, model_params):
        model_params.setdefault('lattice', "Chain")
        CouplingMPOModel.__init__(self, model_params)
