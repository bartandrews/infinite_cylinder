"""Spin valley model based on two-orbital Hubbard model with Hund's coupling.
Hamiltonian based on: "Ferromagnetism and Spin-Valley states in Moire Correlated Insulators"
"""

import numpy as np

from tenpy.models.model import CouplingMPOModel, NearestNeighborModel
from tenpy.tools.params import get_parameter
from tenpy.networks.site import SpinSite, GroupedSite
from tenpy.networks import site


class FermionicTBG5Model(CouplingMPOModel):

    def __init__(self, model_params):
        CouplingMPOModel.__init__(self, model_params)

    def init_sites(self, model_params):

        conserve = get_parameter(model_params, 'conserve', 'Sz', self.name)
        ss = SpinSite(conserve=conserve)

        gs = GroupedSite([ss, ss], labels=['spin', 'valley'], charges='same')
        # gs.add_op('Ntot', gs.Ntotpx + gs.Ntotpy, False)

        site.multi_sites_combine_charges([ss, ss], same_charges=[[(0, 0), (1, 0)]])

        print(sorted(gs.opnames))
        print(gs.state_labels)

        return gs

    def init_terms(self, model_params):

        J = get_parameter(model_params, 'J', 1., self.name, True)
        Js = get_parameter(model_params, 'Js', 1., self.name, True)
        Jv = get_parameter(model_params, 'Jv', 1., self.name, True)

        for u1, u2, dx in self.lat.nearest_neighbors:

            # term 1
            self.add_coupling(1 / 4 * J, u1, 'Spspin Smvalley', u2, 'Smspin Spvalley', dx)
            self.add_coupling(np.conj(1 / 4 * J), u2, 'Spspin Smvalley', u1, 'Smspin Spvalley', -dx)  # H.c.
            self.add_coupling(1 / 4 * J, u1, 'Smspin Smvalley', u2, 'Spspin Spvalley', dx)
            self.add_coupling(np.conj(1 / 4 * J), u2, 'Spspin Spvalley', u1, 'Smspin Smvalley', -dx)  # H.c.
            self.add_coupling(1 / 2 * J, u1, 'Spspin Szvalley', u2, 'Smspin Szvalley', dx)
            self.add_coupling(np.conj(1 / 2 * J), u2, 'Smspin Szvalley', u1, 'Spspin Szvalley', -dx)  # H.c.
            self.add_coupling(1 / 2 * J, u1, 'Szspin Spvalley', u2, 'Szspin Smvalley', dx)
            self.add_coupling(np.conj(1 / 2 * J), u2, 'Szspin Smvalley', u1, 'Szspin Spvalley', -dx)  # H.c.
            self.add_coupling(J, u1, 'Szspin Szvalley', u2, 'Szspin Szvalley', dx)

            factor = 3

            # term 2
            self.add_coupling(factor*(Js+Js)/4., u1, 'Spspin', u2, 'Smspin', dx)
            self.add_coupling(factor*np.conj((Js+Js)/4.), u2, 'Spspin', u1, 'Smspin', -dx)  # h.c.
            # self.add_coupling(factor*(Js-Js)/4., u1, 'Spspin', u2, 'Spspin', dx)
            # self.add_coupling(factor*np.conj((Js-Js)/4.), u2, 'Smspin', u1, 'Smspin', -dx)  # h.c.
            self.add_coupling(factor*Js, u1, 'Szspin', u2, 'Szspin', dx)

            # term 3
            self.add_coupling(factor*(Jv+Jv) / 4., u1, 'Spvalley', u2, 'Smvalley', dx)
            self.add_coupling(factor*np.conj((Jv+Jv) / 4.), u2, 'Spvalley', u1, 'Smvalley', -dx)  # h.c.
            # self.add_coupling(factor*(Jv-Jv) / 4., u1, 'Spvalley', u2, 'Spvalley', dx)
            # self.add_coupling(factor*np.conj((Jv-Jv) / 4.), u2, 'Smvalley', u1, 'Smvalley', -dx)  # h.c.
            self.add_coupling(factor*Jv, u1, 'Szvalley', u2, 'Szvalley', dx)


class FermionicTBG5Chain(FermionicTBG5Model, NearestNeighborModel):

    def __init__(self, model_params):
        model_params.setdefault('lattice', "Chain")
        CouplingMPOModel.__init__(self, model_params)
