"""Five-band model of twisted bilayer graphene.
Hamiltonian based on: "Faithful Tight-binding Models and Fragile Topology of Magic-angle Bilayer Graphene"."""

import numpy as np

from tenpy.models.model import CouplingMPOModel, NearestNeighborModel
from tenpy.tools.params import get_parameter
from tenpy.networks.site import FermionSite, GroupedSite


class FermionicTBG4Model(CouplingMPOModel):

    def __init__(self, model_params):
        CouplingMPOModel.__init__(self, model_params)

    def init_sites(self, model_params):

        conserve = get_parameter(model_params, 'conserve', 'N', self.name)
        fs = FermionSite(conserve=conserve)

        gs = GroupedSite([fs, fs, fs], labels=['pz', 'pp', 'pm'], charges='same')
        gs.add_op('Ntot', gs.Npz + gs.Npp + gs.Npm, False)

        return fs, gs

    def init_terms(self, model_params):

        t = get_parameter(model_params, 't', 1., self.name)
        mu = get_parameter(model_params, 'mu', 0., self.name)

        a = -0.25j * t
        b = 0.2 * t
        c = 0.1 * t
        d = 0.67 * t  # five_band_model lattice assumes d is real

        for u1, u2, dx in self.lat.a1_d:

            self.add_coupling(a*d, u1, 'Cd', u2, 'Cpz', dx, 'JW', True)
            self.add_coupling(np.conj(a)*d, u2, 'Cdpz', u1, 'C', -dx, 'JW', True)  # h.c.

            self.add_coupling(c*d, u1, 'Cd', u2, 'Cpp', dx, 'JW', True)
            self.add_coupling(np.conj(c)*d, u2, 'Cdpp', u1, 'C', -dx, 'JW', True)  # h.c.

            self.add_coupling(b*d, u1, 'Cd', u2, 'Cpm', dx, 'JW', True)
            self.add_coupling(np.conj(b)*d, u2, 'Cdpm', u1, 'C', -dx, 'JW', True)  # h.c.

        for u1, u2, dx in self.lat.a2_d:

            self.add_coupling(np.conj(a)*d, u1, 'Cd', u2, 'Cpz', dx, 'JW', True)
            self.add_coupling(a*d, u2, 'Cdpz', u1, 'C', -dx, 'JW', True)  # h.c.

            self.add_coupling(b*d, u1, 'Cd', u2, 'Cpp', dx, 'JW', True)
            self.add_coupling(np.conj(b)*d, u2, 'Cdpp', u1, 'C', -dx, 'JW', True)  # h.c.

            self.add_coupling(c*d, u1, 'Cd', u2, 'Cpm', dx, 'JW', True)
            self.add_coupling(np.conj(c)*d, u2, 'Cdpm', u1, 'C', -dx, 'JW', True)  # h.c.

        for u1, u2, dx in self.lat.d_d:

            self.add_coupling(d*d, u1, 'Cd', u2, 'C', dx, 'JW', True)
            self.add_coupling(d*d, u2, 'Cd', u1, 'C', -dx, 'JW', True)  # h.c.

        for u1, u2, dx in self.lat.a1_a2:

            # hop from pz orbital
            self.add_coupling(a*np.conj(a), u1, 'Cdpz', u2, 'Cpz', dx, 'JW', True)
            self.add_coupling(a*np.conj(a), u2, 'Cdpz', u1, 'Cpz', -dx, 'JW', True)  # h.c.
            self.add_coupling(a*np.conj(a), u1, 'Cdpp', u2, 'Cpz', dx, 'JW', True)
            self.add_coupling(a*np.conj(a), u2, 'Cdpz', u1, 'Cpp', -dx, 'JW', True)  # h.c.
            self.add_coupling(a*np.conj(a), u1, 'Cdpm', u2, 'Cpz', dx, 'JW', True)
            self.add_coupling(a*np.conj(a), u2, 'Cdpz', u1, 'Cpm', -dx, 'JW', True)  # h.c.

            # hop from pp orbital
            self.add_coupling(a*np.conj(a), u1, 'Cdpz', u2, 'Cpp', dx, 'JW', True)
            self.add_coupling(a*np.conj(a), u2, 'Cdpp', u1, 'Cpz', -dx, 'JW', True)  # h.c.
            self.add_coupling(a*np.conj(a), u1, 'Cdpp', u2, 'Cpp', dx, 'JW', True)
            self.add_coupling(a*np.conj(a), u2, 'Cdpp', u1, 'Cpp', -dx, 'JW', True)  # h.c.
            self.add_coupling(a*np.conj(a), u1, 'Cdpm', u2, 'Cpp', dx, 'JW', True)
            self.add_coupling(a*np.conj(a), u2, 'Cdpp', u1, 'Cpm', -dx, 'JW', True)  # h.c.

            # hop from pm orbital
            self.add_coupling(a*np.conj(a), u1, 'Cdpz', u2, 'Cpm', dx, 'JW', True)
            self.add_coupling(a*np.conj(a), u2, 'Cdpm', u1, 'Cpz', -dx, 'JW', True)  # h.c.
            self.add_coupling(a*np.conj(a), u1, 'Cdpp', u2, 'Cpm', dx, 'JW', True)
            self.add_coupling(a*np.conj(a), u2, 'Cdpm', u1, 'Cpp', -dx, 'JW', True)  # h.c.
            self.add_coupling(a*np.conj(a), u1, 'Cdpm', u2, 'Cpm', dx, 'JW', True)
            self.add_coupling(a*np.conj(a), u2, 'Cdpm', u1, 'Cpm', -dx, 'JW', True)  # h.c.


class FermionicTBG4Chain(FermionicTBG4Model, NearestNeighborModel):

    def __init__(self, model_params):
        model_params.setdefault('lattice', "Chain")
        CouplingMPOModel.__init__(self, model_params)
