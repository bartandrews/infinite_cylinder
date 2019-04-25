"""Five-band model of twisted bilayer graphene.
Hamiltonian based on: "Faithful Tight-binding Models and Fragile Topology of Magic-angle Bilayer Graphene"."""

import numpy as np
import sys

from tenpy.models.model import CouplingMPOModel, NearestNeighborModel
from tenpy.tools.params import get_parameter
from tenpy.networks.site import FermionSite, GroupedSite
from lattices.five_band_model import FiveBandLattice
from tenpy.networks import site


class FermionicTBG4Model(CouplingMPOModel):

    def __init__(self, model_params):
        CouplingMPOModel.__init__(self, model_params)

    def init_sites(self, model_params):

        conserve = get_parameter(model_params, 'conserve', 'N', self.name)
        fs = FermionSite(conserve=conserve)

        gs = GroupedSite([fs, fs, fs], labels=['pz', 'pp', 'pm'], charges='same')
        gs.add_op('Ntot', gs.Npz + gs.Npp + gs.Npm, False)

        site.multi_sites_combine_charges([fs, gs], same_charges=[[(0, 'N'), (1, 'N')]])

        return fs, gs

    def init_lattice(self, model_params):

        choice = get_parameter(model_params, 'lattice', 'FiveBandLattice', self.name)

        if choice != 'FiveBandLattice':
            sys.exit("Error: Please choose the FiveBandLattice for TBG4.")

        Lx = get_parameter(model_params, 'Lx', 1, self.name)
        Ly = get_parameter(model_params, 'Ly', 3, self.name)

        (fs, gs) = self.init_sites(model_params)

        lat = FiveBandLattice(Lx, Ly, gs, fs)

        print(lat.N_sites)

        return lat

    def init_terms(self, model_params):

        t = get_parameter(model_params, 't', 1., self.name)
        mu = get_parameter(model_params, 'mu', 0., self.name)
        U = get_parameter(model_params, 'U', 0., self.name)
        V = get_parameter(model_params, 'V', 0., self.name)

        a = 0.25 * t  # in units of 80 meV
        b = 0.2 * t
        c = 0.1 * t
        d = 0.67 * t  # five_band_model lattice assumes d is real

        field = np.exp(1j*(2*np.pi)*(1/3))  # magnetic field via Peierls substitution (Phi in units of Phi_0)

        for u1, u2, dx in self.lat.a1_d:

            self.add_coupling(field*(-a*1j)*d, u1, 'Cpz', u2, 'Cd', dx, 'JW', True)
            self.add_coupling(np.conj(field*(-a*1j)*d), u2, 'C', u1, 'Cdpz', -dx, 'JW', True)  # h.c.

            self.add_coupling(field*c*d, u1, 'Cpp', u2, 'Cd', dx, 'JW', True)
            self.add_coupling(np.conj(field*c*d), u2, 'C', u1, 'Cdpp', -dx, 'JW', True)  # h.c.

            self.add_coupling(field*b*d, u1, 'Cpm', u2, 'Cd', dx, 'JW', True)
            self.add_coupling(np.conj(field*b*d), u2, 'C', u1, 'Cdpm', -dx, 'JW', True)  # h.c.

        for u1, u2, dx in self.lat.a2_d:

            self.add_coupling(field*(a*1j)*d, u1, 'Cpz', u2, 'Cd', dx, 'JW', True)
            self.add_coupling(np.conj(field*(a*1j)*d), u2, 'C', u1, 'Cdpz', -dx, 'JW', True)  # h.c.

            self.add_coupling(field*b*d, u1, 'Cpp', u2, 'Cd', dx, 'JW', True)
            self.add_coupling(np.conj(field*b*d), u2, 'C', u1, 'Cdpp', -dx, 'JW', True)  # h.c.

            self.add_coupling(field*c*d, u1, 'Cpm', u2, 'Cd', dx, 'JW', True)
            self.add_coupling(np.conj(field*c*d), u2, 'C', u1, 'Cdpm', -dx, 'JW', True)  # h.c.

        for u1, u2, dx in self.lat.d_d:

            self.add_coupling(field*d*d, u1, 'C', u2, 'Cd', dx, 'JW', True)
            self.add_coupling(np.conj(field*d*d), u2, 'C', u1, 'Cd', -dx, 'JW', True)  # h.c.


class FermionicTBG4Chain(FermionicTBG4Model, NearestNeighborModel):

    def __init__(self, model_params):
        model_params.setdefault('lattice', "Chain")
        CouplingMPOModel.__init__(self, model_params)
