"""Spinless fermion chiral-pi-flux model.
Hamiltonian based on: "Fractional Quantum Hall States at Zero Magnetic Field"."""

import numpy as np
import sys

from tenpy.models.model import CouplingMPOModel, NearestNeighborModel
from tenpy.tools.params import get_parameter
from tenpy.networks.site import FermionSite
from lattices.BipartiteSquare import BipartiteSquare


class FermionicPiFluxModel(CouplingMPOModel):

    def __init__(self, model_params):

        CouplingMPOModel.__init__(self, model_params)

    def init_sites(self, model_params):

        conserve = get_parameter(model_params, 'conserve', 'N', self.name)
        site = FermionSite(conserve=conserve)
        return site

    def init_lattice(self, model_params):

        choice = get_parameter(model_params, 'lattice', 'BipartiteSquare', self.name)

        if choice != 'BipartiteSquare':
            sys.exit("Error: Please choose the BipartiteSquare for pi_flux.")

        Lx = get_parameter(model_params, 'Lx', 1, self.name)
        Ly = get_parameter(model_params, 'Ly', 3, self.name)

        fs = self.init_sites(model_params)

        lat = BipartiteSquare(Lx, Ly, fs)

        print(lat.N_sites)

        return lat

    def init_terms(self, model_params):

        t = get_parameter(model_params, 't', -1., self.name, True)
        V = get_parameter(model_params, 'V', 0, self.name, True)
        mu = get_parameter(model_params, 'mu', 0., self.name, True)
        phi_ext = 2*np.pi*get_parameter(model_params, 'phi_ext', 0., self.name)

        t1 = t * np.exp(1j * np.pi/4)
        t2 = t / np.sqrt(2)

        for u in range(len(self.lat.unit_cell)):

            self.add_onsite(mu, 0, 'N', category='mu N')
            self.add_onsite(-mu, 1, 'N', category='mu N')

        for u1, u2, dx in self.lat.NN:

            t1_phi = self.coupling_strength_add_ext_flux(t1, dx, [0, phi_ext])
            self.add_coupling(t1_phi, u1, 'Cd', u2, 'C', dx, 'JW', True, category='t1 Cd_i C_j')
            self.add_coupling(np.conj(t1_phi), u2, 'Cd', u1, 'C', -dx, 'JW', True, category='t1 Cd_i C_j h.c.')

        for u1, u2, dx in self.lat.nNNdashed:

            t2_phi = self.coupling_strength_add_ext_flux(t2, dx, [0, phi_ext])
            self.add_coupling(t2_phi, u1, 'Cd', u2, 'C', dx, 'JW', True, category='t2 Cd_i C_j')
            self.add_coupling(np.conj(t2_phi), u2, 'Cd', u1, 'C', -dx, 'JW', True, category='t2 Cd_i C_j h.c.')

        for u1, u2, dx in self.lat.nNNdotted:

            t2_phi = self.coupling_strength_add_ext_flux(t2, dx, [0, phi_ext])
            self.add_coupling(-t2_phi, u1, 'Cd', u2, 'C', dx, 'JW', True, category='-t2 Cd_i C_j')
            self.add_coupling(-np.conj(t2_phi), u2, 'Cd', u1, 'C', -dx, 'JW', True, category='-t2 Cd_i C_j h.c.')


class FermionicPiFluxChain(FermionicPiFluxModel, NearestNeighborModel):

    def __init__(self, model_params):
        model_params.setdefault('lattice', "Chain")
        CouplingMPOModel.__init__(self, model_params)
