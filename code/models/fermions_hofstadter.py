"""Spinless fermion Hofstadter model.
Hamiltonian based on: "Square Lattice with Magnetic Field", Aidelsburger PhD thesis."""

import numpy as np
import sys

from tenpy.models.model import CouplingMPOModel, NearestNeighborModel
from tenpy.tools.params import get_parameter
from tenpy.networks.site import FermionSite
from lattices.MagneticSquare import MagneticSquare


class FermionicHofstadterModel(CouplingMPOModel):

    def __init__(self, model_params):

        CouplingMPOModel.__init__(self, model_params)

    def init_sites(self, model_params):

        conserve = get_parameter(model_params, 'conserve', 'N', self.name)
        site = FermionSite(conserve=conserve)
        return site

    def init_lattice(self, model_params):

        choice = get_parameter(model_params, 'lattice', 'MagneticSquare', self.name)

        if choice != 'MagneticSquare':
            sys.exit("Error: Please choose the MagneticSquare for hofstadter.")

        Lx = get_parameter(model_params, 'Lx', 1, self.name)
        Ly = get_parameter(model_params, 'Ly', 2, self.name)

        fs = self.init_sites(model_params)

        lat = MagneticSquare(Lx, Ly, fs)

        print(lat.N_sites)

        return lat

    def init_terms(self, model_params):

        t = get_parameter(model_params, 't', -1., self.name, True)
        V = get_parameter(model_params, 'V', 0, self.name, True)
        mu = get_parameter(model_params, 'mu', 0., self.name, True)
        phi_ext = 2*np.pi*get_parameter(model_params, 'phi_ext', 0., self.name)

        alpha = 1/3

        for u1, u2, dx in self.lat.NN_vert:
            t_phi = self.coupling_strength_add_ext_flux(t, dx, [0, phi_ext])
            self.add_coupling(t_phi, u1, 'Cd', u2, 'C', dx, 'JW', True)
            self.add_coupling(np.conj(t_phi), u2, 'Cd', u1, 'C', -dx, 'JW', True)  # h.c.

        for u1, u2, dx in self.lat.NN_h0:
            t_phi = self.coupling_strength_add_ext_flux(t, dx, [0, phi_ext]) * np.exp(-1j*2*np.pi*alpha*0)
            self.add_coupling(t_phi, u1, 'Cd', u2, 'C', dx, 'JW', True)
            self.add_coupling(np.conj(t_phi), u2, 'Cd', u1, 'C', -dx, 'JW', True)  # h.c.

        for u1, u2, dx in self.lat.NN_h1:
            t_phi = self.coupling_strength_add_ext_flux(t, dx, [0, phi_ext]) * np.exp(-1j*2*np.pi*alpha*1)
            self.add_coupling(t_phi, u1, 'Cd', u2, 'C', dx, 'JW', True)
            self.add_coupling(np.conj(t_phi), u2, 'Cd', u1, 'C', -dx, 'JW', True)  # h.c.

        for u1, u2, dx in self.lat.NN_h2:
            t_phi = self.coupling_strength_add_ext_flux(t, dx, [0, phi_ext]) * np.exp(-1j*2*np.pi*alpha*2)
            self.add_coupling(t_phi, u1, 'Cd', u2, 'C', dx, 'JW', True)
            self.add_coupling(np.conj(t_phi), u2, 'Cd', u1, 'C', -dx, 'JW', True)  # h.c.


class FermionicHofstadterChain(FermionicHofstadterModel, NearestNeighborModel):

    def __init__(self, model_params):
        model_params.setdefault('lattice', "Chain")
        CouplingMPOModel.__init__(self, model_params)
