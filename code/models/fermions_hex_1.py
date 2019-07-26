"""Spinless fermion hex_1 model."""

import numpy as np
import sys

from tenpy.models.model import CouplingMPOModel, NearestNeighborModel
from tenpy.tools.params import get_parameter
from tenpy.networks.site import FermionSite
from lattices.MagneticHoneycomb import MagneticHoneycomb


class FermionicHex1Model(CouplingMPOModel):

    def __init__(self, model_params):

        CouplingMPOModel.__init__(self, model_params)

    def init_sites(self, model_params):

        conserve = get_parameter(model_params, 'conserve', 'N', self.name)
        site = FermionSite(conserve=conserve)
        return site

    def init_lattice(self, model_params):

        choice = get_parameter(model_params, 'lattice', 'MagneticHoneycomb', self.name)

        if choice != 'MagneticHoneycomb':
            sys.exit("Error: Please choose the MagneticHoneycomb for hex_1.")

        Lx = get_parameter(model_params, 'Lx', 1, self.name)
        Ly = get_parameter(model_params, 'Ly', 3, self.name)

        fs = self.init_sites(model_params)

        lat = MagneticHoneycomb(Lx, Ly, fs)

        print(lat.N_sites)

        return lat

    def init_terms(self, model_params):

        t = get_parameter(model_params, 't', -1., self.name, True)
        V = get_parameter(model_params, 'V', 0, self.name, True)
        mu = get_parameter(model_params, 'mu', 0., self.name, True)
        phi_ext = 2*np.pi*get_parameter(model_params, 'phi_ext', 0., self.name)

        alpha = 2/5

        # down
        for u1, u2, dx in self.lat.NN0d:
            t_phi = self.coupling_strength_add_ext_flux(t, dx, [0, phi_ext]) * np.exp(-1j*(2*np.pi/3)*alpha*0)
            self.add_coupling(t_phi, u1, 'Cd', u2, 'C', dx, 'JW', True)
            self.add_coupling(np.conj(t_phi), u2, 'Cd', u1, 'C', -dx, 'JW', True)  # h.c.
        for u1, u2, dx in self.lat.NN2d:
            t_phi = self.coupling_strength_add_ext_flux(t, dx, [0, phi_ext]) * np.exp(-1j*(2*np.pi/3)*alpha*2)
            self.add_coupling(t_phi, u1, 'Cd', u2, 'C', dx, 'JW', True)
            self.add_coupling(np.conj(t_phi), u2, 'Cd', u1, 'C', -dx, 'JW', True)  # h.c.
        for u1, u2, dx in self.lat.NN4d:
            t_phi = self.coupling_strength_add_ext_flux(t, dx, [0, phi_ext]) * np.exp(-1j*(2*np.pi/3)*alpha*4)
            self.add_coupling(t_phi, u1, 'Cd', u2, 'C', dx, 'JW', True)
            self.add_coupling(np.conj(t_phi), u2, 'Cd', u1, 'C', -dx, 'JW', True)  # h.c.
        for u1, u2, dx in self.lat.NN6d:
            t_phi = self.coupling_strength_add_ext_flux(t, dx, [0, phi_ext]) * np.exp(-1j*(2*np.pi/3)*alpha*6)
            self.add_coupling(t_phi, u1, 'Cd', u2, 'C', dx, 'JW', True)
            self.add_coupling(np.conj(t_phi), u2, 'Cd', u1, 'C', -dx, 'JW', True)  # h.c.
        for u1, u2, dx in self.lat.NN8d:
            t_phi = self.coupling_strength_add_ext_flux(t, dx, [0, phi_ext]) * np.exp(-1j*(2*np.pi/3)*alpha*8)
            self.add_coupling(t_phi, u1, 'Cd', u2, 'C', dx, 'JW', True)
            self.add_coupling(np.conj(t_phi), u2, 'Cd', u1, 'C', -dx, 'JW', True)  # h.c.

        # upper left
        for u1, u2, dx in self.lat.NN0ul:
            t_phi = self.coupling_strength_add_ext_flux(t, dx, [0, phi_ext]) * np.exp(1j*(np.pi/3)*alpha*(0-1/2))
            self.add_coupling(t_phi, u1, 'Cd', u2, 'C', dx, 'JW', True)
            self.add_coupling(np.conj(t_phi), u2, 'Cd', u1, 'C', -dx, 'JW', True)  # h.c.
        for u1, u2, dx in self.lat.NN2ul:
            t_phi = self.coupling_strength_add_ext_flux(t, dx, [0, phi_ext]) * np.exp(1j*(np.pi/3)*alpha*(2-1/2))
            self.add_coupling(t_phi, u1, 'Cd', u2, 'C', dx, 'JW', True)
            self.add_coupling(np.conj(t_phi), u2, 'Cd', u1, 'C', -dx, 'JW', True)  # h.c.
        for u1, u2, dx in self.lat.NN4ul:
            t_phi = self.coupling_strength_add_ext_flux(t, dx, [0, phi_ext]) * np.exp(1j*(np.pi/3)*alpha*(4-1/2))
            self.add_coupling(t_phi, u1, 'Cd', u2, 'C', dx, 'JW', True)
            self.add_coupling(np.conj(t_phi), u2, 'Cd', u1, 'C', -dx, 'JW', True)  # h.c.
        for u1, u2, dx in self.lat.NN6ul:
            t_phi = self.coupling_strength_add_ext_flux(t, dx, [0, phi_ext]) * np.exp(1j*(np.pi/3)*alpha*(6-1/2))
            self.add_coupling(t_phi, u1, 'Cd', u2, 'C', dx, 'JW', True)
            self.add_coupling(np.conj(t_phi), u2, 'Cd', u1, 'C', -dx, 'JW', True)  # h.c.
        for u1, u2, dx in self.lat.NN8ul:
            t_phi = self.coupling_strength_add_ext_flux(t, dx, [0, phi_ext]) * np.exp(1j*(np.pi/3)*alpha*(8-1/2))
            self.add_coupling(t_phi, u1, 'Cd', u2, 'C', dx, 'JW', True)
            self.add_coupling(np.conj(t_phi), u2, 'Cd', u1, 'C', -dx, 'JW', True)  # h.c.

        # upper right
        for u1, u2, dx in self.lat.NN0ur:
            t_phi = self.coupling_strength_add_ext_flux(t, dx, [0, phi_ext]) * np.exp(1j*(np.pi/3)*alpha*(0+1/2))
            self.add_coupling(t_phi, u1, 'Cd', u2, 'C', dx, 'JW', True)
            self.add_coupling(np.conj(t_phi), u2, 'Cd', u1, 'C', -dx, 'JW', True)  # h.c.
        for u1, u2, dx in self.lat.NN2ur:
            t_phi = self.coupling_strength_add_ext_flux(t, dx, [0, phi_ext]) * np.exp(1j*(np.pi/3)*alpha*(2+1/2))
            self.add_coupling(t_phi, u1, 'Cd', u2, 'C', dx, 'JW', True)
            self.add_coupling(np.conj(t_phi), u2, 'Cd', u1, 'C', -dx, 'JW', True)  # h.c.
        for u1, u2, dx in self.lat.NN4ur:
            t_phi = self.coupling_strength_add_ext_flux(t, dx, [0, phi_ext]) * np.exp(1j*(np.pi/3)*alpha*(4+1/2))
            self.add_coupling(t_phi, u1, 'Cd', u2, 'C', dx, 'JW', True)
            self.add_coupling(np.conj(t_phi), u2, 'Cd', u1, 'C', -dx, 'JW', True)  # h.c.
        for u1, u2, dx in self.lat.NN6ur:
            t_phi = self.coupling_strength_add_ext_flux(t, dx, [0, phi_ext]) * np.exp(1j*(np.pi/3)*alpha*(6+1/2))
            self.add_coupling(t_phi, u1, 'Cd', u2, 'C', dx, 'JW', True)
            self.add_coupling(np.conj(t_phi), u2, 'Cd', u1, 'C', -dx, 'JW', True)  # h.c.
        for u1, u2, dx in self.lat.NN8ur:
            t_phi = self.coupling_strength_add_ext_flux(t, dx, [0, phi_ext]) * np.exp(1j*(np.pi/3)*alpha*(8+1/2))
            self.add_coupling(t_phi, u1, 'Cd', u2, 'C', dx, 'JW', True)
            self.add_coupling(np.conj(t_phi), u2, 'Cd', u1, 'C', -dx, 'JW', True)  # h.c.


class FermionicHex1Chain(FermionicHex1Model, NearestNeighborModel):

    def __init__(self, model_params):
        model_params.setdefault('lattice', "Chain")
        CouplingMPOModel.__init__(self, model_params)
