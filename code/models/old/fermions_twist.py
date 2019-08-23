"""Spinless fermion twist model."""

import numpy as np
import sys

from tenpy.models.model import CouplingMPOModel
from tenpy.tools.params import get_parameter
from tenpy.networks.site import FermionSite
from lattices.MagneticHoneycomb import MagneticHoneycomb


class FermionicTwistModel(CouplingMPOModel):

    def __init__(self, model_params):
        CouplingMPOModel.__init__(self, model_params)

    def init_sites(self, model_params):
        conserve = get_parameter(model_params, 'conserve', 'N', self.name)
        site = FermionSite(conserve=conserve)
        return site

    def init_lattice(self, model_params):
        choice = get_parameter(model_params, 'lattice', 'MagneticTwist', self.name)
        if choice != 'MagneticTwist':
            sys.exit("Error: Please choose the MagneticTwist for twist.")
        Lx = get_parameter(model_params, 'Lx', 1, self.name)
        Ly = get_parameter(model_params, 'Ly', 1, self.name)
        fs = self.init_sites(model_params)
        lat = MagneticHoneycomb(Lx, Ly, fs)
        print(lat.N_sites)
        return lat

    def init_terms(self, model_params):
        t = get_parameter(model_params, 't', -1., self.name, True)
        phi_ext = 2*np.pi*get_parameter(model_params, 'phi_ext', 0., self.name)

        t1 = 0.331 * t
        t2 = -0.01

        p = 2
        q = 7

        alpha = p/q
        numb_sites = 2*q

        # NN ###
        for i in range(0, numb_sites, 2):

            for u1, u2, dx in getattr(self.lat, "NN{}d".format(i)):
                t_phi = self.coupling_strength_add_ext_flux(t1, dx, [0, phi_ext]) \
                        * np.exp(-1j * (2 * np.pi / 3) * alpha * i)
                self.add_coupling(t_phi, u1, 'Cd', u2, 'C', dx, 'JW', True)
                self.add_coupling(np.conj(t_phi), u2, 'Cd', u1, 'C', -dx, 'JW', True)  # h.c.

            for u1, u2, dx in getattr(self.lat, "NN{}ul".format(i)):
                t_phi = self.coupling_strength_add_ext_flux(t1, dx, [0, phi_ext]) \
                        * np.exp(1j * (np.pi / 3) * alpha * (i - 1/2))
                self.add_coupling(t_phi, u1, 'Cd', u2, 'C', dx, 'JW', True)
                self.add_coupling(np.conj(t_phi), u2, 'Cd', u1, 'C', -dx, 'JW', True)  # h.c.

            for u1, u2, dx in getattr(self.lat, "NN{}ur".format(i)):
                t_phi = self.coupling_strength_add_ext_flux(t1, dx, [0, phi_ext]) \
                        * np.exp(1j * (np.pi / 3) * alpha * (i + 1/2))
                self.add_coupling(t_phi, u1, 'Cd', u2, 'C', dx, 'JW', True)
                self.add_coupling(np.conj(t_phi), u2, 'Cd', u1, 'C', -dx, 'JW', True)  # h.c.

        # fifthNN ###
        for i in range(0, numb_sites):

            for u1, u2, dx in getattr(self.lat, "fifthNN{}u".format(i)):
                t_phi = self.coupling_strength_add_ext_flux(t2, dx, [0, phi_ext]) \
                        * np.exp(1j * (2 * np.pi) * alpha * i)
                self.add_coupling(t_phi, u1, 'Cd', u2, 'C', dx, 'JW', True)
                self.add_coupling(np.conj(t_phi), u2, 'Cd', u1, 'C', -dx, 'JW', True)  # h.c.

            for u1, u2, dx in getattr(self.lat, "fifthNN{}bl".format(i)):
                t_phi = self.coupling_strength_add_ext_flux(t2, dx, [0, phi_ext]) \
                        * np.exp(-1j * np.pi * alpha * (i - 3/2))
                self.add_coupling(t_phi, u1, 'Cd', u2, 'C', dx, 'JW', True)
                self.add_coupling(np.conj(t_phi), u2, 'Cd', u1, 'C', -dx, 'JW', True)  # h.c.

            for u1, u2, dx in getattr(self.lat, "fifthNN{}br".format(i)):
                t_phi = self.coupling_strength_add_ext_flux(t2, dx, [0, phi_ext]) \
                        * np.exp(-1j * np.pi * alpha * (i + 3/2))
                self.add_coupling(t_phi, u1, 'Cd', u2, 'C', dx, 'JW', True)
                self.add_coupling(np.conj(t_phi), u2, 'Cd', u1, 'C', -dx, 'JW', True)  # h.c.
