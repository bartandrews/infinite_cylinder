"""Spinless fermion twist model."""

import numpy as np
import sys

from tenpy.models.model import CouplingMPOModel, NearestNeighborModel
from tenpy.tools.params import get_parameter
from tenpy.networks.site import FermionSite
from lattices.MagneticTwist import MagneticTwist


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
        Ly = get_parameter(model_params, 'Ly', 3, self.name)

        fs = self.init_sites(model_params)

        lat = MagneticTwist(Lx, Ly, fs)

        print(lat.N_sites)

        return lat

    def init_terms(self, model_params):

        t = get_parameter(model_params, 't', -1., self.name, True)
        V = get_parameter(model_params, 'V', 0, self.name, True)
        mu = get_parameter(model_params, 'mu', 0., self.name, True)
        phi_ext = 2*np.pi*get_parameter(model_params, 'phi_ext', 0., self.name)

        t1 = 0.331 * t
        t2 = -0.01

        alpha = 2/7

        # NN ###########################################################################################################

        # down
        for u1, u2, dx in self.lat.NN0d:
            t_phi = self.coupling_strength_add_ext_flux(t1, dx, [0, phi_ext]) * np.exp(-1j*(2*np.pi/3)*alpha*0)
            self.add_coupling(t_phi, u1, 'Cd', u2, 'C', dx, 'JW', True)
            self.add_coupling(np.conj(t_phi), u2, 'Cd', u1, 'C', -dx, 'JW', True)  # h.c.
        for u1, u2, dx in self.lat.NN2d:
            t_phi = self.coupling_strength_add_ext_flux(t1, dx, [0, phi_ext]) * np.exp(-1j*(2*np.pi/3)*alpha*2)
            self.add_coupling(t_phi, u1, 'Cd', u2, 'C', dx, 'JW', True)
            self.add_coupling(np.conj(t_phi), u2, 'Cd', u1, 'C', -dx, 'JW', True)  # h.c.
        for u1, u2, dx in self.lat.NN4d:
            t_phi = self.coupling_strength_add_ext_flux(t1, dx, [0, phi_ext]) * np.exp(-1j*(2*np.pi/3)*alpha*4)
            self.add_coupling(t_phi, u1, 'Cd', u2, 'C', dx, 'JW', True)
            self.add_coupling(np.conj(t_phi), u2, 'Cd', u1, 'C', -dx, 'JW', True)  # h.c.
        for u1, u2, dx in self.lat.NN6d:
            t_phi = self.coupling_strength_add_ext_flux(t1, dx, [0, phi_ext]) * np.exp(-1j*(2*np.pi/3)*alpha*6)
            self.add_coupling(t_phi, u1, 'Cd', u2, 'C', dx, 'JW', True)
            self.add_coupling(np.conj(t_phi), u2, 'Cd', u1, 'C', -dx, 'JW', True)  # h.c.
        for u1, u2, dx in self.lat.NN8d:
            t_phi = self.coupling_strength_add_ext_flux(t1, dx, [0, phi_ext]) * np.exp(-1j*(2*np.pi/3)*alpha*8)
            self.add_coupling(t_phi, u1, 'Cd', u2, 'C', dx, 'JW', True)
            self.add_coupling(np.conj(t_phi), u2, 'Cd', u1, 'C', -dx, 'JW', True)  # h.c.
        for u1, u2, dx in self.lat.NN10d:
            t_phi = self.coupling_strength_add_ext_flux(t1, dx, [0, phi_ext]) * np.exp(-1j*(2*np.pi/3)*alpha*10)
            self.add_coupling(t_phi, u1, 'Cd', u2, 'C', dx, 'JW', True)
            self.add_coupling(np.conj(t_phi), u2, 'Cd', u1, 'C', -dx, 'JW', True)  # h.c.
        for u1, u2, dx in self.lat.NN12d:
            t_phi = self.coupling_strength_add_ext_flux(t1, dx, [0, phi_ext]) * np.exp(-1j*(2*np.pi/3)*alpha*12)
            self.add_coupling(t_phi, u1, 'Cd', u2, 'C', dx, 'JW', True)
            self.add_coupling(np.conj(t_phi), u2, 'Cd', u1, 'C', -dx, 'JW', True)  # h.c.

        # upper left
        for u1, u2, dx in self.lat.NN0ul:
            t_phi = self.coupling_strength_add_ext_flux(t1, dx, [0, phi_ext]) * np.exp(1j*(np.pi/3)*alpha*(0-1/2))
            self.add_coupling(t_phi, u1, 'Cd', u2, 'C', dx, 'JW', True)
            self.add_coupling(np.conj(t_phi), u2, 'Cd', u1, 'C', -dx, 'JW', True)  # h.c.
        for u1, u2, dx in self.lat.NN2ul:
            t_phi = self.coupling_strength_add_ext_flux(t1, dx, [0, phi_ext]) * np.exp(1j*(np.pi/3)*alpha*(2-1/2))
            self.add_coupling(t_phi, u1, 'Cd', u2, 'C', dx, 'JW', True)
            self.add_coupling(np.conj(t_phi), u2, 'Cd', u1, 'C', -dx, 'JW', True)  # h.c.
        for u1, u2, dx in self.lat.NN4ul:
            t_phi = self.coupling_strength_add_ext_flux(t1, dx, [0, phi_ext]) * np.exp(1j*(np.pi/3)*alpha*(4-1/2))
            self.add_coupling(t_phi, u1, 'Cd', u2, 'C', dx, 'JW', True)
            self.add_coupling(np.conj(t_phi), u2, 'Cd', u1, 'C', -dx, 'JW', True)  # h.c.
        for u1, u2, dx in self.lat.NN6ul:
            t_phi = self.coupling_strength_add_ext_flux(t1, dx, [0, phi_ext]) * np.exp(1j*(np.pi/3)*alpha*(6-1/2))
            self.add_coupling(t_phi, u1, 'Cd', u2, 'C', dx, 'JW', True)
            self.add_coupling(np.conj(t_phi), u2, 'Cd', u1, 'C', -dx, 'JW', True)  # h.c.
        for u1, u2, dx in self.lat.NN8ul:
            t_phi = self.coupling_strength_add_ext_flux(t1, dx, [0, phi_ext]) * np.exp(1j*(np.pi/3)*alpha*(8-1/2))
            self.add_coupling(t_phi, u1, 'Cd', u2, 'C', dx, 'JW', True)
            self.add_coupling(np.conj(t_phi), u2, 'Cd', u1, 'C', -dx, 'JW', True)  # h.c.
        for u1, u2, dx in self.lat.NN10ul:
            t_phi = self.coupling_strength_add_ext_flux(t1, dx, [0, phi_ext]) * np.exp(1j*(np.pi/3)*alpha*(10-1/2))
            self.add_coupling(t_phi, u1, 'Cd', u2, 'C', dx, 'JW', True)
            self.add_coupling(np.conj(t_phi), u2, 'Cd', u1, 'C', -dx, 'JW', True)  # h.c.
        for u1, u2, dx in self.lat.NN12ul:
            t_phi = self.coupling_strength_add_ext_flux(t1, dx, [0, phi_ext]) * np.exp(1j*(np.pi/3)*alpha*(12-1/2))
            self.add_coupling(t_phi, u1, 'Cd', u2, 'C', dx, 'JW', True)
            self.add_coupling(np.conj(t_phi), u2, 'Cd', u1, 'C', -dx, 'JW', True)  # h.c.

        # upper right
        for u1, u2, dx in self.lat.NN0ur:
            t_phi = self.coupling_strength_add_ext_flux(t1, dx, [0, phi_ext]) * np.exp(1j*(np.pi/3)*alpha*(0+1/2))
            self.add_coupling(t_phi, u1, 'Cd', u2, 'C', dx, 'JW', True)
            self.add_coupling(np.conj(t_phi), u2, 'Cd', u1, 'C', -dx, 'JW', True)  # h.c.
        for u1, u2, dx in self.lat.NN2ur:
            t_phi = self.coupling_strength_add_ext_flux(t1, dx, [0, phi_ext]) * np.exp(1j*(np.pi/3)*alpha*(2+1/2))
            self.add_coupling(t_phi, u1, 'Cd', u2, 'C', dx, 'JW', True)
            self.add_coupling(np.conj(t_phi), u2, 'Cd', u1, 'C', -dx, 'JW', True)  # h.c.
        for u1, u2, dx in self.lat.NN4ur:
            t_phi = self.coupling_strength_add_ext_flux(t1, dx, [0, phi_ext]) * np.exp(1j*(np.pi/3)*alpha*(4+1/2))
            self.add_coupling(t_phi, u1, 'Cd', u2, 'C', dx, 'JW', True)
            self.add_coupling(np.conj(t_phi), u2, 'Cd', u1, 'C', -dx, 'JW', True)  # h.c.
        for u1, u2, dx in self.lat.NN6ur:
            t_phi = self.coupling_strength_add_ext_flux(t1, dx, [0, phi_ext]) * np.exp(1j*(np.pi/3)*alpha*(6+1/2))
            self.add_coupling(t_phi, u1, 'Cd', u2, 'C', dx, 'JW', True)
            self.add_coupling(np.conj(t_phi), u2, 'Cd', u1, 'C', -dx, 'JW', True)  # h.c.
        for u1, u2, dx in self.lat.NN8ur:
            t_phi = self.coupling_strength_add_ext_flux(t1, dx, [0, phi_ext]) * np.exp(1j*(np.pi/3)*alpha*(8+1/2))
            self.add_coupling(t_phi, u1, 'Cd', u2, 'C', dx, 'JW', True)
            self.add_coupling(np.conj(t_phi), u2, 'Cd', u1, 'C', -dx, 'JW', True)  # h.c.
        for u1, u2, dx in self.lat.NN10ur:
            t_phi = self.coupling_strength_add_ext_flux(t1, dx, [0, phi_ext]) * np.exp(1j*(np.pi/3)*alpha*(10+1/2))
            self.add_coupling(t_phi, u1, 'Cd', u2, 'C', dx, 'JW', True)
            self.add_coupling(np.conj(t_phi), u2, 'Cd', u1, 'C', -dx, 'JW', True)  # h.c.
        for u1, u2, dx in self.lat.NN12ur:
            t_phi = self.coupling_strength_add_ext_flux(t1, dx, [0, phi_ext]) * np.exp(1j*(np.pi/3)*alpha*(12+1/2))
            self.add_coupling(t_phi, u1, 'Cd', u2, 'C', dx, 'JW', True)
            self.add_coupling(np.conj(t_phi), u2, 'Cd', u1, 'C', -dx, 'JW', True)  # h.c.

        # fifthNN ######################################################################################################

        # up
        for u1, u2, dx in self.lat.fifthNN0u:
            t_phi = self.coupling_strength_add_ext_flux(t2, dx, [0, phi_ext]) * np.exp(1j * (2 * np.pi) * alpha * 0)
            self.add_coupling(t_phi, u1, 'Cd', u2, 'C', dx, 'JW', True)
            self.add_coupling(np.conj(t_phi), u2, 'Cd', u1, 'C', -dx, 'JW', True)  # h.c.
        for u1, u2, dx in self.lat.fifthNN1u:
            t_phi = self.coupling_strength_add_ext_flux(t2, dx, [0, phi_ext]) * np.exp(1j * (2 * np.pi) * alpha * 1)
            self.add_coupling(t_phi, u1, 'Cd', u2, 'C', dx, 'JW', True)
            self.add_coupling(np.conj(t_phi), u2, 'Cd', u1, 'C', -dx, 'JW', True)  # h.c.
        for u1, u2, dx in self.lat.fifthNN2u:
            t_phi = self.coupling_strength_add_ext_flux(t2, dx, [0, phi_ext]) * np.exp(1j * (2 * np.pi) * alpha * 2)
            self.add_coupling(t_phi, u1, 'Cd', u2, 'C', dx, 'JW', True)
            self.add_coupling(np.conj(t_phi), u2, 'Cd', u1, 'C', -dx, 'JW', True)  # h.c.
        for u1, u2, dx in self.lat.fifthNN3u:
            t_phi = self.coupling_strength_add_ext_flux(t2, dx, [0, phi_ext]) * np.exp(1j * (2 * np.pi) * alpha * 3)
            self.add_coupling(t_phi, u1, 'Cd', u2, 'C', dx, 'JW', True)
            self.add_coupling(np.conj(t_phi), u2, 'Cd', u1, 'C', -dx, 'JW', True)  # h.c.
        for u1, u2, dx in self.lat.fifthNN4u:
            t_phi = self.coupling_strength_add_ext_flux(t2, dx, [0, phi_ext]) * np.exp(1j * (2 * np.pi) * alpha * 4)
            self.add_coupling(t_phi, u1, 'Cd', u2, 'C', dx, 'JW', True)
            self.add_coupling(np.conj(t_phi), u2, 'Cd', u1, 'C', -dx, 'JW', True)  # h.c.
        for u1, u2, dx in self.lat.fifthNN5u:
            t_phi = self.coupling_strength_add_ext_flux(t2, dx, [0, phi_ext]) * np.exp(1j * (2 * np.pi) * alpha * 5)
            self.add_coupling(t_phi, u1, 'Cd', u2, 'C', dx, 'JW', True)
            self.add_coupling(np.conj(t_phi), u2, 'Cd', u1, 'C', -dx, 'JW', True)  # h.c.
        for u1, u2, dx in self.lat.fifthNN6u:
            t_phi = self.coupling_strength_add_ext_flux(t2, dx, [0, phi_ext]) * np.exp(1j * (2 * np.pi) * alpha * 6)
            self.add_coupling(t_phi, u1, 'Cd', u2, 'C', dx, 'JW', True)
            self.add_coupling(np.conj(t_phi), u2, 'Cd', u1, 'C', -dx, 'JW', True)  # h.c.
        for u1, u2, dx in self.lat.fifthNN7u:
            t_phi = self.coupling_strength_add_ext_flux(t2, dx, [0, phi_ext]) * np.exp(1j * (2 * np.pi) * alpha * 7)
            self.add_coupling(t_phi, u1, 'Cd', u2, 'C', dx, 'JW', True)
            self.add_coupling(np.conj(t_phi), u2, 'Cd', u1, 'C', -dx, 'JW', True)  # h.c.
        for u1, u2, dx in self.lat.fifthNN8u:
            t_phi = self.coupling_strength_add_ext_flux(t2, dx, [0, phi_ext]) * np.exp(1j * (2 * np.pi) * alpha * 8)
            self.add_coupling(t_phi, u1, 'Cd', u2, 'C', dx, 'JW', True)
            self.add_coupling(np.conj(t_phi), u2, 'Cd', u1, 'C', -dx, 'JW', True)  # h.c.
        for u1, u2, dx in self.lat.fifthNN9u:
            t_phi = self.coupling_strength_add_ext_flux(t2, dx, [0, phi_ext]) * np.exp(1j * (2 * np.pi) * alpha * 9)
            self.add_coupling(t_phi, u1, 'Cd', u2, 'C', dx, 'JW', True)
            self.add_coupling(np.conj(t_phi), u2, 'Cd', u1, 'C', -dx, 'JW', True)  # h.c.
        for u1, u2, dx in self.lat.fifthNN10u:
            t_phi = self.coupling_strength_add_ext_flux(t2, dx, [0, phi_ext]) * np.exp(1j * (2 * np.pi) * alpha * 10)
            self.add_coupling(t_phi, u1, 'Cd', u2, 'C', dx, 'JW', True)
            self.add_coupling(np.conj(t_phi), u2, 'Cd', u1, 'C', -dx, 'JW', True)  # h.c.
        for u1, u2, dx in self.lat.fifthNN11u:
            t_phi = self.coupling_strength_add_ext_flux(t2, dx, [0, phi_ext]) * np.exp(1j * (2 * np.pi) * alpha * 11)
            self.add_coupling(t_phi, u1, 'Cd', u2, 'C', dx, 'JW', True)
            self.add_coupling(np.conj(t_phi), u2, 'Cd', u1, 'C', -dx, 'JW', True)  # h.c.
        for u1, u2, dx in self.lat.fifthNN12u:
            t_phi = self.coupling_strength_add_ext_flux(t2, dx, [0, phi_ext]) * np.exp(1j * (2 * np.pi) * alpha * 12)
            self.add_coupling(t_phi, u1, 'Cd', u2, 'C', dx, 'JW', True)
            self.add_coupling(np.conj(t_phi), u2, 'Cd', u1, 'C', -dx, 'JW', True)  # h.c.
        for u1, u2, dx in self.lat.fifthNN13u:
            t_phi = self.coupling_strength_add_ext_flux(t2, dx, [0, phi_ext]) * np.exp(1j * (2 * np.pi) * alpha * 13)
            self.add_coupling(t_phi, u1, 'Cd', u2, 'C', dx, 'JW', True)
            self.add_coupling(np.conj(t_phi), u2, 'Cd', u1, 'C', -dx, 'JW', True)  # h.c.

        # bottom left
        for u1, u2, dx in self.lat.fifthNN0bl:
            t_phi = self.coupling_strength_add_ext_flux(t2, dx, [0, phi_ext]) * np.exp(-1j * np.pi * alpha * (0-3/2))
            self.add_coupling(t_phi, u1, 'Cd', u2, 'C', dx, 'JW', True)
            self.add_coupling(np.conj(t_phi), u2, 'Cd', u1, 'C', -dx, 'JW', True)  # h.c.
        for u1, u2, dx in self.lat.fifthNN1bl:
            t_phi = self.coupling_strength_add_ext_flux(t2, dx, [0, phi_ext]) * np.exp(-1j * np.pi * alpha * (1-3/2))
            self.add_coupling(t_phi, u1, 'Cd', u2, 'C', dx, 'JW', True)
            self.add_coupling(np.conj(t_phi), u2, 'Cd', u1, 'C', -dx, 'JW', True)  # h.c.
        for u1, u2, dx in self.lat.fifthNN2bl:
            t_phi = self.coupling_strength_add_ext_flux(t2, dx, [0, phi_ext]) * np.exp(-1j * np.pi * alpha * (2-3/2))
            self.add_coupling(t_phi, u1, 'Cd', u2, 'C', dx, 'JW', True)
            self.add_coupling(np.conj(t_phi), u2, 'Cd', u1, 'C', -dx, 'JW', True)  # h.c.
        for u1, u2, dx in self.lat.fifthNN3bl:
            t_phi = self.coupling_strength_add_ext_flux(t2, dx, [0, phi_ext]) * np.exp(-1j * np.pi * alpha * (3-3/2))
            self.add_coupling(t_phi, u1, 'Cd', u2, 'C', dx, 'JW', True)
            self.add_coupling(np.conj(t_phi), u2, 'Cd', u1, 'C', -dx, 'JW', True)  # h.c.
        for u1, u2, dx in self.lat.fifthNN4bl:
            t_phi = self.coupling_strength_add_ext_flux(t2, dx, [0, phi_ext]) * np.exp(-1j * np.pi * alpha * (4-3/2))
            self.add_coupling(t_phi, u1, 'Cd', u2, 'C', dx, 'JW', True)
            self.add_coupling(np.conj(t_phi), u2, 'Cd', u1, 'C', -dx, 'JW', True)  # h.c.
        for u1, u2, dx in self.lat.fifthNN5bl:
            t_phi = self.coupling_strength_add_ext_flux(t2, dx, [0, phi_ext]) * np.exp(-1j * np.pi * alpha * (5-3/2))
            self.add_coupling(t_phi, u1, 'Cd', u2, 'C', dx, 'JW', True)
            self.add_coupling(np.conj(t_phi), u2, 'Cd', u1, 'C', -dx, 'JW', True)  # h.c.
        for u1, u2, dx in self.lat.fifthNN6bl:
            t_phi = self.coupling_strength_add_ext_flux(t2, dx, [0, phi_ext]) * np.exp(-1j * np.pi * alpha * (6-3/2))
            self.add_coupling(t_phi, u1, 'Cd', u2, 'C', dx, 'JW', True)
            self.add_coupling(np.conj(t_phi), u2, 'Cd', u1, 'C', -dx, 'JW', True)  # h.c.
        for u1, u2, dx in self.lat.fifthNN7bl:
            t_phi = self.coupling_strength_add_ext_flux(t2, dx, [0, phi_ext]) * np.exp(-1j * np.pi * alpha * (7-3/2))
            self.add_coupling(t_phi, u1, 'Cd', u2, 'C', dx, 'JW', True)
            self.add_coupling(np.conj(t_phi), u2, 'Cd', u1, 'C', -dx, 'JW', True)  # h.c.
        for u1, u2, dx in self.lat.fifthNN8bl:
            t_phi = self.coupling_strength_add_ext_flux(t2, dx, [0, phi_ext]) * np.exp(-1j * np.pi * alpha * (8-3/2))
            self.add_coupling(t_phi, u1, 'Cd', u2, 'C', dx, 'JW', True)
            self.add_coupling(np.conj(t_phi), u2, 'Cd', u1, 'C', -dx, 'JW', True)  # h.c.
        for u1, u2, dx in self.lat.fifthNN9bl:
            t_phi = self.coupling_strength_add_ext_flux(t2, dx, [0, phi_ext]) * np.exp(-1j * np.pi * alpha * (9-3/2))
            self.add_coupling(t_phi, u1, 'Cd', u2, 'C', dx, 'JW', True)
            self.add_coupling(np.conj(t_phi), u2, 'Cd', u1, 'C', -dx, 'JW', True)  # h.c.
        for u1, u2, dx in self.lat.fifthNN10bl:
            t_phi = self.coupling_strength_add_ext_flux(t2, dx, [0, phi_ext]) * np.exp(-1j * np.pi * alpha * (10-3/2))
            self.add_coupling(t_phi, u1, 'Cd', u2, 'C', dx, 'JW', True)
            self.add_coupling(np.conj(t_phi), u2, 'Cd', u1, 'C', -dx, 'JW', True)  # h.c.
        for u1, u2, dx in self.lat.fifthNN11bl:
            t_phi = self.coupling_strength_add_ext_flux(t2, dx, [0, phi_ext]) * np.exp(-1j * np.pi * alpha * (11-3/2))
            self.add_coupling(t_phi, u1, 'Cd', u2, 'C', dx, 'JW', True)
            self.add_coupling(np.conj(t_phi), u2, 'Cd', u1, 'C', -dx, 'JW', True)  # h.c.
        for u1, u2, dx in self.lat.fifthNN12bl:
            t_phi = self.coupling_strength_add_ext_flux(t2, dx, [0, phi_ext]) * np.exp(-1j * np.pi * alpha * (12-3/2))
            self.add_coupling(t_phi, u1, 'Cd', u2, 'C', dx, 'JW', True)
            self.add_coupling(np.conj(t_phi), u2, 'Cd', u1, 'C', -dx, 'JW', True)  # h.c.
        for u1, u2, dx in self.lat.fifthNN13bl:
            t_phi = self.coupling_strength_add_ext_flux(t2, dx, [0, phi_ext]) * np.exp(-1j * np.pi * alpha * (13-3/2))
            self.add_coupling(t_phi, u1, 'Cd', u2, 'C', dx, 'JW', True)
            self.add_coupling(np.conj(t_phi), u2, 'Cd', u1, 'C', -dx, 'JW', True)  # h.c.

        # bottom right
        for u1, u2, dx in self.lat.fifthNN0br:
            t_phi = self.coupling_strength_add_ext_flux(t2, dx, [0, phi_ext]) * np.exp(-1j * np.pi * alpha * (0 + 3 / 2))
            self.add_coupling(t_phi, u1, 'Cd', u2, 'C', dx, 'JW', True)
            self.add_coupling(np.conj(t_phi), u2, 'Cd', u1, 'C', -dx, 'JW', True)  # h.c.
        for u1, u2, dx in self.lat.fifthNN1br:
            t_phi = self.coupling_strength_add_ext_flux(t2, dx, [0, phi_ext]) * np.exp(-1j * np.pi * alpha * (1 + 3 / 2))
            self.add_coupling(t_phi, u1, 'Cd', u2, 'C', dx, 'JW', True)
            self.add_coupling(np.conj(t_phi), u2, 'Cd', u1, 'C', -dx, 'JW', True)  # h.c.
        for u1, u2, dx in self.lat.fifthNN2br:
            t_phi = self.coupling_strength_add_ext_flux(t2, dx, [0, phi_ext]) * np.exp(-1j * np.pi * alpha * (2 + 3 / 2))
            self.add_coupling(t_phi, u1, 'Cd', u2, 'C', dx, 'JW', True)
            self.add_coupling(np.conj(t_phi), u2, 'Cd', u1, 'C', -dx, 'JW', True)  # h.c.
        for u1, u2, dx in self.lat.fifthNN3br:
            t_phi = self.coupling_strength_add_ext_flux(t2, dx, [0, phi_ext]) * np.exp(-1j * np.pi * alpha * (3 + 3 / 2))
            self.add_coupling(t_phi, u1, 'Cd', u2, 'C', dx, 'JW', True)
            self.add_coupling(np.conj(t_phi), u2, 'Cd', u1, 'C', -dx, 'JW', True)  # h.c.
        for u1, u2, dx in self.lat.fifthNN4br:
            t_phi = self.coupling_strength_add_ext_flux(t2, dx, [0, phi_ext]) * np.exp(-1j * np.pi * alpha * (4 - 3 / 2))
            self.add_coupling(t_phi, u1, 'Cd', u2, 'C', dx, 'JW', True)
            self.add_coupling(np.conj(t_phi), u2, 'Cd', u1, 'C', -dx, 'JW', True)  # h.c.
        for u1, u2, dx in self.lat.fifthNN5br:
            t_phi = self.coupling_strength_add_ext_flux(t2, dx, [0, phi_ext]) * np.exp(-1j * np.pi * alpha * (5 + 3 / 2))
            self.add_coupling(t_phi, u1, 'Cd', u2, 'C', dx, 'JW', True)
            self.add_coupling(np.conj(t_phi), u2, 'Cd', u1, 'C', -dx, 'JW', True)  # h.c.
        for u1, u2, dx in self.lat.fifthNN6br:
            t_phi = self.coupling_strength_add_ext_flux(t2, dx, [0, phi_ext]) * np.exp(-1j * np.pi * alpha * (6 + 3 / 2))
            self.add_coupling(t_phi, u1, 'Cd', u2, 'C', dx, 'JW', True)
            self.add_coupling(np.conj(t_phi), u2, 'Cd', u1, 'C', -dx, 'JW', True)  # h.c.
        for u1, u2, dx in self.lat.fifthNN7br:
            t_phi = self.coupling_strength_add_ext_flux(t2, dx, [0, phi_ext]) * np.exp(-1j * np.pi * alpha * (7 + 3 / 2))
            self.add_coupling(t_phi, u1, 'Cd', u2, 'C', dx, 'JW', True)
            self.add_coupling(np.conj(t_phi), u2, 'Cd', u1, 'C', -dx, 'JW', True)  # h.c.
        for u1, u2, dx in self.lat.fifthNN8br:
            t_phi = self.coupling_strength_add_ext_flux(t2, dx, [0, phi_ext]) * np.exp(-1j * np.pi * alpha * (8 + 3 / 2))
            self.add_coupling(t_phi, u1, 'Cd', u2, 'C', dx, 'JW', True)
            self.add_coupling(np.conj(t_phi), u2, 'Cd', u1, 'C', -dx, 'JW', True)  # h.c.
        for u1, u2, dx in self.lat.fifthNN9br:
            t_phi = self.coupling_strength_add_ext_flux(t2, dx, [0, phi_ext]) * np.exp(-1j * np.pi * alpha * (9 + 3 / 2))
            self.add_coupling(t_phi, u1, 'Cd', u2, 'C', dx, 'JW', True)
            self.add_coupling(np.conj(t_phi), u2, 'Cd', u1, 'C', -dx, 'JW', True)  # h.c.
        for u1, u2, dx in self.lat.fifthNN10br:
            t_phi = self.coupling_strength_add_ext_flux(t2, dx, [0, phi_ext]) * np.exp(-1j * np.pi * alpha * (10 + 3 / 2))
            self.add_coupling(t_phi, u1, 'Cd', u2, 'C', dx, 'JW', True)
            self.add_coupling(np.conj(t_phi), u2, 'Cd', u1, 'C', -dx, 'JW', True)  # h.c.
        for u1, u2, dx in self.lat.fifthNN11br:
            t_phi = self.coupling_strength_add_ext_flux(t2, dx, [0, phi_ext]) * np.exp(-1j * np.pi * alpha * (11 + 3 / 2))
            self.add_coupling(t_phi, u1, 'Cd', u2, 'C', dx, 'JW', True)
            self.add_coupling(np.conj(t_phi), u2, 'Cd', u1, 'C', -dx, 'JW', True)  # h.c.
        for u1, u2, dx in self.lat.fifthNN12br:
            t_phi = self.coupling_strength_add_ext_flux(t2, dx, [0, phi_ext]) * np.exp(-1j * np.pi * alpha * (12 + 3 / 2))
            self.add_coupling(t_phi, u1, 'Cd', u2, 'C', dx, 'JW', True)
            self.add_coupling(np.conj(t_phi), u2, 'Cd', u1, 'C', -dx, 'JW', True)  # h.c.
        for u1, u2, dx in self.lat.fifthNN13br:
            t_phi = self.coupling_strength_add_ext_flux(t2, dx, [0, phi_ext]) * np.exp(-1j * np.pi * alpha * (13 + 3 / 2))
            self.add_coupling(t_phi, u1, 'Cd', u2, 'C', dx, 'JW', True)
            self.add_coupling(np.conj(t_phi), u2, 'Cd', u1, 'C', -dx, 'JW', True)  # h.c.


class FermionicTwistChain(FermionicTwistModel, NearestNeighborModel):

    def __init__(self, model_params):
        model_params.setdefault('lattice', "Chain")
        CouplingMPOModel.__init__(self, model_params)
