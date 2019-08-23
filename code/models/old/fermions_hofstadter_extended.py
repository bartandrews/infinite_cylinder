"""Spinless fermion Hofstadter model.
Hamiltonian based on: "Square Lattice with Magnetic Field", Aidelsburger PhD thesis."""

import numpy as np
import sys

from tenpy.models.model import CouplingMPOModel
from tenpy.tools.params import get_parameter
from tenpy.networks.site import FermionSite
from lattices.MagneticSquareExtended import MagneticSquareExtended


class FermionicHofstadterExtendedModel(CouplingMPOModel):

    def __init__(self, model_params):
        CouplingMPOModel.__init__(self, model_params)

    def init_sites(self, model_params):
        conserve = get_parameter(model_params, 'conserve', 'N', self.name)
        site = FermionSite(conserve=conserve)
        return site

    def init_lattice(self, model_params):
        choice = get_parameter(model_params, 'lattice', 'MagneticSquareExtended', self.name)
        if choice != 'MagneticSquareExtended':
            sys.exit("Error: Please choose the MagneticSquareExtended for hofstadter_extended.")
        Lx = get_parameter(model_params, 'Lx', 1, self.name)
        Ly = get_parameter(model_params, 'Ly', 2, self.name)
        fs = self.init_sites(model_params)
        lat = MagneticSquareExtended(Lx, Ly, fs)
        print(lat.N_sites)
        return lat

    def init_terms(self, model_params):
        t = get_parameter(model_params, 't', -1., self.name, True)
        V = get_parameter(model_params, 'V', 0, self.name, True)
        phi_ext = 2*np.pi*get_parameter(model_params, 'phi_ext', 0., self.name)

        p = 1
        q = 5

        alpha = p / q
        numb_sites = q

        for u1, u2, dx in self.lat.NN_vert:
            t_phi = self.coupling_strength_add_ext_flux(t, dx, [0, phi_ext])
            self.add_coupling(t_phi, u1, 'Cd', u2, 'C', dx, 'JW', True)
            self.add_coupling(np.conj(t_phi), u2, 'Cd', u1, 'C', -dx, 'JW', True)  # h.c.
            self.add_coupling(V, u1, 'N', u2, 'N', dx)

        for i in range(2*numb_sites):
            for u1, u2, dx in getattr(self.lat, "NN_h{}".format(i)):
                t_phi = self.coupling_strength_add_ext_flux(t, dx, [0, phi_ext]) \
                        * np.exp(-1j * 2 * np.pi * alpha * (i % numb_sites))
                self.add_coupling(t_phi, u1, 'Cd', u2, 'C', dx, 'JW', True)
                self.add_coupling(np.conj(t_phi), u2, 'Cd', u1, 'C', -dx, 'JW', True)  # h.c.
                self.add_coupling(V, u1, 'N', u2, 'N', dx)
