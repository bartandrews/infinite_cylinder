"""Five-band model of twisted bilayer graphene.
Hamiltonian based on: "Faithful Tight-binding Models and Fragile Topology of Magic-angle Bilayer Graphene"."""

import numpy as np
import sys

from tenpy.models.model import CouplingMPOModel, NearestNeighborModel
from tenpy.tools.params import get_parameter
from tenpy.networks.site import FermionSite, GroupedSite
from lattices.FiveBandModel import FiveBandLattice
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

        phi_ext = 2 * np.pi * get_parameter(model_params, 'phi_ext', 0., self.name)

        a = 0.25 * t  # in units of 80 meV
        b = 0.2 * t
        c = 0.1 * t
        d = 0.67 * t  # five_band_model lattice assumes d is real

        field = np.exp(1j*(2*np.pi)*(1/3))  # magnetic field via Peierls substitution (Phi in units of Phi_0)

        for u in range(self.lat.N_cells):

            self.add_onsite(-0.043*t, 0, 'Npz')
            self.add_onsite(0.05*t, 1, 'N')
            self.add_onsite(0.05*t, 2, 'N')
            self.add_onsite(-0.043*t, 3, 'Npz')
            self.add_onsite(0.05*t, 4, 'N')
            self.add_onsite(0.05*t, 5, 'N')

        for u1, u2, dx in self.lat.a1_d:

            coupling_1 = self.coupling_strength_add_ext_flux(field*(-a*1j)*d, dx, [0, phi_ext])
            self.add_coupling(coupling_1, u1, 'Cpz', u2, 'Cd', dx, 'JW', True)
            self.add_coupling(np.conj(coupling_1), u2, 'C', u1, 'Cdpz', -dx, 'JW', True)  # h.c.

            coupling_2 = self.coupling_strength_add_ext_flux(field*c*d, dx, [0, phi_ext])
            self.add_coupling(coupling_2, u1, 'Cpp', u2, 'Cd', dx, 'JW', True)
            self.add_coupling(np.conj(coupling_2), u2, 'C', u1, 'Cdpp', -dx, 'JW', True)  # h.c.

            coupling_3 = self.coupling_strength_add_ext_flux(field*b*d, dx, [0, phi_ext])
            self.add_coupling(coupling_3, u1, 'Cpm', u2, 'Cd', dx, 'JW', True)
            self.add_coupling(np.conj(coupling_3), u2, 'C', u1, 'Cdpm', -dx, 'JW', True)  # h.c.

        for u1, u2, dx in self.lat.a2_d:

            coupling_4 = self.coupling_strength_add_ext_flux(field*(a*1j)*d, dx, [0, phi_ext])
            self.add_coupling(coupling_4, u1, 'Cpz', u2, 'Cd', dx, 'JW', True)
            self.add_coupling(np.conj(coupling_4), u2, 'C', u1, 'Cdpz', -dx, 'JW', True)  # h.c.

            coupling_5 = self.coupling_strength_add_ext_flux(field*b*d, dx, [0, phi_ext])
            self.add_coupling(coupling_5, u1, 'Cpp', u2, 'Cd', dx, 'JW', True)
            self.add_coupling(np.conj(coupling_5), u2, 'C', u1, 'Cdpp', -dx, 'JW', True)  # h.c.

            coupling_6 = self.coupling_strength_add_ext_flux(field*c*d, dx, [0, phi_ext])
            self.add_coupling(coupling_6, u1, 'Cpm', u2, 'Cd', dx, 'JW', True)
            self.add_coupling(np.conj(coupling_6), u2, 'C', u1, 'Cdpm', -dx, 'JW', True)  # h.c.

        for u1, u2, dx in self.lat.d_d:

            coupling_7 = self.coupling_strength_add_ext_flux(field*d*d, dx, [0, phi_ext])
            self.add_coupling(coupling_7, u1, 'C', u2, 'Cd', dx, 'JW', True)
            self.add_coupling(np.conj(coupling_7), u2, 'C', u1, 'Cd', -dx, 'JW', True)  # h.c.

        for u1, u2, dx in self.lat.a_a:

            # from pz orbital

            coupling_8 = self.coupling_strength_add_ext_flux(field*(-1j*a)*np.conj(1j*a), dx, [0, phi_ext])
            coupling_9 = self.coupling_strength_add_ext_flux(field*(-1j*a)*np.conj(b), dx, [0, phi_ext])
            coupling_10 = self.coupling_strength_add_ext_flux(field*(-1j*a)*np.conj(c), dx, [0, phi_ext])

            self.add_coupling(coupling_8, u1, 'Cpz', u2, 'Cdpz', dx, 'JW', True)
            self.add_coupling(np.conj(coupling_8), u2, 'Cpz', u1, 'Cdpz', -dx, 'JW', True)  # h.c.
            self.add_coupling(coupling_9, u1, 'Cpz', u2, 'Cdpp', dx, 'JW', True)
            self.add_coupling(np.conj(coupling_9), u2, 'Cpz', u1, 'Cdpp', -dx, 'JW', True)  # h.c.
            self.add_coupling(coupling_10, u1, 'Cpz', u2, 'Cdpm', dx, 'JW', True)
            self.add_coupling(np.conj(coupling_10), u2, 'Cpz', u1, 'Cdpm', -dx, 'JW', True)  # h.c.

            # from pp orbital

            coupling_11 = self.coupling_strength_add_ext_flux(field*c*np.conj(1j*a), dx, [0, phi_ext])
            coupling_12 = self.coupling_strength_add_ext_flux(field*c*np.conj(b), dx, [0, phi_ext])
            coupling_13 = self.coupling_strength_add_ext_flux(field*c*np.conj(c), dx, [0, phi_ext])

            self.add_coupling(coupling_11, u1, 'Cpp', u2, 'Cdpz', dx, 'JW', True)
            self.add_coupling(np.conj(coupling_11), u2, 'Cpp', u1, 'Cdpz', -dx, 'JW', True)  # h.c.
            self.add_coupling(coupling_12, u1, 'Cpp', u2, 'Cdpp', dx, 'JW', True)
            self.add_coupling(np.conj(coupling_12), u2, 'Cpp', u1, 'Cdpp', -dx, 'JW', True)  # h.c.
            self.add_coupling(coupling_13, u1, 'Cpp', u2, 'Cdpm', dx, 'JW', True)
            self.add_coupling(np.conj(coupling_13), u2, 'Cpp', u1, 'Cdpm', -dx, 'JW', True)  # h.c.

            # from pm orbital

            coupling_14 = self.coupling_strength_add_ext_flux(field*b*np.conj(1j*a), dx, [0, phi_ext])
            coupling_15 = self.coupling_strength_add_ext_flux(field*b*np.conj(b), dx, [0, phi_ext])
            coupling_16 = self.coupling_strength_add_ext_flux(field*b*np.conj(c), dx, [0, phi_ext])

            self.add_coupling(coupling_14, u1, 'Cpm', u2, 'Cdpz', dx, 'JW', True)
            self.add_coupling(np.conj(coupling_14), u2, 'Cpm', u1, 'Cdpz', -dx, 'JW', True)  # h.c.
            self.add_coupling(coupling_15, u1, 'Cpm', u2, 'Cdpp', dx, 'JW', True)
            self.add_coupling(np.conj(coupling_15), u2, 'Cpm', u1, 'Cdpp', -dx, 'JW', True)  # h.c.
            self.add_coupling(coupling_16, u1, 'Cpm', u2, 'Cdpm', dx, 'JW', True)
            self.add_coupling(np.conj(coupling_16), u2, 'Cpm', u1, 'Cdpm', -dx, 'JW', True)  # h.c.


class FermionicTBG4Chain(FermionicTBG4Model, NearestNeighborModel):

    def __init__(self, model_params):
        model_params.setdefault('lattice', "Chain")
        CouplingMPOModel.__init__(self, model_params)
