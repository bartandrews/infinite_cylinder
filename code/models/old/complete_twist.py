"""Complete twist model."""

import numpy as np
import sys

from tenpy.models.model import CouplingMPOModel
from tenpy.tools.params import get_parameter
from tenpy.networks.site import BosonSite, FermionSite, GroupedSite
from lattices.MagneticHoneycomb import MagneticHoneycomb

__all__ = ['BosonicCompleteTwistModel', 'FermionicCompleteTwistModel']


class BosonicCompleteTwistModel(CouplingMPOModel):

    def __init__(self, model_params):
        CouplingMPOModel.__init__(self, model_params)

    def init_sites(self, model_params):
        conserve = get_parameter(model_params, 'conserve', 'N', self.name)

        nmax = 4
        print("Nmax = ", nmax)
        bs = BosonSite(Nmax=nmax, conserve=conserve)

        gs = GroupedSite([bs, bs], labels=['x', 'y'], charges='same')
        gs.add_op('Ntot', gs.Nx + gs.Ny, False)
        # print(sorted(gs.opnames))
        # print(gs.state_labels)
        return gs

    def init_lattice(self, model_params):
        choice = get_parameter(model_params, 'lattice', 'MagneticTwist', self.name)
        if choice != 'MagneticTwist':
            sys.exit("Error: Please choose the MagneticTwist for complete twist.")
        Lx = get_parameter(model_params, 'Lx', 1, self.name)
        Ly = get_parameter(model_params, 'Ly', 1, self.name)
        bs = self.init_sites(model_params)
        lat = MagneticHoneycomb(Lx, Ly, bs)
        print(lat.N_sites)
        return lat

    def init_terms(self, model_params):
        t = get_parameter(model_params, 't', -1., self.name, True)
        V = get_parameter(model_params, 'V', 0, self.name, True)
        phi_ext = 2*np.pi*get_parameter(model_params, 'phi_ext', 0., self.name)

        t1 = 0.331 * t
        t2 = -0.01 * t
        t2dash = 0.097 * t

        p = 2
        q = 9

        alpha = p/q
        numb_sites = 2*q

        # t1 term ###
        for i in range(0, numb_sites, 2):
            for u1, u2, dx in getattr(self.lat, "NN{}d".format(i)):
                t_phi = self.coupling_strength_add_ext_flux(t1, dx, [0, phi_ext]) \
                        * np.exp(-1j * (2 * np.pi / 3) * alpha * i)
                for orbital in ['x', 'y']:
                    self.add_coupling(t_phi, u1, 'Bd'+orbital, u2, 'B'+orbital, dx)
                    self.add_coupling(np.conj(t_phi), u2, 'Bd'+orbital, u1, 'B'+orbital, -dx)  # h.c.
                self.add_coupling(V, u1, 'Ntot', u2, 'Ntot', dx, category='V Ntot_i Ntot_j')
            for u1, u2, dx in getattr(self.lat, "NN{}ul".format(i)):
                t_phi = self.coupling_strength_add_ext_flux(t1, dx, [0, phi_ext]) \
                        * np.exp(1j * (np.pi / 3) * alpha * (i - 1/2))
                for orbital in ['x', 'y']:
                    self.add_coupling(t_phi, u1, 'Bd' + orbital, u2, 'B' + orbital, dx)
                    self.add_coupling(np.conj(t_phi), u2, 'Bd' + orbital, u1, 'B' + orbital, -dx)  # h.c.
                self.add_coupling(V, u1, 'Ntot', u2, 'Ntot', dx, category='V Ntot_i Ntot_j')
            for u1, u2, dx in getattr(self.lat, "NN{}ur".format(i)):
                t_phi = self.coupling_strength_add_ext_flux(t1, dx, [0, phi_ext]) \
                        * np.exp(1j * (np.pi / 3) * alpha * (i + 1/2))
                for orbital in ['x', 'y']:
                    self.add_coupling(t_phi, u1, 'Bd' + orbital, u2, 'B' + orbital, dx)
                    self.add_coupling(np.conj(t_phi), u2, 'Bd' + orbital, u1, 'B' + orbital, -dx)  # h.c.
                self.add_coupling(V, u1, 'Ntot', u2, 'Ntot', dx, category='V Ntot_i Ntot_j')

        # t2 term ###
        for i in range(0, numb_sites):
            for u1, u2, dx in getattr(self.lat, "fifthNN{}u".format(i)):
                t_phi = self.coupling_strength_add_ext_flux(t2, dx, [0, phi_ext]) \
                        * np.exp(1j * (2 * np.pi) * alpha * i)
                for orbital in ['x', 'y']:
                    self.add_coupling(t_phi, u1, 'Bd' + orbital, u2, 'B' + orbital, dx)
                    self.add_coupling(np.conj(t_phi), u2, 'Bd' + orbital, u1, 'B' + orbital, -dx)  # h.c.
            for u1, u2, dx in getattr(self.lat, "fifthNN{}bl".format(i)):
                t_phi = self.coupling_strength_add_ext_flux(t2, dx, [0, phi_ext]) \
                        * np.exp(-1j * np.pi * alpha * (i - 3/2))
                for orbital in ['x', 'y']:
                    self.add_coupling(t_phi, u1, 'Bd' + orbital, u2, 'B' + orbital, dx)
                    self.add_coupling(np.conj(t_phi), u2, 'Bd' + orbital, u1, 'B' + orbital, -dx)  # h.c.
            for u1, u2, dx in getattr(self.lat, "fifthNN{}br".format(i)):
                t_phi = self.coupling_strength_add_ext_flux(t2, dx, [0, phi_ext]) \
                        * np.exp(-1j * np.pi * alpha * (i + 3/2))
                for orbital in ['x', 'y']:
                    self.add_coupling(t_phi, u1, 'Bd' + orbital, u2, 'B' + orbital, dx)
                    self.add_coupling(np.conj(t_phi), u2, 'Bd' + orbital, u1, 'B' + orbital, -dx)  # h.c.

        # t2dash term ###
        for i in range(0, numb_sites):
            # positive term
            for u1, u2, dx in getattr(self.lat, "fifthNN{}u".format(i)):
                t_phi = self.coupling_strength_add_ext_flux(t2dash, dx, [0, phi_ext]) \
                        * np.exp(1j * (2 * np.pi) * alpha * i)
                self.add_coupling(t_phi, u1, 'Bdx', u2, 'By', dx)
                self.add_coupling(np.conj(t_phi), u2, 'Bdy', u1, 'Bx', -dx)  # h.c.
            for u1, u2, dx in getattr(self.lat, "fifthNN{}bl".format(i)):
                t_phi = self.coupling_strength_add_ext_flux(t2dash, dx, [0, phi_ext]) \
                        * np.exp(-1j * np.pi * alpha * (i - 3 / 2))
                self.add_coupling(t_phi, u1, 'Bdx', u2, 'By', dx)
                self.add_coupling(np.conj(t_phi), u2, 'Bdy', u1, 'Bx', -dx)  # h.c.
            for u1, u2, dx in getattr(self.lat, "fifthNN{}br".format(i)):
                t_phi = self.coupling_strength_add_ext_flux(t2dash, dx, [0, phi_ext]) \
                        * np.exp(-1j * np.pi * alpha * (i + 3 / 2))
                self.add_coupling(t_phi, u1, 'Bdx', u2, 'By', dx)
                self.add_coupling(np.conj(t_phi), u2, 'Bdy', u1, 'Bx', -dx)  # h.c.
            # negative term
            for u1, u2, dx in getattr(self.lat, "fifthNN{}u".format(i)):
                t_phi = self.coupling_strength_add_ext_flux(-t2dash, dx, [0, phi_ext]) \
                        * np.exp(1j * (2 * np.pi) * alpha * i)
                self.add_coupling(t_phi, u1, 'Bdy', u2, 'Bx', dx)
                self.add_coupling(np.conj(t_phi), u2, 'Bdx', u1, 'By', -dx)  # h.c.
            for u1, u2, dx in getattr(self.lat, "fifthNN{}bl".format(i)):
                t_phi = self.coupling_strength_add_ext_flux(-t2dash, dx, [0, phi_ext]) \
                        * np.exp(-1j * np.pi * alpha * (i - 3 / 2))
                self.add_coupling(t_phi, u1, 'Bdy', u2, 'Bx', dx)
                self.add_coupling(np.conj(t_phi), u2, 'Bdx', u1, 'By', -dx)  # h.c.
            for u1, u2, dx in getattr(self.lat, "fifthNN{}br".format(i)):
                t_phi = self.coupling_strength_add_ext_flux(-t2dash, dx, [0, phi_ext]) \
                        * np.exp(-1j * np.pi * alpha * (i + 3 / 2))
                self.add_coupling(t_phi, u1, 'Bdy', u2, 'Bx', dx)
                self.add_coupling(np.conj(t_phi), u2, 'Bdx', u1, 'By', -dx)  # h.c.


class FermionicCompleteTwistModel(CouplingMPOModel):

    def __init__(self, model_params):
        CouplingMPOModel.__init__(self, model_params)

    def init_sites(self, model_params):
        conserve = get_parameter(model_params, 'conserve', 'N', self.name)
        fs = FermionSite(conserve=conserve)
        gs = GroupedSite([fs, fs], labels=['x', 'y'], charges='same')
        gs.add_op('Ntot', gs.Nx + gs.Ny, False)
        return gs

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
        V = get_parameter(model_params, 'V', 0, self.name, True)
        phi_ext = 2*np.pi*get_parameter(model_params, 'phi_ext', 0., self.name)

        t1 = 0.331 * t
        t2 = -0.01 * t
        t2dash = 0.097 * t

        p = 2
        q = 9

        alpha = p/q
        numb_sites = 2*q

        # t1 term ###
        for i in range(0, numb_sites, 2):
            for u1, u2, dx in getattr(self.lat, "NN{}d".format(i)):
                t_phi = self.coupling_strength_add_ext_flux(t1, dx, [0, phi_ext]) \
                        * np.exp(-1j * (2 * np.pi / 3) * alpha * i)
                for orbital in ['x', 'y']:
                    self.add_coupling(t_phi, u1, 'Cd'+orbital, u2, 'C'+orbital, dx, 'JW', True)
                    self.add_coupling(np.conj(t_phi), u2, 'Cd'+orbital, u1, 'C'+orbital, -dx, 'JW', True)  # h.c.
                self.add_coupling(V, u1, 'Ntot', u2, 'Ntot', dx, category='V Ntot_i Ntot_j')
            for u1, u2, dx in getattr(self.lat, "NN{}ul".format(i)):
                t_phi = self.coupling_strength_add_ext_flux(t1, dx, [0, phi_ext]) \
                        * np.exp(1j * (np.pi / 3) * alpha * (i - 1/2))
                for orbital in ['x', 'y']:
                    self.add_coupling(t_phi, u1, 'Cd' + orbital, u2, 'C' + orbital, dx, 'JW', True)
                    self.add_coupling(np.conj(t_phi), u2, 'Cd' + orbital, u1, 'C' + orbital, -dx, 'JW', True)  # h.c.
                self.add_coupling(V, u1, 'Ntot', u2, 'Ntot', dx, category='V Ntot_i Ntot_j')
            for u1, u2, dx in getattr(self.lat, "NN{}ur".format(i)):
                t_phi = self.coupling_strength_add_ext_flux(t1, dx, [0, phi_ext]) \
                        * np.exp(1j * (np.pi / 3) * alpha * (i + 1/2))
                for orbital in ['x', 'y']:
                    self.add_coupling(t_phi, u1, 'Cd' + orbital, u2, 'C' + orbital, dx, 'JW', True)
                    self.add_coupling(np.conj(t_phi), u2, 'Cd' + orbital, u1, 'C' + orbital, -dx, 'JW', True)  # h.c.
                self.add_coupling(V, u1, 'Ntot', u2, 'Ntot', dx, category='V Ntot_i Ntot_j')

        # t2 term ###
        for i in range(0, numb_sites):
            for u1, u2, dx in getattr(self.lat, "fifthNN{}u".format(i)):
                t_phi = self.coupling_strength_add_ext_flux(t2, dx, [0, phi_ext]) \
                        * np.exp(1j * (2 * np.pi) * alpha * i)
                for orbital in ['x', 'y']:
                    self.add_coupling(t_phi, u1, 'Cd' + orbital, u2, 'C' + orbital, dx, 'JW', True)
                    self.add_coupling(np.conj(t_phi), u2, 'Cd' + orbital, u1, 'C' + orbital, -dx, 'JW', True)  # h.c.
            for u1, u2, dx in getattr(self.lat, "fifthNN{}bl".format(i)):
                t_phi = self.coupling_strength_add_ext_flux(t2, dx, [0, phi_ext]) \
                        * np.exp(-1j * np.pi * alpha * (i - 3/2))
                for orbital in ['x', 'y']:
                    self.add_coupling(t_phi, u1, 'Cd' + orbital, u2, 'C' + orbital, dx, 'JW', True)
                    self.add_coupling(np.conj(t_phi), u2, 'Cd' + orbital, u1, 'C' + orbital, -dx, 'JW', True)  # h.c.
            for u1, u2, dx in getattr(self.lat, "fifthNN{}br".format(i)):
                t_phi = self.coupling_strength_add_ext_flux(t2, dx, [0, phi_ext]) \
                        * np.exp(-1j * np.pi * alpha * (i + 3/2))
                for orbital in ['x', 'y']:
                    self.add_coupling(t_phi, u1, 'Cd' + orbital, u2, 'C' + orbital, dx, 'JW', True)
                    self.add_coupling(np.conj(t_phi), u2, 'Cd' + orbital, u1, 'C' + orbital, -dx, 'JW', True)  # h.c.

        # t2dash term ###
        for i in range(0, numb_sites):
            # positive term
            for u1, u2, dx in getattr(self.lat, "fifthNN{}u".format(i)):
                t_phi = self.coupling_strength_add_ext_flux(t2dash, dx, [0, phi_ext]) \
                        * np.exp(1j * (2 * np.pi) * alpha * i)
                self.add_coupling(t_phi, u1, 'Cdx', u2, 'Cy', dx, 'JW', True)
                self.add_coupling(np.conj(t_phi), u2, 'Cdy', u1, 'Cx', -dx, 'JW', True)  # h.c.
            for u1, u2, dx in getattr(self.lat, "fifthNN{}bl".format(i)):
                t_phi = self.coupling_strength_add_ext_flux(t2dash, dx, [0, phi_ext]) \
                        * np.exp(-1j * np.pi * alpha * (i - 3 / 2))
                self.add_coupling(t_phi, u1, 'Cdx', u2, 'Cy', dx, 'JW', True)
                self.add_coupling(np.conj(t_phi), u2, 'Cdy', u1, 'Cx', -dx, 'JW', True)  # h.c.
            for u1, u2, dx in getattr(self.lat, "fifthNN{}br".format(i)):
                t_phi = self.coupling_strength_add_ext_flux(t2dash, dx, [0, phi_ext]) \
                        * np.exp(-1j * np.pi * alpha * (i + 3 / 2))
                self.add_coupling(t_phi, u1, 'Cdx', u2, 'Cy', dx, 'JW', True)
                self.add_coupling(np.conj(t_phi), u2, 'Cdy', u1, 'Cx', -dx, 'JW', True)  # h.c.
            # negative term
            for u1, u2, dx in getattr(self.lat, "fifthNN{}u".format(i)):
                t_phi = self.coupling_strength_add_ext_flux(-t2dash, dx, [0, phi_ext]) \
                        * np.exp(1j * (2 * np.pi) * alpha * i)
                self.add_coupling(t_phi, u1, 'Cdy', u2, 'Cx', dx, 'JW', True)
                self.add_coupling(np.conj(t_phi), u2, 'Cdx', u1, 'Cy', -dx, 'JW', True)  # h.c.
            for u1, u2, dx in getattr(self.lat, "fifthNN{}bl".format(i)):
                t_phi = self.coupling_strength_add_ext_flux(-t2dash, dx, [0, phi_ext]) \
                        * np.exp(-1j * np.pi * alpha * (i - 3 / 2))
                self.add_coupling(t_phi, u1, 'Cdy', u2, 'Cx', dx, 'JW', True)
                self.add_coupling(np.conj(t_phi), u2, 'Cdx', u1, 'Cy', -dx, 'JW', True)  # h.c.
            for u1, u2, dx in getattr(self.lat, "fifthNN{}br".format(i)):
                t_phi = self.coupling_strength_add_ext_flux(-t2dash, dx, [0, phi_ext]) \
                        * np.exp(-1j * np.pi * alpha * (i + 3 / 2))
                self.add_coupling(t_phi, u1, 'Cdy', u2, 'Cx', dx, 'JW', True)
                self.add_coupling(np.conj(t_phi), u2, 'Cdx', u1, 'Cy', -dx, 'JW', True)  # h.c.


 # phi=2/7

# # bosons with filled LLL
# product_state = ['1_x 1_y', '1_x 1_y', '0_x 0_y', '0_x 0_y', '0_x 0_y',
#                  '0_x 0_y', '0_x 0_y', '0_x 0_y', '0_x 0_y', '0_x 0_y',
#                  '0_x 0_y', '0_x 0_y', '0_x 0_y', '0_x 0_y']
# # bosons with 1/2-filled LLL
# product_state = ['1_x 1_y', '0_x 0_y', '0_x 0_y', '0_x 0_y', '0_x 0_y',
#                  '0_x 0_y', '0_x 0_y', '0_x 0_y', '0_x 0_y', '0_x 0_y',
#                  '0_x 0_y', '0_x 0_y', '0_x 0_y', '0_x 0_y']
# # fermions with filled LLL
# product_state = ['full_x full_y', 'full_x full_y', 'empty_x empty_y', 'empty_x empty_y', 'empty_x empty_y',
#                  'empty_x empty_y', 'empty_x empty_y', 'empty_x empty_y', 'empty_x empty_y', 'empty_x empty_y',
#                  'empty_x empty_y', 'empty_x empty_y', 'empty_x empty_y', 'empty_x empty_y']
# # fermions with 1/2-filled LLL
# product_state = ['full_x full_y', 'empty_x empty_y', 'empty_x empty_y', 'empty_x empty_y', 'empty_x empty_y',
#                  'empty_x empty_y', 'empty_x empty_y', 'empty_x empty_y', 'empty_x empty_y', 'empty_x empty_y',
#                  'empty_x empty_y', 'empty_x empty_y', 'empty_x empty_y', 'empty_x empty_y']

# phi=2/9

# # bosons with filled LLL
# product_state = ['1_x 1_y', '1_x 1_y', '0_x 0_y', '0_x 0_y',
#                  '0_x 0_y', '0_x 0_y', '0_x 0_y', '0_x 0_y',
#                  '0_x 0_y', '0_x 0_y', '0_x 0_y', '0_x 0_y',
#                  '0_x 0_y', '0_x 0_y', '0_x 0_y', '0_x 0_y',
#                  '0_x 0_y', '0_x 0_y']
# # bosons with 1/2-filled LLL
# product_state = ['1_x 1_y', '0_x 0_y', '0_x 0_y', '0_x 0_y',
#                  '0_x 0_y', '0_x 0_y', '0_x 0_y', '0_x 0_y',
#                  '0_x 0_y', '0_x 0_y', '0_x 0_y', '0_x 0_y',
#                  '0_x 0_y', '0_x 0_y', '0_x 0_y', '0_x 0_y',
#                  '0_x 0_y', '0_x 0_y']
# # fermions with filled LLL
# product_state = ['full_x full_y', 'full_x full_y', 'empty_x empty_y', 'empty_x empty_y',
#                  'empty_x empty_y', 'empty_x empty_y', 'empty_x empty_y', 'empty_x empty_y',
#                  'empty_x empty_y', 'empty_x empty_y', 'empty_x empty_y', 'empty_x empty_y',
#                  'empty_x empty_y', 'empty_x empty_y', 'empty_x empty_y', 'empty_x empty_y',
#                  'empty_x empty_y', 'empty_x empty_y']
# # fermions with filled LLL
# product_state = ['full_x full_y', 'empty_x empty_y', 'empty_x empty_y', 'empty_x empty_y',
#                  'empty_x empty_y', 'empty_x empty_y', 'empty_x empty_y', 'empty_x empty_y',
#                  'empty_x empty_y', 'empty_x empty_y', 'empty_x empty_y', 'empty_x empty_y',
#                  'empty_x empty_y', 'empty_x empty_y', 'empty_x empty_y', 'empty_x empty_y',
#                  'empty_x empty_y', 'empty_x empty_y']

# phi=2/11

# # bosons with filled LLL
# product_state = ['1_x 1_y', '1_x 1_y', '0_x 0_y', '0_x 0_y', '0_x 0_y',
#                  '0_x 0_y', '0_x 0_y', '0_x 0_y', '0_x 0_y', '0_x 0_y',
#                  '0_x 0_y', '0_x 0_y', '0_x 0_y', '0_x 0_y', '0_x 0_y',
#                  '0_x 0_y', '0_x 0_y', '0_x 0_y', '0_x 0_y', '0_x 0_y',
#                  '0_x 0_y', '0_x 0_y']
# # bosons with 1/2-filled LLL
# product_state = ['1_x 1_y', '0_x 0_y', '0_x 0_y', '0_x 0_y', '0_x 0_y',
#                  '0_x 0_y', '0_x 0_y', '0_x 0_y', '0_x 0_y', '0_x 0_y',
#                  '0_x 0_y', '0_x 0_y', '0_x 0_y', '0_x 0_y', '0_x 0_y',
#                  '0_x 0_y', '0_x 0_y', '0_x 0_y', '0_x 0_y', '0_x 0_y',
#                  '0_x 0_y', '0_x 0_y']
# # fermions with filled LLL
# product_state = ['full_x full_y', 'full_x full_y', 'empty_x empty_y', 'empty_x empty_y', 'empty_x empty_y',
#                  'empty_x empty_y', 'empty_x empty_y', 'empty_x empty_y', 'empty_x empty_y', 'empty_x empty_y',
#                  'empty_x empty_y', 'empty_x empty_y', 'empty_x empty_y', 'empty_x empty_y', 'empty_x empty_y',
#                  'empty_x empty_y', 'empty_x empty_y', 'empty_x empty_y', 'empty_x empty_y', 'empty_x empty_y',
#                  'empty_x empty_y', 'empty_x empty_y']
# # fermions with 1/2-filled LLL
# product_state = ['full_x full_y', 'empty_x empty_y', 'empty_x empty_y', 'empty_x empty_y', 'empty_x empty_y',
#                  'empty_x empty_y', 'empty_x empty_y', 'empty_x empty_y', 'empty_x empty_y', 'empty_x empty_y',
#                  'empty_x empty_y', 'empty_x empty_y', 'empty_x empty_y', 'empty_x empty_y', 'empty_x empty_y',
#                  'empty_x empty_y', 'empty_x empty_y', 'empty_x empty_y', 'empty_x empty_y', 'empty_x empty_y',
#                  'empty_x empty_y', 'empty_x empty_y']
# # fermions with double-filled LLL
# product_state = ['full_x full_y', 'full_x full_y', 'empty_x empty_y', 'empty_x empty_y', 'empty_x empty_y',
#                  'empty_x empty_y', 'empty_x empty_y', 'empty_x empty_y', 'empty_x empty_y', 'empty_x empty_y',
#                  'empty_x empty_y', 'full_x full_y', 'full_x full_y', 'empty_x empty_y', 'empty_x empty_y',
#                  'empty_x empty_y', 'empty_x empty_y', 'empty_x empty_y', 'empty_x empty_y', 'empty_x empty_y',
#                  'empty_x empty_y', 'empty_x empty_y']
