"""hex_1_hex_5_orbital model"""

import numpy as np
import sys
import time

from tenpy.models.model import CouplingMPOModel
from tenpy.tools.params import get_parameter
from tenpy.networks.site import BosonSite, FermionSite, GroupedSite
from lattices.old.MagneticHoneycomb import MagneticHoneycomb


class BosonicHex1Hex5OrbitalModel(CouplingMPOModel):

    def __init__(self, model_params):
        CouplingMPOModel.__init__(self, model_params)

    def init_sites(self, model_params):
        Nmax = get_parameter(model_params, 'Nmax', 1, self.name)
        conserve = get_parameter(model_params, 'conserve', 'N', self.name)
        filling = get_parameter(model_params, 'filling', (1, 8), self.name)
        filling = filling[0] / filling[1]
        bs = BosonSite(Nmax=Nmax, conserve=conserve, filling=filling)
        gs = GroupedSite([bs, bs], labels=['x', 'y'], charges='same')
        gs.add_op('Ntot', gs.Nx + gs.Ny, False)
        print(sorted(gs.opnames))
        print(gs.state_labels)
        return gs

    def init_lattice(self, model_params):
        bc_MPS = get_parameter(model_params, 'bc_MPS', 'infinite', self.name)
        order = get_parameter(model_params, 'order', 'default', self.name)
        site = self.init_sites(model_params)
        Lx = get_parameter(model_params, 'Lx', 1, self.name)
        Ly = get_parameter(model_params, 'Ly', 4, self.name)
        qvalue = get_parameter(model_params, 'phi', (1, 4), self.name)[1]
        bc_x = 'periodic' if bc_MPS == 'infinite' else 'open'  # Next line needs default
        bc_x = get_parameter(model_params, 'bc_x', bc_x, self.name)
        bc_y = get_parameter(model_params, 'bc_y', 'cylinder', self.name)
        assert bc_y in ['cylinder', 'ladder']
        bc_y = 'periodic' if bc_y == 'cylinder' else 'open'
        if bc_MPS == 'infinite' and bc_x == 'open':
            raise ValueError("You need to use 'periodic' `bc_x` for infinite systems!")
        lat = MagneticHoneycomb(Lx, Ly, site, order=order, bc=[bc_x, bc_y], bc_MPS=bc_MPS, qvalue=qvalue)
        return lat

    def init_terms(self, model_params):
        t1 = get_parameter(model_params, 't1', 1, self.name)
        t2 = get_parameter(model_params, 't2', 0, self.name)
        t2dash = get_parameter(model_params, 't2dash', 0, self.name)
        # t1, t2, t2dash = 2, 0.05, 0.2  # Yuan parameters
        # t1, t2, t2dash = 0.331, -0.01, 0.097  # Koshino parameters
        # invt2dash = get_parameter(model_params, 'invt2dash', 1, self.name, True)
        # t1, t2, t2dash = 1, -0.025, 1 / invt2dash
        U = get_parameter(model_params, 'U', 0, self.name, True)
        V = get_parameter(model_params, 'V', 0, self.name, True)
        phi_ext = 2 * np.pi * get_parameter(model_params, 'phi_ext', 0., self.name)

        phi_p, phi_q = get_parameter(model_params, 'phi', (1, 4), self.name)
        phi = 2 * np.pi * phi_p / phi_q

        # t1 term ###
        for i in range(0, 2 * phi_q, 2):
            for u1, u2, dx in getattr(self.lat, "NN{}d".format(i)):
                t_phi = self.coupling_strength_add_ext_flux(t1, dx, [0, phi_ext]) \
                        * np.exp(-1j * (phi / 3) * i)
                for orbital in ['x', 'y']:
                    self.add_coupling(t_phi, u1, 'Bd' + orbital, u2, 'B' + orbital, dx)
                    self.add_coupling(np.conj(t_phi), u2, 'Bd' + orbital, u1, 'B' + orbital, -dx)  # h.c.
            for u1, u2, dx in getattr(self.lat, "NN{}ur".format(i)):
                t_phi = self.coupling_strength_add_ext_flux(t1, dx, [0, phi_ext]) \
                        * np.exp(1j * (phi / 6) * (i + 1 / 2))
                for orbital in ['x', 'y']:
                    self.add_coupling(t_phi, u1, 'Bd' + orbital, u2, 'B' + orbital, dx)
                    self.add_coupling(np.conj(t_phi), u2, 'Bd' + orbital, u1, 'B' + orbital, -dx)  # h.c.
            for u1, u2, dx in getattr(self.lat, "NN{}ul".format(i)):
                t_phi = self.coupling_strength_add_ext_flux(t1, dx, [0, phi_ext]) \
                        * np.exp(1j * (phi / 6) * (i - 1 / 2))
                for orbital in ['x', 'y']:
                    self.add_coupling(t_phi, u1, 'Bd' + orbital, u2, 'B' + orbital, dx)
                    self.add_coupling(np.conj(t_phi), u2, 'Bd' + orbital, u1, 'B' + orbital, -dx)  # h.c.

        # t2 term ###
        for i in range(2 * phi_q):
            for u1, u2, dx in getattr(self.lat, "fifthNN{}u".format(i)):
                t_phi = self.coupling_strength_add_ext_flux(t2, dx, [0, phi_ext]) \
                        * np.exp(1j * phi * i)
                t2dash_phi = self.coupling_strength_add_ext_flux(t2dash, dx, [0, phi_ext]) \
                             * np.exp(1j * phi * i)
                for orbital in ['x', 'y']:
                    self.add_coupling(t_phi, u1, 'Bd' + orbital, u2, 'B' + orbital, dx)
                    self.add_coupling(np.conj(t_phi), u2, 'Bd' + orbital, u1, 'B' + orbital, -dx)  # h.c.
                # fifthNN orbital mixing
                self.add_coupling(t2dash_phi, u1, 'Bdx', u2, 'By', dx)
                self.add_coupling(np.conj(t2dash_phi), u2, 'Bdy', u1, 'Bx', -dx)  # h.c.
                self.add_coupling(-t2dash_phi, u1, 'Bdy', u2, 'Bx', dx)
                self.add_coupling(np.conj(-t2dash_phi), u2, 'Bdx', u1, 'By', -dx)  # h.c.
            for u1, u2, dx in getattr(self.lat, "fifthNN{}br".format(i)):
                t_phi = self.coupling_strength_add_ext_flux(t2, dx, [0, phi_ext]) \
                        * np.exp(-1j * (phi / 2) * (i + 3 / 2))
                t2dash_phi = self.coupling_strength_add_ext_flux(t2dash, dx, [0, phi_ext]) \
                             * np.exp(-1j * (phi / 2) * (i + 3 / 2))
                for orbital in ['x', 'y']:
                    self.add_coupling(t_phi, u1, 'Bd' + orbital, u2, 'B' + orbital, dx)
                    self.add_coupling(np.conj(t_phi), u2, 'Bd' + orbital, u1, 'B' + orbital, -dx)  # h.c.
                # fifthNN orbital mixing
                self.add_coupling(t2dash_phi, u1, 'Bdx', u2, 'By', dx)
                self.add_coupling(np.conj(t2dash_phi), u2, 'Bdy', u1, 'Bx', -dx)  # h.c.
                self.add_coupling(-t2dash_phi, u1, 'Bdy', u2, 'Bx', dx)
                self.add_coupling(np.conj(-t2dash_phi), u2, 'Bdx', u1, 'By', -dx)  # h.c.
            for u1, u2, dx in getattr(self.lat, "fifthNN{}bl".format(i)):
                t_phi = self.coupling_strength_add_ext_flux(t2, dx, [0, phi_ext]) \
                        * np.exp(-1j * (phi / 2) * (i - 3 / 2))
                t2dash_phi = self.coupling_strength_add_ext_flux(t2dash, dx, [0, phi_ext]) \
                             * np.exp(-1j * (phi / 2) * (i - 3 / 2))
                for orbital in ['x', 'y']:
                    self.add_coupling(t_phi, u1, 'Bd' + orbital, u2, 'B' + orbital, dx, 'JW', True)
                    self.add_coupling(np.conj(t_phi), u2, 'Bd' + orbital, u1, 'B' + orbital, -dx, 'JW', True)  # h.c.
                # fifthNN orbital mixing
                self.add_coupling(t2dash_phi, u1, 'Bdx', u2, 'By', dx)
                self.add_coupling(np.conj(t2dash_phi), u2, 'Bdy', u1, 'Bx', -dx)  # h.c.
                self.add_coupling(-t2dash_phi, u1, 'Bdy', u2, 'Bx', dx)
                self.add_coupling(np.conj(-t2dash_phi), u2, 'Bdx', u1, 'By', -dx)  # h.c.


class FermionicHex1Hex5OrbitalModel(CouplingMPOModel):

    def __init__(self, model_params):
        CouplingMPOModel.__init__(self, model_params)

    def init_sites(self, model_params):
        conserve = get_parameter(model_params, 'conserve', 'N', self.name)
        filling = get_parameter(model_params, 'n', (1, 9), self.name)
        filling = filling[0] / filling[1]
        fs = FermionSite(conserve=conserve, filling=filling)
        gs = GroupedSite([fs, fs], labels=['x', 'y'], charges='same')
        gs.add_op('Ntot', gs.Nx + gs.Ny, False)
        print(sorted(gs.opnames))
        print(gs.state_labels)
        return gs

    def init_lattice(self, model_params):
        bc_MPS = get_parameter(model_params, 'bc_MPS', 'infinite', self.name)
        order = get_parameter(model_params, 'order', 'default', self.name)
        site = self.init_sites(model_params)
        Lx = get_parameter(model_params, 'LxMUC', 1, self.name)
        Ly = get_parameter(model_params, 'Ly', 6, self.name)
        qvalue = get_parameter(model_params, 'nphi', (1, 3), self.name)[1]
        bc_x = 'periodic' if bc_MPS == 'infinite' else 'open'  # Next line needs default
        bc_x = get_parameter(model_params, 'bc_x', bc_x, self.name)
        bc_y = get_parameter(model_params, 'bc_y', 'cylinder', self.name)
        assert bc_y in ['cylinder', 'ladder']
        bc_y = 'periodic' if bc_y == 'cylinder' else 'open'
        if bc_MPS == 'infinite' and bc_x == 'open':
            raise ValueError("You need to use 'periodic' `bc_x` for infinite systems!")
        lat = MagneticHoneycomb(Lx, Ly, site, order=order, bc=[bc_x, bc_y], bc_MPS=bc_MPS, q=qvalue)
        return lat

    def init_terms(self, model_params):
        t1 = get_parameter(model_params, 't1', 1, self.name)
        t2 = get_parameter(model_params, 't5', 0, self.name)
        t2dash = get_parameter(model_params, 't5dash', 0, self.name)
        # t1, t2, t2dash = 2, 0.05, 0.2  # Yuan parameters
        # t1, t2, t2dash = 0.331, -0.01, 0.097  # Koshino parameters
        # invt2dash = get_parameter(model_params, 'invt2dash', 1, self.name, True)
        # t1, t2, t2dash = 1, -0.025, 1/invt2dash
        U = get_parameter(model_params, 'U', 100, self.name, True)
        V = get_parameter(model_params, 'V', 10, self.name, True)
        phi_ext = 2*np.pi*get_parameter(model_params, 'phi', 0., self.name)

        phi_p, phi_q = get_parameter(model_params, 'nphi', (1, 3), self.name)
        phi = 2 * np.pi * phi_p / phi_q

        # onsite interaction
        for u in range(len(self.lat.unit_cell)):
            print("u in range(len(self.lat.unit_cell)) = ", u)
            self.add_onsite(U, u, 'Nx Ny')

        # t1 term ###
        for i in range(0, 2*phi_q, 2):
            for u1, u2, dx in getattr(self.lat, "NN{}d".format(i)):
                t_phi = self.coupling_strength_add_ext_flux(t1, dx, [0, phi_ext]) \
                        * np.exp(-1j * (phi / 3) * i)
                # t2dash_phi = self.coupling_strength_add_ext_flux(t2dash, dx, [0, phi_ext]) \
                #         * np.exp(-1j * (phi / 3) * i)
                for orbital in ['x', 'y']:
                    self.add_coupling(t_phi, u1, 'Cd'+orbital, u2, 'C'+orbital, dx, 'JW', True)
                    self.add_coupling(np.conj(t_phi), u2, 'Cd'+orbital, u1, 'C'+orbital, -dx, 'JW', True)  # h.c.
                # # NN orbital mixing
                # self.add_coupling(t2dash_phi, u1, 'Cdx', u2, 'Cy', dx)
                # self.add_coupling(np.conj(t2dash_phi), u2, 'Cdy', u1, 'Cx', -dx)  # h.c.
                # self.add_coupling(-t2dash_phi, u1, 'Cdy', u2, 'Cx', dx)
                # self.add_coupling(np.conj(-t2dash_phi), u2, 'Cdx', u1, 'Cy', -dx)  # h.c.
                # NN interaction
                self.add_coupling(V, u1, 'Ntot', u2, 'Ntot', dx)
            for u1, u2, dx in getattr(self.lat, "NN{}ur".format(i)):
                t_phi = self.coupling_strength_add_ext_flux(t1, dx, [0, phi_ext]) \
                        * np.exp(1j * (phi / 6) * (i + 1/2))
                # t2dash_phi = self.coupling_strength_add_ext_flux(t2dash, dx, [0, phi_ext]) \
                #         * np.exp(1j * (phi / 6) * (i + 1 / 2))
                for orbital in ['x', 'y']:
                    self.add_coupling(t_phi, u1, 'Cd' + orbital, u2, 'C' + orbital, dx, 'JW', True)
                    self.add_coupling(np.conj(t_phi), u2, 'Cd' + orbital, u1, 'C' + orbital, -dx, 'JW', True)  # h.c.
                # # NN orbital mixing
                # self.add_coupling(t2dash_phi, u1, 'Cdx', u2, 'Cy', dx)
                # self.add_coupling(np.conj(t2dash_phi), u2, 'Cdy', u1, 'Cx', -dx)  # h.c.
                # self.add_coupling(-t2dash_phi, u1, 'Cdy', u2, 'Cx', dx)
                # self.add_coupling(np.conj(-t2dash_phi), u2, 'Cdx', u1, 'Cy', -dx)  # h.c.
                # NN interaction
                self.add_coupling(V, u1, 'Ntot', u2, 'Ntot', dx)
            for u1, u2, dx in getattr(self.lat, "NN{}ul".format(i)):
                t_phi = self.coupling_strength_add_ext_flux(t1, dx, [0, phi_ext]) \
                        * np.exp(1j * (phi / 6) * (i - 1/2))
                # t2dash_phi = self.coupling_strength_add_ext_flux(t2dash, dx, [0, phi_ext]) \
                #         * np.exp(1j * (phi / 6) * (i - 1 / 2))
                for orbital in ['x', 'y']:
                    self.add_coupling(t_phi, u1, 'Cd' + orbital, u2, 'C' + orbital, dx, 'JW', True)
                    self.add_coupling(np.conj(t_phi), u2, 'Cd' + orbital, u1, 'C' + orbital, -dx, 'JW', True)  # h.c.
                # # NN orbital mixing
                # self.add_coupling(t2dash_phi, u1, 'Cdx', u2, 'Cy', dx)
                # self.add_coupling(np.conj(t2dash_phi), u2, 'Cdy', u1, 'Cx', -dx)  # h.c.
                # self.add_coupling(-t2dash_phi, u1, 'Cdy', u2, 'Cx', dx)
                # self.add_coupling(np.conj(-t2dash_phi), u2, 'Cdx', u1, 'Cy', -dx)  # h.c.
                # NN interaction
                self.add_coupling(V, u1, 'Ntot', u2, 'Ntot', dx)

        # # intermediate interactions ###
        # for i in range(0, 2*phi_q):
        #     for u1, u2, dx in getattr(self.lat, "secondNN{}bl".format(i)):
        #         self.add_coupling(V/np.sqrt(3), u1, 'Ntot', u2, 'Ntot', dx)
        #     for u1, u2, dx in getattr(self.lat, "secondNN{}r".format(i)):
        #         self.add_coupling(V/np.sqrt(3), u1, 'Ntot', u2, 'Ntot', dx)
        #     for u1, u2, dx in getattr(self.lat, "secondNN{}ul".format(i)):
        #         self.add_coupling(V/np.sqrt(3), u1, 'Ntot', u2, 'Ntot', dx)
        # for i in range(0, 2*phi_q, 2):
        #     for u1, u2, dx in getattr(self.lat, "thirdNN{}u".format(i)):
        #         self.add_coupling(V/2, u1, 'Ntot', u2, 'Ntot', dx)
        #     for u1, u2, dx in getattr(self.lat, "thirdNN{}br".format(i)):
        #         self.add_coupling(V/2, u1, 'Ntot', u2, 'Ntot', dx)
        #     for u1, u2, dx in getattr(self.lat, "thirdNN{}bl".format(i)):
        #         self.add_coupling(V/2, u1, 'Ntot', u2, 'Ntot', dx)
        # for i in range(0, 2*phi_q, 2):
        #     for u1, u2, dx in getattr(self.lat, "fourthNN{}bur".format(i)):
        #         self.add_coupling(V/np.sqrt(7), u1, 'Ntot', u2, 'Ntot', dx)
        #     for u1, u2, dx in getattr(self.lat, "fourthNN{}uur".format(i)):
        #         self.add_coupling(V/np.sqrt(7), u1, 'Ntot', u2, 'Ntot', dx)
        #     for u1, u2, dx in getattr(self.lat, "fourthNN{}bul".format(i)):
        #         self.add_coupling(V/np.sqrt(7), u1, 'Ntot', u2, 'Ntot', dx)
        #     for u1, u2, dx in getattr(self.lat, "fourthNN{}uul".format(i)):
        #         self.add_coupling(V/np.sqrt(7), u1, 'Ntot', u2, 'Ntot', dx)
        #     for u1, u2, dx in getattr(self.lat, "fourthNN{}bl".format(i)):
        #         self.add_coupling(V/np.sqrt(7), u1, 'Ntot', u2, 'Ntot', dx)
        #     for u1, u2, dx in getattr(self.lat, "fourthNN{}br".format(i)):
        #         self.add_coupling(V/np.sqrt(7), u1, 'Ntot', u2, 'Ntot', dx)

        # t2 term ###
        for i in range(2*phi_q):
            for u1, u2, dx in getattr(self.lat, "fifthNN{}u".format(i)):
                t_phi = self.coupling_strength_add_ext_flux(t2, dx, [0, phi_ext]) \
                        * np.exp(1j * phi * i)
                t2dash_phi = self.coupling_strength_add_ext_flux(t2dash, dx, [0, phi_ext]) \
                        * np.exp(1j * phi * i)
                for orbital in ['x', 'y']:
                    self.add_coupling(t_phi, u1, 'Cd' + orbital, u2, 'C' + orbital, dx, 'JW', True)
                    self.add_coupling(np.conj(t_phi), u2, 'Cd' + orbital, u1, 'C' + orbital, -dx, 'JW', True)  # h.c.
                # fifthNN orbital mixing
                self.add_coupling(t2dash_phi, u1, 'Cdx', u2, 'Cy', dx)
                self.add_coupling(np.conj(t2dash_phi), u2, 'Cdy', u1, 'Cx', -dx)  # h.c.
                self.add_coupling(-t2dash_phi, u1, 'Cdy', u2, 'Cx', dx)
                self.add_coupling(np.conj(-t2dash_phi), u2, 'Cdx', u1, 'Cy', -dx)  # h.c.
                # # fifthNN interaction
                # self.add_coupling(V/3, u1, 'Ntot', u2, 'Ntot', dx)
            for u1, u2, dx in getattr(self.lat, "fifthNN{}br".format(i)):
                t_phi = self.coupling_strength_add_ext_flux(t2, dx, [0, phi_ext]) \
                        * np.exp(-1j * (phi/2) * (i + 3/2))
                t2dash_phi = self.coupling_strength_add_ext_flux(t2dash, dx, [0, phi_ext]) \
                        * np.exp(-1j * (phi/2) * (i + 3/2))
                for orbital in ['x', 'y']:
                    self.add_coupling(t_phi, u1, 'Cd' + orbital, u2, 'C' + orbital, dx, 'JW', True)
                    self.add_coupling(np.conj(t_phi), u2, 'Cd' + orbital, u1, 'C' + orbital, -dx, 'JW', True)  # h.c.
                # fifthNN orbital mixing
                self.add_coupling(t2dash_phi, u1, 'Cdx', u2, 'Cy', dx)
                self.add_coupling(np.conj(t2dash_phi), u2, 'Cdy', u1, 'Cx', -dx)  # h.c.
                self.add_coupling(-t2dash_phi, u1, 'Cdy', u2, 'Cx', dx)
                self.add_coupling(np.conj(-t2dash_phi), u2, 'Cdx', u1, 'Cy', -dx)  # h.c.
                # # fifthNN interaction
                # self.add_coupling(V/3, u1, 'Ntot', u2, 'Ntot', dx)
            for u1, u2, dx in getattr(self.lat, "fifthNN{}bl".format(i)):
                t_phi = self.coupling_strength_add_ext_flux(t2, dx, [0, phi_ext]) \
                        * np.exp(-1j * (phi/2) * (i - 3/2))
                t2dash_phi = self.coupling_strength_add_ext_flux(t2dash, dx, [0, phi_ext]) \
                        * np.exp(-1j * (phi/2) * (i - 3/2))
                for orbital in ['x', 'y']:
                    self.add_coupling(t_phi, u1, 'Cd' + orbital, u2, 'C' + orbital, dx, 'JW', True)
                    self.add_coupling(np.conj(t_phi), u2, 'Cd' + orbital, u1, 'C' + orbital, -dx, 'JW', True)  # h.c.
                # fifthNN orbital mixing
                self.add_coupling(t2dash_phi, u1, 'Cdx', u2, 'Cy', dx)
                self.add_coupling(np.conj(t2dash_phi), u2, 'Cdy', u1, 'Cx', -dx)  # h.c.
                self.add_coupling(-t2dash_phi, u1, 'Cdy', u2, 'Cx', dx)
                self.add_coupling(np.conj(-t2dash_phi), u2, 'Cdx', u1, 'Cy', -dx)  # h.c.
                # # fifthNN interaction
                # self.add_coupling(V/3, u1, 'Ntot', u2, 'Ntot', dx)

        # # t2dash term ###
        # for i in range(2*phi_q):
        #     # positive term
        #     for u1, u2, dx in getattr(self.lat, "fifthNN{}u".format(i)):
        #         t_phi = self.coupling_strength_add_ext_flux(t2dash, dx, [0, phi_ext]) \
        #                 * np.exp(1j * phi * i)
        #         self.add_coupling(t_phi, u1, 'Cdx', u2, 'Cy', dx, 'JW', True)
        #         self.add_coupling(np.conj(t_phi), u2, 'Cdy', u1, 'Cx', -dx, 'JW', True)  # h.c.
        #     for u1, u2, dx in getattr(self.lat, "fifthNN{}br".format(i)):
        #         t_phi = self.coupling_strength_add_ext_flux(t2dash, dx, [0, phi_ext]) \
        #                 * np.exp(-1j * (phi/2) * (i + 3/2))
        #         self.add_coupling(t_phi, u1, 'Cdx', u2, 'Cy', dx, 'JW', True)
        #         self.add_coupling(np.conj(t_phi), u2, 'Cdy', u1, 'Cx', -dx, 'JW', True)  # h.c.
        #     for u1, u2, dx in getattr(self.lat, "fifthNN{}bl".format(i)):
        #         t_phi = self.coupling_strength_add_ext_flux(t2dash, dx, [0, phi_ext]) \
        #                 * np.exp(-1j * (phi/2) * (i - 3/2))
        #         self.add_coupling(t_phi, u1, 'Cdx', u2, 'Cy', dx, 'JW', True)
        #         self.add_coupling(np.conj(t_phi), u2, 'Cdy', u1, 'Cx', -dx, 'JW', True)  # h.c.
        #     # negative term
        #     for u1, u2, dx in getattr(self.lat, "fifthNN{}u".format(i)):
        #         t_phi = self.coupling_strength_add_ext_flux(-t2dash, dx, [0, phi_ext]) \
        #                 * np.exp(1j * phi * i)
        #         self.add_coupling(t_phi, u1, 'Cdy', u2, 'Cx', dx, 'JW', True)
        #         self.add_coupling(np.conj(t_phi), u2, 'Cdx', u1, 'Cy', -dx, 'JW', True)  # h.c.
        #     for u1, u2, dx in getattr(self.lat, "fifthNN{}br".format(i)):
        #         t_phi = self.coupling_strength_add_ext_flux(-t2dash, dx, [0, phi_ext]) \
        #                 * np.exp(-1j * (phi/2) * (i + 3/2))
        #         self.add_coupling(t_phi, u1, 'Cdy', u2, 'Cx', dx, 'JW', True)
        #         self.add_coupling(np.conj(t_phi), u2, 'Cdx', u1, 'Cy', -dx, 'JW', True)  # h.c.
        #     for u1, u2, dx in getattr(self.lat, "fifthNN{}bl".format(i)):
        #         t_phi = self.coupling_strength_add_ext_flux(-t2dash, dx, [0, phi_ext]) \
        #                 * np.exp(-1j * (phi/2) * (i - 3/2))
        #         self.add_coupling(t_phi, u1, 'Cdy', u2, 'Cx', dx, 'JW', True)
        #         self.add_coupling(np.conj(t_phi), u2, 'Cdx', u1, 'Cy', -dx, 'JW', True)  # h.c.


if __name__ == "__main__":

    t0 = time.time()

    model_params = dict(conserve='N', t1=1, t5=-0.025, t5dash=0.1, n=(int(1), int(9)), nphi=(int(1), int(3)), LxMUC=1, Ly=6, V=10,  # system params
                        bc_MPS='infinite', bc_x='periodic', bc_y='cylinder', order='default',  # MPS params
                        verbose=1, phi_ext=1.)  # utility
    M = FermionicHex1Hex5OrbitalModel(model_params)
    print("max MPO bond dimension = ", max(M.H_MPO.chi))

    print(time.time() - t0)
