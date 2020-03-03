"""hex_1_hex_5_orbital model"""

import numpy as np
import sys
import time

from tenpy.models.model import CouplingMPOModel
from tenpy.tools.params import get_parameter
from tenpy.networks.site import BosonSite, FermionSite, GroupedSite
from tenpy.models.lattice import Honeycomb


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
        qvalue = get_parameter(model_params, 'phi', (1, 4), self.name)[1]
        Lx_MUC = get_parameter(model_params, 'Lx_MUC', 1, self.name)
        Lx = Lx_MUC * qvalue
        Ly = get_parameter(model_params, 'Ly', 4, self.name)
        bc_x = 'periodic' if bc_MPS == 'infinite' else 'open'  # Next line needs default
        bc_x = get_parameter(model_params, 'bc_x', bc_x, self.name)
        bc_y = get_parameter(model_params, 'bc_y', 'cylinder', self.name)
        assert bc_y in ['cylinder', 'ladder']
        bc_y = 'periodic' if bc_y == 'cylinder' else 'open'
        if bc_MPS == 'infinite' and bc_x == 'open':
            raise ValueError("You need to use 'periodic' `bc_x` for infinite systems!")
        lat = Honeycomb(Lx, Ly, site, order=order, bc=[bc_x, bc_y], bc_MPS=bc_MPS)
        return lat

    def init_terms(self, model_params):
        t1 = get_parameter(model_params, 't1', 1, self.name)
        t2 = get_parameter(model_params, 't2', 0, self.name)
        t2dash = get_parameter(model_params, 't2dash', 0, self.name)
        U = get_parameter(model_params, 'U', 0, self.name, True)
        V = get_parameter(model_params, 'V', 0, self.name, True)
        phi_ext = 2 * np.pi * get_parameter(model_params, 'phi_ext', 0., self.name)
        phi_p, phi_q = get_parameter(model_params, 'phi', (1, 4), self.name)
        Lx_MUC = get_parameter(model_params, 'Lx_MUC', 1, self.name)
        phi = 2 * np.pi * phi_p / phi_q

        # onsite interaction
        for u in range(len(self.lat.unit_cell)):
            print("u in range(len(self.lat.unit_cell)) = ", u)
            self.add_onsite(U, u, 'Nx Ny')

        # t1 term ######################################################################################################

        u1, u2, dx = (0, 1, np.array([0, -1]))  # down
        m = np.arange(0, 2 * phi_q * Lx_MUC, 2)
        t_phi = self.coupling_strength_add_ext_flux(t1, dx, [0, phi_ext]) \
                * np.exp(-1j * (phi / 3) * m)[:, np.newaxis]
        for orbital in ['x', 'y']:
            self.add_coupling(t_phi, u1, 'Bd' + orbital, u2, 'B' + orbital, dx)
            self.add_coupling(np.conj(t_phi), u2, 'Bd' + orbital, u1, 'B' + orbital, -dx)  # h.c.

        u1, u2, dx = (0, 1, np.array([-1, 0]))  # upper left
        m = np.roll(np.arange(0, 2 * phi_q * Lx_MUC, 2), -1)  # match convention for strength argument of add_coupling
        t_phi = self.coupling_strength_add_ext_flux(t1, dx, [0, phi_ext]) \
                * np.exp(1j * (phi / 6) * (m - 1 / 2))[:, np.newaxis]
        for orbital in ['x', 'y']:
            self.add_coupling(t_phi, u1, 'Bd' + orbital, u2, 'B' + orbital, dx)
            self.add_coupling(np.conj(t_phi), u2, 'Bd' + orbital, u1, 'B' + orbital, -dx)  # h.c.

        u1, u2, dx = (0, 1, np.array([0, 0]))  # upper right
        m = np.arange(0, 2 * phi_q * Lx_MUC, 2)
        t_phi = self.coupling_strength_add_ext_flux(t1, dx, [0, phi_ext]) \
                * np.exp(1j * (phi / 6) * (m + 1 / 2))[:, np.newaxis]
        for orbital in ['x', 'y']:
            self.add_coupling(t_phi, u1, 'Bd' + orbital, u2, 'B' + orbital, dx)
            self.add_coupling(np.conj(t_phi), u2, 'Bd' + orbital, u1, 'B' + orbital, -dx)  # h.c.

        # t2 term ######################################################################################################

        for sublattice in [0, 1]:  # A and B sublattices

            u1, u2, dx = (sublattice, sublattice, np.array([-1, 2]))  # up
            m = np.roll(np.arange(sublattice, 2 * phi_q * Lx_MUC, 2),
                        -1)  # match convention for strength argument of add_coupling
            t_phi = self.coupling_strength_add_ext_flux(t2, dx, [0, phi_ext]) \
                    * np.exp(1j * phi * m)[:, np.newaxis]
            t2dash_phi = self.coupling_strength_add_ext_flux(t2dash, dx, [0, phi_ext]) \
                         * np.exp(1j * phi * m)[:, np.newaxis]
            for orbital in ['x', 'y']:
                self.add_coupling(t_phi, u1, 'Bd' + orbital, u2, 'B' + orbital, dx)
                self.add_coupling(np.conj(t_phi), u2, 'Bd' + orbital, u1, 'B' + orbital, -dx)  # h.c.
            # fifthNN orbital mixing
            self.add_coupling(t2dash_phi, u1, 'Bdx', u2, 'By', dx)
            self.add_coupling(np.conj(t2dash_phi), u2, 'Bdy', u1, 'Bx', -dx)  # h.c.
            self.add_coupling(-t2dash_phi, u1, 'Bdy', u2, 'Bx', dx)
            self.add_coupling(np.conj(-t2dash_phi), u2, 'Bdx', u1, 'By', -dx)  # h.c.

            u1, u2, dx = (sublattice, sublattice, np.array([-2, 1]))  # bottom right
            m = np.roll(np.arange(sublattice, 2 * phi_q * Lx_MUC, 2),
                        -2)  # match convention for strength argument of add_coupling
            t_phi = self.coupling_strength_add_ext_flux(t2, dx, [0, phi_ext]) \
                    * np.exp(-1j * (phi / 2) * (m + 3 / 2))[:, np.newaxis]
            t2dash_phi = self.coupling_strength_add_ext_flux(t2dash, dx, [0, phi_ext]) \
                         * np.exp(-1j * (phi / 2) * (m + 3 / 2))[:, np.newaxis]
            for orbital in ['x', 'y']:
                self.add_coupling(t_phi, u1, 'Bd' + orbital, u2, 'B' + orbital, dx)
                self.add_coupling(np.conj(t_phi), u2, 'Bd' + orbital, u1, 'B' + orbital, -dx)  # h.c.
            # fifthNN orbital mixing
            self.add_coupling(t2dash_phi, u1, 'Bdx', u2, 'By', dx)
            self.add_coupling(np.conj(t2dash_phi), u2, 'Bdy', u1, 'Bx', -dx)  # h.c.
            self.add_coupling(-t2dash_phi, u1, 'Bdy', u2, 'Bx', dx)
            self.add_coupling(np.conj(-t2dash_phi), u2, 'Bdx', u1, 'By', -dx)  # h.c.

            u1, u2, dx = (sublattice, sublattice, np.array([-1, -1]))  # bottom left
            m = np.roll(np.arange(sublattice, 2 * phi_q * Lx_MUC, 2),
                        -1)  # match convention for strength argument of add_coupling
            t_phi = self.coupling_strength_add_ext_flux(t2, dx, [0, phi_ext]) \
                    * np.exp(-1j * (phi / 2) * (m - 3 / 2))[:, np.newaxis]
            t2dash_phi = self.coupling_strength_add_ext_flux(t2dash, dx, [0, phi_ext]) \
                         * np.exp(-1j * (phi / 2) * (m - 3 / 2))[:, np.newaxis]
            for orbital in ['x', 'y']:
                self.add_coupling(t_phi, u1, 'Bd' + orbital, u2, 'B' + orbital, dx)
                self.add_coupling(np.conj(t_phi), u2, 'Bd' + orbital, u1, 'B' + orbital, -dx)  # h.c.
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
        filling = get_parameter(model_params, 'filling', (1, 9), self.name)
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
        qvalue = get_parameter(model_params, 'phi', (1, 3), self.name)[1]
        Lx_MUC = get_parameter(model_params, 'Lx_MUC', 1, self.name)
        Lx = Lx_MUC * qvalue
        Ly = get_parameter(model_params, 'Ly', 6, self.name)
        bc_x = 'periodic' if bc_MPS == 'infinite' else 'open'  # Next line needs default
        bc_x = get_parameter(model_params, 'bc_x', bc_x, self.name)
        bc_y = get_parameter(model_params, 'bc_y', 'cylinder', self.name)
        assert bc_y in ['cylinder', 'ladder']
        bc_y = 'periodic' if bc_y == 'cylinder' else 'open'
        if bc_MPS == 'infinite' and bc_x == 'open':
            raise ValueError("You need to use 'periodic' `bc_x` for infinite systems!")
        lat = Honeycomb(Lx, Ly, site, order=order, bc=[bc_x, bc_y], bc_MPS=bc_MPS)
        return lat

    def init_terms(self, model_params):
        t1 = get_parameter(model_params, 't1', 1, self.name)
        t2 = get_parameter(model_params, 't2', 0, self.name)
        t2dash = get_parameter(model_params, 't2dash', 0, self.name)
        U = get_parameter(model_params, 'U', 100, self.name, True)
        V = get_parameter(model_params, 'V', 10, self.name, True)
        phi_ext = 2*np.pi*get_parameter(model_params, 'phi_ext', 0., self.name)
        phi_p, phi_q = get_parameter(model_params, 'phi', (1, 3), self.name)
        Lx_MUC = get_parameter(model_params, 'Lx_MUC', 1, self.name)
        phi = 2 * np.pi * phi_p / phi_q

        # onsite interaction
        for u in range(len(self.lat.unit_cell)):
            print("u in range(len(self.lat.unit_cell)) = ", u)
            self.add_onsite(U, u, 'Nx Ny')

        # t1 term ######################################################################################################

        u1, u2, dx = (0, 1, np.array([0, -1]))  # down
        m = np.arange(0, 2*phi_q*Lx_MUC, 2)
        t_phi = self.coupling_strength_add_ext_flux(t1, dx, [0, phi_ext]) \
                * np.exp(-1j * (phi / 3) * m)[:, np.newaxis]
        for orbital in ['x', 'y']:
            self.add_coupling(t_phi, u1, 'Cd'+orbital, u2, 'C'+orbital, dx)
            self.add_coupling(np.conj(t_phi), u2, 'Cd'+orbital, u1, 'C'+orbital, -dx)  # h.c.
        # NN interaction
        self.add_coupling(V, u1, 'Ntot', u2, 'Ntot', dx)

        u1, u2, dx = (0, 1, np.array([-1, 0]))  # upper left
        m = np.roll(np.arange(0, 2 * phi_q * Lx_MUC, 2), -1)  # match convention for strength argument of add_coupling
        t_phi = self.coupling_strength_add_ext_flux(t1, dx, [0, phi_ext]) \
                * np.exp(1j * (phi / 6) * (m - 1 / 2))[:, np.newaxis]
        for orbital in ['x', 'y']:
            self.add_coupling(t_phi, u1, 'Cd' + orbital, u2, 'C' + orbital, dx)
            self.add_coupling(np.conj(t_phi), u2, 'Cd' + orbital, u1, 'C' + orbital, -dx)  # h.c.
        # NN interaction
        self.add_coupling(V, u1, 'Ntot', u2, 'Ntot', dx)

        u1, u2, dx = (0, 1, np.array([0, 0]))  #upper right
        m = np.arange(0, 2*phi_q*Lx_MUC, 2)
        t_phi = self.coupling_strength_add_ext_flux(t1, dx, [0, phi_ext]) \
                * np.exp(1j * (phi / 6) * (m + 1/2))[:, np.newaxis]
        for orbital in ['x', 'y']:
            self.add_coupling(t_phi, u1, 'Cd' + orbital, u2, 'C' + orbital, dx)
            self.add_coupling(np.conj(t_phi), u2, 'Cd' + orbital, u1, 'C' + orbital, -dx)  # h.c.
        # NN interaction
        self.add_coupling(V, u1, 'Ntot', u2, 'Ntot', dx)

        # t2 term ######################################################################################################

        for sublattice in [0, 1]:  #A and B sublattices

            u1, u2, dx = (sublattice, sublattice, np.array([-1, 2]))  # up
            m = np.roll(np.arange(sublattice, 2*phi_q*Lx_MUC, 2), -1)  # match convention for strength argument of add_coupling
            t_phi = self.coupling_strength_add_ext_flux(t2, dx, [0, phi_ext]) \
                    * np.exp(1j * phi * m)[:, np.newaxis]
            t2dash_phi = self.coupling_strength_add_ext_flux(t2dash, dx, [0, phi_ext]) \
                    * np.exp(1j * phi * m)[:, np.newaxis]
            for orbital in ['x', 'y']:
                self.add_coupling(t_phi, u1, 'Cd' + orbital, u2, 'C' + orbital, dx)
                self.add_coupling(np.conj(t_phi), u2, 'Cd' + orbital, u1, 'C' + orbital, -dx)  # h.c.
            # fifthNN orbital mixing
            self.add_coupling(t2dash_phi, u1, 'Cdx', u2, 'Cy', dx)
            self.add_coupling(np.conj(t2dash_phi), u2, 'Cdy', u1, 'Cx', -dx)  # h.c.
            self.add_coupling(-t2dash_phi, u1, 'Cdy', u2, 'Cx', dx)
            self.add_coupling(np.conj(-t2dash_phi), u2, 'Cdx', u1, 'Cy', -dx)  # h.c.

            u1, u2, dx = (sublattice, sublattice, np.array([-2, 1]))  # bottom right
            m = np.roll(np.arange(sublattice, 2*phi_q*Lx_MUC, 2), -2)  # match convention for strength argument of add_coupling
            t_phi = self.coupling_strength_add_ext_flux(t2, dx, [0, phi_ext]) \
                    * np.exp(-1j * (phi/2) * (m + 3/2))[:, np.newaxis]
            t2dash_phi = self.coupling_strength_add_ext_flux(t2dash, dx, [0, phi_ext]) \
                    * np.exp(-1j * (phi/2) * (m + 3/2))[:, np.newaxis]
            for orbital in ['x', 'y']:
                self.add_coupling(t_phi, u1, 'Cd' + orbital, u2, 'C' + orbital, dx)
                self.add_coupling(np.conj(t_phi), u2, 'Cd' + orbital, u1, 'C' + orbital, -dx)  # h.c.
            # fifthNN orbital mixing
            self.add_coupling(t2dash_phi, u1, 'Cdx', u2, 'Cy', dx)
            self.add_coupling(np.conj(t2dash_phi), u2, 'Cdy', u1, 'Cx', -dx)  # h.c.
            self.add_coupling(-t2dash_phi, u1, 'Cdy', u2, 'Cx', dx)
            self.add_coupling(np.conj(-t2dash_phi), u2, 'Cdx', u1, 'Cy', -dx)  # h.c.

            u1, u2, dx = (sublattice, sublattice, np.array([-1, -1]))  # bottom left
            m = np.roll(np.arange(sublattice, 2*phi_q*Lx_MUC, 2), -1)  # match convention for strength argument of add_coupling
            t_phi = self.coupling_strength_add_ext_flux(t2, dx, [0, phi_ext]) \
                    * np.exp(-1j * (phi/2) * (m - 3/2))[:, np.newaxis]
            t2dash_phi = self.coupling_strength_add_ext_flux(t2dash, dx, [0, phi_ext]) \
                    * np.exp(-1j * (phi/2) * (m - 3/2))[:, np.newaxis]
            for orbital in ['x', 'y']:
                self.add_coupling(t_phi, u1, 'Cd' + orbital, u2, 'C' + orbital, dx)
                self.add_coupling(np.conj(t_phi), u2, 'Cd' + orbital, u1, 'C' + orbital, -dx)  # h.c.
            # fifthNN orbital mixing
            self.add_coupling(t2dash_phi, u1, 'Cdx', u2, 'Cy', dx)
            self.add_coupling(np.conj(t2dash_phi), u2, 'Cdy', u1, 'Cx', -dx)  # h.c.
            self.add_coupling(-t2dash_phi, u1, 'Cdy', u2, 'Cx', dx)
            self.add_coupling(np.conj(-t2dash_phi), u2, 'Cdx', u1, 'Cy', -dx)  # h.c.


if __name__ == "__main__":

    t0 = time.time()

    model_params = dict(conserve='N', t1=1, t2=-0.025, t2dash=0.1, filling=(int(1), int(9)), phi=(int(1), int(3)),
                        Lx_MUC=1, Ly=6, V=10,  # system params
                        bc_MPS='infinite', bc_x='periodic', bc_y='cylinder', order='Cstyle',  # MPS params
                        verbose=1, phi_ext=1.)  # utility
    M = FermionicHex1Hex5OrbitalModel(model_params)
    print("max MPO bond dimension = ", max(M.H_MPO.chi))

    print(time.time() - t0)
