"""hex_1 model"""

import numpy as np
import sys
import time

from tenpy.models.model import CouplingMPOModel
from tenpy.tools.params import get_parameter
from tenpy.networks.site import BosonSite, FermionSite
from tenpy.models.lattice import Honeycomb


class BosonicHex1Model(CouplingMPOModel):

    def __init__(self, model_params):
        CouplingMPOModel.__init__(self, model_params)

    def init_sites(self, model_params):
        Nmax = get_parameter(model_params, 'Nmax', 1, self.name)
        conserve = get_parameter(model_params, 'conserve', 'N', self.name)
        filling = get_parameter(model_params, 'filling', (1, 8), self.name)
        filling = filling[0] / filling[1]
        site = BosonSite(Nmax=Nmax, conserve=conserve, filling=filling)
        return site

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
        t = get_parameter(model_params, 't', 1, self.name, True)
        phi_ext = 2 * np.pi * get_parameter(model_params, 'phi_ext', 0., self.name)
        phi_p, phi_q = get_parameter(model_params, 'phi', (1, 4), self.name)
        Lx_MUC = get_parameter(model_params, 'Lx_MUC', 1, self.name)
        phi = 2 * np.pi * phi_p / phi_q

        u1, u2, dx = (0, 1, np.array([0, -1]))  # down
        m = np.arange(0, 2 * phi_q * Lx_MUC, 2)
        t_phi = self.coupling_strength_add_ext_flux(t, dx, [0, phi_ext]) \
                * np.exp(-1j * (phi / 3) * m)[:, np.newaxis]
        self.add_coupling(t_phi, u1, 'Bd', u2, 'B', dx)
        self.add_coupling(np.conj(t_phi), u2, 'Bd', u1, 'B', -dx)  # h.c.

        u1, u2, dx = (0, 1, np.array([-1, 0]))  # upper left
        m = np.roll(np.arange(0, 2 * phi_q * Lx_MUC, 2), -1)  # match convention for strength argument of add_coupling
        t_phi = self.coupling_strength_add_ext_flux(t, dx, [0, phi_ext]) \
                * np.exp(1j * (phi / 6) * (m - 1 / 2))[:, np.newaxis]
        self.add_coupling(t_phi, u1, 'Bd', u2, 'B', dx)
        self.add_coupling(np.conj(t_phi), u2, 'Bd', u1, 'B', -dx)  # h.c.

        u1, u2, dx = (0, 1, np.array([0, 0]))  # upper right
        m = np.arange(0, 2 * phi_q * Lx_MUC, 2)
        t_phi = self.coupling_strength_add_ext_flux(t, dx, [0, phi_ext]) \
                * np.exp(1j * (phi / 6) * (m + 1 / 2))[:, np.newaxis]
        self.add_coupling(t_phi, u1, 'Bd', u2, 'B', dx)
        self.add_coupling(np.conj(t_phi), u2, 'Bd', u1, 'B', -dx)  # h.c.


class FermionicHex1Model(CouplingMPOModel):

    def __init__(self, model_params):
        CouplingMPOModel.__init__(self, model_params)

    def init_sites(self, model_params):
        conserve = get_parameter(model_params, 'conserve', 'N', self.name)
        filling = get_parameter(model_params, 'filling', (1, 9), self.name)
        filling = filling[0] / filling[1]
        site = FermionSite(conserve=conserve, filling=filling)
        return site

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
        t = get_parameter(model_params, 't', 1, self.name, True)
        V = get_parameter(model_params, 'V', 10, self.name, True)
        phi_ext = 2 * np.pi * get_parameter(model_params, 'phi_ext', 0., self.name)
        phi_p, phi_q = get_parameter(model_params, 'phi', (1, 3), self.name)
        Lx_MUC = get_parameter(model_params, 'Lx_MUC', 1, self.name)
        phi = 2 * np.pi * phi_p / phi_q

        u1, u2, dx = (0, 1, np.array([0, -1]))  # down
        m = np.arange(0, 2 * phi_q *Lx_MUC, 2)
        t_phi = self.coupling_strength_add_ext_flux(t, dx, [0, phi_ext]) \
                * np.exp(-1j * (phi / 3) * m)[:, np.newaxis]
        self.add_coupling(t_phi, u1, 'Cd', u2, 'C', dx)
        self.add_coupling(np.conj(t_phi), u2, 'Cd', u1, 'C', -dx)  # h.c.
        # NN interaction
        self.add_coupling(V, u1, 'N', u2, 'N', dx)

        u1, u2, dx = (0, 1, np.array([-1, 0]))  # upper left
        m = np.roll(np.arange(0, 2 * phi_q * Lx_MUC, 2), -1)  # match convention for strength argument of add_coupling
        t_phi = self.coupling_strength_add_ext_flux(t, dx, [0, phi_ext]) \
                * np.exp(1j * (phi / 6) * (m - 1/2))[:, np.newaxis]
        self.add_coupling(t_phi, u1, 'Cd', u2, 'C', dx)
        self.add_coupling(np.conj(t_phi), u2, 'Cd', u1, 'C', -dx)  # h.c.
        # NN interaction
        self.add_coupling(V, u1, 'N', u2, 'N', dx)

        u1, u2, dx = (0, 1, np.array([0, 0]))  # upper right
        m = np.arange(0, 2 * phi_q * Lx_MUC, 2)
        t_phi = self.coupling_strength_add_ext_flux(t, dx, [0, phi_ext]) \
                * np.exp(1j * (phi / 6) * (m + 1 / 2))[:, np.newaxis]
        self.add_coupling(t_phi, u1, 'Cd', u2, 'C', dx)
        self.add_coupling(np.conj(t_phi), u2, 'Cd', u1, 'C', -dx)  # h.c.
        # NN interaction
        self.add_coupling(V, u1, 'N', u2, 'N', dx)


if __name__ == "__main__":

    t0 = time.time()

    model_params = dict(conserve='N', t=1, filling=(int(1), int(9)), phi=(int(1), int(3)),
                        Lx_MUC=1, Ly=6, V=10,  # system params
                        bc_MPS='infinite', bc_x='periodic', bc_y='cylinder', order='default',  # MPS params
                        verbose=1, phi_ext=0)  # utility
    M = FermionicHex1Model(model_params)
    print("max MPO bond dimension = ", max(M.H_MPO.chi))

    print(time.time() - t0)

    # import matplotlib.pyplot as plt
    # ax = plt.gca()
    # M.lat.plot_sites(ax)
    # M.lat.plot_order(ax)
    # M.lat.plot_coupling(ax, M.lat.pairs["nearest_neighbors"])
    # plt.show()
