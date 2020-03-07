# --- python imports
import time
import numpy as np
# --- TeNPy imports
from tenpy.models.model import CouplingMPOModel
from tenpy.tools.params import get_parameter
from tenpy.networks.site import BosonSite, FermionSite
from tenpy.models.lattice import Honeycomb


class HofHex1Hex5Model(CouplingMPOModel):

    def __init__(self, params):
        CouplingMPOModel.__init__(self, params)

    def stats(self, params):
        return get_parameter(params, 'statistics', 'bosons', self.name)

    def init_sites(self, params):
        if self.stats(params) == 'bosons':
            Nmax = get_parameter(params, 'Nmax', 1, self.name)
            conserve = get_parameter(params, 'conserve', 'N', self.name)
            n = get_parameter(params, 'n', (1, 8), self.name)
            n = n[0] / n[1]
            site = BosonSite(Nmax=Nmax, conserve=conserve, filling=n)
        else:
            conserve = get_parameter(params, 'conserve', 'N', self.name)
            n = get_parameter(params, 'n', (1, 9), self.name)
            n = n[0] / n[1]
            site = FermionSite(conserve=conserve, filling=n)
        return site

    def init_lattice(self, params):
        if self.stats(params) == 'bosons':
            nphi_default = (1, 4)
            Ly_default = 4
        else:
            nphi_default = (1, 3)
            Ly_default = 6
        nphi = get_parameter(params, 'nphi', nphi_default, self.name)
        Lx_MUC = get_parameter(params, 'Lx_MUC', 1, self.name)
        Lx = Lx_MUC * nphi[1]
        Ly = get_parameter(params, 'Ly', Ly_default, self.name)
        site = self.init_sites(params)
        order = get_parameter(params, 'order', 'Cstyle', self.name)
        bc_MPS = get_parameter(params, 'bc_MPS', 'infinite', self.name)
        bc_x = 'periodic' if bc_MPS == 'infinite' else 'open'  # Next line needs default
        bc_x = get_parameter(params, 'bc_x', bc_x, self.name)
        bc_y = get_parameter(params, 'bc_y', 'cylinder', self.name)
        assert bc_y in ['cylinder', 'ladder']
        bc_y = 'periodic' if bc_y == 'cylinder' else 'open'
        if bc_MPS == 'infinite' and bc_x == 'open':
            raise ValueError("You need to use 'periodic' `bc_x` for infinite systems!")
        lat = Honeycomb(Lx, Ly, site, order=order, bc_MPS=bc_MPS, bc=[bc_x, bc_y])
        return lat

    def init_terms(self, params):
        if self.stats(params) == 'bosons':
            nphi_default = (1, 4)
            creation = 'Bd'
            annihilation = 'B'
        else:
            nphi_default = (1, 3)
            V = get_parameter(params, 'V', 10, self.name, True)
            creation = 'Cd'
            annihilation = 'C'
        t1 = get_parameter(params, 't1', 1, self.name)
        t2 = get_parameter(params, 't2', 0, self.name)
        mu = get_parameter(params, 'mu', 0., self.name)
        phi_ext_2pi = 2*np.pi*get_parameter(params, 'phi_ext', 0., self.name)
        nphi = get_parameter(params, 'nphi', nphi_default, self.name)
        Lx_MUC = get_parameter(params, 'Lx_MUC', 1, self.name)
        nphi_2pi = 2 * np.pi * nphi[0] / nphi[1]

        for u in range(len(self.lat.unit_cell)):
            self.add_onsite(mu, 0, 'N')

        # 1st-NN term ##################################################################################################

        u1, u2, dx = (0, 1, np.array([0, -1]))  # down
        m = np.arange(0, 2 * nphi[1] * Lx_MUC, 2)
        t_phi = self.coupling_strength_add_ext_flux(t1, dx, [0, phi_ext_2pi]) \
                * np.exp(-1j * (nphi_2pi / 3) * m)[:, np.newaxis]
        self.add_coupling(t_phi, u1, creation, u2, annihilation, dx)
        self.add_coupling(np.conj(t_phi), u2, creation, u1, annihilation, -dx)  # H.c.
        if self.stats(params) == 'fermions':
            self.add_coupling(V, u1, 'N', u2, 'N', dx)

        u1, u2, dx = (0, 1, np.array([-1, 0]))  # upper left
        m = np.roll(np.arange(0, 2 * nphi[1] * Lx_MUC, 2), -1)  # match convention for strength argument of add_coupling
        t_phi = self.coupling_strength_add_ext_flux(t1, dx, [0, phi_ext_2pi]) \
                * np.exp(1j * (nphi_2pi / 6) * (m - 1 / 2))[:, np.newaxis]
        self.add_coupling(t_phi, u1, creation, u2, annihilation, dx)
        self.add_coupling(np.conj(t_phi), u2, creation, u1, annihilation, -dx)  # H.c.
        if self.stats(params) == 'fermions':
            self.add_coupling(V, u1, 'N', u2, 'N', dx)

        u1, u2, dx = (0, 1, np.array([0, 0]))  # upper right
        m = np.arange(0, 2 * nphi[1] * Lx_MUC, 2)
        t_phi = self.coupling_strength_add_ext_flux(t1, dx, [0, phi_ext_2pi]) \
                * np.exp(1j * (nphi_2pi / 6) * (m + 1 / 2))[:, np.newaxis]
        self.add_coupling(t_phi, u1, creation, u2, annihilation, dx)
        self.add_coupling(np.conj(t_phi), u2, creation, u1, annihilation, -dx)  # H.c.
        if self.stats(params) == 'fermions':
            self.add_coupling(V, u1, 'N', u2, 'N', dx)

        # 5th-NN term ##################################################################################################

        for sublattice in [0, 1]:  # A and B sublattices

            u1, u2, dx = (sublattice, sublattice, np.array([-1, 2]))  # up
            m = np.roll(np.arange(sublattice, 2 * nphi[1] * Lx_MUC, 2),
                        -1)  # match convention for strength argument of add_coupling
            t_phi = self.coupling_strength_add_ext_flux(t2, dx, [0, phi_ext_2pi]) \
                    * np.exp(1j * nphi_2pi * m)[:, np.newaxis]
            self.add_coupling(t_phi, u1, creation, u2, annihilation, dx)
            self.add_coupling(np.conj(t_phi), u2, creation, u1, annihilation, -dx)  # H.c.

            u1, u2, dx = (sublattice, sublattice, np.array([-2, 1]))  # bottom right
            m = np.roll(np.arange(sublattice, 2 * nphi[1] * Lx_MUC, 2),
                        -2)  # match convention for strength argument of add_coupling
            t_phi = self.coupling_strength_add_ext_flux(t2, dx, [0, phi_ext_2pi]) \
                    * np.exp(-1j * (nphi_2pi / 2) * (m + 3 / 2))[:, np.newaxis]
            self.add_coupling(t_phi, u1, creation, u2, annihilation, dx)
            self.add_coupling(np.conj(t_phi), u2, creation, u1, annihilation, -dx)  # H.c.

            u1, u2, dx = (sublattice, sublattice, np.array([-1, -1]))  # bottom left
            m = np.roll(np.arange(sublattice, 2 * nphi[1] * Lx_MUC, 2),
                        -1)  # match convention for strength argument of add_coupling
            t_phi = self.coupling_strength_add_ext_flux(t2, dx, [0, phi_ext_2pi]) \
                    * np.exp(-1j * (nphi_2pi / 2) * (m - 3 / 2))[:, np.newaxis]
            self.add_coupling(t_phi, u1, creation, u2, annihilation, dx)
            self.add_coupling(np.conj(t_phi), u2, creation, u1, annihilation, -dx)  # H.c.


# class FerHofHex1Hex5Model(CouplingMPOModel):
#
#     def __init__(self, params):
#         CouplingMPOModel.__init__(self, params)
#
#     def init_sites(self, params):
#         conserve = get_parameter(params, 'conserve', 'N', self.name)
#         filling = get_parameter(params, 'filling', (1, 9), self.name)
#         filling = filling[0] / filling[1]
#         site = FermionSite(conserve=conserve, filling=filling)
#         return site
#
#     def init_lattice(self, params):
#         bc_MPS = get_parameter(params, 'bc_MPS', 'infinite', self.name)
#         order = get_parameter(params, 'order', 'default', self.name)
#         site = self.init_sites(params)
#         qvalue = get_parameter(params, 'phi', (1, 3), self.name)[1]
#         Lx_MUC = get_parameter(params, 'Lx_MUC', 1, self.name)
#         Lx = Lx_MUC * qvalue
#         Ly = get_parameter(params, 'Ly', 6, self.name)
#         bc_x = 'periodic' if bc_MPS == 'infinite' else 'open'  # Next line needs default
#         bc_x = get_parameter(params, 'bc_x', bc_x, self.name)
#         bc_y = get_parameter(params, 'bc_y', 'cylinder', self.name)
#         assert bc_y in ['cylinder', 'ladder']
#         bc_y = 'periodic' if bc_y == 'cylinder' else 'open'
#         if bc_MPS == 'infinite' and bc_x == 'open':
#             raise ValueError("You need to use 'periodic' `bc_x` for infinite systems!")
#         lat = Honeycomb(Lx, Ly, site, order=order, bc=[bc_x, bc_y], bc_MPS=bc_MPS)
#         return lat
#
#     def init_terms(self, params):
#         t1 = get_parameter(params, 't1', 1, self.name)
#         t2 = get_parameter(params, 't2', 0, self.name)
#         # t1, t2 = 1, 0
#         V = get_parameter(params, 'V', 10, self.name, True)
#         phi_ext = 2*np.pi*get_parameter(params, 'phi_ext', 0., self.name)
#         phi_p, phi_q = get_parameter(params, 'phi', (1, 3), self.name)
#         Lx_MUC = get_parameter(params, 'Lx_MUC', 1, self.name)
#         phi = 2 * np.pi * phi_p / phi_q
#
#         # t1 term ######################################################################################################
#
#         u1, u2, dx = (0, 1, np.array([0, -1]))  # down
#         m = np.arange(0, 2 * phi_q * Lx_MUC, 2)
#         t_phi = self.coupling_strength_add_ext_flux(t1, dx, [0, phi_ext]) \
#                 * np.exp(-1j * (phi / 3) * m)[:, np.newaxis]
#         self.add_coupling(t_phi, u1, 'Cd', u2, 'C', dx)
#         self.add_coupling(np.conj(t_phi), u2, 'Cd', u1, 'C', -dx)  # h.c.
#         # NN interaction
#         self.add_coupling(V, u1, 'N', u2, 'N', dx)
#
#         u1, u2, dx = (0, 1, np.array([-1, 0]))  # upper left
#         m = np.roll(np.arange(0, 2 * phi_q * Lx_MUC, 2), -1)  # match convention for strength argument of add_coupling
#         t_phi = self.coupling_strength_add_ext_flux(t1, dx, [0, phi_ext]) \
#                 * np.exp(1j * (phi / 6) * (m - 1 / 2))[:, np.newaxis]
#         self.add_coupling(t_phi, u1, 'Cd', u2, 'C', dx)
#         self.add_coupling(np.conj(t_phi), u2, 'Cd', u1, 'C', -dx)  # h.c.
#         # NN interaction
#         self.add_coupling(V, u1, 'N', u2, 'N', dx)
#
#         u1, u2, dx = (0, 1, np.array([0, 0]))  # upper right
#         m = np.arange(0, 2 * phi_q * Lx_MUC, 2)
#         t_phi = self.coupling_strength_add_ext_flux(t1, dx, [0, phi_ext]) \
#                 * np.exp(1j * (phi / 6) * (m + 1 / 2))[:, np.newaxis]
#         self.add_coupling(t_phi, u1, 'Cd', u2, 'C', dx)
#         self.add_coupling(np.conj(t_phi), u2, 'Cd', u1, 'C', -dx)  # h.c.
#         # NN interaction
#         self.add_coupling(V, u1, 'N', u2, 'N', dx)
#
#         # t2 term ######################################################################################################
#
#         for sublattice in [0, 1]:  # A and B sublattices
#
#             u1, u2, dx = (sublattice, sublattice, np.array([-1, 2]))  # up
#             m = np.roll(np.arange(sublattice, 2 * phi_q * Lx_MUC, 2), -1)  # match convention for strength argument of add_coupling
#             t_phi = self.coupling_strength_add_ext_flux(t2, dx, [0, phi_ext]) \
#                     * np.exp(1j * phi * m)[:, np.newaxis]
#             self.add_coupling(t_phi, u1, 'Cd', u2, 'C', dx)
#             self.add_coupling(np.conj(t_phi), u2, 'Cd', u1, 'C', -dx)  # h.c.
#
#             u1, u2, dx = (sublattice, sublattice, np.array([-2, 1]))  # bottom right
#             m = np.roll(np.arange(sublattice, 2 * phi_q * Lx_MUC, 2), -2)  # match convention for strength argument of add_coupling
#             t_phi = self.coupling_strength_add_ext_flux(t2, dx, [0, phi_ext]) \
#                     * np.exp(-1j * (phi / 2) * (m + 3 / 2))[:, np.newaxis]
#             self.add_coupling(t_phi, u1, 'Cd', u2, 'C', dx)
#             self.add_coupling(np.conj(t_phi), u2, 'Cd', u1, 'C', -dx)  # h.c.
#
#             u1, u2, dx = (sublattice, sublattice, np.array([-1, -1]))  # bottom left
#             m = np.roll(np.arange(sublattice, 2 * phi_q * Lx_MUC, 2), -1)  # match convention for strength argument of add_coupling
#             t_phi = self.coupling_strength_add_ext_flux(t2, dx, [0, phi_ext]) \
#                     * np.exp(-1j * (phi / 2) * (m - 3 / 2))[:, np.newaxis]
#             self.add_coupling(t_phi, u1, 'Cd', u2, 'C', dx)
#             self.add_coupling(np.conj(t_phi), u2, 'Cd', u1, 'C', -dx)  # h.c.


if __name__ == "__main__":

    model_params = dict(conserve='N', t1=1, t2=-0.025, n=(int(1), int(9)), nphi=(int(1), int(3)),
                        Lx_MUC=1, Ly=6, V=10,
                        bc_MPS='infinite', bc_x='periodic', bc_y='cylinder', order='Cstyle',
                        verbose=1, phi_ext=1.)
    M = HofHex1Hex5Model(model_params)

    t0 = time.time()
    print("max MPO bond dimension = ", max(M.H_MPO.chi))
    print("Total time taken (seconds) = ", time.time() - t0)
