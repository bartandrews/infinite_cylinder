# --- python imports
import numpy as np
# --- TeNPy imports
from tenpy.models.model import CouplingMPOModel
from tenpy.models.model import MultiCouplingModel
from tenpy.tools.params import get_parameter
from tenpy.networks.site import BosonSite, FermionSite
# --- infinite_cylinder imports
import functions.func_int as fi


class HofstadterModel(CouplingMPOModel, MultiCouplingModel):

    def __init__(self, params):
        CouplingMPOModel.__init__(self, params)

    def stats(self, params):
        return get_parameter(params, 'statistics', 'bosons', self.name)

    def init_sites(self, params):
        conserve = get_parameter(params, 'conserve', 'N', self.name)
        if self.stats(params) == 'bosons':
            Nmax = get_parameter(params, 'Nmax', 1, self.name)
            n = get_parameter(params, 'n', (1, 8), self.name)
            n = n[0] / n[1]
            site = BosonSite(Nmax=Nmax, conserve=conserve, filling=n)
        else:
            n = get_parameter(params, 'n', (1, 9), self.name)
            n = n[0] / n[1]
            site = FermionSite(conserve=conserve, filling=n)
        return site

    def init_lattice_params(self, params):
        if self.stats(params) == 'bosons':
            nphi_default = (1, 4)
            Ly_default = 4
        else:
            nphi_default = (1, 3)
            Ly_default = 6
        nphi = get_parameter(params, 'nphi', nphi_default, self.name)
        LxMUC = get_parameter(params, 'LxMUC', 1, self.name)
        Lx = LxMUC * nphi[1]
        Ly = get_parameter(params, 'Ly', Ly_default, self.name)
        order = get_parameter(params, 'order', 'Cstyle', self.name)
        bc_MPS = get_parameter(params, 'bc_MPS', 'infinite', self.name)
        bc_x = 'periodic' if bc_MPS == 'infinite' else 'open'  # Next line needs default
        bc_x = get_parameter(params, 'bc_x', bc_x, self.name)
        bc_y = get_parameter(params, 'bc_y', 'cylinder', self.name)
        assert bc_y in ['cylinder', 'ladder']
        bc_y = 'periodic' if bc_y == 'cylinder' else 'open'
        if bc_MPS == 'infinite' and bc_x == 'open':
            raise ValueError("You need to use 'periodic' `bc_x` for infinite systems!")
        bc = [bc_x, bc_y]
        return Lx, Ly, order, bc_MPS, bc

    def init_terms(self, params):
        if self.stats(params) == 'bosons':
            creation, annihilation = 'Bd', 'B'
            V_default, Vrange_default = 0, 0
            nphi_default = (1, 4)
        else:
            creation, annihilation = 'Cd', 'C'
            V_default, Vrange_default = 10, 1
            nphi_default = (1, 3)
        Nmax = get_parameter(params, 'Nmax', 1, self.name, True)
        t1 = get_parameter(params, 't1', 1, self.name, True)
        mu = get_parameter(params, 'mu', 0, self.name)
        V = get_parameter(params, 'V', V_default, self.name, True)
        Vtype = get_parameter(params, 'Vtype', 'Coulomb', self.name)
        Vrange = get_parameter(params, 'Vrange', Vrange_default, self.name)
        nphi = get_parameter(params, 'nphi', nphi_default, self.name)
        nphi_2pi = 2 * np.pi * nphi[0] / nphi[1]
        LxMUC = get_parameter(params, 'LxMUC', 1, self.name)
        phi_2pi = 2 * np.pi * get_parameter(params, 'phi', 0, self.name)

        return creation, annihilation, nphi_default, Nmax, t1, mu, V, Vtype, Vrange, nphi, nphi_2pi, LxMUC, phi_2pi

    def chemical_potential(self, mu, extra_dof=False):
        tot_numb_op = 'N' if not extra_dof else 'Ntot'
        for u in range(len(self.lat.unit_cell)):
            self.add_onsite(mu, 0, tot_numb_op)

    def onsite_interaction(self, U, op1, op2):
        op = op1+' '+op2
        for u in range(len(self.lat.unit_cell)):
            print("u in range(len(self.lat.unit_cell)) = ", u)
            self.add_onsite(U, u, op)

    def offsite_interaction(self, lattice, Nmax, V, Vtype, Vrange, extra_dof=False):
        tot_numb_op = 'N' if not extra_dof else 'Ntot'
        for i in range(1, 11):  # offsite interaction only implemented up to 10th-NN
            if Vrange >= i:
                if Nmax == 1:
                    for u1, u2, dx in fi.NN(lattice, i):
                        self.add_coupling(fi.interaction_strength(lattice, V, Vtype, i-1),
                                          u1, tot_numb_op, u2, tot_numb_op, dx)
                elif Nmax == 2:
                    if lattice != "Squ" or Vrange > 1:
                        raise ValueError("3-body interaction is only implemented for Squ1.")
                    # upper-right corner
                    self.add_multi_coupling(fi.interaction_strength(lattice, V, Vtype, i - 1), 0, tot_numb_op,
                                            [(0, tot_numb_op, [0, 1]), (0, tot_numb_op, [1, 0])])
                    # upper-left corner
                    self.add_multi_coupling(fi.interaction_strength(lattice, V, Vtype, i - 1), 0, tot_numb_op,
                                            [(0, tot_numb_op, [0, 1]), (0, tot_numb_op, [-1, 0])])
                    # lower-right corner
                    self.add_multi_coupling(fi.interaction_strength(lattice, V, Vtype, i - 1), 0, tot_numb_op,
                                            [(0, tot_numb_op, [0, -1]), (0, tot_numb_op, [1, 0])])
                    # lower-left corner
                    self.add_multi_coupling(fi.interaction_strength(lattice, V, Vtype, i - 1), 0, tot_numb_op,
                                            [(0, tot_numb_op, [0, -1]), (0, tot_numb_op, [-1, 0])])
                    # horizontal line
                    self.add_multi_coupling(fi.interaction_strength(lattice, V, Vtype, i - 1), 0, tot_numb_op,
                                            [(0, tot_numb_op, [-1, 0]), (0, tot_numb_op, [1, 0])])
                    # vertical line
                    self.add_multi_coupling(fi.interaction_strength(lattice, V, Vtype, i - 1), 0, tot_numb_op,
                                            [(0, tot_numb_op, [0, -1]), (0, tot_numb_op, [0, 1])])
                else:
                    raise ValueError("N-body interaction is not implemented for N>3.")

    def squ_1_hoppings(self, creation, annihilation, t, nphi, nphi_2pi, LxMUC, phi_2pi):
        u1, u2, dx = (0, 0, np.array([1, 0]))  # right
        t_phi = self.coupling_strength_add_ext_flux(t, dx, [0, phi_2pi])
        self.add_coupling(t_phi, u1, creation, u2, annihilation, dx)
        self.add_coupling(np.conj(t_phi), u2, creation, u1, annihilation, -dx)  # H.c.

        u1, u2, dx = (0, 0, np.array([0, 1]))  # up
        m = np.arange(0, nphi[1] * LxMUC)
        t_phi = self.coupling_strength_add_ext_flux(t, dx, [0, phi_2pi]) \
                * np.exp(-1j * nphi_2pi * m)[:, np.newaxis]
        self.add_coupling(t_phi, u1, creation, u2, annihilation, dx)
        self.add_coupling(np.conj(t_phi), u2, creation, u1, annihilation, -dx)  # H.c.

    def hex_1_hoppings(self, creation, annihilation, t, nphi, nphi_2pi, LxMUC, phi_2pi):
        u1, u2, dx = (0, 1, np.array([0, -1]))  # down
        m = np.arange(0, 2 * nphi[1] * LxMUC, 2)
        t_phi = self.coupling_strength_add_ext_flux(t, dx, [0, phi_2pi]) \
                * np.exp(-1j * (nphi_2pi / 3) * m)[:, np.newaxis]
        self.add_coupling(t_phi, u1, creation, u2, annihilation, dx)
        self.add_coupling(np.conj(t_phi), u2, creation, u1, annihilation, -dx)  # H.c.

        u1, u2, dx = (0, 1, np.array([-1, 0]))  # upper left
        m = np.roll(np.arange(0, 2 * nphi[1] * LxMUC, 2), -1)  # match convention for strength argument of add_coupling
        t_phi = self.coupling_strength_add_ext_flux(t, dx, [0, phi_2pi]) \
                * np.exp(1j * (nphi_2pi / 6) * (m - 1 / 2))[:, np.newaxis]
        self.add_coupling(t_phi, u1, creation, u2, annihilation, dx)
        self.add_coupling(np.conj(t_phi), u2, creation, u1, annihilation, -dx)  # H.c.

        u1, u2, dx = (0, 1, np.array([0, 0]))  # upper right
        m = np.arange(0, 2 * nphi[1] * LxMUC, 2)
        t_phi = self.coupling_strength_add_ext_flux(t, dx, [0, phi_2pi]) \
                * np.exp(1j * (nphi_2pi / 6) * (m + 1 / 2))[:, np.newaxis]
        self.add_coupling(t_phi, u1, creation, u2, annihilation, dx)
        self.add_coupling(np.conj(t_phi), u2, creation, u1, annihilation, -dx)  # H.c.

    def hex_5_hoppings(self, creation, annihilation, t, nphi, nphi_2pi, LxMUC, phi_2pi):
        for sublattice in [0, 1]:  # A and B sublattices

            u1, u2, dx = (sublattice, sublattice, np.array([-1, 2]))  # up
            m = np.roll(np.arange(sublattice, 2 * nphi[1] * LxMUC, 2),
                        -1)  # match convention for strength argument of add_coupling
            t_phi = self.coupling_strength_add_ext_flux(t, dx, [0, phi_2pi]) \
                    * np.exp(1j * nphi_2pi * m)[:, np.newaxis]
            self.add_coupling(t_phi, u1, creation, u2, annihilation, dx)
            self.add_coupling(np.conj(t_phi), u2, creation, u1, annihilation, -dx)  # H.c.

            u1, u2, dx = (sublattice, sublattice, np.array([-2, 1]))  # bottom right
            m = np.roll(np.arange(sublattice, 2 * nphi[1] * LxMUC, 2),
                        -2)  # match convention for strength argument of add_coupling
            t_phi = self.coupling_strength_add_ext_flux(t, dx, [0, phi_2pi]) \
                    * np.exp(-1j * (nphi_2pi / 2) * (m + 3 / 2))[:, np.newaxis]
            self.add_coupling(t_phi, u1, creation, u2, annihilation, dx)
            self.add_coupling(np.conj(t_phi), u2, creation, u1, annihilation, -dx)  # H.c.

            u1, u2, dx = (sublattice, sublattice, np.array([-1, -1]))  # bottom left
            m = np.roll(np.arange(sublattice, 2 * nphi[1] * LxMUC, 2),
                        -1)  # match convention for strength argument of add_coupling
            t_phi = self.coupling_strength_add_ext_flux(t, dx, [0, phi_2pi]) \
                    * np.exp(-1j * (nphi_2pi / 2) * (m - 3 / 2))[:, np.newaxis]
            self.add_coupling(t_phi, u1, creation, u2, annihilation, dx)
            self.add_coupling(np.conj(t_phi), u2, creation, u1, annihilation, -dx)  # H.c.
