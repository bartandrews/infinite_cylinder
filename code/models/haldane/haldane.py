# --- python imports
import numpy as np
# --- TeNPy imports
from tenpy.models.model import CouplingMPOModel
from tenpy.networks.site import BosonSite, FermionSite
# --- infinite_cylinder imports
import functions.func_int as fi


def stats(params):
    return params.get('statistics', 'bosons')


class HaldaneModel(CouplingMPOModel):

    def __init__(self, params):
        CouplingMPOModel.__init__(self, params)

    def init_sites(self, params):
        conserve = params.get('conserve', 'N')
        if stats(params) == 'bosons':
            Nmax = params.get('Nmax', 1)
            n = params.get('n', (1, 2))
            n = n[0] / n[1]
            site = BosonSite(Nmax=Nmax, conserve=conserve, filling=n)
        else:
            n = params.get('n', (1, 3))
            n = n[0] / n[1]
            site = FermionSite(conserve=conserve, filling=n)
        return site

    @staticmethod
    def init_lattice_params(params):
        Lx = params.get('LxMUC', 1)
        Ly = params.get('Ly', 3)
        order = params.get('order', 'Cstyle')
        bc_MPS = params.get('bc_MPS', 'infinite')
        bc_x = 'periodic' if bc_MPS == 'infinite' else 'open'  # Next line needs default
        bc_x = params.get('bc_x', bc_x)
        bc_y = params.get('bc_y', 'cylinder')
        assert bc_y in ['cylinder', 'ladder']
        bc_y = 'periodic' if bc_y == 'cylinder' else 'open'
        if bc_MPS == 'infinite' and bc_x == 'open':
            raise ValueError("You need to use 'periodic' `bc_x` for infinite systems!")
        bc = [bc_x, bc_y]
        return Lx, Ly, order, bc_MPS, bc

    def init_terms(self, params):
        if stats(params) == 'bosons':
            creation, annihilation = 'Bd', 'B'
            V_default, Vrange_default = 0, 0
        else:
            creation, annihilation = 'Cd', 'C'
            V_default, Vrange_default = 0, 0
        Nmax = params.get('Nmax', 1)
        t1 = params.get('t1', 1)
        t2 = params.get('t2', None)
        t3 = params.get('t3', None)
        mu = params.get('mu', 0)
        V = params.get('V', V_default)
        Vtype = params.get('Vtype', 'Coulomb')
        Vrange = params.get('Vrange', Vrange_default)
        phi_2pi = 2 * np.pi * params.get('phi', 0)
        return creation, annihilation, Nmax, t1, t2, t3, mu, V, Vtype, Vrange, phi_2pi

    def chemical_potential(self, mu, extra_dof=False):
        tot_numb_op = 'N' if not extra_dof else 'Ntot'
        for u in range(len(self.lat.unit_cell)):
            self.add_onsite(mu, 0, tot_numb_op)
            self.add_onsite(-mu, 1, tot_numb_op)

    def onsite_interaction(self, U, op1, op2):
        op = op1 + ' ' + op2
        for u in range(len(self.lat.unit_cell)):
            print("u in range(len(self.lat.unit_cell)) = ", u)
            self.add_onsite(U, u, op)

    def offsite_interaction(self, lattice, V, Vtype, Vrange, extra_dof=False):
        tot_numb_op = 'N' if not extra_dof else 'Ntot'
        for i in range(1, 11):  # offsite interaction only implemented up to 10th-NN
            if Vrange >= i:
                for u1, u2, dx in fi.NN(lattice, i):
                    self.add_coupling(fi.interaction_strength(lattice, V, Vtype, i - 1),
                                      u1, tot_numb_op, u2, tot_numb_op, dx)
