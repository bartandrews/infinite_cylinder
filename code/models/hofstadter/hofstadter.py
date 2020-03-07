# --- python imports
import numpy as np
# --- TeNPy imports
from tenpy.models.model import CouplingMPOModel
from tenpy.tools.params import get_parameter
from tenpy.networks.site import BosonSite, FermionSite


class HofstadterModel(CouplingMPOModel):

    def __init__(self, params):
        CouplingMPOModel.__init__(self, params)

    def stats(self, params):
        return get_parameter(params, 'statistics', 'bosons', self.name)

    def init_single_sites(self, params):
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

    def init_lattice_params(self, params):
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
        bc = [bc_x, bc_y]
        return Lx, Ly, site, order, bc_MPS, bc

    def init_basic_terms(self, params):
        if self.stats(params) == 'bosons':
            nphi_default = (1, 4)
            creation = 'Bd'
            annihilation = 'B'
        else:
            nphi_default = (1, 3)
            creation = 'Cd'
            annihilation = 'C'
        t1 = get_parameter(params, 't1', 1, self.name, True)
        mu = get_parameter(params, 'mu', 0., self.name)
        phi_ext_2pi = 2 * np.pi * get_parameter(params, 'phi_ext', 0., self.name)
        nphi = get_parameter(params, 'nphi', nphi_default, self.name)
        Lx_MUC = get_parameter(params, 'Lx_MUC', 1, self.name)
        nphi_2pi = 2 * np.pi * nphi[0] / nphi[1]

        return nphi_default, creation, annihilation, t1, mu, phi_ext_2pi, nphi, nphi_2pi, Lx_MUC
