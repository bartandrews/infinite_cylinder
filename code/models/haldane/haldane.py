# --- python imports
import numpy as np
# --- tenpy imports
from tenpy.models.model import CouplingMPOModel
from tenpy.tools.params import get_parameter
from tenpy.networks.site import BosonSite, FermionSite
from tenpy.models.lattice import Honeycomb


class HalModel(CouplingMPOModel):

    def __init__(self, params):
        CouplingMPOModel.__init__(self, params)

    def stats(self, params):
        return get_parameter(params, 'statistics', 'bosons', self.name)

    def init_sites(self, params):
        conserve = get_parameter(params, 'conserve', 'N', self.name)
        if self.stats(params) == 'bosons':
            Nmax = get_parameter(params, 'Nmax', 1, self.name)
            site = BosonSite(Nmax=Nmax, conserve=conserve)
        else:
            site = FermionSite(conserve=conserve)
        return site

    def init_lattice(self, params):
        site = self.init_sites(params)
        Lx = get_parameter(params, 'LxMUC', 1, self.name)
        Ly = get_parameter(params, 'Ly', 6, self.name)
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
        lat = Honeycomb(Lx, Ly, site, order=order, bc_MPS=bc_MPS, bc=bc)
        return lat

    def init_terms(self, params):
        if self.stats(params) == 'bosons':
            creation, annihilation = 'Bd', 'B'
            V_default, Vrange_default = 0, 0
        else:
            creation, annihilation = 'Cd', 'C'
            V_default, Vrange_default = 10, 1
        t1 = get_parameter(params, 't1', 1., self.name, True)
        t2 = get_parameter(params, 't2',
                           (np.sqrt(129) / 36) * t1 * np.exp(1j * np.arccos(3 * np.sqrt(3 / 43))),
                           self.name, True)
        mu = get_parameter(params, 'mu', 0., self.name, True)
        V = get_parameter(params, 'V', V_default, self.name, True)
        phi_2pi = 2 * np.pi * get_parameter(params, 'phi', 0., self.name)

        for u in range(len(self.lat.unit_cell)):
            self.add_onsite(mu, 0, 'N', category='mu N')
            self.add_onsite(-mu, 1, 'N', category='mu N')

        for u1, u2, dx in self.lat.pairs['nearest_neighbors']:
            t1_phi = self.coupling_strength_add_ext_flux(t1, dx, [0, phi_2pi])
            self.add_coupling(t1_phi, u1, creation, u2, annihilation, dx, category='t1 Ad_i A_j')
            self.add_coupling(np.conj(t1_phi), u2, creation, u1, annihilation, -dx, category='t1 Ad_i A_j H.c.')  # H.c.
            self.add_coupling(V, u1, 'N', u2, 'N', dx, category='V N_i N_j')

        for u1, u2, dx in [(0, 0, np.array([-1, 1])), (0, 0, np.array([1, 0])),
                           (0, 0, np.array([0, -1])), (1, 1, np.array([0, 1])),
                           (1, 1, np.array([1, -1])), (1, 1, np.array([-1, 0]))]:
            t2_phi = self.coupling_strength_add_ext_flux(t2, dx, [0, phi_2pi])
            self.add_coupling(t2_phi, u1, creation, u2, annihilation, dx, category='t2 Ad_i A_j')
            self.add_coupling(np.conj(t2_phi), u2, creation, u1, annihilation, -dx, category='t2 Ad_i A_j H.c.')  # H.c.
