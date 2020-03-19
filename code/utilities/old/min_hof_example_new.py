import numpy as np
import sys
import time

from tenpy.models.model import CouplingMPOModel
from tenpy.tools.params import get_parameter
from tenpy.networks.site import BosonSite, FermionSite
from tenpy.models.lattice import Square
from tenpy.networks.mps import MPS
from tenpy.algorithms import dmrg


class HofstadterModel(CouplingMPOModel):

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
        t1 = get_parameter(params, 't1', 1, self.name, True)
        mu = get_parameter(params, 'mu', 0., self.name)
        V = get_parameter(params, 'V', V_default, self.name, True)
        Vtype = get_parameter(params, 'Vtype', 'Coulomb', self.name)
        Vrange = get_parameter(params, 'Vrange', Vrange_default, self.name)
        nphi = get_parameter(params, 'nphi', nphi_default, self.name)
        nphi_2pi = 2 * np.pi * nphi[0] / nphi[1]
        LxMUC = get_parameter(params, 'LxMUC', 1, self.name)
        phi_2pi = 2 * np.pi * get_parameter(params, 'phi', 0., self.name)

        return creation, annihilation, nphi_default, t1, mu, V, Vtype, Vrange, nphi, nphi_2pi, LxMUC, phi_2pi

    def squ_1_hoppings(self, creation, annihilation, t, nphi, nphi_2pi, LxMUC, phi_2pi):
        u1, u2, dx = (0, 0, np.array([1, 0]))  # right
        t_phi = self.coupling_strength_add_ext_flux(t, dx, [0, phi_2pi])

        print("(right) t_phi = ", t_phi)

        self.add_coupling(t_phi, u1, creation, u2, annihilation, dx)
        self.add_coupling(np.conj(t_phi), u2, creation, u1, annihilation, -dx)  # H.c.

        u1, u2, dx = (0, 0, np.array([0, 1]))  # up
        m = np.arange(0, nphi[1] * LxMUC)
        t_phi = self.coupling_strength_add_ext_flux(t, dx, [0, phi_2pi]) \
                * np.exp(-1j * nphi_2pi * m)[:, np.newaxis]

        print("(up) t_phi = ", t_phi)

        self.add_coupling(t_phi, u1, creation, u2, annihilation, dx)
        self.add_coupling(np.conj(t_phi), u2, creation, u1, annihilation, -dx)  # H.c.


class HofSqu1Model(HofstadterModel):

    def __init__(self, params):
        HofstadterModel.__init__(self, params)

    def init_sites(self, params):
        site = HofstadterModel.init_sites(self, params)

        print(site)

        return site

    def init_lattice(self, params):
        (Lx, Ly, order, bc_MPS, bc) = self.init_lattice_params(params)
        site = self.init_sites(params)
        lat = Square(Lx, Ly, site, order=order, bc_MPS=bc_MPS, bc=bc)

        print(Lx, Ly, order, bc_MPS, bc)

        return lat

    def init_terms(self, params):
        (creation, annihilation, nphi_default, t1, mu, V, Vtype, Vrange, nphi, nphi_2pi, LxMUC, phi_2pi) = \
            HofstadterModel.init_terms(self, params)

        print(creation, annihilation, t1, nphi, nphi_2pi, LxMUC, phi_2pi)

        self.squ_1_hoppings(creation, annihilation, t1, nphi, nphi_2pi, LxMUC, phi_2pi)


if __name__ == "__main__":

        model_params = dict(statistics='bosons', conserve='N', t1=1, n=(int(1), int(10)), nphi=(int(1), int(5)),
                            LxMUC=1, Ly=4, V=0, Vtype='Coulomb', Vrange=0,
                            bc_MPS='infinite', bc_x='periodic', bc_y='cylinder', order='Cstyle',
                            verbose=1, phi=0)
        M = HofSqu1Model(model_params)

        import pprint
        pprint.pprint(list(M.coupling_terms['Bd_i B_j'].to_TermList()))

        import matplotlib.pyplot as plt

        product_state = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0]

        psi = MPS.from_product_state(M.lat.mps_sites(), product_state, bc=M.lat.bc_MPS)

        dmrg_params = {
            'mixer': True,  # setting this to True helps to escape local minima
            'mixer_params': {'amplitude': 1.e-5, 'decay': 1.2, 'disable_after': 30},
            'trunc_params': {
                # 'chi_max': chi_max,
                'svd_min': 1.e-10
            },
            # 'lanczos_params': {
            #     'reortho': True,
            #     'N_cache': 40
            # },
            'chi_list': {0: 9, 10: 49, 20: 100, 40: 50},  # may need to change first value ~39
            'max_E_err': 1.e-6,
            'max_S_err': 1.e-6,
            # 'norm_tol': 1.e-6,
            # 'norm_tol_iter': 1000,
            'max_sweeps': 1000,
            'verbose': 5,
            'N_sweeps_check': 10,
            # 'diag_method': 'lanczos'
        }

        # # engine = dmrg.OneSiteDMRGEngine(psi, M, OneSiteH, dmrg_params)
        # engine = dmrg.TwoSiteDMRGEngine(psi, M, dmrg_params)

        info = dmrg.run(psi, M, dmrg_params)
        E = info['E']

        print(psi.entanglement_entropy()[0])
