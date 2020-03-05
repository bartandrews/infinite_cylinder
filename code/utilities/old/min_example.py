import numpy as np
import sys
import time

from tenpy.models.model import CouplingMPOModel
from tenpy.tools.params import get_parameter
from tenpy.networks.site import BosonSite, FermionSite, GroupedSite
from tenpy.models.lattice import Honeycomb
from tenpy.networks.mps import MPS
from tenpy.algorithms import dmrg


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


if __name__ == "__main__":

        model_params = dict(conserve='N', t1=1, t2=0, t2dash=0, filling=(int(1), int(9)), phi=(int(1), int(3)),
                            Lx_MUC=1, Ly=6, U=100, V=10,  # system params
                            bc_MPS='infinite', bc_x='periodic', bc_y='cylinder', order='default',  # MPS params
                            verbose=1, phi_ext=0)  # utility
        M = FermionicHex1Hex5OrbitalModel(model_params)

        # initial_state that works
        # product_state = ['full_x empty_y', 0, 0, 0, 0, 0, 0, 0, 0, 'full_x empty_y', 0, 0, 0, 0, 0, 0, 0, 0,
        #                  'full_x empty_y', 0, 0, 0, 0, 0, 0, 0, 0, 'full_x empty_y', 0, 0, 0, 0, 0, 0, 0, 0]

        product_state = ['full_x full_y', 0, 0, 0, 0, 0, 0, 0, 0, 'full_x full_y', 0, 0, 0, 0, 0, 0, 0, 0,
                         'full_x full_y', 0, 0, 0, 0, 0, 0, 0, 0, 'full_x full_y', 0, 0, 0, 0, 0, 0, 0, 0]

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
            'diag_method': 'lanczos'
        }

        # engine = dmrg.OneSiteDMRGEngine(psi, M, OneSiteH, dmrg_params)
        engine = dmrg.TwoSiteDMRGEngine(psi, M, dmrg_params)
        engine.run()