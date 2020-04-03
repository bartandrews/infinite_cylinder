import numpy as np
import sys
import time

from tenpy.models.hofstadter import HofstadterBosons
from tenpy.tools.params import get_parameter
from tenpy.networks.site import BosonSite, FermionSite
from tenpy.models.lattice import Square
from tenpy.networks.mps import MPS
from tenpy.algorithms import dmrg


if __name__ == "__main__":

        model_params = dict(Lx=4, Ly=6,
                            # mx=,2 my=1,
                            filling=(1, 8), phi=(1, 4),
                            Jx=1, Jy=1, mu=0, U=0, Nmax=1,
                            bc_MPS='infinite', bc_x='periodic', bc_y='cylinder',
                            conserve='N', order='Cstyle', gauge='landau_x')
        M = HofstadterBosons(model_params)

        # import pprint
        # pprint.pprint(list(M.coupling_terms['Bd_i B_j'].to_TermList()))

        import matplotlib.pyplot as plt

        product_state = [1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0]

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
            'chi_list': {0: 9, 10: 49, 20: 100, 40: 250},  # may need to change first value ~39
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
