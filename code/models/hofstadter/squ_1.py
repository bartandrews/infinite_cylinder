# --- python imports
import time
import numpy as np
# --- TeNPy imports
from tenpy.tools.params import get_parameter
from tenpy.models.lattice import Square
# --- infinite_cylinder imports
from models.hofstadter.hofstadter import HofstadterModel


class HofSqu1Model(HofstadterModel):

    def __init__(self, params):
        super().__init__(params)

    def init_sites(self, params):
        site = self.init_single_sites(params)
        return site

    def init_lattice(self, params):
        (Lx, Ly, site, order, bc_MPS, bc) = self.init_lattice_params(params)
        lat = Square(Lx, Ly, site, order=order, bc_MPS=bc_MPS, bc=bc)
        return lat

    def init_terms(self, params):

        (nphi_default, creation, annihilation, t1, mu, phi_ext_2pi, nphi, nphi_2pi, Lx_MUC) = self.init_basic_terms(params)
        if self.stats(params) == 'fermions':
            V = get_parameter(params, 'V', 10, self.name, True)

        for u in range(len(self.lat.unit_cell)):  # chemical potential
            self.add_onsite(mu, 0, 'N')

        u1, u2, dx = (0, 0, np.array([1, 0]))  # right
        t_phi = self.coupling_strength_add_ext_flux(t1, dx, [0, phi_ext_2pi])
        self.add_coupling(t_phi, u1, creation, u2, annihilation, dx)
        self.add_coupling(np.conj(t_phi), u2, creation, u1, annihilation, -dx)  # H.c.
        if self.stats(params) == 'fermions':
            self.add_coupling(V, u1, 'N', u2, 'N', dx)

        u1, u2, dx = (0, 0, np.array([0, 1]))  # up
        m = np.arange(0, nphi[1] * Lx_MUC)
        t_phi = self.coupling_strength_add_ext_flux(t1, dx, [0, phi_ext_2pi]) \
                * np.exp(-1j * nphi_2pi * m)[:, np.newaxis]
        self.add_coupling(t_phi, u1, creation, u2, annihilation, dx)
        self.add_coupling(np.conj(t_phi), u2, creation, u1, annihilation, -dx)  # H.c.
        if self.stats(params) == 'fermions':
            self.add_coupling(V, u1, 'N', u2, 'N', dx)


if __name__ == "__main__":

    model_params = dict(conserve='N', t1=1, n=(int(1), int(9)), nphi=(int(1), int(3)),
                        Lx_MUC=1, Ly=6, V=10,
                        bc_MPS='infinite', bc_x='periodic', bc_y='cylinder', order='Cstyle',
                        verbose=1, phi_ext=0)
    M = HofSqu1Model(model_params)

    t0 = time.time()
    print("max MPO bond dimension = ", max(M.H_MPO.chi))
    print("Total time taken (seconds) = ", time.time() - t0)
