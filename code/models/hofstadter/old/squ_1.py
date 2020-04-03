# --- python imports
import time
import numpy as np
# --- TeNPy imports
from tenpy.models.model import CouplingMPOModel
from tenpy.tools.params import get_parameter
from tenpy.networks.site import BosonSite, FermionSite
from tenpy.models.lattice import Square


class HofSqu1Model(CouplingMPOModel):

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
        lat = Square(Lx, Ly, site, order=order, bc_MPS=bc_MPS, bc=[bc_x, bc_y])
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
        t1 = get_parameter(params, 't1', 1., self.name, True)
        mu = get_parameter(params, 'mu', 0., self.name)
        phi_2pi = 2 * np.pi * get_parameter(params, 'phi', 0., self.name)
        nphi = get_parameter(params, 'nphi', nphi_default, self.name)
        Lx_MUC = get_parameter(params, 'Lx_MUC', 1, self.name)
        nphi_2pi = 2 * np.pi * nphi[0] / nphi[1]

        for u in range(len(self.lat.unit_cell)):
            self.add_onsite(mu, 0, 'N')

        u1, u2, dx = (0, 0, np.array([1, 0]))  # right
        t_phi = self.coupling_strength_add_ext_flux(t1, dx, [0, phi_2pi])
        self.add_coupling(t_phi, u1, creation, u2, annihilation, dx)
        self.add_coupling(np.conj(t_phi), u2, creation, u1, annihilation, -dx)  # H.c.
        if self.stats(params) == 'fermions':
            self.add_coupling(V, u1, 'N', u2, 'N', dx)

        u1, u2, dx = (0, 0, np.array([0, 1]))  # up
        m = np.arange(0, nphi[1] * Lx_MUC)
        t_phi = self.coupling_strength_add_ext_flux(t1, dx, [0, phi_2pi]) \
                * np.exp(-1j * nphi_2pi * m)[:, np.newaxis]
        self.add_coupling(t_phi, u1, creation, u2, annihilation, dx)
        self.add_coupling(np.conj(t_phi), u2, creation, u1, annihilation, -dx)  # H.c.
        if self.stats(params) == 'fermions':
            self.add_coupling(V, u1, 'N', u2, 'N', dx)


if __name__ == "__main__":

    model_params = dict(conserve='N', t1=1, n=(int(1), int(9)), nphi=(int(1), int(3)),
                        Lx_MUC=1, Ly=6, V=10,
                        bc_MPS='infinite', bc_x='periodic', bc_y='cylinder', order='Cstyle',
                        verbose=1, phi=0)
    M = HofSqu1Model(model_params)

    t0 = time.time()
    print("max MPO bond dimension = ", max(M.H_MPO.chi))
    print("Total time taken (seconds) = ", time.time() - t0)
