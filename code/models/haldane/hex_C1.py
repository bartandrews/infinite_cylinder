# --- python imports
import numpy as np
# --- TeNPy imports
from tenpy.models.lattice import Honeycomb
# --- infinite_cylinder imports
from models.haldane.haldane import HaldaneModel
import functions.func_graph as fg


class HalHexC1Model(HaldaneModel):

    def __init__(self, params):
        HaldaneModel.__init__(self, params)

    def init_sites(self, params):
        site = HaldaneModel.init_sites(self, params)
        return site

    def init_lattice(self, params):
        (Lx, Ly, order, bc_MPS, bc) = self.init_lattice_params(params)
        site = self.init_sites(params)
        lat = Honeycomb(Lx, Ly, site, order=order, bc_MPS=bc_MPS, bc=bc)
        return lat

    def init_terms(self, params):
        (creation, annihilation, Nmax, t1, t2, _, mu, V, Vtype, Vrange, phi_2pi) = HaldaneModel.init_terms(self, params)

        self.chemical_potential(mu)
        self.offsite_interaction("Hex", V, Vtype, Vrange)

        if t2 is None:  # not given
            t2 = (np.sqrt(129) / 36) * t1 * np.exp(1j * np.arccos(3 * np.sqrt(3 / 43)))

        for u1, u2, dx in self.lat.pairs['nearest_neighbors']:
            t1_phi = self.coupling_strength_add_ext_flux(t1, dx, [0, phi_2pi])
            self.add_coupling(t1_phi, u1, creation, u2, annihilation, dx)
            self.add_coupling(np.conj(t1_phi), u2, creation, u1, annihilation, -dx)  # H.c.
        for u1, u2, dx in [(0, 0, np.array([-1, 1])), (0, 0, np.array([1, 0])),
                           (0, 0, np.array([0, -1])), (1, 1, np.array([0, 1])),
                           (1, 1, np.array([1, -1])), (1, 1, np.array([-1, 0]))]:
            t2_phi = self.coupling_strength_add_ext_flux(t2, dx, [0, phi_2pi])
            self.add_coupling(t2_phi, u1, creation, u2, annihilation, dx)
            self.add_coupling(np.conj(t2_phi), u2, creation, u1, annihilation, -dx)  # H.c.


if __name__ == "__main__":

    model_params = dict(statistics='fermions', Nmax=1, conserve='N', t1=1, n=(int(1), int(6)),
                        LxMUC=1, Ly=6, V=1, Vtype='Coulomb', Vrange=1,
                        bc_MPS='infinite', bc_x='periodic', bc_y='cylinder', order='Cstyle',
                        verbose=1, phi=0)
    M = HalHexC1Model(model_params)

    print("max MPO bond dimension = ", max(M.H_MPO.chi))
    fg.plot_MPS_unit_cell(M)
