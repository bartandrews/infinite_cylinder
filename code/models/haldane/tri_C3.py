# --- python imports
import numpy as np
# --- TeNPy imports
from tenpy.models.lattice import Triangular
from tenpy.networks.site import GroupedSite
# --- infinite_cylinder imports
from models.haldane.haldane import HaldaneModel
import functions.func_graph as fg


class HalTriC3Model(HaldaneModel):

    def __init__(self, params):
        HaldaneModel.__init__(self, params)

    def init_sites(self, params):
        site = HaldaneModel.init_sites(self, params)
        gs = GroupedSite([site, site], labels=['A', 'B'], charges='same')
        gs.add_op('Ntot', gs.NA + gs.NB, False)
        print(sorted(gs.opnames))
        print(gs.state_labels)
        return gs

    def init_lattice(self, params):
        (Lx, Ly, order, bc_MPS, bc) = self.init_lattice_params(params)
        site = self.init_sites(params)
        lat = Triangular(Lx, Ly, site, order=order, bc_MPS=bc_MPS, bc=bc)
        return lat

    def init_terms(self, params):
        (creation, annihilation, Nmax, t1, t2, t3, mu, V, Vtype, Vrange, phi_2pi) = HaldaneModel.init_terms(self, params)

        self.chemical_potential(mu, extra_dof=True)
        self.offsite_interaction("Tri", V, Vtype, Vrange, extra_dof=True)

        # t1 = -t1  # try +, since - sign didn't behave as expected for chi=100 FCI (even though it should be -)
        if t2 is None:
            t2 = 0.39 * t1 * 1j
        if t3 is None:
            t3 = -0.34 * t1

        # NN (top-left, top-right, bottom)
        for u1, u2, dx in [(0, 0, np.array([-1, 1])), (0, 0, np.array([1, 0])), (0, 0, np.array([0, -1]))]:
            t1_phi = self.coupling_strength_add_ext_flux(t1, dx, [0, phi_2pi])
            self.add_coupling(t1_phi, u1, creation+'A', u2, annihilation+'B', dx)
            self.add_coupling(np.conj(t1_phi), u2, creation+'B', u1, annihilation+'A', -dx)

        # nNNA (top-left, right, bottom-left)
        for u1, u2, dx in [(0, 0, np.array([-1, 2])), (0, 0, np.array([2, -1])), (0, 0, np.array([-1, -1]))]:
            t2_phi = self.coupling_strength_add_ext_flux(t2, dx, [0, phi_2pi])
            self.add_coupling(t2_phi, u1, creation+'A', u2, annihilation+'A', dx)
            self.add_coupling(np.conj(t2_phi), u2, creation+'A', u1, annihilation+'A', -dx)

        # nNNB (top-right, left, bottom-right)
        for u1, u2, dx in [(0, 0, np.array([1, 1])), (0, 0, np.array([-2, 1])), (0, 0, np.array([1, -2]))]:
            t2_phi = self.coupling_strength_add_ext_flux(t2, dx, [0, phi_2pi])
            self.add_coupling(t2_phi, u1, creation + 'B', u2, annihilation + 'B', dx)
            self.add_coupling(np.conj(t2_phi), u2, creation + 'B', u1, annihilation + 'B', -dx)

        # nnNN (top-left, top-right, bottom)
        for u1, u2, dx in [(0, 0, np.array([-2, 2])), (0, 0, np.array([2, 0])), (0, 0, np.array([0, -2]))]:
            t3_phi = self.coupling_strength_add_ext_flux(t3, dx, [0, phi_2pi])
            self.add_coupling(t3_phi, u1, creation+'A', u2, annihilation+'B', dx)
            self.add_coupling(np.conj(t3_phi), u2, creation+'B', u1, annihilation+'A', -dx)


if __name__ == "__main__":

    model_params = dict(statistics='fermions', Nmax=1, conserve='N', t1=1, n=(int(1), int(6)),
                        LxMUC=1, Ly=6, V=1, Vtype='Coulomb', Vrange=1,
                        bc_MPS='infinite', bc_x='periodic', bc_y='cylinder', order='Cstyle',
                        verbose=1, phi=0)
    M = HalTriC3Model(model_params)

    print("max MPO bond dimension = ", max(M.H_MPO.chi))
    fg.plot_MPS_unit_cell(M)
