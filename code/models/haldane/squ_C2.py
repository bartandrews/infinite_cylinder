# --- python imports
import numpy as np
# --- TeNPy imports
from tenpy.models.lattice import Square
from tenpy.networks.site import GroupedSite
# --- infinite_cylinder imports
from models.haldane.haldane import HaldaneModel
import functions.func_graph as fg


class HalSquC2Model(HaldaneModel):

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
        lat = Square(Lx, Ly, site, order=order, bc_MPS=bc_MPS, bc=bc)
        return lat

    def init_terms(self, params):
        (creation, annihilation, Nmax, t1, t2, t3, mu, V, Vtype, Vrange, phi_2pi) = HaldaneModel.init_terms(self, params)

        self.chemical_potential(mu, extra_dof=True)
        self.offsite_interaction("Squ", V, Vtype, Vrange, extra_dof=True)

        t1_phase = t1 * np.exp(1j * np.pi / 4)
        if t2 is None:
            t2 = t1 / (2 + np.sqrt(2))
        if t3 is None:
            t3 = t1 / (2 + 2*np.sqrt(2))

        for u1, u2, dx in [(0, 0, np.array([0, 1]))]:  # NN up
            t1_phi = self.coupling_strength_add_ext_flux(t1_phase, dx, [0, phi_2pi])
            self.add_coupling(t1_phi, u2, creation+'A', u1, annihilation+'B', dx)
            self.add_coupling(np.conj(t1_phi), u1, creation+'B', u2, annihilation+'A', -dx)

        for u1, u2, dx in [(0, 0, np.array([0, -1]))]:  # NN down
            t1_phi = self.coupling_strength_add_ext_flux(t1_phase, dx, [0, phi_2pi])
            self.add_coupling(t1_phi, u2, creation+'A', u1, annihilation+'B', dx)
            self.add_coupling(np.conj(t1_phi), u1, creation+'B', u2, annihilation+'A', -dx)

        for u1, u2, dx in [(0, 0, np.array([-1, 0]))]:  # NN left
            t1_phi = self.coupling_strength_add_ext_flux(t1_phase, dx, [0, phi_2pi])
            self.add_coupling(t1_phi, u2, creation+'B', u1, annihilation+'A', dx)
            self.add_coupling(np.conj(t1_phi), u1, creation+'A', u2, annihilation+'B', -dx)

        for u1, u2, dx in [(0, 0, np.array([1, 0]))]:  # NN right
            t1_phi = self.coupling_strength_add_ext_flux(t1_phase, dx, [0, phi_2pi])
            self.add_coupling(t1_phi, u2, creation+'B', u1, annihilation+'A', dx)
            self.add_coupling(np.conj(t1_phi), u1, creation+'A', u2, annihilation+'B', -dx)

        for u1, u2, dx in [(0, 0, np.array([1, 1]))]:  # nNN upper-right
            t2_phi = self.coupling_strength_add_ext_flux(t2, dx, [0, phi_2pi])
            self.add_coupling(t2_phi, u1, creation+'A', u2, annihilation+'A', dx)
            self.add_coupling(np.conj(t2_phi), u2, creation+'A', u1, annihilation+'A', -dx)
            self.add_coupling(-t2_phi, u1, creation+'B', u2, annihilation+'B', dx)
            self.add_coupling(np.conj(-t2_phi), u2, creation+'B', u1, annihilation+'B', -dx)

        for u1, u2, dx in [(0, 0, np.array([1, -1]))]:  # nNN lower-right
            t2_phi = self.coupling_strength_add_ext_flux(t2, dx, [0, phi_2pi])
            self.add_coupling(-t2_phi, u1, creation+'A', u2, annihilation+'A', dx)
            self.add_coupling(np.conj(-t2_phi), u2, creation+'A', u1, annihilation+'A', -dx)
            self.add_coupling(t2_phi, u1, creation+'B', u2, annihilation+'B', dx)
            self.add_coupling(np.conj(t2_phi), u2, creation+'B', u1, annihilation+'B', -dx)

        for u1, u2, dx in [(0, 0, np.array([0, 2])), (0, 0, np.array([2, 0]))]:  # nnNN (up, right)
            t3_phi = self.coupling_strength_add_ext_flux(t3, dx, [0, phi_2pi])
            for orbital in ['A', 'B']:
                self.add_coupling(t3_phi, u1, creation+orbital, u2, annihilation+orbital, dx)
                self.add_coupling(np.conj(t3_phi), u2, creation+orbital, u1, annihilation+orbital, -dx)


if __name__ == "__main__":

    model_params = dict(statistics='fermions', Nmax=1, conserve='N', t1=1, n=(int(1), int(1)),
                        LxMUC=1, Ly=3, V=0, Vtype='Coulomb', Vrange=1,
                        bc_MPS='infinite', bc_x='periodic', bc_y='cylinder', order='Cstyle',
                        verbose=1, phi=0)
    M = HalSquC2Model(model_params)

    print("max MPO bond dimension = ", max(M.H_MPO.chi))
    fg.plot_MPS_unit_cell(M)
