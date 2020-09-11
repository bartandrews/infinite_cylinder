# --- python imports
import numpy as np
# --- TeNPy imports
from tenpy.models.lattice import Square
from tenpy.networks.site import GroupedSite
# --- infinite_cylinder imports
from models.haldane.haldane import HaldaneModel
import functions.func_graph as fg


class HalSquC3Model(HaldaneModel):

    def __init__(self, params):
        HaldaneModel.__init__(self, params)

    def init_sites(self, params):
        site = HaldaneModel.init_sites(self, params)
        gs = GroupedSite([site, site, site], labels=['A', 'B', 'C'], charges='same')
        gs.add_op('Ntot', gs.NA + gs.NB + gs.NC, False)
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

        orbs = ['A', 'B', 'C']
        phi = np.pi / 3

        if t2 is None:
            t2 = - t1 / np.sqrt(3)

        for u1, u2, dx in [(0, 0, np.array([1, 0]))]:  # NN right
            t1_phi = self.coupling_strength_add_ext_flux(t1, dx, [0, phi_2pi])
            for l in range(len(orbs)):
                self.add_coupling(t1_phi, u2, creation+orbs[(l+1) % 3], u1, annihilation+orbs[l], dx)
                self.add_coupling(np.conj(t1_phi), u1, creation+orbs[l], u2, annihilation+orbs[(l+1) % 3], -dx)

        for u1, u2, dx in [(0, 0, np.array([0, -1]))]:  # NN down
            for l in range(len(orbs)):
                t1_phi = self.coupling_strength_add_ext_flux(t1 * np.exp(1j * 2*(l+1) * phi), dx, [0, phi_2pi])
                self.add_coupling(t1_phi, u2, creation+orbs[(l+1) % 3], u1, annihilation+orbs[l], dx)
                self.add_coupling(np.conj(t1_phi), u1, creation+orbs[l], u2, annihilation+orbs[(l+1) % 3], -dx)

        for u1, u2, dx in [(0, 0, np.array([1, 1]))]:  # nNN top-right
            for l in range(len(orbs)):
                t2_phi = self.coupling_strength_add_ext_flux(t2*np.exp(-1j * (2*(l+1)-1) * phi), dx, [0, phi_2pi])
                self.add_coupling(t2_phi, u2, creation+orbs[l], u1, annihilation+orbs[l], dx)
                self.add_coupling(np.conj(t2_phi), u1, creation+orbs[l], u2, annihilation+orbs[l], -dx)

        for u1, u2, dx in [(0, 0, np.array([-1, -1]))]:  # nNN bottom-left
            for l in range(len(orbs)):
                t2_phi = self.coupling_strength_add_ext_flux(t2*np.exp(1j * (2*(l+1)-1) * phi), dx, [0, phi_2pi])
                self.add_coupling(t2_phi, u2, creation+orbs[l], u1, annihilation+orbs[l], dx)
                self.add_coupling(np.conj(t2_phi), u1, creation+orbs[l], u2, annihilation+orbs[l], -dx)

        for u1, u2, dx in [(0, 0, np.array([1, -1]))]:  # nNN bottom-right
            for l in range(len(orbs)):
                t2_phi = self.coupling_strength_add_ext_flux(t2*np.exp(1j * (2*(l+1)+1) * phi), dx, [0, phi_2pi])
                self.add_coupling(t2_phi, u2, creation+orbs[(l+2) % 3], u1, annihilation+orbs[l], dx)
                self.add_coupling(np.conj(t2_phi), u1, creation+orbs[l], u2, annihilation+orbs[(l+2) % 3], -dx)


if __name__ == "__main__":

    model_params = dict(statistics='fermions', Nmax=1, conserve='N', t1=1, n=(int(1), int(1)),
                        LxMUC=1, Ly=3, V=0, Vtype='Coulomb', Vrange=1,
                        bc_MPS='infinite', bc_x='periodic', bc_y='cylinder', order='Cstyle',
                        verbose=1, phi=0)
    M = HalSquC3Model(model_params)

    print("max MPO bond dimension = ", max(M.H_MPO.chi))
    fg.plot_MPS_unit_cell(M)
