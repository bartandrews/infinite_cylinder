# --- python imports
import numpy as np
import string
# --- TeNPy imports
from tenpy.models.lattice import Square
from tenpy.networks.site import GroupedSite
# --- infinite_cylinder imports
from models.haldane.haldane import HaldaneModel
import functions.func_graph as fg


class HalSquCNModel(HaldaneModel):

    def __init__(self, params):
        HaldaneModel.__init__(self, params)

    @staticmethod
    def init_orbitals(params):
        C = params.get('C', 3)
        orbs = list(string.ascii_uppercase)[:C]
        print('C, orbs = ', C, orbs)
        return C, orbs

    def init_sites(self, params):
        site = HaldaneModel.init_sites(self, params)
        C, orbs = self.init_orbitals(params)
        gs = GroupedSite([site]*len(orbs), labels=orbs, charges='same')
        # define the Ntot attribute
        Ntot = getattr(gs, 'N'+orbs[0])
        for m in range(1, len(orbs)):
            Ntot += getattr(gs, 'N'+orbs[m])
        gs.add_op('Ntot', Ntot, False)
        print(sorted(gs.opnames))
        print(gs.state_labels)
        return gs

    def init_lattice(self, params):
        (Lx, Ly, order, bc_MPS, bc) = self.init_lattice_params(params)
        site = self.init_sites(params)
        lat = Square(Lx, Ly, site, order=order, bc_MPS=bc_MPS, bc=bc)
        return lat

    def init_terms(self, params):
        (creation, annihilation, Nmax, t1, t2, _, mu, V, Vtype, Vrange, phi_2pi) = HaldaneModel.init_terms(self, params)
        C, orbs = self.init_orbitals(params)

        self.chemical_potential(mu, extra_dof=True)
        self.offsite_interaction("Squ", V, Vtype, Vrange, extra_dof=True)

        phi = np.pi / C

        for u1, u2, dx in [(0, 0, np.array([1, 0]))]:  # NN right
            t1_phi = self.coupling_strength_add_ext_flux(t1, dx, [0, phi_2pi])
            for m in range(len(orbs)):
                self.add_coupling(t1_phi, u2, creation+orbs[(m+1) % 3], u1, annihilation+orbs[m], dx)
                self.add_coupling(np.conj(t1_phi), u1, creation+orbs[m], u2, annihilation+orbs[(m+1) % 3], -dx)

        for u1, u2, dx in [(0, 0, np.array([0, -1]))]:  # NN down
            for m in range(len(orbs)):
                t1_phi = self.coupling_strength_add_ext_flux(t1 * np.exp(1j * 2*(m+1) * phi), dx, [0, phi_2pi])
                self.add_coupling(t1_phi, u2, creation+orbs[(m+1) % 3], u1, annihilation+orbs[m], dx)
                self.add_coupling(np.conj(t1_phi), u1, creation+orbs[m], u2, annihilation+orbs[(m+1) % 3], -dx)

        if C < 5:  # include nNN hoppings

            if t2 is None:
                t2 = - t1 / np.sqrt(C)

            for u1, u2, dx in [(0, 0, np.array([1, 1]))]:  # nNN top-right
                for m in range(len(orbs)):
                    t2_phi = self.coupling_strength_add_ext_flux(t2*np.exp(-1j * (2*(m+1)-1) * phi), dx, [0, phi_2pi])
                    self.add_coupling(t2_phi, u2, creation+orbs[m], u1, annihilation+orbs[m], dx)
                    self.add_coupling(np.conj(t2_phi), u1, creation+orbs[m], u2, annihilation+orbs[m], -dx)

            for u1, u2, dx in [(0, 0, np.array([-1, -1]))]:  # nNN bottom-left
                for m in range(len(orbs)):
                    t2_phi = self.coupling_strength_add_ext_flux(t2*np.exp(1j * (2*(m+1)-1) * phi), dx, [0, phi_2pi])
                    self.add_coupling(t2_phi, u2, creation+orbs[m], u1, annihilation+orbs[m], dx)
                    self.add_coupling(np.conj(t2_phi), u1, creation+orbs[m], u2, annihilation+orbs[m], -dx)

            for u1, u2, dx in [(0, 0, np.array([1, -1]))]:  # nNN bottom-right
                for m in range(len(orbs)):
                    t2_phi = self.coupling_strength_add_ext_flux(t2*np.exp(1j * (2*(m+1)+1) * phi), dx, [0, phi_2pi])
                    self.add_coupling(t2_phi, u2, creation+orbs[(m+2) % 3], u1, annihilation+orbs[m], dx)
                    self.add_coupling(np.conj(t2_phi), u1, creation+orbs[m], u2, annihilation+orbs[(m+2) % 3], -dx)


if __name__ == "__main__":

    model_params = dict(statistics='fermions', Nmax=1, conserve='N', C=3, t1=1, n=(int(1), int(1)),
                        LxMUC=1, Ly=3, V=0, Vtype='Coulomb', Vrange=1,
                        bc_MPS='infinite', bc_x='periodic', bc_y='cylinder', order='Cstyle',
                        verbose=1, phi=0)
    M = HalSquCNModel(model_params)

    print("max MPO bond dimension = ", max(M.H_MPO.chi))
    fg.plot_MPS_unit_cell(M)
