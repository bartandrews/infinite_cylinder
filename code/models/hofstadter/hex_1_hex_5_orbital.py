# --- python imports
import numpy as np
# --- TeNPy imports
from tenpy.tools.params import get_parameter
from tenpy.networks.site import GroupedSite
from tenpy.models.lattice import Honeycomb
# --- infinite_cylinder imports
from models.hofstadter.hofstadter import HofstadterModel


class HofHex1Hex5OrbitalModel(HofstadterModel):

    def __init__(self, params):
        HofstadterModel.__init__(self, params)

    def init_sites(self, params):
        site = HofstadterModel.init_sites(self, params)
        gs = GroupedSite([site, site], labels=['x', 'y'], charges='same')
        gs.add_op('Ntot', gs.Nx + gs.Ny, False)
        print(sorted(gs.opnames))
        print(gs.state_labels)
        return gs

    def init_lattice(self, params):
        (Lx, Ly, order, bc_MPS, bc) = self.init_lattice_params(params)
        site = self.init_sites(params)
        lat = Honeycomb(Lx, Ly, site, order=order, bc_MPS=bc_MPS, bc=bc)
        return lat

    def init_terms(self, params):
        (creation, annihilation, nphi_default, t1, mu, V, Vtype, Vrange, nphi, nphi_2pi, LxMUC, phi_2pi) = \
            HofstadterModel.init_terms(self, params)
        t2 = get_parameter(params, 't2', 0, self.name)
        t2dash = get_parameter(params, 't2dash', 0, self.name)
        U = get_parameter(params, 'U', 0, self.name, True)

        self.chemical_potential(mu, extra_dof=True)
        self.onsite_interaction(U, 'Nx', 'Ny')
        self.offsite_interaction("Hex", V, Vtype, Vrange, extra_dof=True)
        for orbital in ['x', 'y']:
            self.hex_1_hoppings(creation+orbital, annihilation+orbital, t1, nphi, nphi_2pi, LxMUC, phi_2pi)
            self.hex_5_hoppings(creation+orbital, annihilation+orbital, t2, nphi, nphi_2pi, LxMUC, phi_2pi)

        # orbital mixing ###############################################################################################

        for sublattice in [0, 1]:  # A and B sublattices

            u1, u2, dx = (sublattice, sublattice, np.array([-1, 2]))  # up
            m = np.roll(np.arange(sublattice, 2 * nphi[1] * LxMUC, 2),
                        -1)  # match convention for strength argument of add_coupling
            t2dash_phi = self.coupling_strength_add_ext_flux(t2dash, dx, [0, phi_2pi]) \
                         * np.exp(1j * nphi_2pi * m)[:, np.newaxis]
            self.add_coupling(t2dash_phi, u1, creation + 'x', u2, annihilation + 'y', dx)
            self.add_coupling(np.conj(t2dash_phi), u2, creation + 'y', u1, annihilation + 'x', -dx)  # H.c.
            self.add_coupling(-t2dash_phi, u1, creation + 'y', u2, annihilation + 'x', dx)
            self.add_coupling(np.conj(-t2dash_phi), u2, creation + 'x', u1, annihilation + 'y', -dx)  # H.c.

            u1, u2, dx = (sublattice, sublattice, np.array([-2, 1]))  # bottom right
            m = np.roll(np.arange(sublattice, 2 * nphi[1] * LxMUC, 2),
                        -2)  # match convention for strength argument of add_coupling
            t2dash_phi = self.coupling_strength_add_ext_flux(t2dash, dx, [0, phi_2pi]) \
                         * np.exp(-1j * (nphi_2pi / 2) * (m + 3 / 2))[:, np.newaxis]
            self.add_coupling(t2dash_phi, u1, creation + 'x', u2, annihilation + 'y', dx)
            self.add_coupling(np.conj(t2dash_phi), u2, creation + 'y', u1, annihilation + 'x', -dx)  # H.c.
            self.add_coupling(-t2dash_phi, u1, creation + 'y', u2, annihilation + 'x', dx)
            self.add_coupling(np.conj(-t2dash_phi), u2, creation + 'x', u1, annihilation + 'y', -dx)  # H.c.

            u1, u2, dx = (sublattice, sublattice, np.array([-1, -1]))  # bottom left
            m = np.roll(np.arange(sublattice, 2 * nphi[1] * LxMUC, 2),
                        -1)  # match convention for strength argument of add_coupling
            t2dash_phi = self.coupling_strength_add_ext_flux(t2dash, dx, [0, phi_2pi]) \
                         * np.exp(-1j * (nphi_2pi / 2) * (m - 3 / 2))[:, np.newaxis]
            self.add_coupling(t2dash_phi, u1, creation + 'x', u2, annihilation + 'y', dx)
            self.add_coupling(np.conj(t2dash_phi), u2, creation + 'y', u1, annihilation + 'x', -dx)  # H.c.
            self.add_coupling(-t2dash_phi, u1, creation + 'y', u2, annihilation + 'x', dx)
            self.add_coupling(np.conj(-t2dash_phi), u2, creation + 'x', u1, annihilation + 'y', -dx)  # H.c.


if __name__ == "__main__":

    model_params = dict(statistics='fermions', conserve='N', t1=1, t2=-0.025, t2dash=0.1,
                        n=(int(1), int(9)), nphi=(int(1), int(3)),
                        LxMUC=1, Ly=6, V=10, Vtype='Coulomb', Vrange=1,
                        bc_MPS='infinite', bc_x='periodic', bc_y='cylinder', order='Cstyle',
                        verbose=1, phi=1.)
    M = HofHex1Hex5OrbitalModel(model_params)

    print("max MPO bond dimension = ", max(M.H_MPO.chi))
