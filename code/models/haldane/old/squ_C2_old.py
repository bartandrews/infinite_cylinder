# --- python imports
import numpy as np
# --- TeNPy imports
from tenpy.models import lattice
from tenpy.networks.site import GroupedSite
# --- infinite_cylinder imports
from models.haldane.haldane import HaldaneModel


class BipartiteSquare(lattice.Lattice):
    def __init__(self, Lx, Ly, site, **kwargs):
        basis = np.array(([2., 0.], [0., 1.]))
        pos = np.array(([0., 0.], [1., 0.]))
        kwargs.setdefault('basis', basis)
        kwargs.setdefault('positions', pos)

        super().__init__([Lx, Ly], [site, site], **kwargs)
        self.site = site

        self.NN_A_to_B_out = [(0, 1, np.array([0, 0])), (1, 0, np.array([1, 0]))]
        self.NN_A_to_B_in = [(0, 0, np.array([0, 1])), (1, 0, np.array([0, 1]))]
        self.NN_B_to_A_out = [(0, 0, np.array([0, 1])), (1, 0, np.array([0, 1]))]
        self.NN_B_to_A_in = [(0, 1, np.array([0, 0])), (1, 0, np.array([1, 0]))]

        self.nNN_A_pos = [(0, 1, np.array([0, 1])), (1, 0, np.array([1, 1]))]  # red solid
        self.nNN_A_neg = [(0, 1, np.array([-1, 1])), (1, 0, np.array([0, 1]))]  # red dashed
        self.nNN_B_pos = [(0, 1, np.array([-1, 1])), (1, 0, np.array([0, 1]))]  # blue solid
        self.nNN_B_neg = [(0, 1, np.array([0, 1])), (1, 0, np.array([1, 1]))]  # blue dashed

        self.nnNN = [(0, 0, np.array([0, 2])), (1, 1, np.array([0, 2])),
                     (0, 0, np.array([2, 0])), (1, 1, np.array([2, 0]))]


    def plot_lattice(self):
        import matplotlib.pyplot as plt
        ax = plt.gca()
        lat = BipartiteSquare(3, 3, self.site)
        lat.plot_sites(ax)
        lat.plot_coupling(ax, lat.nNN_A_pos, linestyle='solid', color='red')
        lat.plot_coupling(ax, lat.nNN_A_neg, linestyle='dashed', color='red')
        lat.plot_coupling(ax, lat.nNN_B_pos, linestyle='solid', color='blue')
        lat.plot_coupling(ax, lat.nNN_B_neg, linestyle='dashed', color='blue')
        lat.plot_coupling(ax, lat.nnNN, linestyle='dotted', color='black')
        ax.set_aspect('equal')
        plt.show()


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
        lat = BipartiteSquare(Lx, Ly, site, order=order, bc_MPS=bc_MPS, bc=bc)
        return lat

    def init_terms(self, params):
        (creation, annihilation, Nmax, t1, mu, V, Vtype, Vrange, phi_2pi) = HaldaneModel.init_terms(self, params)

        self.chemical_potential(mu, extra_dof=False)
        self.offsite_interaction("Squ", V, Vtype, Vrange, extra_dof=False)

        t2 = t1 / (2 + np.sqrt(2))
        t3 = t1 / (2 + 2*np.sqrt(2))

        for u1, u2, dx in self.lat.NN_A_to_B_out:
            t1_phi = self.coupling_strength_add_ext_flux(t1*np.exp(1j*np.pi/4), dx, [0, phi_2pi])
            self.add_coupling(t1_phi, u2, creation+'B', u1, annihilation+'A', dx)
            self.add_coupling(np.conj(t1_phi), u1, creation+'A', u2, annihilation+'B', -dx)

        for u1, u2, dx in self.lat.NN_A_to_B_in:
            t1_phi = self.coupling_strength_add_ext_flux(t1*np.exp(1j*np.pi/4), dx, [0, phi_2pi])
            self.add_coupling(t1_phi, u1, creation+'A', u2, annihilation+'B', dx)
            self.add_coupling(np.conj(t1_phi), u2, creation+'B', u1, annihilation+'A', -dx)

        for u1, u2, dx in self.lat.NN_B_to_A_out:
            t1_phi = self.coupling_strength_add_ext_flux(t1*np.exp(1j*np.pi/4), dx, [0, phi_2pi])
            self.add_coupling(t1_phi, u2, creation+'A', u1, annihilation+'B', dx)
            self.add_coupling(np.conj(t1_phi), u1, creation+'B', u2, annihilation+'A', -dx)

        for u1, u2, dx in self.lat.NN_B_to_A_in:
            t1_phi = self.coupling_strength_add_ext_flux(t1*np.exp(1j*np.pi/4), dx, [0, phi_2pi])
            self.add_coupling(t1_phi, u1, creation+'B', u2, annihilation+'A', dx)
            self.add_coupling(np.conj(t1_phi), u2, creation+'A', u1, annihilation+'B', -dx)

        for u1, u2, dx in self.lat.nNN_A_pos:
            t2_phi = self.coupling_strength_add_ext_flux(t2, dx, [0, phi_2pi])
            self.add_coupling(t2_phi, u1, creation+'A', u2, annihilation+'A', dx)
            self.add_coupling(np.conj(t2_phi), u2, creation+'A', u1, annihilation+'A', -dx)

        for u1, u2, dx in self.lat.nNN_A_neg:
            t2_phi = self.coupling_strength_add_ext_flux(-t2, dx, [0, phi_2pi])
            self.add_coupling(t2_phi, u1, creation+'A', u2, annihilation+'A', dx)
            self.add_coupling(np.conj(t2_phi), u2, creation+'A', u1, annihilation+'A', -dx)

        for u1, u2, dx in self.lat.nNN_B_pos:
            t2_phi = self.coupling_strength_add_ext_flux(t2, dx, [0, phi_2pi])
            self.add_coupling(t2_phi, u1, creation+'B', u2, annihilation+'B', dx)
            self.add_coupling(np.conj(t2_phi), u2, creation+'B', u1, annihilation+'B', -dx)

        for u1, u2, dx in self.lat.nNN_B_neg:
            t2_phi = self.coupling_strength_add_ext_flux(-t2, dx, [0, phi_2pi])
            self.add_coupling(t2_phi, u1, creation+'B', u2, annihilation+'B', dx)
            self.add_coupling(np.conj(t2_phi), u2, creation+'B', u1, annihilation+'B', -dx)

        for u1, u2, dx in self.lat.nnNN:
            t2_phi = self.coupling_strength_add_ext_flux(t3, dx, [0, phi_2pi])
            for orbital in ['A', 'B']:
                self.add_coupling(-t2_phi, u1, creation+orbital, u2, annihilation+orbital, dx)
                self.add_coupling(-np.conj(t2_phi), u2, creation+orbital, u1, annihilation+orbital, -dx)


if __name__ == "__main__":

    model_params = dict(statistics='fermions', Nmax=1, conserve='N', t1=1, n=(int(1), int(1)),
                        LxMUC=1, Ly=3, V=0, Vtype='Coulomb', Vrange=1,
                        bc_MPS='infinite', bc_x='periodic', bc_y='cylinder', order='Cstyle',
                        verbose=1, phi=0)
    M = HalSquC2Model(model_params)

    print("max MPO bond dimension = ", max(M.H_MPO.chi))
    M.lat.plot_lattice()
