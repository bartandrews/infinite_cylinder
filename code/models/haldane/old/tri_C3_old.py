# --- python imports
import numpy as np
# --- TeNPy imports
from tenpy.models import lattice
from tenpy.networks.site import GroupedSite
# --- infinite_cylinder imports
from models.haldane.haldane import HaldaneModel


class TripartiteTriangular(lattice.Lattice):
    def __init__(self, Lx, Ly, site, **kwargs):
        basis = np.array(([3., 0.], [0.5, 0.5 * np.sqrt(3)]))
        pos = np.array(([0., 0.], [1., 0.], [2., 0.]))
        kwargs.setdefault('basis', basis)
        kwargs.setdefault('positions', pos)

        super().__init__([Lx, Ly], [site, site, site], **kwargs)
        self.site = site

        # top-left, top-right, bottom
        self.NN = [(0, 2, np.array([-1, 1])), (0, 1, np.array([0, 0])), (0, 0, np.array([0, -1])),
                   (1, 0, np.array([0, 1])), (1, 2, np.array([0, 0])), (1, 1, np.array([0, -1])),
                   (2, 1, np.array([0, 1])), (2, 0, np.array([1, 0])), (2, 2, np.array([0, -1]))]
        # top-left, right, bottom-left
        self.nNNA = [(0, 2, np.array([-1, 2])), (0, 2, np.array([0, -1])), (0, 2, np.array([-1, -1])),
                     (1, 0, np.array([0, 2])), (1, 0, np.array([1, -1])), (1, 0, np.array([0, -1])),
                     (2, 1, np.array([0, 2])), (2, 1, np.array([1, -1])), (2, 1, np.array([0, -1]))]
        # top-right, left, bottom-right
        self.nNNB = [(0, 1, np.array([0, 1])), (0, 1, np.array([-1, 1])), (0, 1, np.array([0, -2])),
                     (1, 2, np.array([0, 1])), (1, 2, np.array([-1, 1])), (1, 2, np.array([0, -2])),
                     (2, 0, np.array([1, 1])), (2, 0, np.array([0, 1])), (2, 0, np.array([1, -2]))]
        # top-left, top-right, bottom
        self.nnNN = [(0, 1, np.array([-1, 2])), (0, 2, np.array([0, 0])), (0, 0, np.array([0, -2])),
                     (1, 2, np.array([-1, 2])), (1, 0, np.array([1, 0])), (1, 1, np.array([0, -2])),
                     (2, 0, np.array([0, 2])), (2, 1, np.array([1, 0])), (2, 2, np.array([0, -2]))]

    def plot_lattice(self):
        import matplotlib.pyplot as plt
        ax = plt.gca()
        lat = TripartiteTriangular(3, 3, self.site)
        lat.plot_sites(ax)
        lat.plot_coupling(ax, lat.NN, linestyle='--', color='green')
        lat.plot_coupling(ax, lat.nNNA, linestyle='--', color='red')
        lat.plot_coupling(ax, lat.nNNB, linestyle='--', color='blue')
        lat.plot_coupling(ax, lat.nnNN, linestyle='--', color='black')
        ax.set_aspect('equal')
        plt.show()


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
        lat = TripartiteTriangular(Lx, Ly, site, order=order, bc_MPS=bc_MPS, bc=bc)
        return lat

    def init_terms(self, params):
        (creation, annihilation, Nmax, t1, mu, V, Vtype, Vrange, phi_2pi) = HaldaneModel.init_terms(self, params)

        self.chemical_potential(mu, extra_dof=True)
        self.offsite_interaction("Hex", V, Vtype, Vrange, extra_dof=True)

        t2 = 0.39 * t1 * 1j
        t3 = -0.34 * t1
        for u1, u2, dx in self.lat.NN:
            t1_phi = self.coupling_strength_add_ext_flux(t1, dx, [0, phi_2pi])
            self.add_coupling(t1_phi, u1, creation+'A', u2, annihilation+'B', dx)
            self.add_coupling(np.conj(t1_phi), u2, creation+'B', u1, annihilation+'A', -dx)
        for u1, u2, dx in self.lat.nNNA:
            t2_phi = self.coupling_strength_add_ext_flux(t2, dx, [0, phi_2pi])
            self.add_coupling(t2_phi, u1, creation+'A', u2, annihilation+'A', dx)
            self.add_coupling(np.conj(t2_phi), u2, creation+'A', u1, annihilation+'A', -dx)
        for u1, u2, dx in self.lat.nNNB:
            t2_phi = self.coupling_strength_add_ext_flux(t2, dx, [0, phi_2pi])
            self.add_coupling(t2_phi, u1, creation+'B', u2, annihilation+'B', dx)
            self.add_coupling(np.conj(t2_phi), u2, creation+'B', u1, annihilation+'B', -dx)
        for u1, u2, dx in self.lat.nnNN:
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
    M.lat.plot_lattice()
