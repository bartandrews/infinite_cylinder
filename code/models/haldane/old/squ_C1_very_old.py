# --- python imports
import numpy as np
# --- TeNPy imports
from tenpy.models import lattice
# --- infinite_cylinder imports
from models.haldane.haldane import HaldaneModel


class BipartiteSquare(lattice.Lattice):
    def __init__(self, Lx, Ly, site, **kwargs):
        basis = np.array(([2., 0.], [0., 2.]))
        pos = np.array(([0., 0.], [1., 0.], [0., 1.], [1., 1.]))
        kwargs.setdefault('basis', basis)
        kwargs.setdefault('positions', pos)

        super().__init__([Lx, Ly], [site, site, site, site], **kwargs)
        self.site = site

        # counterclockwise loop (bottom, right, top, left) & spikes (top, bottom, left, right)
        self.NN = [(0, 1, np.array([0, 0])), (1, 3, np.array([0, 0])), (3, 2, np.array([0, 0])), (2, 0, np.array([0, 0])),
                   (2, 0, np.array([0, 1])), (1, 3, np.array([0, -1])), (0, 1, np.array([-1, 0])), (3, 2, np.array([1, 0]))]
        # inner cross (forward slash, back slash) & right extension (top, bottom)
        self.nNNdashed = [(0, 3, np.array([0, 0])), (2, 1, np.array([0, 0])),
                          (3, 0, np.array([1, 1])), (1, 2, np.array([1, -1]))]
        # right cross (forward slash, back slash) & above cross (forward slash, back slash)
        self.nNNdotted = [(1, 2, np.array([1, 0])), (3, 0, np.array([1, 0])),
                          (2, 1, np.array([0, 1])), (3, 0, np.array([0, 1]))]

    def plot_lattice(self):
        import matplotlib.pyplot as plt
        ax = plt.gca()
        lat = BipartiteSquare(3, 3, self.site)
        lat.plot_sites(ax)
        lat.plot_coupling(ax, lat.NN, linestyle='-', color='green')
        lat.plot_coupling(ax, lat.nNNdashed, linestyle='--', color='black')
        lat.plot_coupling(ax, lat.nNNdotted, linestyle='--', color='red')
        ax.set_aspect('equal')
        plt.show()


class HalSquC1Model(HaldaneModel):

    def __init__(self, params):
        HaldaneModel.__init__(self, params)

    def init_sites(self, params):
        site = HaldaneModel.init_sites(self, params)
        return site

    def init_lattice(self, params):
        (Lx, Ly, order, bc_MPS, bc) = self.init_lattice_params(params)
        site = self.init_sites(params)
        lat = BipartiteSquare(Lx, Ly, site, order=order, bc_MPS=bc_MPS, bc=bc)
        return lat

    def init_terms(self, params):
        (creation, annihilation, Nmax, t1, t2, t3, mu, V, Vtype, Vrange, phi_2pi) = HaldaneModel.init_terms(self, params)

        self.chemical_potential(mu, extra_dof=False)
        self.offsite_interaction("Squ", V, Vtype, Vrange, extra_dof=False)

        t1 = -t1
        if t2 is None:
            t2 = t1 * np.exp(1j * np.pi / 4)
        if t3 is None:
            t3 = t1 / np.sqrt(2)

        for u1, u2, dx in self.lat.NN:
            t1_phi = self.coupling_strength_add_ext_flux(t2, dx, [0, phi_2pi])
            self.add_coupling(t1_phi, u1, creation, u2, annihilation, dx)
            self.add_coupling(np.conj(t1_phi), u2, creation, u1, annihilation, -dx)

        for u1, u2, dx in self.lat.nNNdashed:
            t2_phi = self.coupling_strength_add_ext_flux(t3, dx, [0, phi_2pi])
            self.add_coupling(t2_phi, u1, creation, u2, annihilation, dx)
            self.add_coupling(np.conj(t2_phi), u2, creation, u1, annihilation, -dx)

        for u1, u2, dx in self.lat.nNNdotted:
            t2_phi = self.coupling_strength_add_ext_flux(t3, dx, [0, phi_2pi])
            self.add_coupling(-t2_phi, u1, creation, u2, annihilation, dx)
            self.add_coupling(-np.conj(t2_phi), u2, creation, u1, annihilation, -dx)


if __name__ == "__main__":

    model_params = dict(statistics='fermions', Nmax=1, conserve='N', t1=1, n=(int(1), int(1)),
                        LxMUC=1, Ly=3, V=0, Vtype='Coulomb', Vrange=1,
                        bc_MPS='infinite', bc_x='periodic', bc_y='cylinder', order='Cstyle',
                        verbose=1, phi=0)
    M = HalSquC1Model(model_params)

    print("max MPO bond dimension = ", max(M.H_MPO.chi))
    M.lat.plot_lattice()
