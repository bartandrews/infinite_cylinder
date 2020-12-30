# --- python imports
import numpy as np
import matplotlib.pyplot as plt
# --- TeNPy imports
from tenpy.models import lattice
# --- infinite_cylinder imports
from models.bbh.bbh import BBHModel


class QuadpartiteSquare(lattice.Lattice):
    def __init__(self, Lx, Ly, site, **kwargs):
        basis = np.array(([2., 0.], [0., 2.]))
        pos = np.array(([0., 0.], [1., 0.], [0., 1.], [1., 1.]))
        kwargs.setdefault('basis', basis)
        kwargs.setdefault('positions', pos)

        super().__init__([Lx, Ly], [site, site, site, site], **kwargs)
        self.site = site

        # gamma_x
        self.gamma_x = [(2, 3, np.array([0, 0])), (1, 0, np.array([0, 0]))]
        # pos_gamma_y
        self.pos_gamma_y = [(1, 3, np.array([0, 0]))]
        # neg_gamma_y
        self.neg_gamma_y = [(2, 0, np.array([0, 0]))]

        # lambda_x
        self.lambda_x = [(0, 1, np.array([-1, 0])), (2, 3, np.array([-1, 0]))]
        # pos_lambda_y
        self.pos_lambda_y = [(1, 3, np.array([0, -1]))]
        # neg_lambda_y
        self.neg_lambda_y = [(0, 2, np.array([0, -1]))]

    def plot_lattice(self):
        import matplotlib.pyplot as plt
        ax = plt.gca()
        lat = QuadpartiteSquare(3, 3, self.site)
        lat.plot_sites(ax)
        lat.plot_coupling(ax, lat.gamma_x, linestyle='-', color='red')
        lat.plot_coupling(ax, lat.pos_gamma_y, linestyle='-', color='red')
        lat.plot_coupling(ax, lat.neg_gamma_y, linestyle='-', color='red')
        lat.plot_coupling(ax, lat.lambda_x, linestyle='--', color='blue')
        lat.plot_coupling(ax, lat.pos_lambda_y, linestyle='--', color='blue')
        lat.plot_coupling(ax, lat.neg_lambda_y, linestyle='--', color='blue')
        ax.set_aspect('equal')
        lat.plot_order(ax, linestyle='dotted', color='k')
        plt.show()


class BBHSqu2dModel(BBHModel):

    def __init__(self, params):
        BBHModel.__init__(self, params)

    def init_sites(self, params):
        site = BBHModel.init_sites(self, params)
        return site

    def init_lattice(self, params):
        (Lx, Ly, order, bc_MPS, bc) = self.init_lattice_params(params)
        site = self.init_sites(params)
        lat = QuadpartiteSquare(Lx, Ly, site, order=order, bc_MPS=bc_MPS, bc=bc)
        return lat

    def init_terms(self, params):
        (creation, annihilation, Nmax, t1, t2, t3, t4) = BBHModel.init_terms(self, params)

        for u1, u2, dx in self.lat.gamma_x:
            self.add_coupling(t1, u1, creation, u2, annihilation, dx)
            self.add_coupling(np.conj(t1), u2, creation, u1, annihilation, -dx)
        for u1, u2, dx in self.lat.pos_gamma_y:
            self.add_coupling(t2, u1, creation, u2, annihilation, dx)
            self.add_coupling(np.conj(t2), u2, creation, u1, annihilation, -dx)
        for u1, u2, dx in self.lat.neg_gamma_y:
            self.add_coupling(-t2, u1, creation, u2, annihilation, dx)
            self.add_coupling(np.conj(-t2), u2, creation, u1, annihilation, -dx)

        for u1, u2, dx in self.lat.lambda_x:
            self.add_coupling(t3, u1, creation, u2, annihilation, dx)
            self.add_coupling(np.conj(t3), u2, creation, u1, annihilation, -dx)
        for u1, u2, dx in self.lat.pos_lambda_y:
            self.add_coupling(t4, u1, creation, u2, annihilation, dx)
            self.add_coupling(np.conj(t4), u2, creation, u1, annihilation, -dx)
        for u1, u2, dx in self.lat.neg_lambda_y:
            self.add_coupling(-t4, u1, creation, u2, annihilation, dx)
            self.add_coupling(np.conj(-t4), u2, creation, u1, annihilation, -dx)


if __name__ == "__main__":

    model_params = dict(statistics='fermions', Nmax=1, conserve='N', t1=0.5, t2=1, t3=0.5, t4=1, n=(int(1), int(2)),
                        LxMUC=3, Ly=3,
                        bc_MPS='infinite', bc_x='periodic', bc_y='cylinder', order='Cstyle',
                        verbose=1)
    M = BBHSqu2dModel(model_params)

    print("max MPO bond dimension = ", max(M.H_MPO.chi))
    M.lat.plot_lattice()
