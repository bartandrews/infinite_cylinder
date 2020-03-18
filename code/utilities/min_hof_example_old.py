import numpy as np

from tenpy.models import lattice
from tenpy.networks import site

from tenpy.models.model import CouplingMPOModel
from tenpy.tools.params import get_parameter
from tenpy.networks.site import BosonSite
from tenpy.networks.mps import MPS
from tenpy.algorithms import dmrg


class MagneticSquare(lattice.Lattice):

    def __init__(self, Lx, Ly, siteA, q, reorder, **kwargs):

        numb_sites = q

        basis = np.array(([numb_sites, 0.], [0, 1]))

        pos_list = []
        for i in range(numb_sites):
            pos_list.append([i, 0])
        pos = tuple(pos_list)

        kwargs.setdefault('order', 'default')
        kwargs.setdefault('bc', 'periodic')
        kwargs.setdefault('bc_MPS', 'infinite')
        kwargs.setdefault('basis', basis)
        kwargs.setdefault('positions', pos)

        super().__init__([Lx, Ly], [siteA] * numb_sites, **kwargs)

        # redefine order of the MPS sites
        if reorder:
            self.order = lattice.get_order_grouped([Lx, Ly, numb_sites], [(i,) for i in range(0, numb_sites)])

        NN_horiz_list = []
        for i in range(numb_sites):
            if i == numb_sites - 1:
                NN_horiz_list.append((i, 0, np.array([1, 0])))
            else:
                NN_horiz_list.append((i, i + 1, np.array([0, 0])))
        self.NN_horiz = NN_horiz_list

        for i in range(numb_sites):
            setattr(self, "NN_v{}".format(i), [(i, i, np.array([0, 1]))])


def plot_lattice(q, reorder):
    numb_sites = q

    import matplotlib.pyplot as plt

    ax = plt.gca()
    fs = site.FermionSite()
    lat = MagneticSquare(1, 4, fs, q, reorder, basis=[[numb_sites, 0], [0, 1]])
    lat.plot_sites(ax)
    lat.plot_coupling(ax, lat.NN_horiz, linestyle='-', color='black')

    for i in range(numb_sites):
        lat.plot_coupling(ax, getattr(lat, "NN_v{}".format(i)), linestyle='-', color='C{}'.format(i))

    lat.plot_order(ax)
    ax.set_aspect('equal')
    plt.show()


class BosonicHofstadterModel(CouplingMPOModel):

    def __init__(self, model_params):
        CouplingMPOModel.__init__(self, model_params)

    def init_sites(self, model_params):
        Nmax = get_parameter(model_params, 'Nmax', 1, self.name)
        conserve = get_parameter(model_params, 'conserve', 'N', self.name)
        filling = get_parameter(model_params, 'filling', (1, 8), self.name)
        filling = filling[0] / filling[1]
        site = BosonSite(Nmax=Nmax, conserve=conserve, filling=filling)

        print(site)

        return site

    def init_lattice(self, model_params):
        bc_MPS = get_parameter(model_params, 'bc_MPS', 'infinite', self.name)
        order = get_parameter(model_params, 'order', 'default', self.name)
        site = self.init_sites(model_params)
        Lx = get_parameter(model_params, 'Lx', 1, self.name)
        Ly = get_parameter(model_params, 'Ly', 4, self.name)
        phi = get_parameter(model_params, 'phi', 4, self.name)
        bc_x = 'periodic' if bc_MPS == 'infinite' else 'open'  # Next line needs default
        bc_x = get_parameter(model_params, 'bc_x', bc_x, self.name)
        bc_y = get_parameter(model_params, 'bc_y', 'cylinder', self.name)
        assert bc_y in ['cylinder', 'ladder']
        bc_y = 'periodic' if bc_y == 'cylinder' else 'open'
        if bc_MPS == 'infinite' and bc_x == 'open':
            raise ValueError("You need to use 'periodic' `bc_x` for infinite systems!")

        print(Lx, Ly, site, order, [bc_x, bc_y], bc_MPS)

        reorder = get_parameter(model_params, 'reorder', True, self.name)
        lat = MagneticSquare(Lx, Ly, site, phi[1], reorder, order=order, bc=[bc_x, bc_y], bc_MPS=bc_MPS)
        return lat

    def init_terms(self, model_params):

        t = get_parameter(model_params, 't', 1., self.name, True)
        phi_ext = 2 * np.pi * get_parameter(model_params, 'phi_ext', 0., self.name)

        phi_p, phi_q = get_parameter(model_params, 'phi', (1, 4), self.name)
        phi = 2 * np.pi * phi_p / phi_q

        print(t, (phi_p, phi_q), phi, phi_ext)

        for u1, u2, dx in self.lat.NN_horiz:
            print("lat_NN_horiz = ", u1, u2, dx)
            t_phi = self.coupling_strength_add_ext_flux(t, dx, [0, phi_ext])

            print("(right) t_phi = ", t_phi)

            self.add_coupling(t_phi, u1, 'Bd', u2, 'B', dx)
            self.add_coupling(np.conj(t_phi), u2, 'Bd', u1, 'B', -dx)  # h.c.

        for i in range(phi_q):
            for u1, u2, dx in getattr(self.lat, "NN_v{}".format(i)):
                print("lat_NN_v = ", u1, u2, dx)
                t_phi = self.coupling_strength_add_ext_flux(t, dx, [0, phi_ext]) * np.exp(-1j * phi * i)

                print("(up) t_phi = ", t_phi)

                self.add_coupling(t_phi, u1, 'Bd', u2, 'B', dx)
                self.add_coupling(np.conj(t_phi), u2, 'Bd', u1, 'B', -dx)  # h.c.


if __name__ == "__main__":

    reorder_param = True

    plot_lattice(5, reorder=reorder_param)
    model_params = dict(conserve='N', t=1, filling=(int(1), int(10)), phi=(int(1), int(5)), Lx=1,
                        Ly=4, Nmax=1,  # system params
                        bc_MPS='infinite', bc_x='periodic', bc_y='cylinder', order='Cstyle',  # MPS params
                        verbose=1, phi_ext=0, reorder=reorder_param)  # utility
    M = BosonicHofstadterModel(model_params)
    # model_params2 = dict(conserve='N', t=1, filling=(int(1), int(10)), phi=(int(1), int(5)), Lx=1,
    #                     Ly=4, Nmax=1,  # system params
    #                     bc_MPS='infinite', bc_x='periodic', bc_y='cylinder', order='Cstyle',  # MPS params
    #                     verbose=1, phi_ext=0, reorder=True)  # utility
    # M2 = BosonicHofstadterModel(model_params2)
    #
    # assert M.H_MPO.L == M2.H_MPO.L
    # assert M.H_MPO.is_hermitian()
    # assert M2.H_MPO.is_hermitian()
    # assert M.H_MPO.is_equal(M2.H_MPO)

    import pprint
    pprint.pprint(list(M.coupling_terms['Bd_i B_j'].to_TermList()))

    product_state = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0]

    psi = MPS.from_product_state(M.lat.mps_sites(), product_state, bc=M.lat.bc_MPS)

    dmrg_params = {
        'mixer': True,  # setting this to True helps to escape local minima
        'mixer_params': {'amplitude': 1.e-5, 'decay': 1.2, 'disable_after': 30},
        'trunc_params': {
            # 'chi_max': chi_max,
            'svd_min': 1.e-10
        },
        # 'lanczos_params': {
        #     'reortho': True,
        #     'N_cache': 40
        # },
        'chi_list': {0: 9, 10: 49, 20: 100, 40: 50},  # may need to change first value ~39
        'max_E_err': 1.e-6,
        'max_S_err': 1.e-6,
        # 'norm_tol': 1.e-6,
        # 'norm_tol_iter': 1000,
        'max_sweeps': 1000,
        'verbose': 5,
        'N_sweeps_check': 10,
        'diag_method': 'lanczos'
    }

    info = dmrg.run(psi, M, dmrg_params)
    E = info['E']

    print(psi.entanglement_entropy()[0])
