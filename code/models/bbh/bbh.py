# --- python imports
import numpy as np
# --- TeNPy imports
from tenpy.models.model import CouplingMPOModel
from tenpy.networks.site import BosonSite, FermionSite
from tenpy.models.lattice import Square
# --- infinite_cylinder imports
import functions.func_graph as fg


def stats(params):
    return params.get('statistics', 'bosons')


class BBHModel(CouplingMPOModel):

    def __init__(self, params):
        CouplingMPOModel.__init__(self, params)

    def init_sites(self, params):
        conserve = params.get('conserve', 'N')
        if stats(params) == 'bosons':
            Nmax = params.get('Nmax', 1)
            n = params.get('n', (1, 2))
            n = n[0] / n[1]
            site = BosonSite(Nmax=Nmax, conserve=conserve, filling=n)
        else:
            n = params.get('n', (1, 2))
            n = n[0] / n[1]
            site = FermionSite(conserve=conserve, filling=n)
        return site

    def init_lattice(self, params):
        Lx = params.get('LxMUC', 3)
        Ly = params.get('Ly', 3)
        order = params.get('order', 'Cstyle')
        bc_MPS = params.get('bc_MPS', 'infinite')
        bc_x = 'periodic' if bc_MPS == 'infinite' else 'open'  # Next line needs default
        bc_x = params.get('bc_x', bc_x)
        bc_y = params.get('bc_y', 'cylinder')
        assert bc_y in ['cylinder', 'ladder']
        bc_y = 'periodic' if bc_y == 'cylinder' else 'open'
        if bc_MPS == 'infinite' and bc_x == 'open':
            raise ValueError("You need to use 'periodic' `bc_x` for infinite systems!")
        bc = [bc_x, bc_y]
        site = self.init_sites(params)
        lat = Square(Lx, Ly, site, order=order, bc_MPS=bc_MPS, bc=bc)
        return lat

    def init_terms(self, params):
        if stats(params) == 'bosons':
            creation, annihilation = 'Bd', 'B'
        else:
            creation, annihilation = 'Cd', 'C'
        Nmax = params.get('Nmax', 1)
        t1 = params.get('t1', 1)
        t2 = params.get('t2', 1)
        t3 = params.get('t3', 1)
        t4 = params.get('t4', 1)

        # gamma_x
        t1_coupling = np.array([[np.conj(t1), t1, np.conj(t1), t1, np.conj(t1), t1],
                               [0, 0, 0, 0, 0, 0],
                               [np.conj(t1), t1, np.conj(t1), t1, np.conj(t1), t1],
                               [0, 0, 0, 0, 0, 0],
                               [np.conj(t1), t1, np.conj(t1), t1, np.conj(t1), t1],
                               [0, 0, 0, 0, 0, 0]])
        # gamma_y
        t2_coupling = np.array([[np.conj(+t2), 0, np.conj(+t2), 0, np.conj(+t2), 0],
                                [t2, 0, t2, 0, t2, 0],
                                [np.conj(+t2), 0, np.conj(+t2), 0, np.conj(+t2), 0],
                                [t2, 0, t2, 0, t2, 0],
                                [np.conj(+t2), 0, np.conj(+t2), 0, np.conj(+t2), 0],
                                [t2, 0, t2, 0, t2, 0]])
        # lambda_x
        t3_coupling = np.array([[0, 0, 0, 0, 0, 0],
                                [np.conj(t3), np.conj(t3), np.conj(t3), np.conj(t3), np.conj(t3), np.conj(t3)],
                                [0, 0, 0, 0, 0, 0],
                                [np.conj(t3), np.conj(t3), np.conj(t3), np.conj(t3), np.conj(t3), np.conj(t3)],
                                [0, 0, 0, 0, 0, 0],
                                [np.conj(t3), np.conj(t3), np.conj(t3), np.conj(t3), np.conj(t3), np.conj(t3)]])
        # lambda_y
        t4_coupling = np.array([[0, np.conj(+t4), 0, np.conj(+t4), 0, np.conj(+t4)],
                                [0, t4, 0, t4, 0, t4],
                                [0, np.conj(+t4), 0, np.conj(+t4), 0, np.conj(+t4)],
                                [0, t4, 0, t4, 0, t4],
                                [0, np.conj(+t4), 0, np.conj(+t4), 0, np.conj(+t4)],
                                [0, t4, 0, t4, 0, t4]])

        self.add_coupling(t1_coupling, 0, creation, 0, annihilation, np.array([1, 0]), plus_hc=True, category='t1')
        self.add_coupling(t2_coupling, 0, creation, 0, annihilation, np.array([0, 1]), plus_hc=True, category='t2')
        self.add_coupling(t3_coupling, 0, creation, 0, annihilation, np.array([1, 0]), plus_hc=True, category='t3')
        self.add_coupling(t4_coupling, 0, creation, 0, annihilation, np.array([0, 1]), plus_hc=True, category='t4')


if __name__ == "__main__":

    model_params = dict(statistics='fermions', conserve='N', t1=1, t2=1, t3=0, t4=0, n=(int(1), int(4)),
                        LxMUC=6, Ly=6, bc_MPS='infinite', bc_x='periodic', bc_y='cylinder', order='Cstyle', verbose=1)
    M = BBHModel(model_params)

    print("max MPO bond dimension = ", max(M.H_MPO.chi))
    fg.plot_MPS_unit_cell(M)
