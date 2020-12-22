# --- TeNPy imports
from tenpy.models.model import CouplingMPOModel
from tenpy.networks.site import BosonSite, FermionSite


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

    @staticmethod
    def init_lattice_params(params):
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
        return Lx, Ly, order, bc_MPS, bc

    def init_terms(self, params):
        if stats(params) == 'bosons':
            creation, annihilation = 'Bd', 'B'
        else:
            creation, annihilation = 'Cd', 'C'
        Nmax = params.get('Nmax', 1)
        t1 = params.get('t1', 0.5)
        t2 = params.get('t2', 1)
        t3 = params.get('t3', 0.5)
        t4 = params.get('t4', 1)
        return creation, annihilation, Nmax, t1, t2, t3, t4
