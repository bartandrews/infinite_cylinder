# --- TeNPy imports
from tenpy.models.lattice import Honeycomb
# --- infinite_cylinder imports
from models.hofstadter.hofstadter import HofstadterModel
import functions.func_graph as fg


class HofHex1Model(HofstadterModel):

    def __init__(self, params):
        HofstadterModel.__init__(self, params)

    def init_sites(self, params):
        site = HofstadterModel.init_sites(self, params)
        return site

    def init_lattice(self, params):
        (Lx, Ly, order, bc_MPS, bc) = self.init_lattice_params(params)
        site = self.init_sites(params)
        lat = Honeycomb(Lx, Ly, site, order=order, bc_MPS=bc_MPS, bc=bc)
        return lat

    def init_terms(self, params):
        (creation, annihilation, nphi_default, t1, mu, V, Vtype, Vrange, nphi, nphi_2pi, LxMUC, phi_2pi) = \
            HofstadterModel.init_terms(self, params)

        self.chemical_potential(mu)
        self.hex_1_hoppings(creation, annihilation, t1, nphi, nphi_2pi, LxMUC, phi_2pi)
        self.offsite_interaction("Hex", V, Vtype, Vrange)


if __name__ == "__main__":

    model_params = dict(statistics='fermions', conserve='N', t1=1, n=(int(1), int(9)), nphi=(int(1), int(3)),
                        LxMUC=1, Ly=6, V=10, Vtype='Coulomb', Vrange=1,
                        bc_MPS='infinite', bc_x='periodic', bc_y='cylinder', order='Cstyle',
                        verbose=1, phi=0)
    M = HofHex1Model(model_params)

    print("max MPO bond dimension = ", max(M.H_MPO.chi))
    fg.plot_MPS_unit_cell(M)
