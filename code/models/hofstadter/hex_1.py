# --- TeNPy imports
from tenpy.models.lattice import Honeycomb
# --- infinite_cylinder imports
from models.hofstadter.hofstadter import HofstadterModel


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
        (creation, annihilation, nphi_default, t1, mu, V, nphi, nphi_2pi, Lx_MUC, phi_2pi) = \
            HofstadterModel.init_terms(self, params)

        self.chemical_potential(mu)

        hex1_int_flag = True if self.stats(params) == 'fermions' else False
        self.hex_1_hoppings(creation, annihilation, t1, V, nphi, nphi_2pi, Lx_MUC, phi_2pi,
                            interaction=hex1_int_flag)


if __name__ == "__main__":

    model_params = dict(statistics='fermions', conserve='N', t1=1, n=(int(1), int(9)), nphi=(int(1), int(3)),
                        Lx_MUC=1, Ly=6, V=10,
                        bc_MPS='infinite', bc_x='periodic', bc_y='cylinder', order='Cstyle',
                        verbose=1, phi=0)
    M = HofHex1Model(model_params)

    print("max MPO bond dimension = ", max(M.H_MPO.chi))

    # import matplotlib.pyplot as plt
    # ax = plt.gca()
    # M.lat.plot_sites(ax)
    # M.lat.plot_order(ax)
    # M.lat.plot_coupling(ax, M.lat.pairs["nearest_neighbors"])
    # plt.show()
