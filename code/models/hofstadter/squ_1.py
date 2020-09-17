# --- python imports
import csv
import numpy as np
import os
# --- TeNPy imports
from tenpy.models.lattice import Square
# --- infinite_cylinder imports
from models.hofstadter.hofstadter import HofstadterModel
import functions.func_graph as fg


class HofSqu1Model(HofstadterModel):

    def __init__(self, params):
        HofstadterModel.__init__(self, params)

    def init_sites(self, params):
        site = HofstadterModel.init_sites(self, params)
        return site

    def init_lattice(self, params):
        (Lx, Ly, order, bc_MPS, bc) = self.init_lattice_params(params)
        site = self.init_sites(params)
        lat = Square(Lx, Ly, site, order=order, bc_MPS=bc_MPS, bc=bc)
        return lat

    def init_terms(self, params):
        (creation, annihilation, nphi_default, Nmax, t1, mu, V, Vtype, Vrange, nphi, nphi_2pi, LxMUC, phi_2pi) = \
            HofstadterModel.init_terms(self, params)

        ###########
        # t_aniso #
        ###########

        t_aniso = 1
        r = params.get('r', None)
        if r is not None:
            # read in squ_1_hoppings_ratio.dat
            path = 'code/models/hofstadter'
            with open(os.path.join(path, 'squ_1_hoppings_ratio.dat'), 'r') as csvfile:
                converter = csv.reader(csvfile, delimiter='\t')
                tx_ty = []
                ratios = [[] for i in range(5)]
                for row in converter:
                    tx_ty.append(float(row[0]))
                    for i in range(5):
                        ratios[i].append(float(row[i + 1]))
            # determine chern number
            if nphi == [1, 4]:
                C = 1
            elif nphi == [4, 7]:
                C = 2
            elif nphi == [4, 11]:
                C = 3
            elif nphi == [4, 15]:
                C = 4
            elif nphi == [4, 19]:
                C = 5
            else:
                raise ValueError("Custom gap-to-width ratio is not implemented for this value of nphi.")
            # find closest corresponding hoppings
            lst = np.asarray(ratios[C - 1])
            idx = (np.abs(lst - r)).argmin()
            t_aniso = tx_ty[idx]
            print(f"Closest tx-to-ty ratio corresponding to desired gap-to-width ratio {r} is {tx_ty[idx]}.")

        self.chemical_potential(mu)
        self.squ_1_hoppings(creation, annihilation, t1, nphi, nphi_2pi, LxMUC, phi_2pi, t_anisotropy=t_aniso)
        self.offsite_interaction("Squ", Nmax, V, Vtype, Vrange)


if __name__ == "__main__":

    model_params = dict(statistics='fermions', Nmax=1, conserve='N', t1=1, n=(int(1), int(12)), nphi=(int(1), int(4)),
                        LxMUC=1, Ly=12, V=10, Vtype='Coulomb', Vrange=1,
                        bc_MPS='infinite', bc_x='periodic', bc_y='cylinder', order='Cstyle',
                        verbose=1, phi=0)
    M = HofSqu1Model(model_params)

    print("max MPO bond dimension = ", max(M.H_MPO.chi))
    fg.plot_MPS_unit_cell(M)
